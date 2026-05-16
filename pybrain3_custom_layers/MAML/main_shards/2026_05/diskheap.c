/* # Dylan Kenneth Eliot 
 *
 * diskheap.c - LD_PRELOAD disk-backed heap
 * All heap allocations are written to disk on free.
 * No locks, no TLS, no pthread. Re-entrancy via atomic flag.
 *
 * Usage:
 *   DISKHEAP_FILE=/path/to/file PYTHONMALLOC=malloc \
 *   LD_PRELOAD=/root/diskheap.so python3.10 script.py
 */
#define _GNU_SOURCE
#include <stdint.h>
#include <stddef.h>
#include <string.h>
#include <unistd.h>
#include <fcntl.h>
#include <sys/mman.h>
#include <stdlib.h>
#include <stdio.h>
#include <dlfcn.h>
#include <sys/syscall.h>

/* ── bootstrap bump allocator (no malloc needed during dlsym) ─────────────── */
static char   boot[131072];   /* 128KB static buffer */
static size_t boot_pos;

static void *boot_alloc(size_t n) {
    n = (n + 15) & ~(size_t)15;
    if (boot_pos + n > sizeof(boot)) return (void*)0;
    void *p = boot + boot_pos;
    boot_pos += n;
    return p;
}

static int in_boot(void *p) {
    return (char*)p >= boot && (char*)p < boot + sizeof(boot);
}

/* ── real allocator pointers ─────────────────────────────────────────────── */
static void *(*real_malloc)(size_t);
static void  (*real_free)(void*);
static void *(*real_realloc)(void*, size_t);
static void *(*real_calloc)(size_t, size_t);

/* ── disk file ───────────────────────────────────────────────────────────── */
static int g_fd = -1;

/* ── re-entrancy: per-thread via __thread, initialized to 0 by default ───── */
/* We avoid TLS malloc by using a static array indexed by tid % NTHREADS.    */
#define NTHREADS 64
static volatile int g_reent[NTHREADS];

static int reent_get(void) {
    pid_t tid = (pid_t)syscall(186); /* SYS_gettid aarch64 */
    return g_reent[(unsigned)tid % NTHREADS];
}
static void reent_set(int v) {
    pid_t tid = (pid_t)syscall(186);
    g_reent[(unsigned)tid % NTHREADS] = v;
}

/* ── init ────────────────────────────────────────────────────────────────── */
static int g_ready;

static void diskheap_init(void) __attribute__((constructor));
static void diskheap_init(void) {
    /* resolve real functions first using boot allocator for any
       internal mallocs dlsym may trigger */
    real_malloc  = dlsym(RTLD_NEXT, "malloc");
    real_free    = dlsym(RTLD_NEXT, "free");
    real_realloc = dlsym(RTLD_NEXT, "realloc");
    real_calloc  = dlsym(RTLD_NEXT, "calloc");

    const char *path = getenv("DISKHEAP_FILE");
    if (!path) path = "/tmp/diskheap.bin";

    g_fd = open(path, O_RDWR|O_CREAT|O_TRUNC, 0600);
    /* g_fd < 0 means disk backing unavailable; we'll fall through to real_malloc */

    g_ready = 1;
}

static void diskheap_fini(void) __attribute__((destructor));
static void diskheap_fini(void) {
    if (g_fd >= 0) {
        close(g_fd);
        const char *path = getenv("DISKHEAP_FILE");
        if (!path) path = "/tmp/diskheap.bin";
        unlink(path);
        g_fd = -1;
    }
}

/* magic to identify our allocations */
#define MAGIC 0xD15CD15CU

typedef struct {
    uint32_t magic;
    uint32_t pad;
    size_t   size;
    off_t    disk_off;
} Hdr;

#define HSIZ  sizeof(Hdr)

/* ── disk helpers ────────────────────────────────────────────────────────── */

/* append `size` bytes from `buf` to disk file, return offset */
static off_t disk_append(const void *buf, size_t size) {
    off_t off = lseek(g_fd, 0, SEEK_END);
    if (off < 0) return -1;
    size_t rem = size;
    const char *p = buf;
    while (rem) {
        ssize_t w = write(g_fd, p, rem);
        if (w <= 0) return -1;
        p += w; rem -= (size_t)w;
    }
    return off;
}

/* read `size` bytes from disk at `off` into `buf` */
static void disk_read(off_t off, void *buf, size_t size) {
    pread(g_fd, buf, size, off);
}

/* ── interceptors ────────────────────────────────────────────────────────── */

void *malloc(size_t size) {
    if (!real_malloc)
        return boot_alloc(size);
    if (reent_get())
        return real_malloc(size);
    reent_set(1);

    void *shadow = (void*)0;
    if (g_ready && g_fd >= 0) {
        shadow = real_malloc(HSIZ + size);
        if (shadow) {
            Hdr *h = shadow;
            h->magic    = MAGIC;
            h->size     = size;
            h->disk_off = -1;
            reent_set(0);
            return (char*)shadow + HSIZ;
        }
    }

    reent_set(0);
    return real_malloc(size);
}

void free(void *ptr) {
    if (!ptr) return;
    if (in_boot(ptr)) return;
    if (!real_free) return;
    if (!g_ready || g_fd < 0 || reent_get()) {
        real_free(ptr);
        return;
    }
    reent_set(1);
    Hdr *h = (Hdr*)((char*)ptr - HSIZ);
    if (h->magic != MAGIC) { reent_set(0); real_free(ptr); return; }
    disk_append(ptr, h->size);
    real_free(h);
    reent_set(0);
}

void *calloc(size_t nmemb, size_t size) {
    if (!real_calloc)
        return boot_alloc(nmemb * size);
    /* always bypass hook for calloc to avoid glibc arena deadlock */
    return real_calloc(nmemb, size);
}

void *realloc(void *ptr, size_t size) {
    if (!real_realloc) return malloc(size);
    if (!ptr) return malloc(size);
    if (!size) { free(ptr); return (void*)0; }
    if (in_boot(ptr)) return real_realloc(ptr, size);
    if (!g_ready || g_fd < 0 || reent_get())
        return real_realloc(ptr, size);

    Hdr *h = (Hdr*)((char*)ptr - HSIZ);
    if (h->magic != MAGIC) return real_realloc(ptr, size);
    size_t old = h->size;
    void *np = malloc(size);
    if (!np) return (void*)0;
    memcpy(np, ptr, old < size ? old : size);
    free(ptr);
    return np;
}
