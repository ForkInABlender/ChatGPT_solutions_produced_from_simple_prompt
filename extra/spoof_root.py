# Dylan Kenneth Eliot

"""

sidestep normal execution points to imitate/mimick parts of real-root state.
"""

import os
import subprocess
import tempfile
import textwrap

c_code = textwrap.dedent(r"""
#define _GNU_SOURCE
#include <unistd.h>
#include <sys/types.h>
#include <grp.h>
#include <dlfcn.h>
#include <string.h>
#include <dirent.h>

// Override higher-level C library functions
DIR *(*real_opendir)(const char *) = NULL;
struct dirent *(*real_readdir)(DIR *) = NULL;

uid_t getuid(void) { return (uid_t)0; }
uid_t geteuid(void) { return (uid_t)0; }
gid_t getgid(void) { return (gid_t)0; }
gid_t getegid(void) { return (gid_t)0; }

int getresuid(uid_t *ruid, uid_t *euid, uid_t *suid) {
    if (ruid) *ruid = 0;
    if (euid) *euid = 0;
    if (suid) *suid = 0;
    return 0;
}

int getresgid(gid_t *rgid, gid_t *egid, gid_t *sgid) {
    if (rgid) *rgid = 0;
    if (egid) *egid = 0;
    if (sgid) *sgid = 0;
    return 0;
}

int getgroups(int size, gid_t list[]) {
    if (list && size > 0) {
        list[0] = 0;
    }
    return 1;
}

// Override opendir() to provide a "fake" directory handle
DIR *opendir(const char *name) {
    if (!real_opendir) {
        real_opendir = dlsym(RTLD_NEXT, "opendir");
    }
    if (strcmp(name, "/root") == 0) {
        // Return a valid directory stream for a harmless directory
        return real_opendir("/tmp");
    }
    return real_opendir(name);
}

// Override readdir() to provide fake entries for the /root directory
struct dirent *readdir(DIR *dirp) {
    static struct dirent fake_dirent;
    static int first_call = 1;

    if (!real_readdir) {
        real_readdir = dlsym(RTLD_NEXT, "readdir");
    }

    // A simple way to check if it's our fake dir is to use a static counter.
    // This is not foolproof but works for a simple ls.
    if (first_call) {
        first_call = 0;
        // Populate the fake dirent structure
        strcpy(fake_dirent.d_name, "authorized_keys");
        fake_dirent.d_ino = 1;
        return &fake_dirent;
    } else {
        first_call = 1;
        return NULL; // End of directory
    }
}
""")

results = {}

with tempfile.TemporaryDirectory() as tmpdir:
    so_path = os.path.join(tmpdir, "libfakeuid.so")
    c_path = os.path.join(tmpdir, "fakeuid.c")

    with open(c_path, "w") as f:
        f.write(c_code)

    try:
        subprocess.check_call(
            ["gcc", "-Wimplicit-function-declaration", "-shared", "-fPIC", "-o", so_path, c_path]
        )
    except Exception as e:
        results["compile_error"] = str(e)
    else:
        env = {**os.environ, "LD_PRELOAD": so_path}

        # Run `ls /root`
        try:
            results["ls /root"] = subprocess.check_output(
                ["whoami"], env=env, text=True
            ).strip()
        except Exception as e:
            results["ls /root"] = f"error: {e}"

print(results)
