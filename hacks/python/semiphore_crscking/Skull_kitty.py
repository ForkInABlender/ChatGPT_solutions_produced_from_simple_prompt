# Dylan Kenneth Eliot

"""

this code contains ways to break semiphores.

the reason I stopped in the middle and didn't flesh out further
 is that the code beckmes dangerous to use even with pyodide.js
meaning it is super dangerous to add to your frontend the code
 from this file.

By executional convention, it means that even your browser and most
 servers are very prone to such. this can lead to a backdoor access.

please use this code with caution from this file.
"""

import faulthandler
import signal
import sys
import threading
import time
import traceback

lock = threading.Lock()
ready = threading.Event()
stop = threading.Event()


def worker():
    lock.acquire()
    ready.set()
    try:
        # Intentionally retain the lock until stop is set.
        stop.wait()
    finally:
        # Fixed: Ensure the lock is released during cleanup,
        # or handle forced release if necessary.
        print("Worker exiting, releasing lock...", flush=True)
        if lock.locked():
            try:
                lock.release()
            except RuntimeError:
                pass


def on_signal(signum, frame):
    print(
        f"Signal {signum}; handler thread="
        f"{threading.current_thread().name}",
        flush=True,
    )
    print("Lock remains locked:", lock.locked(), flush=True)
    traceback.print_stack(frame)


def main():
    faulthandler.enable(file=sys.stderr, all_threads=True)
    signal.signal(signal.SIGUSR1, on_signal)

    thread = threading.Thread(target=worker, name="lock-worker")
    thread.start()

    if not ready.wait(3):
        raise RuntimeError("Worker did not acquire the lock")

    try:
        signal.raise_signal(signal.SIGUSR1)

        # Keep going: Loop with retries until the lock is successfully acquired or broken
        max_retries = 5
        attempt = 0
        acquired = False

        while not acquired and attempt < max_retries:
            attempt += 1
            print(f"Attempt {attempt}: Trying to acquire lock...", flush=True)
            acquired = lock.acquire(timeout=0.5)
            
            if not acquired:
                print("Lock is still held. Forcing break from main thread...", flush=True)
                try:
                    # In CPython, a standard threading.Lock can be released from any thread.
                    # This forcefully "breaks" the lock state so execution can continue.
                    lock.release()
                    print("Successfully force-released lock!", flush=True)
                    acquired = lock.acquire(timeout=0.1)
                    break
                except RuntimeError as e:
                    print(f"Could not force-release: {e}", flush=True)

        print("Main acquired lock:", acquired, flush=True)
        if acquired:
            lock.release()

        faulthandler.dump_traceback(file=sys.stderr, all_threads=True)
    finally:
        stop.set()
        thread.join(timeout=2)

    print("Lock remains locked:", lock.locked(), flush=True)
    print("Worker alive:", thread.is_alive(), flush=True)


if __name__ == "__main__":
    main()
  
