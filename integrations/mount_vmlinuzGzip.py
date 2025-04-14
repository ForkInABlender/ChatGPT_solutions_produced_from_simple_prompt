# Dylan Kenneth Eliot

"""
This is the gzipped linux kernel to access your vmlinuz kernel itself. 

"""

import os
import gzip
import errno
from fuse import FUSE, Operations, FuseOSError

class KernelFS(Operations):
    def __init__(self, kernel_image_path):
        self.kernel_image_path = kernel_image_path
        print(f"[INIT] Mounting {self.kernel_image_path}...")

    def getattr(self, path, fh=None):
        print(f"[DEBUG] getattr called for: {path}")
        if path == '/':
            return {
                'st_mode': 0o755 | 0o040000,  # directory
                'st_nlink': 2
            }
        elif path == '/vmlinuz':
            file_size = os.path.getsize(self.kernel_image_path)
            return {
                'st_mode': 0o444 | 0o100000,  # read-only regular file
                'st_nlink': 1,
                'st_size': file_size
            }
        else:
            raise FuseOSError(errno.ENOENT)

    def readdir(self, path, fh):
        print(f"[DEBUG] readdir called for: {path}")
        if path == '/':
            return ['.', '..', 'vmlinuz']
        else:
            raise FuseOSError(errno.ENOENT)

    def open(self, path, flags):
        print(f"[DEBUG] open called for: {path}")
        if path != '/vmlinuz':
            raise FuseOSError(errno.ENOENT)
        return 0

    def read(self, path, size, offset, fh):
        print(f"[DEBUG] read called for: {path}, size: {size}, offset: {offset}")
        if path != '/vmlinuz':
            raise FuseOSError(errno.ENOENT)

        with open(self.kernel_image_path, 'rb') as f:
            f.seek(offset)
            return f.read(size)

def mount_kernel_image(kernel_image_path, mount_point):
    os.makedirs(mount_point, exist_ok=True)
    fuse = FUSE(KernelFS(kernel_image_path), mount_point, foreground=True, ro=True)

# Replace with your actual paths
if __name__ == "__main__":
    kernel_image_path = './vmlinuz-5.11.0-051100-generic'  # Update this
    mount_point = './kernel'           # Update this
    mount_kernel_image(kernel_image_path, mount_point)
