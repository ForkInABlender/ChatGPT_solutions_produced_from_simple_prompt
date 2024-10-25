# Dylan Kenneth Eliot

"""
This is technically a usb device, but in doing so, it also occasionally grants root privileges and has other issues. 

Other than that, under proper usage, all file operations are completely normal.



"""

from fuse import FUSE, Operations
from collections import defaultdict
import errno
import stat
import time

class VirtualUSB(Operations):
    def __init__(self):
        self.files = defaultdict(lambda: {
            'st_mode': (stat.S_IFREG | 0o644),
            'st_nlink': 1,
            'st_size': 0,
            'st_ctime': time.time(),
            'st_mtime': time.time(),
            'st_atime': time.time(),
            'data': b''
        })
        self.files['/'] = {
            'st_mode': (stat.S_IFDIR | 0o755),
            'st_nlink': 2,
            'st_ctime': time.time(),
            'st_mtime': time.time(),
            'st_atime': time.time()
        }

    def getattr(self, path, fh=None):
        """Retrieve the attributes of a path."""
        if path not in self.files:
            raise OSError(errno.ENOENT, f"{path} not found")
        return self.files[path]

    def readdir(self, path, fh):
        """Read directory contents."""
        if path not in self.files or not stat.S_ISDIR(self.files[path]['st_mode']):
            raise OSError(errno.ENOENT, f"{path} not found or not a directory")
        return ['.', '..'] + [name[1:] for name in self.files if name != '/' and name.startswith(path)]

    def create(self, path, mode):
        """Create a new file."""
        self.files[path] = {
            'st_mode': (stat.S_IFREG | mode),
            'st_nlink': 1,
            'st_size': 0,
            'st_ctime': time.time(),
            'st_mtime': time.time(),
            'st_atime': time.time(),
            'data': b''
        }
        return 0

    def mkdir(self, path, mode):
        """Create a new directory."""
        self.files[path] = {
            'st_mode': (stat.S_IFDIR | mode),
            'st_nlink': 2,
            'st_ctime': time.time(),
            'st_mtime': time.time(),
            'st_atime': time.time()
        }
        parent = '/' if path.count('/') == 1 else path[:path.rfind('/')]
        if parent in self.files:
            self.files[parent]['st_nlink'] += 1

    def open(self, path, flags):
        """Open a file."""
        if path not in self.files:
            raise OSError(errno.ENOENT, f"{path} not found")
        return 0

    def write(self, path, data, offset, fh):
        """Write data to a file."""
        if path not in self.files:
            raise OSError(errno.ENOENT, f"{path} not found")

        file_data = self.files[path]['data']
        new_data = file_data[:offset] + data + file_data[offset + len(data):]
        self.files[path]['data'] = new_data
        self.files[path]['st_size'] = len(new_data)
        self.files[path]['st_mtime'] = time.time()
        return len(data)

    def read(self, path, size, offset, fh):
        """Read data from a file."""
        if path not in self.files:
            raise OSError(errno.ENOENT, f"{path} not found")
        file_data = self.files[path]['data']
        return file_data[offset:offset + size]

    def unlink(self, path):
        """Remove a file."""
        if path in self.files:
            del self.files[path]
        else:
            raise OSError(errno.ENOENT, f"{path} not found")

    def rmdir(self, path):
        """Remove a directory."""
        if path in self.files and stat.S_ISDIR(self.files[path]['st_mode']):
            if any(p.startswith(path) and p != path for p in self.files):
                raise OSError(errno.ENOTEMPTY, f"{path} is not empty")
            del self.files[path]
        else:
            raise OSError(errno.ENOENT, f"{path} not found or not a directory")

    def setattr(self, path, attr, fh=None):
        """Set attributes (e.g., permissions, size)."""
        if path not in self.files:
            raise OSError(errno.ENOENT, f"{path} not found")
        for key, value in attr.items():
            if key in self.files[path]:
                self.files[path][key] = value
        return 0

    def truncate(self, path, length, fh=None):
        """Truncate a file to a specific length."""
        if path not in self.files:
            raise OSError(errno.ENOENT, f"{path} not found")
        data = self.files[path]['data']
        if length < len(data):
            self.files[path]['data'] = data[:length]
        else:
            self.files[path]['data'] = data + b'\x00' * (length - len(data))
        self.files[path]['st_size'] = length
        self.files[path]['st_mtime'] = time.time()

    def flush(self, path, fh):
        """Flush any cached information to storage."""
        # No actual storage backend, so just a noop.
        return 0

def main(mountpoint):
    """Mount the virtual USB filesystem."""
    FUSE(
        VirtualUSB(),
        mountpoint,
        foreground=True,
        nothreads=True,
        debug=True,
        ro=False  # Ensure read-write mode
    )

if __name__ == '__main__':
    mount_point = 'mount_usb'
    main(mount_point)
