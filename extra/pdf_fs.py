# Dylan Kenneth Eliot

"""
Break down all pages into component images, and pages separately as file data, then treat as a singular mount point from a PDF document.

Now all PDFs are fuse'd files too :)

"""

import sys, os, errno, stat, fitz # (pip install fusepy PyMuPDF )
from fuse import FUSE, Operations, FuseOSError

class PDFFilesystem(Operations):
    def __init__(self, pdf_path):
        self.doc = fitz.open(pdf_path)
        self.files = {}
        self.images = {}
        self._build_file_structure()
    def _build_file_structure(self):
        for page_num in range(len(self.doc)):
            page_path = f"/page_{page_num + 1}.txt"
            self.files[page_path] = self.doc[page_num].get_text()
        image_index = 1
        for i in range(1, self.doc.xref_length()):
            info = self.doc.xref_object(i)
            if isinstance(info, bytes):
                info = info.decode('utf-8', errors='ignore')
            if info.startswith('<<') and '/Image' in info:
                try:
                    base_image = self.doc.extract_image(i)
                    image_bytes = base_image["image"]
                    ext = base_image["ext"]
                    image_path = f"/images/image_{image_index}.{ext}"
                    self.images[image_path] = image_bytes
                    image_index += 1
                except Exception:
                    continue  # Ignore broken images silently

    def readdir(self, path, fh):
        if path == '/':
            return ['.', '..'] + [name[1:] for name in self.files.keys()] + ['images']
        elif path == '/images':
            return ['.', '..'] + [os.path.basename(name) for name in self.images.keys()]
        else:
            return ['.', '..']

    def getattr(self, path, fh=None):
        if path == '/':
            return dict(st_mode=(stat.S_IFDIR | 0o755), st_nlink=2)
        if path == '/images':
            return dict(st_mode=(stat.S_IFDIR | 0o755), st_nlink=2)
        if path in self.files:
            return dict(st_mode=(stat.S_IFREG | 0o644), st_nlink=1, st_size=len(self.files[path].encode('utf-8')))
        if path.startswith('/images/'):
            if path in self.images:
                return dict(st_mode=(stat.S_IFREG | 0o644), st_nlink=1, st_size=len(self.images[path]))
        raise FuseOSError(errno.ENOENT)

    def open(self, path, flags):
        if path in self.files or path in self.images:
            return 0
        raise FuseOSError(errno.ENOENT)

    def read(self, path, size, offset, fh):
        if path in self.files:
            data = self.files[path].encode('utf-8')
            return data[offset:offset + size]
        if path in self.images:
            data = self.images[path]
            return data[offset:offset + size]
        raise FuseOSError(errno.ENOENT)

if len(sys.argv) != 3:
    print(f"Usage: {sys.argv[0]} <pdf_file> <mount_point>")
    exit(1)
pdf_path = sys.argv[1]
mount_point = sys.argv[2]
FUSE(PDFFilesystem(pdf_path), mount_point, nothreads=True, foreground=True)
