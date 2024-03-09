# Dylan Kenneth Eliot & GPT-4-Plugins ( Beta Edition )

"""
This is a simplified way to structure a fuse filesystem in a spreadsheet.

Now it can be used for ip based filesystem remoting without the ssh security model one normally would have to fight.

``soffice --accept="socket,host=localhost,port=2002;urp;" --norestore --nologo --nodefault`` to start the remote server. From here it should be configure for ip.



"""


#!/usr/bin/env python
import os
import stat
import errno
import fuse
import base64
import pyoo
from time import time

# Ensure FUSE API compatibility
fuse.fuse_python_api = (0, 2)

# Connect to LibreOffice
desktop = pyoo.Desktop('localhost', 2002)

# Attempt to open an existing spreadsheet or create a new one
try:
    doc = desktop.open_spreadsheet('filesystem.ods')
except FileNotFoundError:
    doc = desktop.create_spreadsheet()
    doc.save_as('filesystem.ods')

# Select the first sheet in the spreadsheet
sheet = doc.sheets[0]

# Spreadsheet utility functions
def ss_read_row(row):
    return [cell.value for cell in sheet[row, :2]]  # Assuming name and content are in columns A and B

def ss_write_row(row, name, content):
    sheet[row, 0].value, sheet[row, 1].value = name, content
    doc.save()

def ss_clear_row(row):
    sheet[row, 0].value, sheet[row, 1].value = "", ""
    doc.save()

# PyooFS class extending fuse.Fuse
class PyooFS(fuse.Fuse):
    def __init__(self, *args, **kwargs):
        super(PyooFS, self).__init__(*args, **kwargs)
        self.root = Directory("/")

    def getattr(self, path):
        node = self.root.find(path)
        if node is None:
            return -errno.ENOENT
        return node.stat

    def readdir(self, path, offset):
        node = self.root.find(path)
        if not isinstance(node, Directory):
            return -errno.ENOTDIR
        for item in node.contents:
            yield fuse.Direntry(item.name)

    def open(self, path, flags):
        node = self.root.find(path)
        if node is None:
            return -errno.ENOENT
        if isinstance(node, Directory):
            return -errno.EISDIR
        return 0  # Success

    def read(self, path, size, offset):
        node = self.root.find(path)
        if node is None:
            return -errno.ENOENT
        return node.content[offset:offset + size]

    def write(self, path, buf, offset):
        node = self.root.find(path)
        if node is None:
            return -errno.ENOENT
        node.content = node.content[:offset] + buf + node.content[offset + len(buf):]
        node.sync()
        return len(buf)

    def create(self, path, mode):
        parent_path, name = os.path.split(path)
        parent = self.root.find(parent_path)
        if parent is None or not isinstance(parent, Directory):
            return -errno.ENOENT
        if parent.find(name):
            return -errno.EEXIST
        new_file = parent.add(File(name))
        new_file.sync()
        return 0

    def unlink(self, path):
        parent_path, name = os.path.split(path)
        parent = self.root.find(parent_path)
        if parent and parent.remove(name):
            return 0
        return -errno.ENOENT

    def load_from_sheet(self):
        global sheet
        row = 0
        while True:
            data = ss_read_row(row)
            if not data[0]:  # Assuming the first column is the filename
                break
            path, content = data
            self.create(path, None)
            node = self.root.find(path)
            if node:
                node.content = base64.b64decode(content)
            row += 1

class Node:
    def __init__(self, name):
        self.name = name
        self.stat = self.default_stat()

    def default_stat(self):
        st = fuse.Stat()
        st.st_mode = stat.S_IFREG | 0o644
        st.st_nlink = 1
        st.st_size = 0
        st.st_mtime = st.st_atime = st.st_ctime = time()
        return st

    def sync(self):
        pass

class File(Node):
    def __init__(self, name, row=None):
        super(File, self).__init__(name)
        self.content = b''
        self.row = row  # Spreadsheet row where this file's data is stored
        self.stat.st_size = len(self.content)

    def sync(self):
        global sheet, doc
        if self.row is None:
            # Find the first empty row to store new file data
            self.row = self.find_empty_row()
        encoded_content = base64.b64encode(self.content).decode('utf-8')
        ss_write_row(self.row, self.name, encoded_content)
        doc.save()

    @staticmethod
    def find_empty_row():
        global sheet
        row = 0
        while True:
            # Assuming the first column contains the file name. Adjust as needed.
            if not sheet[row, 0].value:
                break
            row += 1
        return row

class Directory(Node):
    def __init__(self, name):
        super(Directory, self).__init__(name)
        self.stat.st_mode = stat.S_IFDIR | 0o755
        self.contents = []

    def add(self, node):
        self.contents.append(node)
        return node

    def remove(self, name):
        for i, item in enumerate(self.contents):
            if item.name == name:
                self.contents.pop(i)
                return True
        return False

    def find(self, path):
        if path == '/':
            return self
        path = path.lstrip('/')
        for content in self.contents:
            if content.name == path.split('/')[0]:
                if isinstance(content, Directory):
                    return content.find('/'.join(path.split('/')[1:]))
                return content
        return None

if __name__ == '__main__':
    fs = PyooFS()
    fs.load_from_sheet()
    fs.parse(errex=1)
    fs.main()
    fuse.Fuse(fs, "_borsh/", nothreads=True, foreground=True)
  
