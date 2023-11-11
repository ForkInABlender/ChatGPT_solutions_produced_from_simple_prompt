# Dylan Kenneth Eliot & GPT-4-plugins (Alpha Edition)

"""
soffice --accept="pipe,name=my_pipe;urp;" --norestore --nologo --nodefault

use the above command to start a linux named pipe.

Then use this script. 

This can be updated so it uses an IP instead of a named pipe as well.

What this does is encode the file content and data into a savable format (base64) to compress and decompress data.

This is useful for when content changes or file accesses are addressed. Similar to inodes on NTFS/Fat/Fat32/XFat/EFI partitions are laid out.


"""


import pyoo
import os, stat, base64
import fuse

# Connect to LibreOffice Calc
desktop = pyoo.Desktop(pipe='my_pipe')
doc = desktop.open_spreadsheet('path_to_your_file.ods')
sheet = doc.sheets[0]  # Assuming we are working with the first sheet

# Functions to interact with LibreOffice Calc
def calc_get_cell(sheet, col, row):
    return sheet[row, col].value

def calc_write_cell(sheet, col, row, data):
    sheet[row, col].value = data

class DefaultStat(fuse.Stat):
    def __init__(self):
        self.st_mode = 0
        self.st_ino = 0
        self.st_dev = 0
        self.st_nlink = 0
        self.st_uid = 0
        self.st_gid = 0
        self.st_size = 0
        self.st_atime = 0
        self.st_mtime = 0
        self.st_ctime = 0

class File:
    def __init__(self, path, row):
        self.name = path
        self.stat = DefaultStat()
        self.row = row

class Directory:
    def __init__(self, name, row):
        self.contents = []
        self.stat = DefaultStat()
        self.name = name
        self.isdir = True
        self.row = row

    def add(self, item):
        self.contents.append(item)

    def get_item_at_path(self, path):
        target = os.path.basename(path)

        if target == ".." or target == ".":
            return Directory(path)

        for item in self.contents:
            if item.name == target:
                return item

        return None

class SpreadsheetFS(fuse.Fuse):
    def getattr(self, path):
        st = DefaultStat()

        if path == "/":
            st.st_mode = stat.S_IFDIR | 0o755
            st.st_nlink = 2
        else:
            target = top_level_dir.get_item_at_path(path)
            if target is not None:
                st.st_mode = stat.S_IFREG | 0o444
                st.st_nlink = 1
                st.st_size = len(self.read_file_contents(path))
            else:
                return -fuse.ENOENT
        return st

    def readdir(self, path, offset):
        for item in top_level_dir.contents:
            yield fuse.Direntry(item.name)

    def open(self, path, flags):
        target = top_level_dir.get_item_at_path(path)
        if target is None:
            return -fuse.ENOENT
        return 0

    def write(self, path, buf, offset):
        self.update_file_contents(path, buf, offset)
        return len(buf)

    def read(self, path, size, offset):
        target = top_level_dir.get_item_at_path(path)
        if target is None:
            return -fuse.ENOENT

        file_contents = self.read_file_contents(path)
        size = min(size, len(file_contents))
        offset = min(offset, len(file_contents))
        return file_contents[offset:offset + size]

    def unlink(self, path):
        self.delete_file(path)
        return 0

    def create(self, path, mode):
        self.create_file(path)
        return 0

    def truncate(self, path, length):
        target = top_level_dir.get_item_at_path(path)
        if target is None:
            return -fuse.ENOENT

        file_contents = self.read_file_contents(path)
        file_contents = file_contents[:length]

        self.update_file_contents(path, file_contents)
        return 0

    def read_file_contents(self, path):
        target = top_level_dir.get_item_at_path(path)
        contents = calc_get_cell(sheet, 'A', target.row)
        return base64.b64decode(contents)

    def update_file_contents(self, path, buf, offset=0):
        target = top_level_dir.get_item_at_path(path)

        # Get / unpack file contents
        curr_file_contents = self.read_file_contents(path)

        new_file_contents = curr_file_contents[:offset] + buf + curr_file_contents[offset + len(buf):]
        new_file_contents = base64.b64encode(new_file_contents)

        calc_write_cell(sheet, 'A', target.row, new_file_contents.decode('ascii'))

    def create_file(self, path):
        increment_num_files()
        calc_write_cell(sheet, 'B', tot_num_files, os.path.basename(path))
        top_level_dir.contents.append(File(os.path.basename(path), tot_num_files))

    def delete_file(self, path):
        target = top_level_dir.get_item_at_path(path)
        if target is None:
            return

        calc_write_cell(sheet, 'B', target.row, "")
        top_level_dir.contents.remove(target)

# Additional functions and global variables (like increment_num_files, top_level_dir, etc.) are needed for full functionality.

# Initialize the filesystem
top_level_dir = Directory("", 0)
tot_num_files = 0

def main():
    usage = """
    LibreOffice Calc Spreadsheet Filesystem

    """ + fuse.Fuse.fusage
    server = SpreadsheetFS(version="%prog " + fuse.__version__,
                           usage=usage,
                           dash_s_do="setsingle")
    server.parse(errex=1)
    server.main()

if __name__ == "__main__":
    main()
