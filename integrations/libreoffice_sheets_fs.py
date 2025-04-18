# Dylan Kenneth Eliot
"""
first, you run:

`soffice "--accept=socket,host=localhost,port=2002;urp;StarOffice.ComponentContext" --norestore --nologo --headless &`

then you run this script. Make sure libreoffice is installed before doing so.

This module assumes you use python3.8 with fusepy & pyoo installed as well.

Each cell from A1 to Z100 has 4 files: ['value', 'formula', 'type', 'formatted']

Each "file" is read, write accessible of "Sheet1/{cell-id}". This means you have effectively an efficient way to store information and retrieve it if used
 [correctly].


"""


import os
import stat
import pyoo
from errno import ENOENT
from fuse import FUSE, Operations, FuseOSError

class SpreadsheetFS(Operations):
    def __init__(self, file_path):
        self.file_path = file_path
        ctx = pyoo.Desktop('localhost', 2002)
        self.doc = ctx.open_spreadsheet(file_path)
        self.sheet_names = {sheet.name for sheet in self.doc.sheets}
        self.buffers = {}  # Cache for in-progress writes

    def readdir(self, path, fh):
        parts = path.strip("/").split("/")
        if path == "/":
            return ['.', '..'] + list(self.sheet_names)
        elif len(parts) == 1 and parts[0] in self.sheet_names:
            sheet = self.doc.sheets[parts[0]]
            entries = ['.', '..']
            for row in range(100):
                for col in range(26):
                    if sheet[row, col].value is not None:
                        cell_name = f"{chr(65 + col)}{row + 1}"
                        entries.append(cell_name)
            return entries
        elif len(parts) == 2 and parts[0] in self.sheet_names:
            return ['.', '..', 'value', 'formula', 'type', 'formatted']
        raise FuseOSError(ENOENT)

    def getattr(self, path, fh=None):
        parts = path.strip("/").split("/")
        try:
            if path == "/":
                return {'st_mode': (stat.S_IFDIR | 0o755), 'st_nlink': 2}
            if len(parts) == 1 and parts[0] in self.sheet_names:
                return {'st_mode': (stat.S_IFDIR | 0o755), 'st_nlink': 2}
            if len(parts) == 2:
                # It's a cell directory like /Sheet1/A1
                return {'st_mode': (stat.S_IFDIR | 0o755), 'st_nlink': 2}
            if len(parts) == 3:
                sheet_name, cell_name, attr = parts
                if sheet_name not in self.sheet_names:
                    raise FuseOSError(ENOENT)
                if not cell_name or len(cell_name) < 2:
                    raise FuseOSError(ENOENT)
                col = ord(cell_name[0].upper()) - 65
                row = int(cell_name[1:]) - 1
                sheet = self.doc.sheets[sheet_name]
                if not (0 <= row < 100 and 0 <= col < 26):
                    raise FuseOSError(ENOENT)
                if attr not in {"value", "formula", "type", "formatted"}:
                    raise FuseOSError(ENOENT)
                value = str(self._get_cell_attr(sheet[row, col], attr))
                # Make cell attributes writable (0o666)
                return {
                    'st_mode': (stat.S_IFREG | 0o666),
                    'st_nlink': 1,
                    'st_size': len(value.encode())
                }
            raise FuseOSError(ENOENT)
        except Exception as e:
            print(f"getattr error on path={path}: {e}")
            raise FuseOSError(ENOENT)

    def open(self, path, flags):
        return 0

    def read(self, path, size, offset, fh):
        parts = path.strip("/").split("/")
        try:
            if len(parts) == 3:
                sheet_name, cell_name, attr = parts
                col = ord(cell_name[0].upper()) - 65
                row = int(cell_name[1:]) - 1
                sheet = self.doc.sheets[sheet_name]
                content = str(self._get_cell_attr(sheet[row, col], attr)).encode()
                return content[offset:offset + size]
        except Exception as e:
            print(f"read error: {e}")
            return b''

    def write(self, path, data, offset, fh):
        parts = path.strip("/").split("/")
        try:
            if len(parts) == 3:
                sheet_name, cell_name, attr = parts
                if sheet_name not in self.sheet_names:
                    raise FuseOSError(ENOENT)
                if not cell_name or len(cell_name) < 2:
                    raise FuseOSError(ENOENT)
                col = ord(cell_name[0].upper()) - 65
                row = int(cell_name[1:]) - 1
                sheet = self.doc.sheets[sheet_name]
                if not (0 <= row < 100 and 0 <= col < 26):
                    raise FuseOSError(ENOENT)
                cell = sheet[row, col]
                
                # Use the buffer logic from t8.py
                key = path
                if key not in self.buffers:
                    self.buffers[key] = bytearray()
                buf = self.buffers[key]

                # Expand buffer if needed
                while len(buf) < offset:
                    buf += b'\x00'
                
                buf[offset:offset + len(data)] = data
                self.flush(path, fh)  # Immediately flush after writing
                return len(data)
            else:
                raise FuseOSError(ENOENT)
        except Exception as e:
            print(f"write error: {e}")
            raise FuseOSError(ENOENT)

    def flush(self, path, fh):
        if path not in self.buffers:
            return 0

        parts = path.strip("/").split("/")
        sheet_name, cell_name, attr = parts
        col = ord(cell_name[0].upper()) - 65
        row = int(cell_name[1:]) - 1

        try:
            value = self.buffers[path].decode().strip()
            sheet = self.doc.sheets[sheet_name]
            cell = sheet[row, col]

            if attr == "value":
                cell.value = value
            elif attr == "formula":
                cell.formula = value
            else:
                raise FuseOSError(ENOENT)

            self.doc.store()  # Save to file
            print(f"Updated {sheet_name}!{cell_name} ({attr}) = {value}")

        except Exception as e:
            print(f"flush() error: {e}")
        finally:
            del self.buffers[path]
        return 0

    def truncate(self, path, length):
        self.buffers[path] = bytearray(length)
        return 0

    def create(self, path, mode, fi=None):
        return 0  # Allow overwrite on open

    def _get_cell_attr(self, cell, attr):
        try:
            if attr == "value":
                return cell.value if cell.value is not None else ""
            elif attr == "formula":
                return cell.formula if cell.formula is not None else ""
            elif attr == "type":
                return getattr(cell, "Type", "") or ""
            elif attr == "formatted":
                return cell.string if hasattr(cell, "string") else ""
        except Exception as e:
            print(f"Error fetching attribute '{attr}' from cell: {e}")
            return ""

if __name__ == '__main__':
    import sys
    if len(sys.argv) != 3:
        print("Usage: python3 spreadsheet_fuse.py ")
        exit(1)
    spreadsheet_path = sys.argv[1]
    mount_point = sys.argv[2]
    os.makedirs(mount_point, exist_ok=True)
    FUSE(SpreadsheetFS(spreadsheet_path), mount_point, nothreads=True, foreground=True, rw=True, default_permissions=True)
