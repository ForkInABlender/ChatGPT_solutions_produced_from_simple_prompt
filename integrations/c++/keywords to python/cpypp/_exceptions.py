#ignore
__all__ = [ "PytppError", "DirectiveError", "ExpandError", "ExpressionError", 
            "ManyEndIfsError", "UnclosedIfBlocksError", "MoreThanOneFileError",
            "OverwriteError", "NotFileOrPathWarn", "FileNotFoundError" ]

import os

def _fmt_message(message, e, linenum=None, filename=None):
     result = "Pytpp: " + e.__class__.__name__ + ": " + message
     if linenum: linenum = str(linenum)
     if linenum or filename:
        if filename:
            # path, name = os.path.split(filename)
            name = filename
            result += " (" + name
            if linenum:
                if "," in linenum: result += ", lines " + linenum
                else             : result += ", line " + linenum 
            result += ")"
        else:
            if "," in linenum: result += " (lines " + str(linenum) + ")"
            else             : result += " (line " + str(linenum) + ")"

     return result

class PytppError(Exception): pass

class DirectiveError(PytppError):
    def __init__(self, directive, linenum=None, filename=None):
        message = "invalid '{}'".format(directive)        
        PytppError.__init__(self, _fmt_message(message, self, linenum, filename))

class ExpandError(PytppError):
    def __init__(self, text, linenum=None, filename=None):
        message = text
        PytppError.__init__(self, _fmt_message(message, self, linenum, filename))

class ExpressionError(PytppError):
    def __init__(self, expr, exc_info, linenum=None, filename=None):
        message = "invalid '{}' expression ({})".format(expr, exc_info)
        PytppError.__init__(self, _fmt_message(message, self, linenum, filename))

class ManyEndIfsError(PytppError):
    def __init__(self, linenum=None, filename=None):
        message = "trying to remove more blocks than exists"
        PytppError.__init__(self, _fmt_message(message, self, linenum, filename))

class MoreThanOneFileError(PytppError):
    def __init__(self, linenum=None, filename=None):
        message = "only one file name must be used with -o option"
        PytppError.__init__(self, _fmt_message(message, self, linenum, filename))

class NotFileOrPathWarn(PytppError):
    def __init__(self, value, linenum=None, filename=None):
        message = "'{}' is not a file or a path".format(value)
        PytppError.__init__(self, _fmt_message(message, self, linenum, filename))

class FileNotFoundError(PytppError):
    def __init__(self, value, linenum=None, filename=None):
        message = "File {} not found".format(value)
        PytppError.__init__(self, _fmt_message(message, self, linenum, filename))

class UnclosedIfBlocksError(PytppError):
    def __init__(self, linenum=None, filename=None):
        message = "unclosed Ifblocks found"
        PytppError.__init__(self, _fmt_message(message, self, linenum, filename))

class OverwriteError(PytppError):
    def __init__(self, input_file, linenum=None, filename=None):
        message = "file '{}' will be overwritten. If you want this use -f flag".format(input_file)
        PytppError.__init__(self, _fmt_message(message, self, linenum, filename))
