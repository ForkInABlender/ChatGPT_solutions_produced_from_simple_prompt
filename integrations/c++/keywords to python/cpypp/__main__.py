#ignore
import os
import sys
import errno
import argparse
from ._cpypp import py_preprocessor, _PROG, __version__
from ._exceptions import *

class cpypp_command_line(object):

    def __init__(self):

        self._defines_dict = {}

    # eval expressions
    def _eval(self, expr, locals=None):
        if "__import__" in expr: raise ExpressionError(expr, "invalid __import__  usage")
        if "lambda"     in expr: raise ExpressionError(expr, "invalid lambda usage")
        try                  : return eval(expr, {}, {})
        except Exception as e: raise ExpressionError(expr, str(e))

    # -----------------------------------------------------------------------------------------

    def _find_files(self, dirname, extensions=".py", maxlevels=10):

        ext_list = extensions.split("|")

        if dirname.endswith(os.sep): dirname = dirname[:-1]
        # if dirname.startswith("." + os.sep) or dirname == ".": dirname = dirname[2:]
        minlevels = len(dirname.split(os.sep)) if dirname else 0
        last_cur = None

        for cur, _paths, files in os.walk(dirname or "."):
            if last_cur != cur:
               if self._quiet < 1: print("Listing '{}' ...".format(cur))
               last_cur = cur
            for file in files:
                valid_extension = any([ext for ext in ext_list if file.endswith(ext)])
                if valid_extension:
                   # if cur.startswith("." + os.sep) or cur == ".": cur = cur[2:]
                   curlevels = len(cur.split(os.sep)) if cur else 0
                   if  (curlevels - minlevels) <= maxlevels:
                       yield os.path.join(cur, file)

    # -----------------------------------------------------------------------------------------

    # process a batch of files
    def process_all(self):

        try:
            if self._defines:
               for item in self._defines:
                   name, value = item.split("=",1) if "=" in item else (item, True)
                   if not isinstance(value, bool): value = self._eval(value, {})
                   self._defines_dict.update({name:value})

            if len(self._files) > 1 and self._output_file:
                raise MoreThanOneFileError()

            if  self._compile: 
                self._remove_meta = self._force = self._output_file = None

            for file_or_path in self._files:

                if os.path.isfile(file_or_path):
                    self.process_file(file_or_path, self._output_file)
                elif os.path.isdir(file_or_path):
                    self.process_path(file_or_path)
                else: 
                    e = NotFileOrPathWarn(file_or_path)
                    if self._quiet < 2: print("Warning {}".format(str(e)))

        except PytppError as e:
            if self._quiet < 2: print(str(e))
        except Exception as e:
            if self._quiet < 2: print("Pytpp: {}: {}".format(e.__class__.__name__, str(e)))

    # -----------------------------------------------------------------------------------------

    def process_file(self, input_file, output_file=None):

        def cache_from_source(name, legacy):
            if sys.version_info < (3,0): legacy=True
            if legacy: return name + "c"
            parts = name.split(os.sep)
            path_parts, name_parts = parts[:-1], parts[-1].split(".")
            path_parts.append("__pycache__")
            if len(name_parts) < 2: name_parts.append("py")
            name_parts.insert(-1, "cpython-{}{}".format(sys.version_info.major, sys.version_info.minor))
            return os.path.join(os.path.join(*path_parts), ".".join(name_parts)) + "c"

        try:
            if  not self._compile:
                output_file = output_file or input_file[0:-len(input_file.split('.')[-1])-1] + ".py"
            else:
                output_file = cache_from_source(input_file, self._legacy)

            if self._path:

               if  output_file.startswith("." + os.sep): 
                   output_file = output_file[len(os.sep)+1:]
               output_file = os.path.join(self._path, output_file)
               path, file  = os.path.split(output_file)

               if  not os.path.isdir(path):
                  try : os.makedirs(path)
                  except OSError as exc :  # Python ≥ 2.5
                      if exc.errno != errno.EEXIST or not os.path.isdir(path): raise exc

            if  input_file == output_file and not self._force:
                raise OverwriteError(input_file)

            if output_file != "-":
               if self._compile:
                   if self._quiet < 1: print("Preprocessing and compiling '{}' ...".format(input_file))
               else:
                   if self._quiet < 1: print("Preprocessing '{}' to '{}' ...".format(input_file, output_file))

            py_preprocessor(input_file=input_file,
                                output_file=output_file, 
                                defines=self._defines_dict.copy(), 
                                remove_meta=self._remove_meta,
                                run=False, 
                                compile=self._compile,
                                legacy=self._legacy).parse()

        except PytppError as e:
            if self._quiet < 2: print(str(e))
        except Exception as e:
            if self._quiet < 2: print("Pytpp: {}: {}".format(e.__class__.__name__, str(e)))

    # -----------------------------------------------------------------------------------------

    def process_path(self, input_path):

        for input_file in self._find_files(input_path):
            self.process_file(input_file)

    # -----------------------------------------------------------------------------------------

    def parse_args(self):

        formatter_class=lambda prog: argparse.HelpFormatter (
                                        "<python2.7+> -m cpypp", 
                                     )

        parser = argparse.ArgumentParser(
                      prog=_PROG, 
                      description="c-style Python preprocessor",
                      formatter_class=formatter_class
                 )

        parser.add_argument("-v", "--version", dest="version", action="store_true", 
                            help=("print Pytpp version"))

        parser.add_argument("-d", dest="defines", action='append', metavar="name",
                            help='same as #define. Ex. -d debug or -d "var=2+2" (eval 4)')

        parser.add_argument("-e", dest="extensions", metavar="EXT",
                            help=("include files with only these extensions. default is '.py' and "
                                  "extensions must be separated with '|' char. "
                                  "Ex. '.py|.pypp'"))

        parser.add_argument("-l", type=int, default=999, dest='maxlevels',
                            help=("levels to recurse into subdirectories. Use '0' to don't recurse. "
                                  "Default is no limit"))

        parser.add_argument("-p", metavar='PATH',  dest='path', default=None,
                            help=("directory to prepend to file names and paths "
                                  "before save processed files. The full path will be "
                                  "created if it does not exists"))

        parser.add_argument("-f", dest="force", action="store_true", 
                            help=("force overwrite of files when output file name has the "
                                  "same name of input file name"))

        parser.add_argument("-r", dest="remove_meta", action="store_true",
                            help="remove meta tags and commented lines from final code")

        parser.add_argument("-o", dest="output_file", metavar="FILE",
                            help=("output file name when you are preprocessing just one file at once. "
                                  "Use '-o -' to stdout"))

        parser.add_argument("-q", action='count', dest="quiet", default=0,
                            help=("output only error messages; '-qq' will suppress "
                                  "the error messages as well"))

        parser.add_argument("-c", "--compileall" ,  dest="compile", action="store_true",
                            help=("compile each file after preprocessing. "
                                  "When this option is used, no preprocessed source file will be saved to disk "
                                  "and options '-o', '-r' and '-f' are discarded"))

        parser.add_argument("-b", dest="legacy", action="store_true",
                            help=("use legacy (pre-PEP3147) compiled file locations. Valid only when "
                                  "'-c' is used"))

        parser.add_argument("files", nargs="*", metavar="FILE|PATH",
                            help="one or more files or paths names")

        args = parser.parse_args()
        for k, v in args.__dict__.items():
            setattr(self, "_" + k, v)

        if self._version: print("{} vr {}, python vr {}".format( _PROG, __version__, 
                                                                 sys.version.split()[0]))
        return self

if __name__ == "__main__" : 
   try:
      import codecs
      if sys.version_info > (3,0): stderr, stdout = sys.stderr.buffer, sys.stdout.buffer
      else                       : stderr, stdout = sys.stderr, sys.stdout

      if sys.stdout.encoding != 'UTF-8': sys.stdout = codecs.getwriter('utf-8')(stdout, 'strict')
      if sys.stderr.encoding != 'UTF-8': sys.stderr = codecs.getwriter('utf-8')(stderr, 'strict')    
  
      cpypp_command_line().parse_args().process_all()
   except Exception as e:
      print(str(e))


