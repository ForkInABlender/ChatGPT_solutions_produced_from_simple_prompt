#ignore
import io
import os
import ast
import imp
import sys
import traceback
import py_compile

from ._exceptions import *
from . import *

__version__="1.0.2"
__author__="Wellington Rats"
__email__="wellrats@gmail.com"

_PROG = "cpypp"

class py_preprocessor(object):

    _reserved_words = ["False", "None", "True", "and", "as", "assert", "async", "await", "break", "class", 
                       "continue", "def", "del", "elif", "else", "except", "finally", "for", "from", "global", 
                       "if", "import", "in", "is", "lambda", "nonlocal", "not", "or", "pass", "raise", "return", 
                       "try", "while", "with", "yield"]

    _directives = ["define", "undef", "exclude", "endexclude", "ignore", "endignore",
                   "if", "ifndef", "ifdef", "else", "elif", "elifdef", 
                   "endif", "endififdef", "endifif", "endifall",
                   "expand", "include", "includeident"]

    _parse_params = [ "defines", "remove_meta",  "run", "compile", "legacy", "expand_open", "expand_close", 
                      "expand_escape", "escape", "comment" ]

    _pypp_dont_call_sys_exit = "__PYPP_DONT_CALL_SYS_EXIT__"
    _pypp_dont_expand_names  = "__PYPP_DONT_EXPAND_NAMES__"
    _pypp_dont_expand_exprs  = "__PYPP_DONT_EXPAND_EXPRS__"

    _parse_defaults = dict(defines={}, remove_meta=True, run=True, expand_open="#{", expand_close="}#", 
                           escape="#", comment="# ", expand_escape="!")

    def __init__(self, input_file=sys.argv[0],
                       module_name="__main__",
                       output_file=None, **kwargs): 


        for arg in py_preprocessor._parse_params:
            setattr(self, arg, kwargs.get(arg, py_preprocessor._parse_defaults.get(arg)))

        # public variables
        self.input_file = input_file
        if self.input_file.endswith(".pyc"): self.input_file = self.input_file[:-1]
        self.output_file = output_file
        self.module_name = module_name or "__main__"
        self.read_encoding = sys.stdin.encoding
        self.write_encoding = sys.stdout.encoding

        self._reset_internal()

        # declaring directive consts
        for directive in py_preprocessor._directives:
            setattr(self, "_" + directive.upper(), self.escape + directive)

    # reseting internal things to parse a second file

    def _reset_internal(self):

        self._batch_mode = not self.run
        self._run_mode   = self.run
        self._linenum = 0
        self._excludeblock = False
        self._ignoreblock = False
        self._ifblocks = []
        self._output_buffer = u""
        self._expand_open_len = len(self.expand_open)
        self._expand_close_len = len(self.expand_close)
        self._already_parsed = False
        self._ident_step = 0
        self._did_something = False


        self.defines.update(__VERSION__=sys.version.split()[0])
        if sys.version_info >= (3,0): self.defines.update(dict(__PYTHON3__=True))
        else                        : self.defines.update(dict(__PYTHON2__=True))

        if  self._batch_mode: self.defines.update(dict(__MODEBAT__=True, __MODE__ = "batch"))
        else                : self.defines.update(dict(__MODERUN__ = True, __MODE__="run"))

    # the #define directive
    def define(self, name, value=True): 
        self.defines[name.strip()] = value

    # the #undef directive
    def undef(self, name): self.defines.pop(name.strip(), None)
        
    # search: if define is defined
    def defined(self, name): return name.strip() in self.defines.keys()

    # faced to avoid from executing in first step of compiler
    def parsed(self): return False

    # parsing/processing
    def parse(self, input_file=None, module_name=None, output_file=None, **kwargs):

        if self._already_parsed: return False
        self._already_parsed = True

        for arg in py_preprocessor._parse_params:
            if arg in kwargs: setattr(self, arg, kwargs.get(arg))
       
        if input_file  : self.input_file  = input_file
        if module_name : self.module_name = module_name
        if output_file : self.output_file = output_file

        if not self.input_file: 
            self.input_file = self.module_name.replace(".", os.sep) + ".py"
        if self.input_file.endswith(".pyc") : 
            self.input_file = self.input_file[:-1]

        self._reset_internal()
        ident_adjust = ident_saved = 0

        # Let's do magic
        if not os.path.exists(self.input_file): return False

        try:
            with io.open(os.path.join(self.input_file), "r", encoding=self.read_encoding) as input_file:
                self._lines = input_file.read().splitlines()

            last_line_was_empty = False
            # process the input file
            for line in self._lines:
                self._linenum += 1

                # to supress or not to supress
                result = self._lexer(line)
                # while len(result) < 4: result += (None,)
                comment, ignore, replace, ident = result
                if comment or ignore or replace: self._did_something = True
                if replace is not None: line = replace
                line = line.strip()

                # check if needs to reident line
                if line and self._ident_step > 1:

                    if self._ident_step in (2,3) and ident_saved >= ident:
                       self._ident_step = 0

                    if self._ident_step == 2: 
                       self._ident_step = 3
                       ident_adjust = (ident - ident_saved)

                    if self._ident_step == 3: 
                       ident -= ident_adjust

                spcs = " " * ident

                # process and output
                if ignore:

                    pass

                elif self.remove_meta and comment: 

                    pass

                elif not line and self.remove_meta: 

                    if not last_line_was_empty:
                        self._output_buffer += "\n"
                        last_line_was_empty = True

                elif comment:

                    last_line_was_empty = False
                    if line.startswith("# "): self._output_buffer += spcs + line + "\n"
                    else                    : self._output_buffer += spcs + self.comment + line + "\n"

                else:
                    
                    last_line_was_empty = False
                    if not self._ignoreblock: line = self._expander(line)
                    self._output_buffer += spcs + line + "\n"

                if self._ident_step == 1: ident_saved, self._ident_step = ident, 2

        finally:                
            #Warnings for unclosed #ifdef blocks
            if self._ifblocks:
                linenums = ",".join([ str(block.get("linenum")) for block in self._ifblocks ])
                raise UnclosedIfBlocksError(linenums, self.input_file)

        if self.compile : self._post_process_compile()
        elif self.run   : self._post_process_run()
        else            : self._post_process_norun()

    # evaluate
    def _lexer(self, original_line):
    # return values are (comment, ignore, replace, ident) (bool, bool, str)
        # strip and count spaces before first non space char

        ident, len_line = 0, len(original_line)
        if  original_line.startswith(" "):
            while (ident < len_line) and (original_line[ident] == " "): ident += 1
        line = (original_line[ident:] if ident else original_line).strip()

        # optimization
        # if not (self._ifblocks or self._excludeblock or self._ignoreblock) and not line.startswith(self.escape):
        #      return False, False, None, ident

        directive, compl = line.split(None, 1) if " " in line else (line, "")

        # handle #endignore directives

        # if self._ignoreblock is True then leave while not parse a #endignore
        if self._ignoreblock and directive != self._ENDIGNORE:
            return False, False, None, ident

        # ignore parse calls in final code
        elif "py_preprocessor" in line: 
            return True, False, None, ident

        # match PYPP.parse( and PYPP.parsed(
        elif "PYPP.parse" in line and not line.startswith("#"):
            if line.startswith("if ") : self._ident_step = 1
            return True, False, None, ident

        elif "False and False" in line and line.startswith("if "):
            self._ident_step = 1
            return True, False, None, ident

        # handle #define directives
        elif directive == self._ENDIGNORE:

            self._test_is_empty(directive, compl)
            self._ignoreblock = False
            return True, False, None, ident

        elif directive == self._DEFINE:

            name, expr = compl.split(None, 1) if " " in compl else (compl, None)
            self._test_is_alnum(directive, name)            
            if expr : value = self._eval(expr)
            else    : value = self.defined(name)

            self.define(name, value)
            return True, False, None, ident

        # handle #undef directives
        elif directive == self._UNDEF:

            self._test_is_alnum(directive, compl)
            self.undef(compl)
            return True, False, None, ident

        # handle #exclude directives
        elif directive == self._EXCLUDE:

            self._test_is_empty(directive, compl)
            self._excludeblock = True
            return True, False, None, ident

        # handle #endexclude directives
        elif directive == self._ENDEXCLUDE:

            self._test_is_empty(directive, compl)
            self._excludeblock = False
            return True, False, None, ident

        elif directive == self._IGNORE:

            self._test_is_empty(directive, compl)
            self._ignoreblock = True
            return True, False, None, ident

        # handle #if directives 
        elif directive == self._IF:

            self._test_not_empty(directive, compl)
            name="block_if_{}".format(self._linenum)
            self._ifblocks_append(name, value=self._eval(compl))
            return True, False, None, ident

        # handle #ifnotdef directives (is the same as: #ifdef X #else)
        elif directive == self._IFNDEF:

            self._test_is_alnum(directive, compl)
            self._ifblocks_append(compl, value=not self.defined(compl))
            return True, False, None, ident

        # handle #ifdef directives
        elif directive == self._IFDEF:

            self._test_is_alnum(directive, compl)
            self._ifblocks_append(compl)
            return True, False, None, ident

        # handle #elseifdef directives
        elif directive == self._ELIFDEF:

            self._test_is_alnum(directive, compl)
            self._test_ifblocks(directive)
            self._ifblocks_replace_last(compl, value=self.defined(compl))
            return True, False, None, ident

        # handle #elseifdef directives
        elif directive == self._ELIF:

            self._test_is_not_empty(directive, compl)
            self._test_ifblocks(directive)
            self._ifblocks_replace_last(None, value=self._eval(compl))
            return True, False, None, ident

        # handle #else directives
        elif directive == self._ELSE:

            self._test_is_empty(directive, compl)
            self._test_ifblocks(directive)
            self._ifblocks_invert_last()
            return True, False, None, ident

        # handle #endififdef directives
        elif directive == self._ENDIFIFDEF:

            self._test_is_alnum(directive, compl)
            if len(self._ifblocks) : self._ifblocks.pop(-1)
            else                    : raise ManyEndIfsError(self._linenum, self.input_file)
            self._ifblocks_append(compl)
            return True, False, None, ident

        #handle #endifif directives
        elif directive == self._ENDIFIF:

            self._test_not_empty(directive, compl)
            if len(self._ifblocks) : self._ifblocks.pop(-1)
            else                    : raise ManyEndIfsError(self._linenum, self.input_file)
            name="block_if_{}".format(self._linenum)
            self._ifblocks_append(name, value=self._eval(compl))
            return True, False, None, ident

        # handle #endifall directives
        elif directive == self._ENDIFALL:

            self._test_is_empty(directive, compl)
            self._ifblocks = []
            return True, False, None, ident

        # handle #endif and #endif numb directives
        elif directive == self._ENDIF:

            if compl: self._test_is_digit(directive, compl)
            try     : number = int(compl)
            except  : number = 1

            for i in range(0, number): 
                if  len(self._ifblocks): self._ifblocks.pop(-1)
                else                    : raise ManyEndIfsError(self._linenum, self.input_file)

            return True, False, None, ident

        elif directive == self._INCLUDE:

            self._test_not_empty(directive, compl)
            file_name, i = self._find_file(compl), self._linenum
            if not file_name: raise FileNotFoundError(compl, self._linenum, self.input_file)
            self._lines[i:i] = open(file_name).read().splitlines()

            return True, False, None, ident

        elif directive == self._INCLUDEIDENT:

            self._test_not_empty(directive, compl)
            file_name, spcs, i = self._find_file(compl), " " * ident, self._linenum
            if not file_name: raise FileNotFoundError(compl, self._linenum, self.input_file)
            self._lines[i:i] = [ (spcs + line) for line in open(file_name).read().splitlines() ]

            return True, False, None, ident

        # handle #expand directive
        # protects from compile errors
        elif directive == self._EXPAND:

             return self._ifblocks_not_true(), False, compl, ident

        else: #No directive --> execute

            # process the excludeblock
            if self._excludeblock is True:
                return True, False, None, ident

            # process the ifblock
            elif self._ifblocks: # is True:
                return self._ifblocks_not_true(), False, None, ident

            #here can add other stuff for deleting comnments eg
            else:
                return False, False, None, ident

    # expand #{ ..expr ...}# to values
    def _expander(self, line):

        def find_token(line, start):
          
            pos, len_line, collecting, escape, lastchar, token = start, len(line), False, False, None, ""
            while (pos < len_line):
                c = line[pos]
                if  collecting:
                    if (c.isalnum() or c == "_") : token += c
                    else                         : return (start - (1 if escape else 0), pos)
                else:
                    if (c.isalpha() or c == "_"): 
                        if lastchar == self.expand_escape: escape = True
                        collecting, start, token = True, pos, c
                lastchar=c
                pos += 1

            return -1, -1

        # Expand #{ .. expressions }#
        if not bool(self.defines.get(py_preprocessor._pypp_dont_expand_exprs)):

            start, end, interations = -1, 0, 0
            while True:

                start = line.find(self.expand_open, end)
                end   = line.find(self.expand_close, end)
                if start < 0 and end < 0:
                    break
                elif start >= 0 and end < 0:
                    raise ExpandError("missing end expander '{}'".format(self.expand_close), 
                                      self._linenum, self.input_file)
                if start < 0 and end >= 0: 
                    raise ExpandError("missing begin expander '{}'".format(self.expand_open), 
                                      self._linenum, self.input_file)

                self._did_something = True

                expr     = line[ start + self._expand_open_len : end]  ## start is before KEY_OPEN
                end     += self._expand_close_len                      ##end now is after KEY_CLOSE
                removed  = (end - start)
                func = str
                if   expr.endswith(",s"): func = str;  expr = expr[:-2]
                elif expr.endswith(",r"): func = repr; expr = expr[:-2]

                value    = func(self._eval(expr))
                inserted = len(value)

                line     = line[:start] + value + line[end:]
                end     += (inserted - removed)


        # Expand names
        if not bool(self.defines.get(py_preprocessor._pypp_dont_expand_names)):

            start, end, interations = -1, 0, 0
            while True:

                start, end = find_token(line, end)
                if start < 0 and end < 0: break
                token = line[start:end]
                if token.startswith(self.expand_escape): escape = True; token = token[1:]
                else                                   : escape = False

                if token in self.defines.keys() and not token in py_preprocessor._reserved_words:

                    self._did_something = True

                    if escape: value = token; removed = 1
                    else     : value = repr(self.defines.get(token)); removed  = end - start
                    inserted = len(value)

                    line = line[:start] + str(value) + line[end:]
                    end = start + (inserted - removed)

            return line

    # eval expressions
    def _eval(self, expr, locals=None):
        if "__import__" in expr: 
            raise ExpressionError(expr, "invalid __import__ usage", self._linenum, self.input_file)
        if "lambda" in expr: 
            raise ExpressionError(expr, "invalid lambda usage", self._linenum, self.input_file)
        try : 
            if locals is None: locals = self.defines
            return eval(expr, {}, locals)
        except Exception as e:
            raise ExpressionError(expr, str(e), self._linenum, self.input_file)

    #returning: validness of #ifdef #else block
    # if any nested condition is false so supress line
    def _ifblocks_not_true(self): 
        for cond in self._ifblocks:
            if not cond.get("value"): return True
        return False

    #append a new name to an if_block
    def _ifblocks_append(self, name, value=None):
        if value is None: value = self.defined(name)
        attrs = dict(value=value, name=name, linenum=self._linenum, already_true=value)
        self._ifblocks.append(attrs)

    # Replace the condition of last block (if true only if was not true before)
    # for work with #elif directives
    def _ifblocks_replace_last(self, name=None, value=None):

        attrs = self._ifblocks[-1]
        if name:
            attrs["name"] = name
            if value == None: value = self.defined(name)

        if value and attrs["already_true"]: value = False
        attrs["value"] = value

    # invert condition of last block. For #else directives
    def _ifblocks_invert_last(self):
        return self._ifblocks_replace_last(value=not self._ifblocks[-1]["value"])

    # quick test functions
    def _test_is_alnum(self, directive, name):
        if name in py_preprocessor._reserved_words: 
            self._directive_error(directive)
        else:
            try    : ast.parse("{}=None".format(name))
            except : self._directive_error(directive)

    _test_is_digit  = lambda self, d, name : self._directive_error(d) if not name.isdigit() else None
    _test_is_empty  = lambda self, d, name : self._directive_error(d) if name else None
    _test_not_empty = lambda self, d, name : self._directive_error(d) if not name else None
    _test_ifblocks  = lambda self, d    : self._directive_error(s) if not self._ifblocks else None

    # look for a file name using sys.path or not
    def _find_file(self, compl):
        if compl.startswith("<") and compl.endswith(">"):
            file_name = compl[1:-1]; use_sys_path=True
        else: 
            file_name = self._eval(compl); use_sys_path=False

        if use_sys_path:
            for path in sys.path:
                full_file_name = os.path.join(path, file_name)
                if os.path.isfile(full_file_name): return full_file_name

        if os.path.isfile(file_name): return file_name
        # if os.path.isfile(os.path.join(".", file_name)): return os.path.join(".", file_name)

        return None

    # error handling
    def _directive_error(self, directive):
        raise DirectiveError(directive, self._linenum, self.input_file)

    # redo traceback
    def _rewrite_traceback(self):
        trace = traceback.format_exc().splitlines()
        index = 0
        for line in trace:
            if index == (len(trace) - 2):
                print(line.replace("<string>", self.input_file))
            else:
                print(line)
            index += 1

    # post-processor compile
    def _post_process_compile(self):

        try:

            parts = self.input_file.rsplit(".", 1)
            if self._did_something: 
                tmp_output_file = parts[0] + "__." + parts[1]
                with io.open(tmp_output_file, "w", encoding=self.write_encoding) as f:
                     f.write(self._output_buffer)
            else: 
                tmp_output_file = self.input_file

            py_compile.compile(tmp_output_file, self.output_file, self.input_file)

        except Exception as e:
            print(str(e))
            exc_type, exc_value, exc_traceback = sys.exc_info()
            traceback.print_tb(exc_traceback)

        finally:
            if tmp_output_file != self.input_file and os.path.isfile(tmp_output_file): 
                os.remove(tmp_output_file)

    # post-processor with run
    def _post_process_run(self):
        parts = self.input_file.rsplit(".", 1)
        self.output_file = parts[0] + "__." + parts[1]

        with io.open(self.output_file, "w", encoding=self.write_encoding) as f:
             f.write(self._output_buffer)

        # if this module is loaded as a library override the import
        if self.module_name != "__main__" : self._override_import()  
        else                              : self._on_the_fly()

    # post-processor with norun
    def _post_process_norun(self):

        if not self.output_file or self.output_file == "-": sys.stdout.write(self._output_buffer); return

        with io.open(self.output_file, "w", encoding=self.write_encoding) as f:
            f.write(self._output_buffer)

    # postprocessor - override an import
    def _override_import(self):
        try:

            module_name = self.module_name
            tmp_module_name = module_name + "__"

            pyc_output_file = self.output_file + "c"
            pyc_input_file  = self.input_file + "c"

            if module_name in sys.modules: del sys.modules[module_name]
            sys.path_importer_cache.clear()
            sys.modules[module_name] = __import__(tmp_module_name, globals() , locals(), ["*"])
            if tmp_module_name in sys.modules: del sys.modules[tmp_module_name]

        except:
            self._rewrite_traceback()

        finally:
            # remove tmp (.py & .pyc) files
            if os.path.exists(pyc_input_file)   : os.remove(pyc_input_file)
            if os.path.exists(pyc_output_file)  : os.remove(pyc_output_file)            
            if os.path.exists(self.output_file) : os.remove(self.output_file)
            

    # postprocessor - on-the-fly execution
    def _on_the_fly(self):
        try:            
            with io.open(self.output_file, "r", encoding=self.read_encoding) as f:
                 exec(f.read(), dict(__name__="__main__"))
        except:
            self._rewrite_traceback()
        finally:
            # remove tmp (.py & .pyc) files
            pyc_output_file = self.output_file + "c"
            if os.path.exists(self.output_file): os.remove(self.output_file)
            if os.path.exists(pyc_output_file): os.remove(pyc_output_file)
            if  not bool(self.defines.get(py_preprocessor._pypp_dont_call_sys_exit)):
                sys.exit(0)
