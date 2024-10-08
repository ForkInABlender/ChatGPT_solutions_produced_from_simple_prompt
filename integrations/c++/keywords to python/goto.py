# Dylan Kenneth Eliot & GPT-4o ( Alpha Edition )

"""

In this version of goto.py, it uses exceptions and labels have no arguments. 


This version is python2.7 to current compatible, with current stable version of python being [python3.10].

"""
class Goto(Exception):
    pass


def goto(label, globals_=None):
    if globals_ is None:
        # Automatically grab the global scope of the calling module if not provided
        f_globs = __import__("inspect").currentframe().f_back.f_globals
    if label in f_globs:
        func = f_globs[label]
        func()
    else:
        raise ValueError(f"Label {label} not defined")

def label(func):
    def wrapper(*args, **kwargs):
        try:
            return func(*args, **kwargs)
        except Goto as jump:
            if jump.args[0] == func.__name__:
                return wrapper(*args, **kwargs)
            else:
                raise
    return wrapper


"""
# Usage example
@label
def start():
    print("This is the start.")
    goto("end")
    goto("middle")

@label
def middle():
    print("This is the middle.")

@label
def end():
    print("This is the end.")

start()
#"""
