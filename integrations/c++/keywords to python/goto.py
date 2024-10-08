# Dylan Kenneth Eliot & GPT-4o ( Alpha Edition )

"""

In this version of goto.py, it uses exceptions and labels have no arguments. 


This version is python2.7 to current compatible, with current stable version of python being [python3.10].

"""
class Goto(Exception):
    pass

def goto(label, globals_=None):
    if globals_ is None:
        globals_ = __import__("inspect").currentframe().f_back.f_globals
    if label in globals_:
        func = globals_[label]
        if hasattr(func, 'is_label') and func.is_label:
            func()
        else:
            raise ValueError(f"Invalid Label {label}; decorate with `@label` to treat it as a label")
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
    wrapper.is_label = True
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
