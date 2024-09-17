# Dylan Kenneth Eliot & GPT-4o ( Alpha Edition )

"""

In this version of goto.py, it uses exceptions and labels have no arguments. 


"""
class Goto(Exception):
    pass

def goto(label, globals_=globals()):
    if label in globals_:
        func = globals_[label]
        func()
    else:
        raise ValueError("Label %s not defined" % label)

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
