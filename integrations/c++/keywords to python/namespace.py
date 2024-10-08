# Dylan Kenneth Eliot & GPT-4o ( Alpha Edition )

"""
This namespace module is aimed for imitation of the namespace definition normally found in c++ when using "using namespace".

This version is "with {name-ofnamepsace-1} as __builtins__" compatible. Meaning you can defer multiple namespace objects to __builtins__ under a with statement, and still use python functions as per normal.

"""

class Namespace:
    def __init__(self):
        self.functions = {}
        self.active = False
    def __call__(self, func):
        def wrapper(*args, **kwargs):
            if self.active:
                return func(*args, **kwargs)
            else:
                return None
        self.functions[func.__name__] = wrapper
        return wrapper
    def __getattr__(self, name):
        return self.functions.get(name)
    def __enter__(self):
        self.active = True
        return self
    def __exit__(self, exc_type, exc_value, traceback):
        self.active = False
