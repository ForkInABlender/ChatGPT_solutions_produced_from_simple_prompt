# Dylan Kenneth Eliot

"""
Mind you, this is the slimmed down version.

It minimally functions as a brython.js version of the cpython module 'hypothesis'.

"""


import itertools, sys

class ExhaustiveFail(Exception): pass
assume=lambda condition:eval('if not condition: raise ExhaustiveFail("assumption failed")')
class Strategies:
    @staticmethod
    def integers(min_value=-2, max_value=2): return list(range(min_value, max_value + 1))

def given(**kwargs_values):
    names, pools = list(kwargs_values.keys()), [list(v) for v in kwargs_values.values()]
    def decorator(fn):
        def wrapper(*args, **kwargs):
            tested_inputs = []
            for values in itertools.product(*pools):
                examples = dict(zip(names, values))
                try: tested_inputs.append(fn(*args, **examples))
                except ExhaustiveFail: continue
                except Exception as e: tested_inputs.append(examples); print("Falsifying example:", examples,"\b\nError:", repr(e)); raise
            return tested_inputs
        return wrapper
    return decorator
