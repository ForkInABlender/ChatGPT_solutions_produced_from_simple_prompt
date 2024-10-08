# Dylan Kenneth Eliot

"""

This file is not my own works, but a replica of Solomon Ucko's answer that works.

Example usage as noted by the stack-overflow user:
'''
with Switch(2) as case:
    while case(0):
        print('foo')
        break
    while case(1, 2, 3):
        print('bar')
        break
    while case(4, 5):
        print('baz')
        break
    while case.default:
        print('default')
        break
'''

"""

class Switch:
    def __init__(self, value):
        self.value = value
        self._entered = False
        self._broken = False
        self._prev = None

    def __enter__(self):
        return self

    def __exit__(self, type, value, traceback):
        return False # Allows a traceback to occur

    def __call__(self, *values):
        if self._broken:
            return False
        
        if not self._entered:
            if values and self.value not in values:
                return False
            self._entered, self._prev = True, values
            return True
        
        if self._prev is None:
            self._prev = values
            return True
        
        if self._prev != values:
            self._broken = True
            return False
        
        if self._prev == values:
            self._prev = None
            return False
    
    @property
    def default(self):
        return self()
