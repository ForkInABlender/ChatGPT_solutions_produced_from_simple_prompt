# Dylan Kenneth Eliot & GPT-4 ( Alpha Edition )

"""

now it can remember and pass its cache between layers. or transformers....


"""

from pybrain3.structure import LinearLayer

class MemoizationLayer(LinearLayer):
    def __init__(self, *args, **kwargs):
        super(MemoizationLayer, self).__init__(*args, **kwargs)
        self.cache = {}

    def _forwardImplementation(self, inbuf, outbuf):
        input_tuple = tuple(inbuf)  # Ensure the input is hashable
        if input_tuple in self.cache:
            outbuf[:] = self.cache[input_tuple]
        else:
            super(MemoizationLayer, self)._forwardImplementation(inbuf, outbuf)
            self.cache[input_tuple] = outbuf.copy()
