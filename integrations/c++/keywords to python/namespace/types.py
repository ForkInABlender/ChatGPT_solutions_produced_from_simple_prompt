import types
import typing

class Namespace(types.SimpleNamespace):
    __setitem__ = types.SimpleNamespace.__setattr__

    def __call__(self, entity: typing.Union[types.FunctionType, typing.Type]) -> None:
        self[entity.__name__] = entity
