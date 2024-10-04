import sys
import types
import typing

from . import base

def modcall(module: typing.Union[str, types.ModuleType], func: callable) -> None:
    class CallableModule(base.AbstractCallableModule):
        __call__ = staticmethod(func)

    module: types.ModuleType = \
    (
        module
        if isinstance(module, types.ModuleType)
        else sys.modules[module]
    )

    module.__class__: typing.Type[CallableModule] = CallableModule
