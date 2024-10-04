import abc
import types
import typing

class AbstractCallableModule(types.ModuleType, abc.ABC):
    @staticmethod
    @abc.abstractmethod
    def __call__(*args: typing.Any, **kwargs: typing.Any) -> None:
        raise NotImplementedError
