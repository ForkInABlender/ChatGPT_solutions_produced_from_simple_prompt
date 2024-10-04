from .exceptions import InterfaceException


class interface:
    def __init__(self, interface_klass, *args, **kwargs):
        self.required_methods = self.get_method_list(interface_klass)

    def __call__(self, Klass, *args, **kwargs):
        missing_methods = [method for method in self.required_methods if method not in self.get_method_list(Klass)]
        if not len(missing_methods) == 0:
            method_list = ','.join(missing_methods)
            raise InterfaceException(f'Missing interface methods: {method_list}')

        class NewKlass(*Klass.__bases__):
            def __init__(self, *args, **kwargs):
                self.og = Klass(*args, **kwargs)

            def __getattribute__(self, attr):
                try:
                    x = super().__getattribute__(attr)
                except AttributeError:
                    pass
                else:
                    return x

                return self.og.__getattribute__(attr)
        return NewKlass

    @staticmethod
    def get_method_list(klass):
        return [func for func in dir(klass) if callable(getattr(klass, func)) and not func.startswith("__")]
