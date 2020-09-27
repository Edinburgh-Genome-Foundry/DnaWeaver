import functools


def memoize(obj):
    """Memoize the function results"""
    cache = obj.cache = {}

    @functools.wraps(obj)
    def memoizer(*args, **kwargs):
        key = str(args) + str(kwargs)
        if key not in cache:
            cache[key] = obj(*args, **kwargs)
        return cache[key]

    return memoizer


def functions_list_to_string(functions):
    def name_function(fun):
        name = str(fun)
        return "unnamed function" if name.startswith("<function") else name

    return " ; ".join([name_function(fun) for fun in functions])
