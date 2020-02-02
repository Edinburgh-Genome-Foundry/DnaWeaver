class ObjectDict(dict):
    def __getattr__(self, key):
        if key in self.__dict__:
            return self.__dict__[key]
        dict.__getattribute__(self, key)

    def __setattr__(self, key, value):
        self[key] = value
        self.__dict__[key] = value

    @staticmethod
    def from_dict(d):
        obj = ObjectDict(
            {
                key: (
                    [
                        ObjectDict.from_dict(e) if isinstance(e, dict) else e
                        for e in value
                    ]
                    if isinstance(value, (list, tuple))
                    else (
                        ObjectDict.from_dict(value)
                        if isinstance(value, dict)
                        else value
                    )
                )
                for (key, value) in d.items()
            }
        )
        for key, value in obj.items():
            sanitized_key = key.replace(" ", "_").replace(".", "_")
            obj.__dict__[sanitized_key] = value
        return obj