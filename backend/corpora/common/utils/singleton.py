class Singleton(type):
    """
    Used to turn any class into a singleton-like object. All classes that inherit this class will refer to the same
    instance of the class.

    https://python-3-patterns-idioms-test.readthedocs.io/en/latest/Singleton.html
    https://stackoverflow.com/questions/6760685/creating-a-singleton-in-python
    """

    __instances = {}

    def __call__(cls, *args, **kwargs):
        if cls not in cls.__instances:
            cls.__instances[cls] = super(Singleton, cls).__call__(*args, **kwargs)
        return cls.__instances[cls]
