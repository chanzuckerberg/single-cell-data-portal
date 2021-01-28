class ProcessingCancelled(Exception):
    def __init__(self, status, *args, **kwargs) -> None:
        super().__init__(status, *args, **kwargs)


class ProcessingFailed(Exception):
    def __init__(self, status, *args, **kwargs) -> None:
        super().__init__(status, *args, **kwargs)
