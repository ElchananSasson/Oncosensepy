class InvalidDataSetException(Exception):
    def __init__(self, message):
        super().__init__(message)


class InvalidPathException(Exception):
    def __init__(self, message):
        super().__init__(message)


class InvalidDirectoryPathException(Exception):
    def __init__(self, message):
        super().__init__(message)
