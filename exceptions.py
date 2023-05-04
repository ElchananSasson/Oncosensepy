class InvalidDataSetException(Exception):
    def __init__(self, message):
        super().__init__(message)


class InvalidPathException(Exception):
    def __init__(self, message):
        super().__init__(message)


class InvalidDirectoryPathException(Exception):
    def __init__(self, message):
        super().__init__(message)


class InvalidColumnsException(Exception):
    def __init__(self, message):
        super().__init__(message)


class NegativeNumberException(Exception):
    def __init__(self, message):
        super().__init__(message)


class InvalidCellLineException(Exception):
    def __init__(self, message):
        super().__init__(message)


class InvalidUIDException(Exception):
    def __init__(self, message):
        super().__init__(message)
