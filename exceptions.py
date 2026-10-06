class BioForgeError(Exception):
    pass


class FastaFormatError(BioForgeError):
    pass


class InvalidSequenceError(BioForgeError):
    pass


class DataFileError(BioForgeError):
    def __init__(self, message, line_number=None):
        self.line_number = line_number
        super().__init__(message)


BioforgeError = BioForgeError
InvalidsequenceError = InvalidSequenceError
