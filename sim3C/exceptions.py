class Sim3CException(Exception):
    """Module base exception class"""
    def __init__(self, message):
        super(Sim3CException, self).__init__(message)


class NoCutSitesException(Sim3CException):
    """Occurs when a target template contains no cutsites for a specified restriction enzyme"""
    def __init__(self, seq_name, enz_name):
        super(NoCutSitesException, self).__init__(
            'sequence [{}] had no cutsites for enzyme [{}]'.format(seq_name, enz_name))


class FastaException(Sim3CException):
    """Occurs when Bio.SeqIO read calls throw an exception"""
    def __init__(self, file_name):
        super(FastaException, self).__init__(
            'Failed to read from sequence file [{}]. Is it FASTA formatted?'.format(file_name))


class OutOfBoundsException(Sim3CException):
    """Raised when coordinates lie out of range of replicon"""
    def __init__(self, pos, maxpos):
        super(OutOfBoundsException, self).__init__(
            "exceeded maximum template length {} > {}".format(pos, maxpos))


class EmptyRegistryException(Sim3CException):
    """No registry was empty when attempting to act upon its contents"""
    pass


class MonochromosomalException(Sim3CException):
    """A method require more than one chromosome was invoked on a monochromosomal cell"""
    pass
