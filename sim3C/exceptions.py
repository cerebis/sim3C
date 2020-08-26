class Sim3CException(Exception):
    """Module base exception class"""
    def __init__(self, message):
        super().__init__(message)


class NoCutSitesException(Sim3CException):
    """Occurs when a target template contains no cutsites for a specified restriction enzyme"""
    def __init__(self, enzyme_name):
        super(NoCutSitesException, self).__init__("no restriction sites were found for {}".format(enzyme_name))


class NoRepliconsException(Sim3CException):
    """Occurs when a target cell contains no replicons. This can be caused by supply references
    which contain no cut-sites """
    def __init__(self, cell_name):
        super(NoRepliconsException, self).__init__("cell {} contained no replicons".format(cell_name))


class FastaException(Sim3CException):
    """Occurs when Bio.SeqIO read calls throw an exception"""
    def __init__(self, file_name):
        super().__init__('Failed to read from sequence file [{}]. Is it FASTA formatted?'.format(file_name))


class OutOfBoundsException(Sim3CException):
    """Raised when coordinates lie out of range of replicon"""
    def __init__(self, pos, maxpos):
        super().__init__("exceeded maximum template length {} > {}".format(pos, maxpos))


class EmptyRegistryException(Sim3CException):
    """No registry was empty when attempting to act upon its contents"""
    pass


class MonochromosomalException(Sim3CException):
    """A method require more than one chromosome was invoked on a monochromosomal cell"""
    pass
