import Bio.SeqIO
import Bio.bgzf
import argparse
import io_utils
import tqdm
from Bio.Restriction import Restriction


def get_enzyme_instance(enz_name):
    """
    Fetch an instance of a given restriction enzyme by its name.
    :param enz_name: the case-sensitive name of the enzyme
    :return: RestrictionType the enzyme instance
    """
    return getattr(Restriction, enz_name)


def open_file(fname):
    """
    Open either a bgzip or uncompressed file. Bgzip status is determined
    by a simplistic suffix (.gz) check.

    :param fname: the file to open
    :return: open stream
    """
    if args.input.endswith('gz'):
        hndl = Bio.bgzf.open(args.input, 'rb')
    else:
        hndl = open(args.input, 'r')
    return hndl


class SiteCounter(object):

    def __init__(self, enzyme_names):
        """
        Simple class to count the total number of enzymatic cut sites for the given
        list if enzymes.
        :param enzyme_names: a list of enzyme names (proper case sensitive spelling a la NEB)
        """
        self.enzymes = [get_enzyme_instance(en) for en in enzyme_names]

    def count_sites(self, seq, linear=True):
        """
        Count the number of sites found in the given sequence, where all enzymes defined
        at initialization are included.
        :param seq: Bio.Seq object
        :param linear: True, sequence is treated as linear.
        :return: the total number of sites .
        """
        return sum(len(en.search(seq, linear)) for en in self.enzymes)


parser = argparse.ArgumentParser(description='Collect various information about contigs/scaffolds output from SPAdes')
parser.add_argument('-e', '--enzymes', required=True, action='append',
                    help='A comma delimited set of case-sensitive enzyme names (NEB)')
parser.add_argument('--tip-size', type=int, help='Restrict analysis to tips of this length [None]')
parser.add_argument('--min-len', type=int, default=0, help='Minimum sequence length to analyze')
parser.add_argument('input', metavar='FASTA', help='Input FASTA to analyze')
parser.add_argument('output', help='Output file name')
args = parser.parse_args()

counter = SiteCounter(args.enzymes)
seq_info = []

# estimate total sequences before beginning.
nseq = 0
with open_file(args.input) as in_hndl:
    for l in in_hndl:
        if l[0] == '>':
            nseq += 1

# collect coverage, length and number of sites for each sequence
for si in tqdm.tqdm(Bio.SeqIO.parse(open_file(args.input), 'fasta'), total=nseq):
    li = len(si.seq)

    if li < args.min_len:
        continue

    row_i = {'seqid': si.id, 'length': li}

    # assume SPAdes continues to record coverage as the last element of the seq id.
    try:
        t = si.id.split('_')
        if len(t) <= 1:
            raise IOError("Sequence ID \"{}\" doesn't conform to SPAdes format.".format(si.id))
        row_i['coverage'] = float(t[-1])
    except TypeError as e:
        raise IOError('Error while parsing sequence id for coverage info. Is this output from SPAdes?')

    if args.tip_size:
        if li < 2*args.tip_size:
            # small contigs simply divide their extent in half
            l_sites = counter.count_sites(si.seq[:li/2])
            r_sites = counter.count_sites(si.seq[-li/2:])
        else:
            l_sites = counter.count_sites(si.seq[:args.tip_size])
            r_sites = counter.count_sites(si.seq[-args.tip_size:])

        row_i['l_sites'] = l_sites
        row_i['r_sites'] = r_sites
    else:
        row_i['sites'] = counter.count_sites(si.seq)

    seq_info.append(row_i)

# for reference, record input parameters as well as the collected information
report = {
    'enzymes': args.enzymes,
    'tip_size': args.tip_size,
    'min_len': args.min_len,
    'source': args.input,
    'seq_info': seq_info
}

# persist as a YAML file
with open(args.output, 'w') as out_hndl:
    io_utils.write_to_stream(out_hndl, report, fmt='yaml')
