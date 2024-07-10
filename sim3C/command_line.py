import logging
import time
import os
import dnaio
import itertools

from Bio import SeqIO

import sim3C.random as random
from sim3C._version import version_stamp, runtime_info
from sim3C.abundance import generate_profile
from sim3C.art import ILLUMINA_PROFILES
from sim3C.exceptions import Sim3CException, FastaException
from sim3C.simulator import SequencingStrategy


__log_name__ = 'sim3C.log'


def init_log(verbose):
    """
    Initialise the runtime logger for both console and file output.

    :param verbose: set console verbosity level.
    :return: logger
    """
    logging.captureWarnings(True)
    logger = logging.getLogger('main')

    # root log listens to everything
    root = logging.getLogger('')
    root.setLevel(logging.DEBUG)

    # log message format
    formatter = logging.Formatter(fmt='%(levelname)-8s | %(asctime)s | %(name)7s | %(message)s')

    # Runtime console listens to INFO by default
    ch = logging.StreamHandler()
    if verbose:
        ch.setLevel(logging.DEBUG)
    else:
        ch.setLevel(logging.INFO)
    ch.setFormatter(formatter)
    root.addHandler(ch)

    fh = logging.FileHandler(__log_name__, mode='a')
    fh.setLevel(logging.DEBUG)
    fh.setFormatter(formatter)
    root.addHandler(fh)

    return logger


def check_suffixes(file_name1: str, file_name2: str = None) -> None:
    """
    Check that the suffixes on specified output files match and are understood by dnaio
    :param file_name1: first output file
    :param file_name2: optional second output file
    """
    def check_allowed(_fname: str) -> str:
        """
        Iterate over allowed suffixes and return the first match, or raise an exception
        :param _fname: file name to check
        :return: found suffix
        """
        for allowed in [''.join(li) for li in itertools.product(['.fq', '.fa'], ['', '.gz', '.bz2'])]:
            if _fname.endswith(allowed):
                return allowed
        raise Sim3CException('Output files must end in  [.fq, .fa, .fq.gz, fq.bz2, .fa.gz, or fa.bz2]')

    suffix1 = check_allowed(file_name1)
    if file_name2 is not None:
        suffix2 = check_allowed(file_name2)
        if suffix1 != suffix2:
            raise Sim3CException(f'Paired output files should have matching extensions: {suffix1}, {suffix2}')


def main():
    import argparse
    import sys

    #
    # Commandline interface
    #
    parser = argparse.ArgumentParser(description='Simulate HiC read pairs')
    parser.add_argument('-V', '--version', default=False, action='version',
                        version=version_stamp(False), help='Version')
    parser.add_argument('-v', '--verbose', default=False, action='store_true', help='Verbose output')

    parser.add_argument('--convert', dest='convert_symbols', default=False, action='store_true',
                        help='Convert unsupported symbols in sequence to N (required by Art)')

    parser.add_argument('-r', '--seed', metavar='INT', type=int, default=int(time.time()),
                        help="Random seed for initialising number generator")
    parser.add_argument('-m', '--method', default='hic', choices=['hic', 'meta3c', 'dnase'],
                        help='Library preparation method [hic]')
    parser.add_argument('-e', '--enzyme', metavar='NEB_NAME', action='append',
                        help='Case-sensitive NEB enzyme name. Use multiple times for multiple enzymes')
    parser.add_argument('-n', '--num-pairs', metavar='INT', type=int, required=True,
                        help='Number of read-pairs generate')
    parser.add_argument('-l', '--read-length', metavar='INT', type=int, required=True,
                        help='Length of reads from Hi-C fragments')
    parser.add_argument('--prefix', default='SIM3C', help='Prefix for read names [SIM3C]')
    parser.add_argument('--insert-mean', metavar='INT', type=int, default=400,
                        help='Mean insert size [400]')
    parser.add_argument('--insert-sd', metavar='INT', type=int, default=50,
                        help='Standard deviation of insert sizes [50]')
    parser.add_argument('--insert-min', metavar='INT', type=int, default=100,
                        help='Minimum allowed insert size [100]')
    parser.add_argument('--insert-max', metavar='INT', type=int, default=None,
                        help='Maximum allowed insert size [None]')

    parser.add_argument('--bridge-adapter',
                        help='Bridge adapter sequence (for dnase mode)')

    parser.add_argument('--linear', default=False, action='store_true',
                        help='Treat all replicons as linear molecules')
    parser.add_argument('--efficiency', metavar='FLOAT', type=float,
                        help='HiC/Meta3C efficiency factor [hic, dnase: 0.5 or meta3c: 0.02]')
    parser.add_argument('--anti-rate', metavar='FLOAT', type=float, default=0.2,
                        help='Rate of anti-diagonal fragments [0.2]')
    parser.add_argument('--trans-rate', metavar='FLOAT', type=float, default=0.1,
                        help='Rate of inter-replicon (trans) fragment formation [0.1]')
    parser.add_argument('--spurious-rate', metavar='FLOAT', type=float, default=0.02,
                        help='Rate of spurious fragment formation [0.02]')

    parser.add_argument('-P', '--profile', dest='profile_in', metavar='FILE',
                        help='Community abundance profile')
    parser.add_argument('--profile-name', metavar='FILE', default='profile.tsv',
                        help='Output file name for a procedural community profile')
    parser.add_argument('--profile-format', metavar='STRING', default='table',
                        choices=['table', 'toml'], help='Profile definition format (table, toml) [table]')

    parser.add_argument('--dist', metavar='DISTNAME', choices=['equal', 'uniform', 'lognormal'],
                        help='Abundance profile distribution choices: equal, uniform, lognormal')
    parser.add_argument('--lognorm-mu', metavar='FLOAT', type=float, default='1',
                        help='Log-normal relative abundance mu parameter')
    parser.add_argument('--lognorm-sigma', metavar='FLOAT', type=float, default='1',
                        help='Log-normal relative abundance sigma parameter')

    parser.add_argument('--simple-reads', default=False, action='store_true', help='Do not simulate sequencing errors')
    parser.add_argument('--machine-profile', help='An ART sequencing machine profile [EmpMiSeq250]',
                        default='EmpMiSeq250', choices=[*ILLUMINA_PROFILES])
    parser.add_argument('--ins-rate', type=float, default=9.e-5, help='Insert rate [9e-5]')
    parser.add_argument('--del-rate', type=float, default=1.1e-4, help='Deletion rate [1.1e-4]')

    parser.add_argument(dest='genome_seq', metavar='FASTA', help='Genome sequences for the community')
    parser.add_argument(dest='output_file_1', metavar='OUTPUT1', help='Output Hi-C reads file (interleaved or R1)')
    parser.add_argument(dest='output_file_2', metavar='OUTPUT2', help='Output Hi-C reads (R2)', nargs='?')
    args = parser.parse_args()

    logger = init_log(args.verbose)

    logger.debug(runtime_info())
    logger.debug(sys.version.replace('\n', ' '))

    try:

        if 'community_table' in args and args.dist:
            raise Sim3CException('Cannot define abundance both explicitly as a table (-t) and a distribution (--dist).')

        if args.enzyme:
            if args.method == 'dnase':
                raise Sim3CException('No enzyme can be specified for the dnase method')
            if len(args.enzyme) == 1:
                args.enzyme.append(None)
            elif len(args.enzyme) > 2:
                raise Sim3CException('A maximum of two enzymes can be specified')
        elif args.method != 'dnase':
            raise Sim3CException('At least one enzyme must be specified')

        if not args.bridge_adapter:
            if args.method == 'dnase':
                logger.warning('No bridge adapter sequence was specified for the dnase method')
        elif args.method != 'dnase':
            raise Sim3CException('Bridge adapter sequence is only valid for the dnase method')

        #
        # Prepare community abundance profile, either procedurally or from a file
        #
        #   Note: currently, all sequences for a single taxon are
        #   treated equally.
        #
        if not args.profile_in and not args.dist:
            logger.error('An abundance profile must be supplied either as a file or procedurally')
            sys.exit(1)

        if args.dist:
            # generate a procedural profile.ebug
            # the number of taxa is defined by number of sequences. i.e. monochromosomal organisms

            if os.path.basename(args.profile_name) != args.profile_name:
                logger.error('Arguments to profile-name should not contain path information')
                sys.exit(1)

            profile_path = os.path.join(os.path.dirname(args.output_file_1), args.profile_name)
            if os.path.exists(profile_path):
                logger.error('Delete or move previous procedural abundance profile: {}'.format(profile_path))
                sys.exit(1)

            seq_index = None
            try:
                seq_index = SeqIO.index(args.genome_seq, 'fasta')
                seq_names = list(seq_index)
                if len(seq_names) == 0:
                    logger.error('No valid records were found while scanning supplied fasta file {}'
                                 .format(args.genome_seq))
            except Exception:
                raise FastaException(args.genome_seq)
            finally:
                if seq_index:
                    seq_index.close()

            profile = generate_profile(seq_names, mode=args.dist, lognorm_mu=args.lognorm_mu,
                                       lognorm_sigma=args.lognorm_sigma)

            # present result to console
            profile.write_table(sys.stdout)

            # save result to file
            with open(profile_path, 'wt') as h_out:
                profile.write_table(h_out)

            # generated profile will be used downstream
            args.profile_in = profile_path

        if not args.efficiency:
            if args.method == 'hic':
                args.efficiency = 0.5
            elif args.method == 'meta3c':
                args.efficiency = 0.02
            elif args.method == 'dnase':
                args.efficiency = 0.5
            else:
                raise Sim3CException('BUG: missing efficiency for method {}'.format(args.method))

        # list of CLI arguments to pass as parameters to the simulation
        kw_names = ['prefix', 'machine_profile', 'insert_mean', 'insert_sd', 'insert_min', 'insert_max',
                    'anti_rate', 'spurious_rate', 'trans_rate',
                    'efficiency',
                    'ins_rate', 'del_rate',
                    'simple_reads', 'linear', 'convert_symbols', 'profile_format', 'bridge_adapter']

        # extract these parameters from the parsed arguments
        kw_args = {k: v for k, v in vars(args).items() if k in kw_names}

        # When none, seed will be drawn randomly
        random.init_state(args.seed)

        # initialise a sequencing strategy for this community
        # and the given experimental parameters
        strategy = SequencingStrategy(args.profile_in, args.genome_seq, args.enzyme,
                                      args.num_pairs, args.method, args.read_length, **kw_args)

        # Run the simulation
        if args.output_file_2 is not None:
            # split R1/R2 output
            check_suffixes(args.output_file_1, args.output_file_2)
            logger.info('Output read-pairs will be split into R1/R2 files')
            with dnaio.open(args.output_file_1, args.output_file_2, interleaved=False, mode='w') as writer:
                strategy.run(writer)
        else:
            # interleaved output
            check_suffixes(args.output_file_1)
            logger.info('Output read-pairs will be interleaved')
            with dnaio.open(args.output_file_1, mode='w', interleaved=True) as writer:
                strategy.run(writer)

    except Sim3CException as ex:
        logger.error(str(ex))
        sys.exit(1)

    except Exception as ex:
        logger.exception(ex)
        sys.exit(1)
