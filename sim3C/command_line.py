import sys
if not (sys.version_info.major == 2 and sys.version_info.minor == 7):
    print('Python {0.major}.{0.minor} is incompatible with sim3C, please use Python 2.7'.format(sys.version_info))
    sys.exit(1)

import logging
import time
import os

from Bio import SeqIO

from sim3C._version import version_stamp, runtime_info
from sim3C.abundance import generate_profile
from sim3C.art import ILLUMINA_PROFILES
from sim3C.exceptions import Sim3CException, FastaException
from sim3C.io_utils import open_output
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

    parser.add_argument('-C', '--compress', choices=['gzip', 'bzip2'], default=None,
                        help='Compress output files')

    parser.add_argument('-r', '--seed', metavar='INT', type=int, default=int(time.time()),
                        help="Random seed for initialising number generator")
    parser.add_argument('-m', '--method', default='hic', choices=['hic', 'meta3c', 'dnase'],
                        help='Library preparation method [hic]')
    parser.add_argument('-e', '--enzyme', dest='enzyme_name', default=None,
                        help='Restriction enzyme (case-sensitive) [NlaIII]')

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

    parser.add_argument('--create-cids', default=False, action='store_true',
                        help='Simulate chromosome interacting domains')
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

    parser.add_argument('--dist', metavar='DISTNAME', choices=['equal', 'uniform', 'lognormal'],
                        help='Abundance profile distribution choices: equal, uniform, lognormal')
    parser.add_argument('--lognorm-mu', metavar='FLOAT', type=float, default='1',
                        help='Log-normal relative abundance mu parameter')
    parser.add_argument('--lognorm-sigma', metavar='FLOAT', type=float, default='1',
                        help='Log-normal relative abundance sigma parameter')

    parser.add_argument('--simple-reads', default=False, action='store_true', help='Do not simulate sequencing errors')
    parser.add_argument('--machine-profile', help='An ART sequencing machine profile [EmpMiSeq250]',
                        default='EmpMiSeq250', choices=ILLUMINA_PROFILES.keys())
    parser.add_argument('--ins-rate', type=float, default=9.e-5, help='Insert rate [9e-5]')
    parser.add_argument('--del-rate', type=float, default=1.1e-4, help='Deletion rate [1.1e-4]')

    parser.add_argument('--bridge-adapter', dest='bridge_seq', default='',
                        help='bridge adapter sequence (for dnase mode)')

    parser.add_argument('--sam', dest='sam_out', default='',
                        help='Output file name for sam file. No output sam file if not specified.')

    parser.add_argument(dest='genome_seq', metavar='FASTA',
                        help='Genome sequences for the community')
    parser.add_argument(dest='output_file', metavar='OUTPUT',
                        help='Output Hi-C reads file')
    args = parser.parse_args()

    logger = init_log(args.verbose)

    logger.debug(runtime_info())
    logger.debug(sys.version.replace('\n', ' '))

    try:

        if 'community_table' in args and args.dist:
            raise RuntimeError('Cannot define abundance both explicitly as a table (-t) and a distribution (--dist).')

        if args.method == 'dnase':
            if args.enzyme_name:
                raise RuntimeError('The dnase method does not accept an enyzme specification.')
        elif not args.enzyme_name:
            raise RuntimeError('No enzyme was specified')

        # bridge adapter settings
        if args.method == 'dnase':
            if args.bridge_seq == '':
                logger.warning('Bridge adapter seq is not set. sim3C creates reads without bridge adapter seq.')
        else:
            if args.bridge_seq != '':
                raise RuntimeError('Bridge adapter seq is valid only in the dnase method.')

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

            profile_path = os.path.join(os.path.dirname(args.output_file), args.profile_name)
            if os.path.exists(profile_path):
                logger.error('Delete or move previous procedural abundance profile: {}'.format(profile_path))
                sys.exit(1)

            seq_index = None
            try:
                seq_index = SeqIO.index(args.genome_seq, 'fasta')
                seq_names = list(seq_index)
            except Exception:
                raise FastaException(args.genome_seq)
            finally:
                if seq_index:
                    seq_index.close()

            profile = generate_profile(args.seed, seq_names, mode=args.dist,
                                       lognorm_mu=args.lognorm_mu, lognorm_sigma=args.lognorm_sigma)

            # present result to console
            profile.write_table(sys.stdout)

            # save result to file
            with open(profile_path, 'w') as h_out:
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

        # list of CLI arguments to pass as parameters to the simulation
        kw_names = ['prefix', 'machine_profile', 'insert_mean', 'insert_sd', 'insert_min', 'insert_max',
                    'anti_rate', 'spurious_rate', 'trans_rate',
                    'efficiency',
                    'ins_rate', 'del_rate',
                    'create_cids', 'simple_reads', 'linear', 'convert_symbols', 'bridge_seq', 'sam_out']

        # extract these parameters from the parsed arguments
        kw_args = {k: v for k, v in vars(args).items() if k in kw_names}

        # initialise a sequencing strategy for this community
        # and the given experimental parameters
        strategy = SequencingStrategy(args.seed, args.profile_in, args.genome_seq, args.enzyme_name,
                                      args.num_pairs, args.method, args.read_length, **kw_args)

        # Run the simulation
        with open_output(args.output_file, mode='w', compress=args.compress) as out_stream:
            strategy.run(out_stream)


    except Sim3CException as ex:
        logger.error(str(ex))
        sys.exit(1)

    except Exception as ex:
        logger.exception(ex)
        sys.exit(1)
