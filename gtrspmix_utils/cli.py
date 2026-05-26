
import sys
import subprocess
import logging
from pathlib import Path
import argparse
import datetime
from gtrspmix_utils.constants import __version__,__ascii__

class MyArgFormatter(argparse.RawTextHelpFormatter):
    def _format_action_invocation(self, action):
        if not action.option_strings or action.nargs == 0:
            return super()._format_action_invocation(action)
        return ', '.join(action.option_strings)


def get_options(start_time):
    description = f'{__ascii__}{__version__} \n# \
                \nDescription: \
                \n  Optimize multiple exchangeability matrices (GTRs) with \
                \n  Profile mixture models \
                  '
    parser = argparse.ArgumentParser(description=description, add_help=True, \
        formatter_class=MyArgFormatter)

    input_group = parser.add_argument_group('Input data')
    input_group.add_argument('-s', '--fasta', type=Path, action='store', required=True, \
        help='[Required] Train fasta file')
    input_group.add_argument('-t', '-te', '--tree', type=Path, action='store', required=True, \
        help='[Required] Tree file in newick format')

    general_group = parser.add_argument_group('General flags')
    general_group.add_argument('--opt-gtr', action='store_true', \
        help='[Optional] Enable GTR Optimization ')
    outdir = './GTRspmix_out_' + start_time.strftime('%Y%m%d')
    general_group.add_argument('-o', '--outdir', type=Path, action='store', default=outdir, \
        help=f'[Optional] Outdir \n(default={outdir})')
    general_group.add_argument('-iqtree', '--iqtree', type=str, action='store', default='iqtree3', \
        help='[Optional] IQ-TREE bin file PATH (default=iqtree3)')
    general_group.add_argument('-gotree', '--gotree', type=str, action='store', default='gotree', \
        help='[Optional] gotree bin file PATH (default=gotree)')
    general_group.add_argument('-nt', '-T', type=str, action='store', default=1, \
        help='[Optional] Num threads for IQ-TREE \n(default=1)')
    general_group.add_argument('-mem', '--mem', type=str, action='store', default='90%', \
        help='[Optional] Max RAM for IQ-TREE \n(default=90%%)')
    general_group.add_argument('-seed', '--seed', type=str, action='store', \
        default=start_time.strftime('%S%f'), \
        help='[Optional] Ramdom seed for IQ-TREE ')
    general_group.add_argument('-v', '--version', action='version', \
        version = f'%(prog)s {__version__}', \
        help=f'Show version and exit \nVersion: {__version__}')
    general_group.add_argument('--keep-tmp', action='store_true', default=False,
        help='[Optional] Keep intermediate alignment files.\
\nWARNING - These files can be extremely large depending on scale values\
\n(default: False)')
    general_group.add_argument('--overwrite', action='store_true', default=False,
        help='[Optional] Overwrite outdir \n(default: False)')

    model_group = parser.add_argument_group('Model settings')
    model_group.add_argument('-m', '--model', type=str, action='store', default=None, \
        help='[Optional] Use pre-defined GTRspmix model as starting model \n(e.g. S10pfamC60)')
    model_group.add_argument('-m-gtr20', '--m-gtr20', type=str, action='store', default='POISSON', \
        help='[Optional] Initial Exchangeability matrix \n(default=POISSON)')
    model_group.add_argument('-m-rate', '-mrate', '--mrate', '--m-rate', type=str, action='store', default='G4', \
        help='[Optional] Rate model name \nGamma or FreeRate is supported \n(default=G4)')
    model_group.add_argument('-nf', '--nexus-few', type=Path, action='store', default=None, \
        help='[Optional] NEXUS file of fewer profile model \n(e.g. MEOW10.nex)')
    model_group.add_argument('-nm', '--nexus-many', type=Path, action='store', default=None, \
        help='[Optional] NEXUS file of more profile model \n(e.g. MEOW60.nex)')
    model_group.add_argument('-km', '--kmeans', type=int, action='store', default=None, \
        help='[Optional] Use Kmeans++ to make profile clusters\
\nSpecify num of clusters (exchangeability matrices) \n(e.g. 10)')
    model_group.add_argument('-n', '--nexus', type=Path, action='store', default=None, \
        help='[Optional] NEXUS file for Kmeans or re-starting')
    model_group.add_argument('-j', '--json', type=Path, action='store', default=None, \
        help='[Optional] JSON file with profile cluster from previous run')

    epsilon_group = parser.add_argument_group('Convergence threshold')
    epsilon_group.add_argument('-me', '--epsilon', type=float, action='store', default=10, \
        help='[Optional] Epsilon value for total LL \n(default=10)')
    epsilon_group.add_argument('-me-theta', '--epsilon-theta', type=float, action='store', default=0.01, \
        help='[Optional] Epsilon value for IQ-TREE theta estimation \n(default=0.01)')
    epsilon_group.add_argument('-me-gtr', '--epsilon-gtr', type=float, action='store', default=0.99, \
        help='[Optional] Epsilon value for IQ-TREE GTR optimization \n(default=0.99)')
    # epsilon_group.add_argument('-me-pro', '--epsilon-profile', type=float, action='store', default=0.99, \
    #     help='[Optional] Epsilon value for IQ-TREE Profile optimization \n(default=0.01)')
    epsilon_group.add_argument('--scale-gtr', type=float, action='store', default=10, \
        help='[Optional] Scaling factor for GTR soft-partitioning \n(default=10)')
    epsilon_group.add_argument('--forever', action='store_true', default=False,
        help='[Optional] Run until manually terminated \n(i.e., epsilon=0.0)')
    # epsilon_group.add_argument('--scale-profile', type=float, action='store', default=100, \
    #     help='[Experimental] Scaling factor for Profile soft-partitioning \n(default=100)')

    # restart_group = parser.add_argument_group('Flags for re-starting')
    # restart_group.add_argument('-j', '--json', type=Path, action='store', default=None, \
    #     help='[Optional] JSON file with frequency cluster')
#     restart_group.add_argument('-it', '--num-it', type=int, action='store', default=0, \
#         help='[Optional] Number of iteration of NEXUS file \
# \n(e.g. if you re-start from model_3.nex, specify \"3\")')

    experimental_group = parser.add_argument_group('Profile optimization (Experimental)')
    experimental_group.add_argument('--opt-profile-FO', action='store_true', \
        help='[Experimental] Enable Profile Optimization using ML estimation (FO)')
    experimental_group.add_argument('--opt-profile-F', action='store_true', \
        help='[Experimental] Enable Profile Optimization by observed frequency (F)')
    experimental_group.add_argument('-me-pro', '--epsilon-profile', type=float, action='store', default=0.99, \
        help='[Experimental] Epsilon value for IQ-TREE Profile optimization \n(default=0.01)')
    experimental_group.add_argument('--scale-profile', type=float, action='store', default=1000, \
        help='[Experimental] Scaling factor for Profile soft-partitioning \n(default=1000)')

    return parser.parse_args()



def check_args(args):
    logger = logging.getLogger('gtrspmix')
    # mode check
    if not (args.opt_gtr or args.opt_profile_FO or args.opt_profile_F):
        logger.error(f'Argument Error: Please specify "--opt-gtr" and/or "--opt-profile-FO" or "--opt-profile-F".')
        sys.exit(1)
    elif args.opt_profile_FO and args.opt_profile_F:
        logger.error(f'Argument Error: Do not specify "--opt-profile-FO" and "--opt-profile-F".')
        sys.exit(1)
    elif (args.opt_profile_FO or args.opt_profile_F) and not args.opt_gtr and args.nexus != None:
        mode = 'ProfileOnly'
    elif args.opt_gtr and args.model != None:
        mode = 'PreDefined'
    elif args.opt_gtr and args.nexus != None and args.json != None:
        mode = 'ReStart'
    elif args.opt_gtr and args.nexus_few != None and args.nexus_many != None:
        mode = 'FromScratch SPPC'
    elif args.opt_gtr and args.nexus != None and args.kmeans != None:
        mode = 'FromScratch Kmeans'
    else:
        logger.error(f'Argument Error: Please specify the correct combination of input MODELs.')
        logger.error(f'  - GTRspmix modes')
        logger.error('    - PreDefined: --opt-gtr AND --model')
        logger.error('    - SPPC: --opt-gtr AND --nexus_few AND --nexus_many')
        logger.error('    - Kmeans: --opt-gtr AND --nexus AND --kmeans')
        logger.error('    - ReStart: --opt-gtr AND --nexus AND --json')
        logger.error(f'  - Profile  Only mode')
        logger.error('    - ProfileOnly: --opt-profile AND --nexus')
        sys.exit(1)

    # dependencies check
    r = subprocess.run(args.iqtree+['-V'], stdout=subprocess.PIPE)
    if r.returncode != 0:
        logger.error(f'IQ-TREE check Failed: {args.iqtree} is not found in PATH.')
        sys.exit(1)
    r = subprocess.run(args.gotree+['version'], stdout=subprocess.PIPE)
    if r.returncode != 0:
        logger.error(f'gotree check Failed: {args.gotree} is not found in PATH.')
        sys.exit(1)

    # Profile Optimization
    if args.opt_profile_FO:
        logger.warning(f'Profile Optimization (FO) is ENABLED.')
        logger.warning(f'  - This flag is Experimental and may be Unstable.')

    if args.opt_profile_F:
        logger.warning(f'Profile Optimization (F) is ENABLED.')
        logger.warning(f'  - This flag is Experimental and may be Unstable.')

    if args.overwrite:
        logger.warning(f'"--overwrite" is ENABLED.')
        logger.warning(f'  - previous output will be removed.')

    # keep tmp files
    if args.keep_tmp:
        logger.warning(f'"--keep-tmp" is ENABLED.')
        logger.warning(f'  - This flag is storage-intensive.')



    return mode

