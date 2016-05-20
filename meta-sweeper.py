#!/usr/bin/env python
"""
meta-sweeper - for performing parametric sweeps of simulated
metagenomic sequencing experiments.
Copyright (C) 2016 "Matthew Z DeMaere"

This program is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published
by the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful, but
WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program.  If not, see <http://www.gnu.org/licenses/>.
"""

import subprocess
import argparse
import time
import sys
import os


class CustomFormatter(argparse.HelpFormatter):

    def __init__(self, prog):
        super(CustomFormatter, self).__init__(prog)
        self._max_help_position = 50

parser = argparse.ArgumentParser(description='Launch meta-sweeper', formatter_class=CustomFormatter)

parser.add_argument('-s', '--seed', type=int, default=int(time.time()), help='Random seed')
parser.add_argument('-n', '--num-jobs', default=1, help='Number of concurrent jobs [1]')
parser.add_argument('-r', '--resume', action='store_true', default=False, help='Resume previous run')
parser.add_argument('-w', '--work-dir', default='work', help='Nextflow working directory')
parser.add_argument('--profile', default='standard', choices=['standard', 'sge', 'pbs'],
                    help='Execution profile [standard]')
parser.add_argument('--work-flow', required=True, choices=['sweep', 'other'],
                    help='Workflow to run')
parser.add_argument('output', help='Output directory')

args = parser.parse_args()

cmd_args = ['nextflow',
            '-C', 'test.config',
            'run', '{0}.nf'.format(args.work_flow),
            '-w', args.work_dir,
            '-qs', str(args.num_jobs),
            '--seed', str(args.seed),
            '-profile', args.profile,
            '--output', args.output]

if args.resume:
    if not os.exists(args.work):
        print 'Cannot resume, the working directory {0} does not exist'.format(args.work_dir)
        sys.exit(1)
    cmd_args += ['-resume']

subprocess.call(cmd_args)

