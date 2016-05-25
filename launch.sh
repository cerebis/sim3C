#!/bin/bash
# meta-sweeper - for performing parametric sweeps of simulated
# metagenomic sequencing experiments.
# Copyright (C) 2016 "Matthew Z DeMaere"
# 
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published
# by the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
# 
# This program is distributed in the hope that it will be useful, but
# WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
# 
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <http://www.gnu.org/licenses/>.

function test_exists () {
    command -v $1 >/dev/null 2>&1 || { echo >&2 "Required tool \"$1\" was not on path.  Aborting."; exit 1; }
}

function canonical_path () {
    if [ ! -e $1 ]
    then
        echo "$1 does not exist"
        exit 1
    else
        echo `/usr/bin/env python -c "import os, sys; print os.path.realpath(sys.argv[1])" $1`
    fi
}

#
# Extract the location of this script to determine the full path
# for meta-sweeper home. This will be set as METASWEEPER_HOME.
#
# Set the location of meta-sweeper
script_path=`canonical_path $0`
export METASWEEPER_HOME=${script_path%/*}

# Nextflow should be on path
test_exists "nextflow"

# Invoke nextflow, passing additional arguments from cmdline
nextflow $@
