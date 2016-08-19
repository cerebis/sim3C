#!/bin/bash
#
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
#

if [ $# -ne 1 ]
then
	echo "Usage: [file path]"
	exit 1
fi

if [ ! -e $1 ]
then
	echo "File path \"$1\" does not exist"
	exit 1
fi

WAIT_TIME=3
MAX_TIME=12
N=0
while true
do
    ISOPEN=`lsof $1 2>/dev/null | awk '!/^COMMAND/ {print $1,$2,$3,$9}'`
    if [ -z "$ISOPEN" ]
    then
        break
    fi

    echo "$1 is still open by a process, waiting..."

    sleep $WAIT_TIME
    ((N++))

    if (( N*WAIT_TIME >= MAX_TIME ))
    then
	echo "Exceeded maximum wait of $MAX_TIME seconds"
	exit 1
    fi
done


