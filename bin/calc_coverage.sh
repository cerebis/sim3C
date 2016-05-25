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
    echo "Usage: [bam file]"
    exit 1
fi

bedtools genomecov -ibam $1 | \
awk '
BEGIN{n=0} 
{
    # ignore whole genome records
    if ($1 != "genome") {
        # store names as they appear in repository
        # we use this to preserve file order 
        if (!($1 in seq_cov)) {
            name_repo[n++]=$1
        }
        # sum uses relative weights from histogram
        seq_cov[$1]+=$2*$3/$4
    }
}
END{
    for (i=0; i<n; i++) {
        print i+1, name_repo[i], seq_cov[name_repo[i]]
    }
}'

