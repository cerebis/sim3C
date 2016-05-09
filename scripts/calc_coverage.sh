#!/bin/bash

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

