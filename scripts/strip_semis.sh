#!/bin/bash

if [ $# -ne 1 ]
then
    echo "Usage: [fasta file]"
    exit 1
fi

sed -i -r 's/^(>[^;]+); $/\1/' $1

