#!/bin/bash

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


