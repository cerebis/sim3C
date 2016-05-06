#!/bin/bash 

SOURCE=${BASH_SOURCE[0]}
while [ -h "$SOURCE" ]; do
	SOURCE=`readlink $SOURCE`
done
DIR=`dirname $SOURCE`
base="test.phiX.a5"
echo "[test_a5] $DIR/bin/a5_pipeline.pl $DIR/example/phiX_p1.fastq $DIR/example/phiX_p2.fastq $base"
$DIR/bin/a5_pipeline.pl $DIR/example/phiX_p1.fastq $DIR/example/phiX_p2.fastq $base > $base.out 2> $base.err 

if [ ! -s $base.contigs.fasta ] ; then 
	echo "Test run of A5 unsuccessful. No contigs produced."
	exit
elif [ ! -s $base.final.scaffolds.fasta ] ; then
	echo "Test run of A5 unsuccessful. No scaffolds produced."
	exit
else 
	echo "A5 successfully produced contigs and scaffolds."
fi

NCHAR_REF=`cat $DIR/example/phiX.a5.final.scaffolds.fasta | wc -c | sed -e 's/\ //g'`
NCHAR_FINAL=`cat test.phiX.a5.final.scaffolds.fasta | wc -c | sed -e 's/\ //g'`

NCHAR_DIFF=`expr $NCHAR_REF - $NCHAR_FINAL`
if [ $NCHAR_DIFF -lt 0 ]; then
	NCHAR_DIFF=`expr $NCHAR_FINAL - $NCHAR_REF`
fi

if [ $NCHAR_DIFF -gt 100 ]; then
	echo "Test run of A5 unsuccessful."	
else 
	echo "Test run of A5 successful. Removing temporary files."
	rm -rf $base.*
fi
