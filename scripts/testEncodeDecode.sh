#!/usr/bin/env bash

#============================
# USAGE:
# ./testEncodeDecode.sh PATH_TO_MINCE file1.fastq [file2.fastq]
#
# If two files are given, they will be treated as paired-end reads
#
# This will report whether mince was able to encode and decode correctly
#===========================

MINCE="$1"

if [ ${#} -eq 2 ] ; then
    IN="$2"
    BN=`basename "$IN" .fastq`
    gawk 'BEGIN {n=0} {if(n++ % 4 == 1) print $1}' $IN | sort > $BN-sorted.reads
    echo "Compressing $IN to $BN"
    $MINCE -e -r $IN -o $BN-ENC
    echo "Decompressing $BN-ENC to $BN-DEC"
    $MINCE -d -i $BN-ENC -o $BN-DEC
    echo ""
    echo "Comparing $BN"
    grep -v ">" $BN-DEC.fa | sort > $BN-DEC-sorted.reads
    cmp $BN-sorted.reads $BN-DEC-sorted.reads
    if [ ${?} -eq 0 ] ; then
        echo "Looks ok; deleting files"
        rm $BN-sorted.reads $BN-DEC-sorted.reads
        rm $BN-ENC*.lz
        rm $BN-DEC.fa
    fi
else
    IN1="$2"
    IN2="$3"
    BN=`basename "$IN1" .fastq`
    paste $IN1 $IN2 \
        | gawk 'BEGIN {n=0} {if(n++ % 4 == 1) print $0}' \
        | sort > $BN-sorted.reads

    echo "Compressing $IN to $BN"
    $MINCE -e -l IU -1 $IN1 -2 $IN2 -o $BN-ENC
    echo "Decompressing $BN-ENC to $BN-DEC"
    $MINCE -d -i $BN-ENC -o $BN-DEC
    echo ""
    echo "Comparing $BN"
    paste $BN-DEC1.fa $BN-DEC2.fa | grep -v ">" | sort > $BN-DEC-sorted.reads
    cmp $BN-sorted.reads $BN-DEC-sorted.reads
    if [ ${?} -eq 0 ] ; then
        echo "Looks ok; deleting files"
        rm $BN-sorted.reads $BN-DEC-sorted.reads
        rm $BN-DEC?.fa
        rm $BN-ENC*.lz
    fi
fi
