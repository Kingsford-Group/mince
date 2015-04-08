#!/usr/bin/env bash
MINCE="$1"

if [ ${#} -eq 2 ] ; then
    IN="$2"
    BN=`basename "$IN" .fastq`
    gawk 'BEGIN {n=0} {if(n++ % 4 == 1) print $1}' $IN | sort > $BN-sorted.reads
    echo "Compressing $IN to $BN"
    $MINCE -e -r $IN -o $BN-ENC
    echo "Decompressing $BN-ENC to $BN-DEC"
    $MINCE -d -i $BN-ENC -o $BN-DEC
    echo "Comparing $BN"
    grep -v ">" $BN-DEC.fa | sort > $BN-DEC-sorted.reads
    cmp $BN-sorted.reads $BN-DEC-sorted.reads
else
    IN1="$2"
    IN2="$3"
    BN=`basename "$IN1" .fastq`
    paste $IN1 $IN2 \
        | gawk 'BEGIN {n=0} {if(n++ % 4 == 1) print $1}' \
        | sort > $BN-sorted.reads

    echo "Compressing $IN to $BN"
    $MINCE -e -l IU -1 $IN1 -2 $IN2 -o $BN-ENC
    echo "Decompressing $BN-ENC to $BN-DEC"
    $MINCE -d -i $BN-ENC -o $BN-DEC
    echo "Comparing $BN"
    paste $BN-DEC_1.fa $BN-DEC_2.fa | grep -v ">" | sort > $BN-DEC-sorted.reads
    cmp $BN-sorted.reads $BN-DEC-sorted.reads
fi
