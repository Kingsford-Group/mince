#!/bin/bash

input=$1
output=$2

awk '{OFS="\t"; getline seq; \
                getline sep; \
                getline qual; \
                print $0,seq,sep,qual}' $input | \
shuf | \
awk '{OFS="\n"; print $1,$2,$3,$4}' \
> $output
