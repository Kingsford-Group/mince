#!/bin/bash

OPTIND=1 # reset OPTIND

showHelp () {
echo "usage: $0 -i <input_dir>"
}

inputDir=""
outputDir=""
method=""
mince="../build/src/mince"
scalce="/home/robp/scalceOrig/scalce-v2.7/scalce"
quip="/home/robp/SoftwareStaging/quip-1.1.8/src/quip"
fastqz="/home/robp/SoftwareStaging/fastqz/fastqz"

while getopts "h?i:o:m:" opt; do
    case "$opt" in
        h|\?)
            showHelp
            exit 0
            ;;
        i)
            inputDir=$OPTARG
            ;;
        o)
            outputDir=$OPTARG
            ;;
        m)
            method=$OPTARG
            ;;
    esac
done

mkdir -p $outputDir

if [ "$method" = "quip" ]
then
    for fn in $inputDir/*.fastq; do
        echo $fn
        bn=`basename $fn`
        outputName=${bn%%\.fastq}.qp
        echo "$quip -a -v $fn"
        $quip -a -v $fn
        echo "mv ${fn}.qp $outputDir/$outputName"
        mv ${fn}.qp $outputDir/$outputName
        #$mince -e -i $fn -o $outputDir/$outputName
        #echo "plzip -c $outputDir/$outputName > $outputDir/$outputName.lz"
        #plzip -c $outputDir/$outputName > $outputDir/$outputName.lz
    done
fi

if [ "$method" = "scalce" ]
then
    echo "hi"
fi

if [ "$method" = "fastqz" ]
then
    for fn in $inputDir/*_merged.canonical.fastq; do
        echo $fn
        bn=`basename $fn`
        outputName=${bn%%_merged.canonical.fastq}
        echo "$fastqz c $fn $outputDir/$outputName"
        $fastqz c $fn $outputDir/$outputName
    done
fi


if [ "$method" = "mince" ]
then
    for fn in $inputDir/*_merged.canonical.fastq; do
        echo $fn
        bn=`basename $fn`
        outputName=${bn%%_merged.canonical.fastq}.mince
        echo "$mince -e -i $fn -o $outputDir/$outputName"
        $mince -e -i $fn -o $outputDir/$outputName
        echo "plzip -c $outputDir/$outputName.seqs > $outputDir/$outputName.seqs.lz"
        plzip -c $outputDir/$outputName.seqs > $outputDir/$outputName.seqs.lz
        echo "plzip -c $outputDir/$outputName.offs > $outputDir/$outputName.offs.lz"
        plzip -c $outputDir/$outputName.offs > $outputDir/$outputName.offs.lz
    done
fi

