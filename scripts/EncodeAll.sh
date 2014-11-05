#!/bin/bash

OPTIND=1 # reset OPTIND

showHelp () {
echo "usage: $0 -i <input_dir> -o <output_dir> -m {scalce, mince, quip, fastqz}"
}

inputDir=""
outputDir=""
method=""
pe=0
mince="/home/robp/mince/build/src/mince"
scalce="/home/robp/scalce-v2.7/scalce"
#quip="/home/robp/SoftwareStaging/quip-1.1.8/src/quip"
fastqz="/home/robp/fastqz/fastqz"

while getopts "h?pi:o:m:" opt; do
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
        p)
            pe=1
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
    #parallel 'mkdir job{#} && cd job{#} && ~/scalce-v2.7/scalce {} -T 6 -o /data/mince/scalce_results/{/.} && cd .. && rm -fr job{#}' ::: $inputDir/*fastq
    parallel 'echo processing {}; var={/.}; bn=${var%%_1}; mkdir job{#} && cd job{#} && ~/scalce-v2.7/scalce {} -r -T 6 -o /data/mince/scalce_results/$bn && cd .. && rm -fr job{#}' ::: $inputDir/*_1.fastq
fi

if [ "$method" = "fastqz" ]
then
    # W00t; GNU Parallel is *awesome*
    parallel -j 8 $fastqz c {} $outputDir/{/.} ::: `find $inputDir -mindepth 1 ! -name '*_[12].fastq'`
    #for fn in $inputDir/*.fastq; do
        #echo $fn
        #bn=`basename $fn`
        #outputName=${bn%%\.fastq}
        #echo "$fastqz c $fn $outputDir/$outputName"
        #$fastqz c $fn $outputDir/$outputName
    #done
fi


if [ "$method" = "mince" ]
then
    #variantArgs=("" "-n");   
    #variants=("" "_norc");
    variantArgs=("");
    variants=("");
    vaLen=${#variants[@]};
    mkdir -p $outputDir
    for (( i=1; i<${vaLen}+1; i++ ));
    do 
        suffix=${variants[$i-1]};
        mkdir -p $outputDir${suffix}
    done

    files=`find $inputDir -mindepth 1 -name *.fastq ! -name '*_[12].fastq'`
    echo $files
    for fn in $files; do
        echo file: $fn;
        bn=`basename $fn`;
        echo "**" ${vaLen} ${variantArgs[1]} "**"
        for (( i=1; i<${vaLen}+1; i++ ));
        do 
            arg=${variantArgs[$i-1]};
            suffix=${variants[$i-1]};
            thisOutputDir=${outputDir}${suffix};
            outputName=${bn%%\.fastq}.mince
            echo "$mince -b 12 -e ${arg} -i $fn -o $thisOutputDir/$outputName"
            $mince -b 12 -e ${arg} -i $fn -o $thisOutputDir/$outputName
            echo "plzip -c $thisOutputDir/$outputName.seqs > $thisOutputDir/$outputName.seqs.lz"
            plzip -c $thisOutputDir/$outputName.seqs > $thisOutputDir/$outputName.seqs.lz
            pigz -c $thisOutputDir/$outputName.seqs > $thisOutputDir/$outputName.seqs.gz
            echo "plzip -c $thisOutputDir/$outputName.offs > $thisOutputDir/$outputName.offs.lz"
            plzip -c $thisOutputDir/$outputName.offs > $thisOutputDir/$outputName.offs.lz
            pigz -c $thisOutputDir/$outputName.offs > $thisOutputDir/$outputName.offs.gz
            echo "plzip -c $thisOutputDir/$outputName.flips > $thisOutputDir/$outputName.flips.lz"
            plzip -c $thisOutputDir/$outputName.flips > $thisOutputDir/$outputName.flips.lz
            pigz -c $thisOutputDir/$outputName.flips > $thisOutputDir/$outputName.flips.gz
            if [ "$arg" == "" ]
            then
                decodeName=${outputName%%\.mince}.decoded.fa
                echo "$mince -d -i ${thisOutputDir}/${outputName} -o ${thisOutputDir}/../mince_decoded/${decodeName}"
                #$mince -d -i ${thisOutputDir}/${outputName} -o ${thisOutputDir}/../mince_decoded/${decodeName}
            fi	
        done
    done
fi

