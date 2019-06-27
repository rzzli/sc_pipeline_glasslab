#!/bin/bash

#filenames are tumor1 ... tumorX, each should have its own tumor1.sample sheet
WKDIR=$1 #directory with .sample files
OUTDIR=$2 #dir of output
JSONPATH=$3 #path for json file parameter
rfile=$4 #full path to rfile
genome=$5 #full path to ref genome
core=$6 # num of core, default 32
mem=$7 #memory, default 128

[ -z "$genome" ] && genome=/data/scratch/ouyang/sc10X/refdata-cellranger-GRCh38-3.0.0 || echo $genome
[ -z "$core" ] && core=32 || echo $core
[ -z "$mem" ] && mem=128 || echo $mem
#core=
#mem=
cd ${OUTDIR}
for sample in $(ls ${WKDIR}/*.sample)
do
  bname=$(basename ${sample} .sample)
  mkdir ${bname}
  cd ${bname}
  logname="${bname}_QM.csv"
  touch ${logname}

#filename=$1

    {
    read
    while IFS=$'\t' read -r sID    location    index10x    prefix
    do
    /data/scratch/ouyang/sc10X/cellranger-3.0.1/cellranger count --id=${sID} \
                    --sample=${prefix} \
                    --transcriptome=${genome} \
                    --fastqs=${location} \
                    --localcores=${core} \
                    --localmem=${mem} \
                    --nosecondary \
                    --indices=${index10x}
    done
    } < $sample # run each row of the .sample file



## concat csv files from sub directory
    for dir in */ # loop through microglia1 2 3 ....
        do
        if [ -s ${logname} ]
        then
        tail -n +2 ${dir}/outs/metrics_summary.csv | cat>>${logname}  # if exist, add after line two
        else
        cat ${dir}/outs/metrics_summary.csv>>${logname} # if empty, add header and content starting line two
        fi

        done
    cd ..
done
###########
echo "done with fastq to count file"


## part 2

#WKDIR should contain .sample files $1
#OUTDIR is where results are stored
#basename(cur_sample) should be tumor1, tumor2 etc
# change directory to outdir
cd ${OUTDIR}

mkdir -p Results
mkdir -p intermediate
for sample in $(ls ${WKDIR}/*.sample)
do
  bname=$(basename ${sample} .sample)  #basename
  sample_base=$(basename ${sample})  #control.sample rather than full path to control.sample
  mkdir Results/${bname}
  mkdir intermediate/${bname}
  Rscript ${rfile} ${sample} ${OUTDIR} ${JSONPATH}
  echo "finished ${bname} script"


done
