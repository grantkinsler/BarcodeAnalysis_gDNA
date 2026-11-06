#!/bin/bash
#BSUB -J cellranger
#BSUB -o cellranger.%J.out
#BSUB -e cellranger.%J.error
#BSUB -n 32
#BSUB -M 131072
#BSUB -R "span[hosts=1] rusage [mem=131072]"


## USAGE information
export PATH=/home/gharm/cellranger3/cellranger-3.1.0:$PATH
module load bcl2fastq2/v2.20.0.422
cd /project/shafferslab/Guillaume/10X_exp1_reanalysis/
cellranger mkfastq --run=/project/shafferslab/Guillaume/10X_exp1_reanalysis/20190611_10X1_r1_seq1 --samplesheet=/project/shafferslab/Guillaume/10X_exp1_reanalysis/20190611_10X1_r1_seq1_ProcessingFiles/mkfastq/R1_10X1_samplesheet.csv --ignore-dual-index --jobmode=local --localcores=32 --localmem=128 --qc --id=20190611_10X1_r1_seq1_fastq
