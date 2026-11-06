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
cellranger mkfastq --run=/project/shafferslab/Guillaume/10X_exp1_reanalysis/20190808_10X1_BC_r1_r2_seq1 --samplesheet=/project/shafferslab/Guillaume/10X_exp1_reanalysis/20190808_10X1_BC_r1_r2_seq1_ProcessingFiles/mkfastq/BC_samplesheet.csv  --jobmode=local --localcores=32 --localmem=128 --qc --use-bases-mask=Y28,I8,I8,Y123 --id=20190808_10X1_BC_r1_r2_seq1_fastq
