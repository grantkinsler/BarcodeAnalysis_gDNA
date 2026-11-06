#!/bin/bash
#BSUB -J cellranger
#BSUB -o cellranger.%J.out
#BSUB -e cellranger.%J.error
#BSUB -n 32
#BSUB -M 131072
#BSUB -R "span[hosts=1] rusage [mem=131072]"

## USAGE information
export PATH=/home/gharm/cellranger-6.1.2:$PATH
module load bcl2fastq2/v2.20.0.422
cd /project/shafferslab/Guillaume/10X_exp3_analysis/

cellranger mkfastq --run=/project/shafferslab/Guillaume/10X_exp3_analysis/20211215_10X3_All_cDNA_Barcodes_Seq1 --csv=/project/shafferslab/Guillaume/10X_exp3_analysis/20211215_10X3_All_cDNA_Barcodes_Seq1_ProcessingFiles/BC_samplesheet.csv  --jobmode=local --localcores=32 --localmem=128 --use-bases-mask=Y28,I8,I8,Y123 --id=20211215_10X3_All_cDNA_Barcodes_Seq1_FASTQ
