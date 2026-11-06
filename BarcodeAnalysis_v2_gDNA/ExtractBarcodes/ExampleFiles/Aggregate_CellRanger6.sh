#!/bin/bash
#BSUB -J cellranger
#BSUB -o cellranger.%J.out
#BSUB -e cellranger.%J.error
#BSUB -n 32
#BSUB -M 131072
#BSUB -R "span[hosts=1] rusage [mem=131072]"

export PATH=/home/gharm/cellranger-6.1.2:$PATH
cd /project/shafferslab/Guillaume/10X_exp3_analysis
cellranger aggr --id=20211203_10X3_aggr_unormalized \
                 --csv=/project/shafferslab/Guillaume/10X_exp3_analysis/Aggregate_ProcessingFiles_unormalized/Aggregate.csv \
                 --normalize=none \
                 --jobmode=local \
                 --localcores=32 \
                 --localmem=128 \
                 --nosecondary
