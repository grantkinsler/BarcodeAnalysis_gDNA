#!/bin/bash
#BSUB -J cellranger
#BSUB -o cellranger.%J.out
#BSUB -e cellranger.%J.error
#BSUB -n 32
#BSUB -M 131072
#BSUB -R "span[hosts=1] rusage [mem=131072]"

export PATH=/home/gharm/cellranger3/cellranger-3.1.0:$PATH
cd /project/shafferslab/Guillaume/10X_exp1_reanalysis
cellranger aggr --id=20210305_10X1_aggr_unormalized \
                 --csv=/project/shafferslab/Guillaume/10X_exp1_reanalysis/Aggregate_ProcessingFiles/Aggregate.csv \
                 --normalize=none \
                 --jobmode=local \
                 --localcores=32 \
                 --localmem=128
