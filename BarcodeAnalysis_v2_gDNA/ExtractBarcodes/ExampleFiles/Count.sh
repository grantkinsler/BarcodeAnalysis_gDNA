#!/bin/bash
#BSUB -J cellranger
#BSUB -o cellranger.%J.out
#BSUB -e cellranger.%J.error
#BSUB -n 32
#BSUB -M 131072
#BSUB -R "span[hosts=1] rusage [mem=131072]"

export PATH=/home/gharm/cellranger3/cellranger-3.1.0:$PATH
cd /project/shafferslab/Guillaume/10X_exp1_reanalysis
cellranger count --id=R1Enriched1_count \
                 --libraries=/project/shafferslab/Guillaume/10X_exp1_reanalysis/R1Enriched1_count_ProcessingFiles/Library.csv \
                 --transcriptome=/home/gharm/cellranger3/refdata-cellranger-GRCh38-3.0.0 \
                 --feature-ref=/project/shafferslab/Guillaume/10X_exp1_reanalysis/BarcodeProcessing/Barcode_output/CellRanger_inputs/FeatureReference.csv \
                 --jobmode=local \
                 --localcores=32 \
                 --localmem=128
