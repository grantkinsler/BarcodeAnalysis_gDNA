# BarcodeAnalysis_v2_gDNA

Barcode Analysis for barcode_v2 gDNA, based on the original `SydShafferLab/BarcodeAnalysis`. A general description on how to run this code can be found in [this wiki](https://github.com/SydShafferLab/BarcodeAnalysis/wiki). 

At this point, this script has been adapted for gDNA analysis and NOT for 10X barcodes. Please use `https://github.com/SydShafferLab/BarcodeAnalysis_v2_10X` for 10X analysis. 

Please refer to the `SydShafferLab/BarcodeAnalysis` wiki for how to edit the `paths_and_variables.json` file. These are the main changes from the original .json file: 
1. `barcodeSource`: keep as "gDNA" only. 
1. `strtseq` and `strtseq_revcomp` have been edited for barcode_v2. 
1. `bclen`: 20
1. `sc_mm`: 4
1. `starseqMatch`: 70
1. `spike_in_seqs` have been edited for new spike-ins that were made for barcode_v2.

Notes:
1. This code was edited to only process one lane file per sample (i.e. only "_L001"). If you have multiple lane .fastq files per sample, please rename the filenames such that all sample filenames have a unique string before "_L001". 
