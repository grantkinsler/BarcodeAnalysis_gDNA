
import Bio
from Bio import SeqIO
from Bio.Seq import Seq
from rapidfuzz import fuzz
import matplotlib as plt
from matplotlib import pyplot
import numpy as np
import glob
import os
import subprocess
from time import sleep
import statistics
import shutil
import re
import jstyleson
import json


##### Users should modulate these parametes to run this script

barcode_max_length = 90 #max length that can be inputed in starcode

#### END OF USER DEFINED PARAMETERS ##### NOTE: there are many more parameters that can be changed but this required you to go into the code below and understand all the different funtions

# Find paths_and_variables.json file
path_to_script = os.path.abspath(os.getcwd())

path = path = os.path.expanduser(path_to_script + "/paths_and_variables.json")

# read paths_and_variables.json file
with open(path, 'r') as myfile:
    data=myfile.read()

result_dict =  jstyleson.loads(data) # Raise Exception

scripts=result_dict['scripts']         #path to all barcode analyssi scripts
Fastqfolder10x=result_dict['Fastqfolder10x'] #Folder that contains all folders containing FASTQ files generated from sequencing the barcodes
FastqfoldergDNA=result_dict['FastqfoldergDNA']#Folder that contains all folders containing gDNA FASTQ files generated from sequencing the barcodes.
Outfolder= result_dict['Outfolder']    #folder you want outputs go go into (dont make this folder, this scipt will make it)
strtseq= result_dict['strtseq']        #common sequence right before starcode starts
barcodeSource = result_dict['barcodeSource']    #determine whether the data has barcodes from 10x and gDNA "both" or only from 10x "10x"
GSAMP= result_dict['GSAMP']            #Define which samples should be run together in starcode
bclen = result_dict['bclen']           #length to keep from sequenced barcode
strtseq =  result_dict['strtseq']      #common sequence right before starcode starts
strtseq_revcomp =  result_dict['strtseq_revcomp'] #rev_comp common sequence right before starcode starts
startseqMatch =  result_dict['startseqMatch']# The percentage match you for startseq to be called as correct in a barcode
sc_mm =  result_dict['sc_mm']




#define funtion to determine file has something in the first line (consider files with nothing in the first line as empty)
def empty(fname):
    with open(fname) as f:
       return f.readline() == ""

def any_gDNA(split_name):
    check = True
    if len(split_name.strip("gDNA")) < 2:
        check = False
    return check

#define funtion to determine if folder exist
def does_folder_exist(path_to_folder):
    if not os.path.exists(path_to_folder):
        os.mkdir(path_to_folder)
    else:
        raise Exception("folder {} already exists".format(path_to_folder))

#define any new paths
raw_seq = Outfolder + "/raw_sequences/"
raw_seq_comb = raw_seq + "combined/"
raw_seq_not_comb = raw_seq + "not_combined/"
qScore_path = Outfolder  + "/qScore/"
qScore_comb = qScore_path + "combined"
qScore_not_comb = qScore_path + "not_combined"
# Add starcode to PATH
os.environ["PATH"] += os.pathsep + scripts + '/starcode/'

# Make any necessary folders
path_to_folders = [Outfolder,raw_seq,raw_seq_comb,raw_seq_not_comb,qScore_path,qScore_comb,qScore_not_comb]

# checking whether folder/directory exists
for path_to_folder in path_to_folders:
    does_folder_exist(path_to_folder)

print("startseq:" + strtseq)


#-------------------------------------10x--------------------------------------------
if barcodeSource == 'both' or barcodeSource == '10x':

    #unzip all files created b 10x for barcode runs
    gunzipCommand = ['gunzip', '-r', Fastqfolder10x]
    subprocess.call(gunzipCommand)
    print("unzipped cDNA")


    #get all Read2 fastq file paths
    all_R2_10x  = []
    for grp in GSAMP:
            for smp in grp:
                all_R2_10x.append(glob.glob(Fastqfolder10x + "/**/"  + smp + "*_R2*.fastq", recursive = True))

    all_R2_10x = [item for sublist in all_R2_10x for item in sublist]

    #define samples
    samples = []
    for paths in all_R2_10x:
        samples.append(paths.split("/")[-1].split("_L0")[0])
    samples = list(set(samples))
    samples.sort()

    #loop through all the Read 2 fastq files (which contain barcode sequences)
    for sample in samples:
        s_fastq = []
        for path in all_R2_10x:
            if "/"+ sample in path:
                s_fastq.append(path)
        #get the sequences in all the fastqs
        #conatains all sequences for a sample
        seqs = []
        qscore = []
        for fsmp in s_fastq:
            for record in SeqIO.parse(fsmp, "fastq"):
                    seqs.append(str(record.seq))
                    qscore.append(str(record.letter_annotations["phred_quality"]))


        #determine which seqeunces are actual barcodes by determining which start with the constant sequence before the barcode

        # variables: the constant sequence before the barcode, percent match to call it a barcode (here it is 70)

        # since we stagger priming sight, here we determine what the correct stagger is for a given fastq file
        strt_ind = []
        for lines in seqs:
            c_ind = len(lines.split(strtseq)[0])
            if c_ind < len(lines):
                strt_ind.append(c_ind)

        bc_strt = statistics.mode(strt_ind)
        print("bc_strt ", bc_strt)

        # modify sequences so that those that have a start sequence (where there can be erros) get replaces with a pefect start sequence, and get rid of stagger so that all sequence start with the perfect start sequences
        modseq1 = seqs
        # for i in seqs:
        #     modseq1.append(i)

        #trim barcodes
        # variable: how long do you want your barcode
        sc_input = []
        qscore_input = []
        for i in modseq1:
            sc_input.append(i[0:barcode_max_length+len(strtseq)])

        #write files with these edited barcodes ( these are used as the input into starcode)
        f = open(raw_seq_not_comb + "/" "raw_seq" + "_" + sample +".txt","w")
        f.write('\n'.join(sc_input))
        f.close()

        qscore_input = []
        for i in qscore:
            qscore_i = str(json.loads(i)[0:barcode_max_length+len(strtseq)])
            qscore_input.append(qscore_i)

        #write files qscore
        f = open(qScore_not_comb + "/" "qScore" + "_" + sample +".txt","w")
        f.write('\n'.join(qscore_input))
        f.close()
        sleep(20)



#-------------------------------------gDNA--------------------------------------------

if barcodeSource == 'both' or barcodeSource == 'gDNA':

    #unzip all files
    gunzipCommand = ['gunzip', '-r', FastqfoldergDNA]
    subprocess.call(gunzipCommand)
    print("unzipped gDNA")

    #get all gDNA Read1 fastq file paths
    all_R1_gDNA_unfilt = glob.glob(FastqfoldergDNA + "/**/*_R1*.fastq", recursive = True)
    all_R1_gDNA_unfilt.sort()
    print(all_R1_gDNA_unfilt[0])

    # Remove any Read1 fastq files that you dont care about
    all_R1_gDNA_temp = []
    for grp in GSAMP:
        for smp in grp:
            if str(smp) in str(all_R1_gDNA_unfilt):
                all_R1_gDNA_temp.append([s for s in all_R1_gDNA_unfilt if str(smp) in s])

    all_R1_gDNA = [item for sublist in all_R1_gDNA_temp for item in sublist]

    #define samples Read1
    samples_R1_gDNA = []
    for paths in all_R1_gDNA:
        samples_R1_gDNA.append(paths.split("/")[-1].split("_L0")[0])
    samples_R1_gDNA = list(set(samples_R1_gDNA))
    samples_R1_gDNA.sort()
    print(samples_R1_gDNA[0])


    #loop through all the Read 1 fastq files from gDNA (which contain barcode sequences) and write the textfile for starcode inputs
    readsKeep_gDNA = []
    counter = 0
    file_endpoints_gDNA = []
    for sample in samples_R1_gDNA:
        file_endpoints_gDNA.append([])
        file_endpoints_gDNA[counter].insert(0,0)
        s_fastq = []
        for path in all_R1_gDNA:
            if "/"+ sample in path:
                s_fastq.append(path)
        # get JUST the sequences in all the fastqs contains all sequences for a sample
        seqs = []
        qscore = []
        for fsmp in s_fastq:
            for record in SeqIO.parse(fsmp, "fastq"):
                seqs.append(str(record.seq.reverse_complement()))
                hold_q = str(record.letter_annotations["phred_quality"])
                hold_q = hold_q[len(hold_q)::-1]

                hold_q = hold_q.replace(']','')
                hold_q = hold_q.replace('[','')
                hold_q = ''.join(['[',hold_q,']'])

                qscore.append(hold_q)
            file_endpoints_gDNA[counter].append(len(seqs))

        #determine which seqeunces are actual barcodes by determining which start with the constant sequence before the barcode

        # variables: the constant sequence before the barcode, percent match to call it a barcode (here it is 70)

        # since we stagger priming sight, here we determine what the correct stagger is for a given fastq file
        strt_ind = []
        for lines in seqs:
            c_ind = len(lines.split(strtseq)[0])
            if c_ind < len(lines):
                strt_ind.append(c_ind)

        print("seq ", len(seqs))
        print("strt_ind ", len(strt_ind))
        bc_strt = statistics.mode(strt_ind)
        print("bc_strt ", bc_strt)

        # modify sequences so that those that have a start sequence (where there can be erros) get replaces with a pefect start sequence, and get rid of stagger so that all sequence start with the perfect start sequences
        modseq1 = seqs

        #trim barcodes
        # variable: how long do you want your barcode
        sc_input = []
        for i in modseq1:
            sc_input.append(i[0:barcode_max_length+len(strtseq)])

        #write files with these edited barcodes ( these are used as the input into starcode)
        f = open(raw_seq_not_comb + "/" "raw_seq" + "_" + sample +".txt","w")
        f.write('\n'.join(sc_input))
        f.close()

        qscore_input = []
        for i in qscore:
            qscore_i = str(json.loads(i)[0:barcode_max_length+len(strtseq)])
            qscore_input.append(qscore_i)

        #write files qscore
        f = open(qScore_not_comb + "/" "qScore" + "_" + sample +".txt","w")
        f.write('\n'.join(qscore_input))
        f.close()
        sleep(20)










count_path = 0
# Get all of the files in the raw_seq folder and concatenate them
for grp in GSAMP:

    count_path = count_path+1
    count = 0
    with open(raw_seq_comb + '/raw_seq_group'+ str(count_path) + '_comb.txt', 'w') as outfile:
        for smp in grp:
            smpf = (raw_seq_not_comb + "raw_seq_" + smp +".txt")
            with open(smpf) as infile:
                for line in infile:
                    outfile.write(line)
                count = count + 1
                if count < len(grp):
                    outfile.write("\n") # only add line if it is not the last sample

count_path = 0
# Get all of the files in the qScore folder and concatenate them
for grp in GSAMP:

    count_path = count_path+1
    count = 0
    with open(qScore_comb + '/qScore_group'+ str(count_path) + '_comb.txt', 'w') as outfile:
        for smp in grp:
            smpf = (qScore_not_comb + "/" "qScore" + "_" + smp +".txt")
            with open(smpf) as infile:
                for line in infile:
                    outfile.write(line)
                count = count + 1
                if count < len(grp):
                    outfile.write("\n") # only add line if it is not the last sample

print("")

print(' Done with Step 0 :D')

print("")
