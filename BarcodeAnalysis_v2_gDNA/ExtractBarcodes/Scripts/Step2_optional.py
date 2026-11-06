import Bio
from Bio import SeqIO
from Bio.Seq import Seq
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
import numpy as np
import glob
import os

import subprocess
from time import sleep
import statistics
import shutil

from rapidfuzz import fuzz
import jstyleson

from Functions.LV_distance import LV_distance
from Functions.check_barcodes import check_barcodes

print("Running...")
print(" ")

# Find paths_and_variables.json file
path_to_script = os.path.abspath(os.getcwd())

path = path = os.path.expanduser(path_to_script + "/paths_and_variables.json")

# read .json file
with open(path, 'r') as myfile:
    data=myfile.read()

result_dict =  jstyleson.loads(data) # Raise Exception

Outfolder= result_dict['Outfolder']    #folder you want outputs go go into (dont make this folder, this scipt will make it)
strtseq= result_dict['strtseq']        #common sequence right before starcode starts
num_of_barcodes_to_use= result_dict['num_of_barcodes_to_use']   # number of barcodes to use to calculate LV (2 random samplings are taken)

print("startseq:" + strtseq)

# #define any new paths
raw_seq_path = glob.glob(Outfolder + "/raw_sequences/**/*.txt")

check_out = Outfolder + "/check_barcodes/"

# Make necessary folder
os.mkdir(check_out)

# # ### For each of the samples run LV_distance
LV_distance(num_of_barcodes_to_use, raw_seq_path, check_out)

# # ### Ven
check_barcodes(check_out,raw_seq_path,strtseq)

print("")

print(' Done with PreBarcodeQuant :D')

print("")



