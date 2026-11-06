from polyleven import levenshtein
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
import random
import os

def LV_distance(num_of_barcodes_to_use, raw_seq_path, Outfolder):

# Function takes in list of barcodes (subsets it by num_of_bases_to_use) and runs LV distance on a number of random barcodes.

# Variables:
# cell_barcodes = list of barcodes
# num_of_bases_to_use = number of barcodes to use from cell_barcodes
# output_filename = name of output file for pdf plot of the LV distance

    print(" ")
    print("     LV_distance() ...")
    print("                   ...")

    #define new paths
    LV_out = Outfolder + "LV_distance/"

    # Make necessary folder
    os.mkdir(LV_out)

    for grp in raw_seq_path:

        with open(grp, 'r') as infile:

             cell_barcodes = infile.readlines()

        output_filename = LV_out +  grp.split('/')[-1].split('.txt')[0] + '_LV_distance.pdf'
        print("     Working on " + output_filename)
        print(" ")
        with PdfPages(output_filename) as pdf:
            for barcode_length in [5,10,20,30,50,70]:
                print("                 Length of barcodes = " + str(barcode_length))

                cell_barcodes_chop = [x[0:barcode_length] for x in cell_barcodes]

                # Initialize list
                seq_rand_list = []
                levenshtein_distances = []

                # set a differnt random seed
                for x in [1,23]:

                    random.seed(x)

                    ### Choose num_of_barcodes_to_use random barcodes to check LV distance
                    seq_rand_list = random.sample(cell_barcodes_chop, num_of_barcodes_to_use)

                    for j,seq in enumerate(seq_rand_list):
                        for i in range(len(seq_rand_list)):
                            if j != i:
                                levenshtein_distances.append(levenshtein(seq,seq_rand_list[i]))

                # plot
                plot_f = plt.figure()

                weights = np.ones_like(levenshtein_distances)/float(len(levenshtein_distances))
                y_plot, _, _ = plt.hist(levenshtein_distances,150, weights=weights)

                plt.hist(levenshtein_distances,150, weights=weights);

                # aesthetics
                plt.ylim([0,np.max(y_plot) + 0.05])
                plt.xlim([-1,71])
                plt.xlabel('Levenshtein Distance')
                plt.ylabel('Frequency (Normalized)')
                plt.title("num_of_barcodes=" + str(num_of_barcodes_to_use) + " barcode_len=" + str(barcode_length))
                plt.suptitle(output_filename.split("/")[-1], y=0.99, fontsize=13)

                # Save/close
                pdf.savefig(plot_f)
                plt.close()
