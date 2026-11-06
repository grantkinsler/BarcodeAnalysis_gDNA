# Summary of barcodes

# Run fuzzywuzzy to get %

from Bio import SeqIO
from rapidfuzz import fuzz
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
import numpy as np
import glob
import statistics
import shutil
import os
import re
from matplotlib_venn import venn3
import json
import warnings

def check_barcodes(check_out,raw_seq_path,strtseq):
	warnings.filterwarnings( "ignore", module = "matplotlib_venn/*" ) # for plotting

	print(" ")
	print("     check_barcodes() ...")
	print(" ")

	#define new paths
	check_out_fuzz = check_out + "/Match_to_StartSeq/"
	check_out_venn = check_out + "/Venn_diagram/"
	check_out_sum = check_out + "/summary/"

	# Make necessary folder
	os.mkdir(check_out_fuzz)
	os.mkdir(check_out_venn)
	os.mkdir(check_out_sum)


	for grp in raw_seq_path:
		with open(grp, 'r') as infile:
		     seqs = infile.readlines()

		#Change path so it matches the qScore
		qscore_path = grp.replace('raw_seq','qScore')
		qscore_path = qscore_path.replace('uences','')

		with open(qscore_path, 'r') as infile:
		     qscore = infile.readlines()

		sample = grp.split('/')[-1]
		#determine which seqeunces are actual barcodes by determining which start with the constant sequence before the barcode

		# variables: the constant sequence before the barcode, percent match to call it a barcode (here it is 70)

		# since we stagger priming sight, here we determine what the correct stagger is for a given fastq file
		strt_ind = []
		for lines in seqs:
		    c_ind = len(lines.split(strtseq)[0])
		    if c_ind < len(lines):
		        strt_ind.append(c_ind)
		print("seqs", len(seqs))
		print("strt_ind", len(strt_ind))
		bc_strt = statistics.mode(strt_ind)

		#------------------- Sequence run -------------------
		seq_base_run_no = []
		no_good_strings = ["AAAAA","TTTTT","GGGGG","CCCCC","NN"] # Any sequence we don't want 
		for i,lines in enumerate(seqs):
		    seq_i = lines[bc_strt:(len(strtseq)+bc_strt)]
		    if len([1 for x in no_good_strings if len(re.findall(x, seq_i) ) ]) > 0:
		        seq_base_run_no.append(i)
		        

		#------------------- Q-score -------------------
		qscore_no = []
		for i,lines in enumerate(qscore):

			qscore_i = json.loads(lines);
			qscore_i = np.array(qscore_i);

			# print(qscore_i)
			# print(type(qscore_i))
			# qscore_i = [y.strip().split(']')[0].split(' ') for y in qscore_i .split('[') if y!=''][0]
			# qscore_i = np.array([int(i) for i in qscore_i])
			# qscore_i = qscore_i[bc_strt:(len(strtseq)+bc_strt)]

			if len(np.where(qscore_i <= 14)[0]) > 5: #reads with Qscore <14 more than 5 times will get removed
				qscore_no.append(i)


		#------------------- Fuzz -------------------
		fuzz_range = [10,20,30,40,50]+list(range(60,101,2))
		fuzzy_values = []

		cnt_fuzz = 0

		for fuzzy_percent in fuzz_range:
			# modify sequences so that those that have a start sequence (where there can be erros) get replaces with a pefect start sequence, and get rid of stagger so that all sequence start with the perfect start sequences
			cnt_good = 0
			fuzzy_no_hold = []

			for i,lines in enumerate(seqs):
			    strt = lines[bc_strt:(len(strtseq)+bc_strt)]
			    pctmatch = (fuzz.ratio(strtseq,strt))

			    if pctmatch >= fuzzy_percent:
			    	cnt_good += 1

			fuzzy_values.append(cnt_good/(len(seqs)+1) ) 


		#------------------- Fuzz Venn plot -------------------
		fuzz_plot_val = [50,70,90,100]
		fuzzy_no = []

		cnt_fuzz = 0

		for fuzzy_percent in fuzz_plot_val:
			# modify sequences so that those that have a start sequence (where there can be erros) get replaces with a pefect start sequence, and get rid of stagger so that all sequence start with the perfect start sequences
			cnt_good = 0
			fuzzy_no_hold = []

			for i,lines in enumerate(seqs):
			    strt = lines[bc_strt:(len(strtseq)+bc_strt)]
			    pctmatch = (fuzz.ratio(strtseq,strt))

			    if pctmatch < fuzzy_percent:
			    	fuzzy_no_hold.append(i)

			fuzzy_no.append(fuzzy_no_hold)

			if len(fuzzy_no_hold) == 0:
				print("No sequence has less than " + str(fuzz_plot_val[cnt_fuzz]) + "percent overlap with the startseq") # make this a break


		with PdfPages(check_out_fuzz + sample + "_check_barcode_fuzN.pdf") as pdf: # startseq_specified_error
			plot_f = plt.figure()
			plt.plot(fuzz_range ,fuzzy_values)

			# asthetics
			plt.xlabel('percent match to startseq')
			plt.ylabel('Percent of sequences that contain the startseq')
			plt.title(" number_of_sequences=" + str(len(seqs)) + "  min_value=" + "{:.2f}".format(np.min(fuzzy_values))  )
			plt.suptitle(sample + "_Match_to_StartSeq", y=0.99, fontsize=13)
			pdf.savefig(plot_f)
			plt.close()

		colors = ['darkviolet','deepskyblue','pink']
		with PdfPages(check_out_venn + sample + "_check_barcode_venn.pdf") as pdf: # startseq_specified_error
			for cnt_i,fuzz_i in enumerate(fuzz_plot_val):
				plot_f = plt.figure()
				set1 = set(seq_base_run_no)
				set2 = set(fuzzy_no[cnt_i])
				set3 = set(qscore_no)

				ax = plt.gca() 

				total = len(set1.union(set2,set3))

				venn3([set1, set2,set3], ('Base_repeat', 'Match_to_StartSeq','Qscore'),subset_label_formatter=lambda x: f"{(x/total):1.0%}");

				plt.suptitle(sample + "_venn" +  "   pct_startseq_match=" + str(fuzz_i), y=0.99);
				plt.title("Percent of total reads: Qscore=" + str("{:3.0f}".format(100*len(qscore_no)/len(seqs))) + " Base_repeat=" + str("{:3.0f}".format(100*len(seq_base_run_no)/len(seqs))) + " Match_to_StartSeq=" + str("{:3.0f}".format(100*len(fuzzy_no[cnt_i])/len(seqs))) );
				pdf.savefig(plot_f);
				plt.close()

		with open(check_out_sum + sample +"_startseq_summary.txt","w") as write_file:
			write_file.write("Quality check on startsequence:" + "\n")
			write_file.write("Total_sequences = %d" % len(seqs) + "\n")
			write_file.write("           num reads removed by Qscore = %d" % len(qscore_no) + "\n")
			write_file.write("Percent_of_total_Qscore = " + str("{:3.0f}".format(100*len(qscore_no)/len(seqs))) + "\n")
			write_file.write("           num reads removed by Base_repeat = %d" % len(seq_base_run_no) + "\n")
			write_file.write("Percent_of_total_Base_repeat = " + str("{:3.0f}".format(100*len(seq_base_run_no)/len(seqs))) + "\n")
			write_file.write("..." + "\n")
			write_file.write("Varying Match_to_StartSeq:" + "\n")
			for cnt_i,fuzz_i in enumerate(fuzz_plot_val):
				write_file.write("           num reads removed by Match_to_StartSeq_"+str(fuzz_i) + " = " + str(len(fuzzy_no[cnt_i])) + "\n")
				write_file.write("Percent_of_total_Match_to_StartSeq_"+str(fuzz_i) + " = " + str("{:3.0f}".format(100*len(fuzzy_no[cnt_i])/len(seqs))) + "\n")

	print(" ")
	print("Done :D")
	print(" ")