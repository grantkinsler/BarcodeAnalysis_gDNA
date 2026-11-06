#
from glob import glob
from Bio import SeqIO
from Bio.Seq import Seq
import statistics
from rapidfuzz import fuzz
from time import sleep

def cDNA_10x_after_starcode(mod_R2,GSAMP,filt_haveStart_10x,filt_WSN_10x):
    #get all Read2 fastq file paths after starcode
    all_R2_10x_mod_unfilt = glob(mod_R2 + "/**/*_R2*.fastq", recursive = True)
    all_R2_10x_mod_unfilt.sort()

    # Remove any Read2 fastq files that you dont care about
    all_R2_10x_mod_temp = [] 
    for grp in GSAMP:
        for smp in grp:
            if str(smp) in str(all_R2_10x_mod_unfilt):
                all_R2_10x_mod_temp.append([s for s in all_R2_10x_mod_unfilt if str(smp) in s])

    all_R2_10x_mod = [item for sublist in all_R2_10x_mod_temp for item in sublist]

    #define samples
    samples_R2_10x_mod = []
    for paths in all_R2_10x_mod:

        samples_R2_10x_mod.append(paths.split("/")[-1].split("_L0")[0])
    samples_R2_10x_mod = list(set(samples_R2_10x_mod))
    samples_R2_10x_mod.sort()


    # for 10x samples, build a boolean of sequences that follow WSN nucleotide motif
    file_endpoints_10x_postsc = []
    bool_WSN_10x = []
    counter = 0
    for sample in samples_R2_10x_mod:

        s_fastq = []
        for path in all_R2_10x_mod:
            if "/"+ sample in path:
                s_fastq.append(path)
                

        
        # Build boolean that sequences with runs of 4 or containing N are not used
        bool_WSN_10x.append([])
        file_endpoints_10x_postsc.append([])
        seqs = []
        for fsmp in s_fastq:
            for record in SeqIO.parse(fsmp, "fastq"):
                if 'AAAA' in record.seq:
                    bool_WSN_10x[counter].append(0)
                elif 'TTTT' in record.seq:
                    bool_WSN_10x[counter].append(0)
                elif 'GGGG' in record.seq:
                    bool_WSN_10x[counter].append(0)
                elif 'CCCC' in record.seq:
                    bool_WSN_10x[counter].append(0)
                elif 'NN' in record.seq:
                    bool_WSN_10x[counter].append(0)
                else:
                    bool_WSN_10x[counter].append(1)
                seqs.append(record.seq)

            file_endpoints_10x_postsc[counter].append(len(seqs))
            
        file_endpoints_10x_postsc[counter].insert(0,0)
        counter = counter + 1

    # Use the booleans to write new files only using only good read2s in 10x samples
    counter = 0
    for sample in samples_R2_10x_mod:

        s_fastq = []
        for path in all_R2_10x_mod:
            if "/"+ sample in path:
                s_fastq.append(path)
        
        # get all of the fastq lines from the modified fastqs
        lines = [] # advance counter as I read lines and only append those where the boolean is 1
        for fsmp in s_fastq:
            with open(fsmp) as f:
                for line in f:
                    lines.append(line)

        lane_counter = 0
        for fsamp in s_fastq: 

            #write files with these edited barcodes ( these are used as the input into starcode)
            f = open(filt_WSN_10x + "/" + fsamp.split('/')[-1],"w")
            # print('Beginning count: ' + str(file_endpoints_10x_postsc[counter][lane_counter]))
            # print('End count: '+ str(file_endpoints_10x_postsc[counter][lane_counter+1]))
            for j in range(file_endpoints_10x_postsc[counter][lane_counter],file_endpoints_10x_postsc[counter][lane_counter+1]):
                #add print statements
                if bool_WSN_10x[counter][j] == 1:
                    f.write(lines[(4*j)])
                    f.write(lines[(4*j)+1])
                    f.write(lines[(4*j)+2])
                    f.write(lines[(4*j)+3])
            f.close()
            
            lane_counter = lane_counter + 1
            
            
            
        counter = counter + 1

    #get all Read1 fastq file paths after starcode
    all_R1_10x_mod_unfilt = glob(filt_haveStart_10x + "/**/*_R1*.fastq", recursive = True)
    all_R1_10x_mod_unfilt.sort()

    # Remove any Read1 fastq files that you dont care about
    all_R1_10x_mod_temp = [] 
    for grp in GSAMP:
        for smp in grp:
            if str(smp) in str(all_R1_10x_mod_unfilt):
                all_R1_10x_mod_temp.append([s for s in all_R1_10x_mod_unfilt if str(smp) in s])

    all_R1_10x_mod = [item for sublist in all_R1_10x_mod_temp for item in sublist]

    #define samples
    samples_R1_10x_mod = []
    for paths in all_R1_10x_mod:

        samples_R1_10x_mod.append(paths.split("/")[-1].split("_L0")[0])
    samples_R1_10x_mod = list(set(samples_R1_10x_mod))
    samples_R1_10x_mod.sort()


    # Use the booleans to write new files only using only good read2s in 10x samples
    counter = 0
    for sample in samples_R1_10x_mod:

        s_fastq = []
        for path in all_R1_10x_mod:
            if "/"+ sample in path:
                s_fastq.append(path)
        
        # get all of the fastq lines from the modified fastqs
        lines = [] # advance counter as I read lines and only append those where the boolean is 1
        for fsmp in s_fastq:
            with open(fsmp) as f:
                for line in f:
                    lines.append(line)

        lane_counter = 0
        for fsamp in s_fastq: 

            #write files with these edited barcodes ( these are used as the input into starcode)
            f = open(filt_WSN_10x + "/" + fsamp.split('/')[-1],"w")
            # print('Beginning count: ' + str(file_endpoints_10x_postsc[counter][lane_counter]))
            # print('End count: '+ str(file_endpoints_10x_postsc[counter][lane_counter+1]))
            for j in range(file_endpoints_10x_postsc[counter][lane_counter],file_endpoints_10x_postsc[counter][lane_counter+1]):
                #add print statements
                if bool_WSN_10x[counter][j] == 1:
                    f.write(lines[(4*j)])
                    f.write(lines[(4*j)+1])
                    f.write(lines[(4*j)+2])
                    f.write(lines[(4*j)+3])
            f.close()
            
            lane_counter = lane_counter + 1
            
            
            
        counter = counter + 1

    #get all Read1 fastq file paths after starcode
    all_I1_10x_mod_unfilt = glob(filt_haveStart_10x + "/**/*_I1*.fastq", recursive = True)
    all_I1_10x_mod_unfilt.sort()

    # Remove any Read1 fastq files that you dont care about
    all_I1_10x_mod_temp = [] 
    for grp in GSAMP:
        for smp in grp:
            if str(smp) in str(all_I1_10x_mod_unfilt):
                all_I1_10x_mod_temp.append([s for s in all_I1_10x_mod_unfilt if str(smp) in s])

    all_I1_10x_mod = [item for sublist in all_I1_10x_mod_temp for item in sublist]

    #define samples
    samples_I1_10x_mod = []
    for paths in all_I1_10x_mod:

        samples_I1_10x_mod.append(paths.split("/")[-1].split("_L0")[0])
    samples_I1_10x_mod = list(set(samples_I1_10x_mod))
    samples_I1_10x_mod.sort()


    # Use the booleans to write new files only using only good read2s in 10x samples
    counter = 0
    for sample in samples_I1_10x_mod:

        s_fastq = []
        for path in all_I1_10x_mod:
            if "/"+ sample in path:
                s_fastq.append(path)
        
        # get all of the fastq lines from the modified fastqs
        lines = [] # advance counter as I read lines and only append those where the boolean is 1
        for fsmp in s_fastq:
            with open(fsmp) as f:
                for line in f:
                    lines.append(line)

        lane_counter = 0
        for fsamp in s_fastq: 

            #write files with these edited barcodes ( these are used as the input into starcode)
            f = open(filt_WSN_10x + "/" + fsamp.split('/')[-1],"w")
            # print('Beginning count: ' + str(file_endpoints_10x_postsc[counter][lane_counter]))
            # print('End count: '+ str(file_endpoints_10x_postsc[counter][lane_counter+1]))
            for j in range(file_endpoints_10x_postsc[counter][lane_counter],file_endpoints_10x_postsc[counter][lane_counter+1]):
                #add print statements
                if bool_WSN_10x[counter][j] == 1:
                    f.write(lines[(4*j)])
                    f.write(lines[(4*j)+1])
                    f.write(lines[(4*j)+2])
                    f.write(lines[(4*j)+3])
            f.close()
            
            lane_counter = lane_counter + 1
            
            
            
        counter = counter + 1

    #get all Read1 fastq file paths after starcode
    all_I2_10x_mod_unfilt = glob(filt_haveStart_10x + "/**/*_I2*.fastq", recursive = True)
    all_I2_10x_mod_unfilt.sort()

    # Remove any Read1 fastq files that you dont care about
    all_I2_10x_mod_temp = [] 
    for grp in GSAMP:
        for smp in grp:
            if str(smp) in str(all_I2_10x_mod_unfilt):
                all_I2_10x_mod_temp.append([s for s in all_I2_10x_mod_unfilt if str(smp) in s])

    all_I2_10x_mod = [item for sublist in all_I2_10x_mod_temp for item in sublist]

    #define samples
    samples_I2_10x_mod = []
    for paths in all_I2_10x_mod:

        samples_I2_10x_mod.append(paths.split("/")[-1].split("_L0")[0])
    samples_I2_10x_mod = list(set(samples_I2_10x_mod))
    samples_I2_10x_mod.sort()


    # Use the booleans to write new files only using only good read2s in 10x samples
    counter = 0
    for sample in samples_I2_10x_mod:

        s_fastq = []
        for path in all_I2_10x_mod:
            if "/"+ sample in path:
                s_fastq.append(path)
        
        # get all of the fastq lines from the modified fastqs
        lines = [] # advance counter as I read lines and only append those where the boolean is 1
        for fsmp in s_fastq:
            with open(fsmp) as f:
                for line in f:
                    lines.append(line)

        lane_counter = 0
        for fsamp in s_fastq: 

            #write files with these edited barcodes ( these are used as the input into starcode)
            f = open(filt_WSN_10x + "/" + fsamp.split('/')[-1],"w")
            # print('Beginning count: ' + str(file_endpoints_10x_postsc[counter][lane_counter]))
            # print('End count: '+ str(file_endpoints_10x_postsc[counter][lane_counter+1]))
            for j in range(file_endpoints_10x_postsc[counter][lane_counter],file_endpoints_10x_postsc[counter][lane_counter+1]):
                #add print statements
                if bool_WSN_10x[counter][j] == 1:
                    f.write(lines[(4*j)])
                    f.write(lines[(4*j)+1])
                    f.write(lines[(4*j)+2])
                    f.write(lines[(4*j)+3])
            f.close()
            
            lane_counter = lane_counter + 1
            
        counter = counter + 1
