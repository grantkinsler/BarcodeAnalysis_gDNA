# This function takes in Fastq files taken from 10x cDNA barcodes and trims the sequences
# so that starcode can run properly
from glob import glob
from Bio import SeqIO
import statistics
from rapidfuzz import fuzz
from time import sleep

def cDNA_10x_starcode_prep(Fastqfolder10x,sc_in,GSAMP,bclen,filt_haveStart_10x,strtseq,startseqMatch,sc_mm):

    #get all Read2 fastq file paths
    all_R2_10x_unfilt = glob(Fastqfolder10x + "/**/*_R2*.fastq", recursive = True)
    all_R2_10x_unfilt.sort()

    # Remove any Read2 fastq files that you dont care about
    all_R2_10x_temp = []
    for grp in GSAMP:
        for smp in grp:
            if str(smp) in str(all_R2_10x_unfilt):
                all_R2_10x_temp.append([s for s in all_R2_10x_unfilt if str(smp) in s])

    all_R2_10x = [item for sublist in all_R2_10x_temp for item in sublist]

    #define samples
    samples_R2_10x = []
    for paths in all_R2_10x:

        samples_R2_10x.append(paths.split("/")[-1].split("_L0")[0])
    samples_R2_10x = list(set(samples_R2_10x))
    samples_R2_10x.sort()

    #loop through all the Read 2 fastq files from 10x (which contain barcode sequences) and write the textfile for starcode inputs
    readsKeep_10x = []
    counter = 0
    file_endpoints_10x = []
    for sample in samples_R2_10x:
        file_endpoints_10x.append([])
        file_endpoints_10x[counter].insert(0,0)

        s_fastq = []
        for path in all_R2_10x:
            if "/"+ sample in path:
                s_fastq.append(path)
        # get JUST the sequences in all the fastqs contains all sequences for a sample
        seqs = []
        for fsmp in s_fastq:
            for record in SeqIO.parse(fsmp, "fastq"):
                seqs.append(str(record.seq))
            file_endpoints_10x[counter].append(len(seqs))

        #determine which seqeunces are actual barcodes by determining which start with the constant sequence before the barcode

        # variables: the constant sequence before the barcode, percent match to call it a barcode (here it is 70)

        # since we stagger priming sight, here we determine what the correct stagger is for a given fastq file
        strt_ind = []
        for lines in seqs:
            c_ind = len(lines.split(strtseq)[0])
            if c_ind < len(lines):
                strt_ind.append(c_ind)

        bc_strt = statistics.mode(strt_ind)

        # modify sequences so that those that have a start sequence (where there can be erros) get replaces with a pefect start sequence, and get rid of stagger so that all sequence start with the perfect start sequences
        modseq1 = [] # For writing just the good sequences
        modseq2 = [] # For writing good and bad sequnces
        readsKeep_10x.append([])
        for i in seqs:
            strt = i[bc_strt:(len(strtseq)+bc_strt)]
            pctmatch = (fuzz.ratio(strtseq,strt))

            if pctmatch >= startseqMatch:
                trim = i[bc_strt + len(strtseq) :]
                modseq1.append(strtseq + trim)
                modseq2.append(strtseq + trim)
                readsKeep_10x[counter].append(True)
            else :
                readsKeep_10x[counter].append(False)
                modseq2.append(i)

        #trim barcodes
        # variable: how long do you want your barcode
        sc_input = []
        for i in modseq1:
            sc_input.append(i[0:bclen+len(strtseq)])
            
        #write files with these edited barcodes ( these are used as the input into starcode)
        f = open(sc_in + "/" "sc_input" + "_" + sample +".txt","w")
        f.write('\n'.join(sc_input))
        f.close()
        sleep(20)
        
        
        # get all of the fastq lines from the fastqs
        lines = []
        for fsmp in s_fastq:
            with open(fsmp) as f:
                for line in f:
                    lines.append(line)
                    
        # Write the edited sequence and trimmed file  to a fastq folder
        lane_counter = 0
        for fsamp in s_fastq:

            #write files with these edited barcodes ( these are used as the input into starcode)
            with open(filt_haveStart_10x + "/" + fsamp.split('/')[-1],"w") as f:
                # print('Beginning count: ' + str(file_endpoints_10x[counter][lane_counter]))
                # print('End count: '+ str(file_endpoints_10x[counter][lane_counter+1]))
                for j in range(file_endpoints_10x[counter][lane_counter],file_endpoints_10x[counter][lane_counter+1]):
                    if readsKeep_10x[counter][j] == 1:
                        f.write(lines[(4*j)])
                        f.write(modseq2[j][0:bclen+len(strtseq)]+'\n') # pulls the trimmed sequence
                        f.write(lines[(4*j)+2])
                        f.write(lines[(4*j)+3][bc_strt:(bclen+len(strtseq)+bc_strt)]+'\n')
            
            lane_counter = lane_counter + 1

        counter = counter+1
      
    #get all 10x Read1 fastq file paths
    all_R1_10x_unfilt = glob(Fastqfolder10x + "/**/*_R1*.fastq", recursive = True)
    all_R1_10x_unfilt.sort()

    # Remove any Read1 fastq files that you dont care about
    all_R1_10x_temp = []
    for grp in GSAMP:
        for smp in grp:
            if str(smp) in str(all_R1_10x_unfilt):
                all_R1_10x_temp.append([s for s in all_R1_10x_unfilt if str(smp) in s])

    all_R1_10x = [item for sublist in all_R1_10x_temp for item in sublist]

    #define samples Read1
    samples_R1_10x = []
    for paths in all_R1_10x:

        samples_R1_10x.append(paths.split("/")[-1].split("_L0")[0])
    samples_R1_10x = list(set(samples_R1_10x))
    samples_R1_10x.sort()

    # Loop through the corresponding Read1 files and remove the bad reads
    counter = 0
    for sample in samples_R1_10x:
        s_fastq = []
        for path in all_R1_10x:
            if "/"+ sample in path:
                s_fastq.append(path)

        # get all of the fastq lines from the fastqs
        lines = []
        for fsmp in s_fastq:
            with open(fsmp) as f:
                for line in f:
                    lines.append(line)
                    
        # Write the edited sequence and trimmed file  to a fastq folder
        lane_counter = 0
        for fsamp in s_fastq:

            #write files with these edited barcodes ( these are used as the input into starcode)
            with open(filt_haveStart_10x + "/" + fsamp.split('/')[-1],"w") as f:
                # print('Beginning count: ' + str(file_endpoints_10x[counter][lane_counter]))
                # print('End count: '+ str(file_endpoints_10x[counter][lane_counter+1]))
                for j in range(file_endpoints_10x[counter][lane_counter],file_endpoints_10x[counter][lane_counter+1]):
                    if readsKeep_10x[counter][j] == 1:
                        f.write(lines[(4*j)])
                        f.write(lines[(4*j)+1]) # pulls the trimmed sequence
                        f.write(lines[(4*j)+2])
                        f.write(lines[(4*j)+3])
            
            lane_counter = lane_counter + 1

        
        counter = counter+1

    #get all 10x Index1 fastq file paths
    all_I1_10x_unfilt = glob(Fastqfolder10x + "/**/*_I1*.fastq", recursive = True)
    all_I1_10x_unfilt.sort()
    # Remove any Read1 fastq files that you dont care about
    all_I1_10x_temp = []
    for grp in GSAMP:
        for smp in grp:
            if str(smp) in str(all_I1_10x_unfilt):
                all_I1_10x_temp.append([s for s in all_I1_10x_unfilt if str(smp) in s])

    all_I1_10x = [item for sublist in all_I1_10x_temp for item in sublist]

    #define samples Read1
    samples_I1_10x = []
    for paths in all_I1_10x:

        samples_I1_10x.append(paths.split("/")[-1].split("_L0")[0])
    samples_I1_10x = list(set(samples_I1_10x))
    samples_I1_10x.sort()

    # Loop through the corresponding index1 files and remove the bad reads
    counter = 0
    for sample in samples_I1_10x:
        s_fastq = []
        for path in all_I1_10x:
            if "/"+ sample in path:
                s_fastq.append(path)

        # get all of the fastq lines from the fastqs
        lines = []
        for fsmp in s_fastq:
            with open(fsmp) as f:
                for line in f:
                    lines.append(line)
                    
        # Write the edited sequence and trimmed file  to a fastq folder
        lane_counter = 0
        for fsamp in s_fastq:

            #write files with these edited barcodes ( these are used as the input into starcode)
            with open(filt_haveStart_10x + "/" + fsamp.split('/')[-1],"w") as f:
                # print('Beginning count: ' + str(file_endpoints_10x[counter][lane_counter]))
                # print('End count: '+ str(file_endpoints_10x[counter][lane_counter+1]))
                for j in range(file_endpoints_10x[counter][lane_counter],file_endpoints_10x[counter][lane_counter+1]):
                    if readsKeep_10x[counter][j] == 1:
                        f.write(lines[(4*j)])
                        f.write(lines[(4*j)+1]) # pulls the trimmed sequence
                        f.write(lines[(4*j)+2])
                        f.write(lines[(4*j)+3])
            
            lane_counter = lane_counter + 1

        
        counter = counter+1
      
    #get all 10x Index2 fastq file paths
    all_I2_10x_unfilt = glob(Fastqfolder10x + "/**/*_I2*.fastq", recursive = True)
    all_I2_10x_unfilt.sort()

    # Remove any Read1 fastq files that you dont care about
    all_I2_10x_temp = []
    for grp in GSAMP:
        for smp in grp:
            if str(smp) in str(all_I2_10x_unfilt):
                all_I2_10x_temp.append([s for s in all_I2_10x_unfilt if str(smp) in s])

    all_I2_10x = [item for sublist in all_I2_10x_temp for item in sublist]

    #define samples Read1
    samples_I2_10x = []
    for paths in all_I2_10x:

        samples_I2_10x.append(paths.split("/")[-1].split("_L0")[0])
    samples_I2_10x = list(set(samples_I2_10x))
    samples_I2_10x.sort()

    # Loop through the corresponding index2 files and remove the bad reads
    counter = 0
    for sample in samples_I2_10x:
        s_fastq = []
        for path in all_I2_10x:
            if "/"+ sample in path:
                s_fastq.append(path)

        # get all of the fastq lines from the fastqs
        lines = []
        for fsmp in s_fastq:
            with open(fsmp) as f:
                for line in f:
                    lines.append(line)
                    
        # Write the edited sequence and trimmed file  to a fastq folder
        lane_counter = 0
        for fsamp in s_fastq:

            #write files with these edited barcodes ( these are used as the input into starcode)
            with open(filt_haveStart_10x + "/" + fsamp.split('/')[-1],"w") as f:
                # print('Beginning count: ' + str(file_endpoints_10x[counter][lane_counter]))
                # print('End count: '+ str(file_endpoints_10x[counter][lane_counter+1]))
                for j in range(file_endpoints_10x[counter][lane_counter],file_endpoints_10x[counter][lane_counter+1]):
                    if readsKeep_10x[counter][j] == 1:
                        f.write(lines[(4*j)])
                        f.write(lines[(4*j)+1]) # pulls the trimmed sequence
                        f.write(lines[(4*j)+2])
                        f.write(lines[(4*j)+3])
            
            lane_counter = lane_counter + 1

        
        counter = counter+1

    return samples_R2_10x

