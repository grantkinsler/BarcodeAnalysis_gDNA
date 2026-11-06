#Script will take in gDNA outputs from sc_outputs and analyze them


import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
import numpy as np
import jstyleson
import glob
import json
import warnings
import pandas as pd
import itertools
import os



#----------------load variables from paths_and_variables.json file----------------------
#---------------------------------------------------------------------------------------

#define funtion to determine if folder exist
def does_folder_exist(path_to_folder):
    if not os.path.exists(path_to_folder):
        os.mkdir(path_to_folder)


# Find paths_and_variables.json file
path_to_script = os.path.abspath(os.getcwd())

path = path = os.path.expanduser(path_to_script + "/paths_and_variables.json")

# read paths_and_variables.json file
with open(path, 'r') as myfile:
    data=myfile.read()

result_dict =  jstyleson.loads(data) # Raise Exception

Outfolder= result_dict['Outfolder']    #folder you want outputs go go into (dont make this folder, this scipt will make it)
GSAMP= result_dict['GSAMP']            #Define which samples should be run together in starcode

spike_in_added= result_dict['spike_in_added']#"yes" = Spike-ins added "no" = Spike-ins not added
spike_in_color= result_dict['spike_in_color']#List of colors for each spike-in
spike_in_value= result_dict['spike_in_value']#List of number of cells for each spike-in
spike_in_seqs=result_dict["spike_in_seqs"]   #Sequence of known spike-ins without startseq 
strtseq=result_dict["strtseq"]               #Common sequence right before starcode starts

#define any new paths
barcode_quantify_folder =  Outfolder + "/barcode_quantify/"
path_to_feature_ref =      Outfolder + "/CellRanger_inputs/" + "FeatureReference_filtered.csv" 
path_to_sc_output_files = Outfolder + "/starcode_outputs/sc_output_counts_*.txt"

# Make any necessary folders
path_to_folders = [barcode_quantify_folder]

# checking whether folder/directory exists
for path_to_folder in path_to_folders:
    does_folder_exist(path_to_folder)


#---------------------------------------------------------------------------------------
#---------------------------------------------------------------------------------------



print("Running ")
print(" ")


#get all gDNA Read1 fastq file paths
all_sc_output = glob.glob(path_to_sc_output_files, recursive = True)

counter = 0
for grp in GSAMP:

    # Remove any files not in grp
    subset_sc_output = []
    for smp in grp:

        hold_path = []
        for path in all_sc_output:
            #Is the path one of the samples in group grp
            if smp in path:
                subset_sc_output.append(path)


    #------ Lets make the custom FeatureRef ------
    df_ref = pd.read_csv(path_to_feature_ref)
    ind_dict_ref = dict((k,i) for i,k in enumerate(df_ref['sequence']))

    for path in subset_sc_output:

        # Load in sc_output
    	df_output = pd.read_csv(path,header = None, sep='\t', engine='python')

        # Remove startseq
    	df_output[0] = df_output[0].str[len(strtseq):]

        # Use a dictionary to find intersections faster
    	ind_dict_output = dict((k,i) for i,k in enumerate(df_output[0]))
    	inter = set(df_ref['sequence']).intersection(df_output[0])
    	indices_ref = [ ind_dict_ref[x] for x in inter ]
    	indices_out = [ ind_dict_output[x] for x in inter ]

        # Determine what values belong to the current sc_output
    	hold_zero = [0]*len(df_ref['sequence'])
    	for x in inter:
    		hold_zero[ind_dict_ref[x]] = df_output[1][ind_dict_output[x] ] 
    	df_ref[path.split('/')[-1]] = hold_zero

    # Save the values in a column with apropriate name
    counter += 1
    filepath = barcode_quantify_folder +  'FeatureReference_filtered_group'+ str(counter) + '_counts.csv'
    df_ref.to_csv(filepath)  




#---------------------------------------------------------------------------------------
#---------------------------------------------------------------------------------------



print("     FeatureReference created ")
print(" ")



# Reads per million
def rpm_norm(reads):
    scaling_factor = np.nansum(reads)/1000000
    return [i/scaling_factor for i in reads]

def rpm_to_cells(m,b,rpm):
    return [k*m +b for k in rpm]

def get_r2_numpy(x, y):
    slope, intercept = np.polyfit(x, y, 1)
    r_squared = 1 - (sum((y - (slope * x + intercept))**2) / ((len(y) - 1) * np.var(y, ddof=1)))
    return slope, intercept, r_squared

#---------------------------------------Plots---------------------------------------------------

counter = 0
for grp in GSAMP:

    counter += 1
    print("     group"+ str(counter) +' :')
    print(" ")
    path_to_group_folder = barcode_quantify_folder + "plots_group"+ str(counter) + "/"
    os.mkdir(path_to_group_folder)


    #Load in Data
    filepath = barcode_quantify_folder  + '/FeatureReference_filtered_group'+ str(counter) + '_counts.csv'
    df_ref = pd.read_csv(filepath)

    #----------------------------Find and label spike ins------------------
    df_ref['spike_in_c'] = ['k']*len(df_ref.iloc[:,0])

    if spike_in_added == 'yes':
        for sp,spike_in_seq in enumerate(spike_in_seqs):
            df_ref.loc[df_ref['sequence'].str.match(spike_in_seq),'spike_in_c'] = spike_in_color[sp]


    #----------------------------Remove 0 and 1----------------------------
    print("         Removing zeros and ones")
    print(" ")

    df_plot = df_ref.iloc[:,-len(grp)-1:]

    df_plot = df_plot.replace(0, np.nan, regex=True)
    df_plot = df_plot.replace(1, np.nan, regex=True)

    #----------------------------Reads per million--------------------------
    for i in range(len(grp)):
        df_plot.iloc[:,i] = rpm_norm(df_plot.iloc[:,i])

    #-------------------------Get all combinations--------------------------
    all_combinations = itertools.combinations([g for g in range(len(grp))],2)

    print("             Plot:")
    #-----------------------------------------------------------------------
    #----------------------------Spike-ins ---------------------------------
    if spike_in_added== 'yes':
        #Keep Spike-ins color
        all_spike_in_color_in_grp = []
        all_spike_in_value_in_grp = []
        all_line_fit = []
        summary_spike_in = []

        with PdfPages(path_to_group_folder + 'Spike_in_group'+ str(counter)+'.pdf') as pdf:
            for i in range(len(grp)):

                print("                   Spike - "  + grp[i])
                plot_f = plt.figure()

                y = df_plot.iloc[:,i] 
                col_y = df_plot.iloc[:,-1]

                new_y =[]
                new_col_y =[]

                for n,m in enumerate(col_y):
                    if m in spike_in_color:
                        #Clean up any possibly lost barcodes
                        if np.isnan(y[n]) == False:
                            new_y.append(y[n])
                            new_col_y.append(m)
                summary_spike_in.append(len(new_col_y))

                #Make the order of spike-ins found the same as the spike_in_values
                x = new_col_y
                test_list = list(zip( x,range(len(x))))
                sort_order = spike_in_color # initializing sort order
                test_list.sort(key = lambda i: sort_order.index(i[0])) 

                new_y = [new_y[a] for _,a in test_list]

                #Remove colors that had a nan
                new_spike_in_color = []
                new_spike_in_value = []
                remove_new_col_y = new_col_y
                for c,col in enumerate(spike_in_color):
                    if col in remove_new_col_y:
                        remove_new_col_y.remove(col)

                        new_spike_in_color.append(col)
                        new_spike_in_value.append(spike_in_value[c])

                print("new_y:", new_y)
                print("new_spike_in_value:", new_spike_in_value)
                
                #Fit line
                m, b = np.polyfit(new_y,new_spike_in_value, 1)

                all_line_fit.append([m,b])

                x_line = np.linspace(0, max(new_y), 2)
                y_line = [k*m+b for k in x_line]


                #Plot
                scatt = plt.scatter(new_y,new_spike_in_value, s=3, c=new_spike_in_color)
                
                plt.suptitle('Group'+ str(counter), y=0.99);
                plt.title(grp[i])
                plt.xlabel('[RPM]')
                plt.ylabel('[Cells]')
                plt.plot(x_line,y_line,'-',color = 'k', alpha=0.2)


                for c,cop in enumerate(new_spike_in_color):
                    plt.plot([], [], '.'+cop, label=new_spike_in_value[c])
                plt.legend(title="Spike-ins",bbox_to_anchor=(1.05, 1.0), loc='upper left')

                plt.tight_layout()

                #Save
                pdf.savefig(plot_f);
                plt.close()

                all_spike_in_color_in_grp.append(new_spike_in_color)
                all_spike_in_value_in_grp.append(new_spike_in_value)


    #------------------------------------------------------------------------------------
    #----------------------------Barplot-------------------------------------------------             
    summary_filter_in = []
    print("             ")
    with PdfPages(path_to_group_folder + 'barplot_group'+ str(counter)+'.pdf') as pdf:

        
        for i in range(len(grp)):

            print("                   Bar - " + grp[i])
            plot_f = plt.figure()

            y = df_plot.iloc[:,i]
            col_y = df_plot.iloc[:,-1]

            new_y =[]
            new_col_y =[]
            for x,m in enumerate(y):
                if np.isnan(m) == False:
                    new_y.append(m)
                    new_col_y.append(col_y[x])
            summary_filter_in.append(len(new_y))

            if spike_in_added == 'yes':
                #Convert from RPM to Cells
                m,b = all_line_fit[i]
                new_y = [k*m for k in new_y]

            #Order
            index_y = sorted(range(len(new_y)), reverse=True,key=new_y.__getitem__)
            new_col_y = [new_col_y[k] for k in index_y]
            new_y = sorted(new_y, reverse=True)

            #Plot
            plt.bar(range(len(new_y)),new_y, color=new_col_y)
            plt.suptitle('Group'+ str(counter), y=0.99);
            plt.title(grp[i])
            plt.xlabel('Barcodes')
            
            if spike_in_added == 'yes':
                plt.ylabel('Cells')
                new_spike_in_color = all_spike_in_color_in_grp[i]
                new_spike_in_value = all_spike_in_value_in_grp[i]

                for c,cop in enumerate(new_spike_in_color):
                    plt.plot([], [], '.'+cop, label=new_spike_in_value[c] )
                plt.legend(title="Spike-ins",bbox_to_anchor=(1.05, 1.0), loc='upper left')
            else:
                plt.ylabel('RPM')

            plt.tight_layout()

            #Save
            pdf.savefig(plot_f);
            plt.close()

 
    
    #---------------------------------------------------------------------
    #----------------------------Scatter plot ----------------------------
            
    print("             ")
    with PdfPages(path_to_group_folder+ 'scatter_group'+ str(counter)+'.pdf') as pdf:

        
        for combination in all_combinations:

            print("                   Scatt - " + grp[combination[0]] + " and " + grp[combination[1]])
            plot_f = plt.figure()

            x = df_plot.iloc[:,combination[0]] 
            y = df_plot.iloc[:,combination[1]] 
            col_y = df_plot.iloc[:,-1]

            if spike_in_added == 'yes':
                #Convert from RPM to Cells
                m_x,b_x = all_line_fit[combination[0]]
                x = rpm_to_cells(m_x,0,x)

                m_y,b_y = all_line_fit[combination[0]]
                y = rpm_to_cells(m_y,0,y)

            
            new_x =[]
            new_y =[]
            no_pike_x = []
            no_pike_y = []
            new_col_y =[]
            for i,m in enumerate(x):
                if np.isnan(m) == False:
                    if np.isnan(y[i]) == False:


                        new_x.append(m)
                        new_y.append(y[i])
                        new_col_y.append(col_y[i])

                        if col_y[i] == 'k':    # without spike in
                            no_pike_x.append(m)
                            no_pike_y.append(y[i])

            # fit line and get r^2
            slope, intercept, r_squared = get_r2_numpy(np.array(no_pike_x), np.array(no_pike_y) ) 

            x_fit = [0,max(no_pike_x)]
            y_fit =[k*slope +intercept for k in x_fit]

            #Plot
            plt.scatter(new_x,new_y, s=3, c=new_col_y)
            plt.plot(x_fit,y_fit,'-',color = 'k', alpha=0.2)
            plt.suptitle('Group'+ str(counter), y=0.99);
            plt.title(grp[combination[0]] + " and " + grp[combination[1]])

            fit_name = "R^2=" +"{:.2f}".format(r_squared) + "     slope=" + "{:.2f}".format(slope)+ "     y-intercept=" + "{:.2f}".format(intercept)

            plt.text(0,np.max(new_y),fit_name)

            if spike_in_added == 'yes':
                plt.xlabel(grp[combination[0]] + ' [Cells]')
                plt.ylabel(grp[combination[1]] + ' [Cells]')
                for c,cop in enumerate(spike_in_color ):
                    plt.plot([], [], '.'+cop, label=spike_in_value[c])
                plt.legend(title="Spike-ins",bbox_to_anchor=(1.05, 1.0), loc='upper left')
            else:
                plt.xlabel(grp[combination[0]] + ' [RPM]')
                plt.ylabel(grp[combination[1]] + ' [RPM]')

            plt.tight_layout()

            #Save
            pdf.savefig(plot_f);
            plt.close()

    with open(path_to_group_folder+ 'summary_group'+ str(counter)+".txt","w") as write_file:
        write_file.write("Summary for group"+ str(counter) + "\n")
        write_file.write(" " + "\n")

        if spike_in_added == 'yes':
            write_file.write("Number of spike-ins added: ="+ "{:.2f}".format(len(spike_in_value)) + "\n")
            for i,g in enumerate(grp):
                write_file.write(g+ "\n")
                write_file.write("      Number of spike-ins captured: = "  + "{:.2f}".format(summary_spike_in[i]) + "\n")
                write_file.write("      Percentage of spike-ins captured: = %"  + "{:.2f}".format(100*summary_spike_in[i]/len(spike_in_value)) + "\n")
            write_file.write(" " + "\n")

        write_file.write("Total barcodes captured = %d" % len(df_ref.iloc[:,0]) + "\n")
        for i,g in enumerate(grp):
            write_file.write(g+ "\n")
            write_file.write("      Barcodes with counts greater than 1: = "  + "{:.2f}".format(summary_filter_in[i]) + "\n")
            write_file.write("      Percentage of barcodes with counts greater than 1: = %"  + "{:.2f}".format(100*summary_filter_in[i]/len(df_ref.iloc[:,0])) + "\n")
        write_file.write(" " + "\n")


print("     Done :D ")
print(" ")







