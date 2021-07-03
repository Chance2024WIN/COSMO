## script name: dict_no_strains_pred_operon_TP_FP_FN.py

#####################################################################################################

## Description:
## script to create a dictionary from files containing operon predictions,
## per operon and per strain, that is output of script: predicted_operons_counts_per_lineage.py

## Requirements:
# path to files ending in *calls_&_cov_all_lineages.txt

## Output
## TWO output files:
# 1. "_calls_per_lineage_&_condition.txt"
# AND
# 2. All files ending in: *CONT_calls_&_cov_all_lineages.txt

## keys = operon name
## values = isolate/family name...
## with count of TP, TN or FP calls per isolate/family

####################################################################################################

import os
import string
import glob
import collections
import pandas as pd

print("Please add the path to the directory containing all files with extension: \n e.g. home/user/Desktop/*D2*calls_&_cov_all_lineages.txt\n")
call_files = input()

# call_files = "/home/tracey/Desktop/PhD_Project/data/software/Operon_algorithms/Python_new/tb_operon_detection/RESULTS/test/*D2*calls_&_cov_all_lineages.txt"
call_files_path = glob.glob(call_files)
outfile_all_name =  "Intersection_calls_ALL_lineages.txt"
# outfile_control = "All_control_calls.txt"
# outfile_exper = "All_exper_calls.txt"

outfile_all_strains = open(outfile_all_name, "w")
# outfile_all_control = open(outfile_control, "w")
# outfile_all_exp = open(outfile_exper,"w")


TP_list_ALL, FP_list_ALL, FN_list_ALL = [],[],[]
TP_ALL_strains = []
counts_dict_TP, counts_dict_FP, counts_dict_FN = dict(), dict(), dict()

CONTROL_dict = {}
EXPERIMENTAL_dict = {}

for files in call_files_path:
    TP_list, FP_list, FN_list = [],[],[]
    family = os.path.basename(files).rstrip("D2_calls_&_cov_all_lineages.txt").strip()
    # print(family)
    outfile_name = family + "_calls_per_lineage_&_condition.txt"

    open_call_files = open(files)
    TP_list, FP_list, FN_list = [],[],[]
    for line_files in open_call_files:
        if line_files.startswith("Call"):
            continue
        column_all = line_files.strip().split("\t")
        strain = column_all[6].split("_")[0] + "_" + column_all[6].split("_")[1]#.strip("_D2_d1_F5_f5")
        # print(strain)

        call = column_all[0]
        # print(isolate)
        if call == "TP":
            TP_operon = column_all[1] #+ "\t" + column_all[2]
            TP_list.append(TP_operon)
            TP_list_ALL.append(TP_operon)
            if TP_operon not in counts_dict_TP: #
                counts_dict_TP[TP_operon]= [strain] #{family} for making a set
            else:
                counts_dict_TP[TP_operon].append(strain) #.add(family)



            if strain.split("_")[0].startswith("C4") or strain.split("_")[0].startswith(("M20","M23","M26","M29","M32","M35")) or strain.split("_")[1].startswith("CONT"):
                if TP_operon not in CONTROL_dict:
                    CONTROL_dict[TP_operon]= [strain]
                else:
                    CONTROL_dict[TP_operon].append(strain)


            elif strain.split("_")[0].startswith(("M21","M24","M27","M30","M33","M36")) or strain.split("_")[1].startswith("RIF"):
                if TP_operon not in EXPERIMENTAL_dict:
                    EXPERIMENTAL_dict[TP_operon]= [strain]
                else:
                    EXPERIMENTAL_dict[TP_operon].append(strain)
                # Experimental.append(line_files)
                # exp_line = line_files
                # outfile_all_exp.write(exp_line)

        elif call == "FP":
            FP_operon = column_all[1] #+ "\t" + column_all[2]
            FP_predicted_operon = column_all[2]
            FP_list.append(FP_operon)
            FP_list_ALL.append(FP_operon)
            # print(FP_list_ALL)
            # break

            if FP_operon not in counts_dict_FP: #
                counts_dict_FP[FP_operon]= [strain]  #if list is empty, open a new list with the first value          #{family} #this can be used in case I want a set
            else:
                counts_dict_FP[FP_operon].append(strain) # append to list if list already exists                      #add(family) # this is to add to a set

        elif "FN" in call:
            FN_operon =  column_all[1] #+ "\t" + column_all[2]
            FN_predicted_operon = column_all[2]
            FN_list.append(FN_operon)
            FN_list_ALL.append(FN_operon)

            if FN_operon not in counts_dict_FN: #
                counts_dict_FN[FN_operon]= [strain]#{family}
            else:
                counts_dict_FN[FN_operon].append(strain) #.add(family)

##############################################################################################

## add operon as key to dictionary
## add each MTB family as a set value to the operon (key) in the dictionary

## which TP, FP and FN calls are common and unique across ALL lineages and conditions

##############################################################################################

ouline_dict_heading_TP = "Which families called the operon below as TP?\n"
outfile_all_strains.write(ouline_dict_heading_TP)
for key in counts_dict_TP:
    TP_collection = collections.Counter (counts_dict_TP[key])
    ouline_dict_TP = f"{key}   {TP_collection}   {len(counts_dict_TP[key])}" + "\n" # print out the key and value from dictionary, as well as the length of the set. In this case, the value is a set
    outfile_all_strains.write(ouline_dict_TP)
outfile_all_strains.write("\n\n")


heading_counts_TP_ALL = "Number of isolates which called the operons as TP\n"
outfile_all_strains.write(heading_counts_TP_ALL)
counts_TP_ALL = pd.Series(TP_list_ALL).value_counts()#+ "\t" + "TP"
outfile_all_strains.write(str(counts_TP_ALL))
outfile_all_strains.write("\n\n")



ouline_dict_heading_FP = "Which families called the operon below as FP?\n"
outfile_all_strains.write(ouline_dict_heading_FP)
for key in counts_dict_FP:
    FP_collection = collections.Counter (counts_dict_FP[key])
    ouline_dict_FP = f"{key}    {FP_collection}   {len(counts_dict_FP[key])}" + "\n"
    outfile_all_strains.write(ouline_dict_FP)
outfile_all_strains.write("\n\n")

heading_counts_FP_ALL = "Number of isolates which called the operons as FP\n"
outfile_all_strains.write(heading_counts_FP_ALL)
counts_FP_ALL = pd.Series(FP_list_ALL).value_counts()
outfile_all_strains.write(str(counts_FP_ALL))
outfile_all_strains.write("\n\n")


ouline_dict_heading_FN = "Which families called the operon below as FN?\n"
outfile_all_strains.write(ouline_dict_heading_FN)
for key in counts_dict_FN:
    FN_collection = collections.Counter (counts_dict_FN[key])
    ouline_dict_FN = f"{key}   {FN_collection}   {len(counts_dict_FN[key])}" + "\n"
    outfile_all_strains.write(ouline_dict_FN)
outfile_all_strains.write("\n\n")


heading_counts_FN_ALL = "Number of isolates which called the operons as FN\n"
outfile_all_strains.write(heading_counts_FN_ALL)
counts_FN_ALL = pd.Series(FN_list_ALL).value_counts()
print("The file: Intersection_calls_ALL_lineages.txt, has been successfully printed to the directory where this script is stored")
outfile_all_strains.write(str(counts_FN_ALL))
outfile_all_strains.write("\n\n")
