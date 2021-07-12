## script name: calc_mean_operons.py

#################################################################

## Description:
## script to calculate the Differential expression
## Between control versus RIF-stress operons
## creates files with columns containing the operon name,
## plus the strain name under control versus exp condition
# eg. M2_CONT, M2_EXP
## Prints out the Fold change (FC)

## NB!! This file is hardcoded to only work on files produced by
# COSMO with format strain-name_D*_d*_F*_f*.csv, eg.
# Strain8_D4_D2_F5_f10_csv, where D, d, F and f,
# correspond to the expression values and fold differences,
# as per the README/usage of COSMO on Github:
# https://github.com/SANBI-SA/COSMO

##################################################################

import os.path
import glob
import string
import numpy as np
from statistics import mean
import math
from operator import itemgetter

operon_file_path = input("Please add the name of the path to the operon file 'Combined_operon_list.txt'\n eg. /home/user/Desktop/Combined_operon_list.txt:\n")
operon_file = open(operon_file_path)
file_dir = input("Please add the name of the path to the csv files produced by COSMO: \n eg. /home/user/Desktop/Strain29_Strain36/*.csv:\n")

#Test file paths
# /home/me/Desktop/PhD_Project/data/Python_scripts/NEW_correct_combined_operon_bins.txt
# operon_file_path  = "/home/me/Desktop/PhD_Project/data/Python_scripts/NEW_correct_combined_operon_bins.txt"
# operon_file_path  = "/home/me/Desktop/PhD_Project/data/software/Operon_algorithms/Python_new/tb_operon_detection/2_operons_no_genes.txt"
# file_dir = "/home/me/Desktop/PhD_Project/data/software/Operon_algorithms/Python_new/tb_operon_detection/output/C5/*.csv"



file_paths = glob.glob(file_dir)
outfile_name = f"{file_dir.split('/')[-2]}_FC_operons_control_vs_RIF.txt"
outfile = open(outfile_name, 'w')

operon_names = []
for lines in operon_file:
    operons = lines.strip().split("\t")[0]
    if operons.startswith("Rv"):
        operon_names.append(operons)

dict_key_list = []
list_keys, dict_name_before, dict_name_after = [], {},{} # the dict_name_before will contain the operon name and CONTROL expr value, while the after will contain the RIF

for operon in operon_names:
    for filenames in file_paths:
        file_name = os.path.basename(filenames)
        strain_original = file_name.rstrip(".csv")
        dict_keys = f'{operon}_{strain_original}' # make the dictionary key the "operon name_strain_name, eg. Rv0046c-Rv0047c_M24_L2_D2_d1_F5_f5"
        dict_key_list.append(dict_keys)
# print(dict_key_list)

print("\nYour files have been successfully read in. \nPlease be patient. The next step will take a while.\n\n")

###########################################################
# Loop through the dictionary keys and csv files
# containing all the genes and their expression levels
###########################################################
for dict_key in dict_key_list:
    op = dict_key.split("_")[0]
    # strain_key = dict_key
    strain = ("_").join(dict_key.split("_")[1:]) #have to change this to 1 for real files. It's 2 for test files
    # print(strain)
    # print(dict_key)
    for filenames in file_paths:
        file_name = filenames.split("/")[-1].strip(".csv") # have to remove #strip("test_"). for real files!!!!!!!!!!
        # print (file_name)
        # break
        opened_files = open(filenames)
        for file_line in opened_files:
            if len(file_line.split(',')) > 7: #csv files with more than 7 columns contain operons
                if not file_line.split(',')[7].startswith("Rv"):
                #file_line.split(',')[7] == "" or file_line.split(',')[7].startswith("EBG") of : #none of our operons have these genes
                    continue
                gene = file_line.split(',')[7] # the genes are in the 7th column
                # gene_list.append(gene)
                strand = file_line.split(',')[4]
                expres_level= float(file_line.split(',')[8])

                gene_int = int(gene.strip(string.ascii_letters)) #convert the gene name to an integer
                # print (gene, gene_int, expres_level, strain)
                # break

#                 ########################################################################################################
#                 # add key (strain_name + operon_name) +  value (expression levels of genes) to the relevant dictionary
#                 #########################################################################################################
#
                #operons and genes on reverse strand have a "c" in them
                #check if gene lies within or on the boundary of an operon
                # strain name comes from dictionary key and file name from the actual file name
                if "c" in op and "c" in gene and \
                file_name == strain and \
                (gene_int == int(op.split("-")[0].strip(string.ascii_letters)) \
                or gene_int == int(op.split("-")[1].strip(string.ascii_letters)) \
                or int(op.split("-")[0].strip(string.ascii_letters)) < gene_int < int(op.split("-")[1].strip(string.ascii_letters))):
                    # print(gene, op, expres_level, strain, file_name)
#
                    dict_name_before.setdefault(dict_key, []) #Set the default value type of the dictionary key, to a list
                    dict_name_after.setdefault(dict_key, [])

#                   #only start adding to the dict if gene name is part of the operon...."CONT"
                    control = ["C4","M20","M23","M26","M29","M32","M35","C5_CONT","C6_CONT","R4370_CONT"]
                    rif = ["M21","M24","M27","M30","M33","M36","C5_RIF","C6_RIF","R4370_RIF"]
                    for element in control:
                        if element in dict_key and dict_key in dict_name_before:
                            dict_name_before[dict_key].append(expres_level)
                            # print(gene, op, expres_level, strain, file_name, "cont")
                    for element in rif:
                        # print (element)
                        #if the element, eg. M21, is in the dict_key and....
                        if element in dict_key and dict_key in dict_name_after:  #only start adding to the dict if gene name is part of the operon...."RIF"
                            dict_name_after[dict_key].append(expres_level)
                            # print(gene, op, expres_level, strain, file_name, "rif")

#
                elif not "c" in op and not "c" in gene and \
                file_name == strain and \
                (gene_int == int(op.split("-")[0].strip(string.ascii_letters)) \
                or gene_int == int(op.split("-")[1].strip(string.ascii_letters)) \
                or int(op.split("-")[0].strip(string.ascii_letters)) < gene_int < int(op.split("-")[1].strip(string.ascii_letters))):
                    # print(gene_int, op, expres_level, strain)

                    dict_name_before.setdefault(dict_key, []) #Set default value of key as a list
                    dict_name_after.setdefault(dict_key, [])
            #
                    control = ["C4","M20","M23","M26","M29","M32","M35", "C5_CONT","C6_CONT","R4370_CONT"]
                    rif = ["M21","M24","M27","M30","M33","M36","C5_RIF","C6_RIF","R4370_RIF"]
                    for element in control:
                        if element in dict_key and dict_key in dict_name_before: # in other files like C6, I just write "if 'CONT' in dict_key and dict_key in dict_name_before"
                            dict_name_before[dict_key].append(expres_level)
                            # print(gene, op, expres_level, strain)
                    for element in rif:
                        if element in dict_key and dict_key in dict_name_after:  #only start adding to the dict if gene name is part of the operon...."RIF"
                            dict_name_after[dict_key].append(expres_level)
                            # print(gene, op, expres_level, strain)

#remove key and value if dict is empty
dict_name_before = dict(filter(itemgetter(1),dict_name_before.items()))
dict_name_after = dict(filter(itemgetter(1),dict_name_after.items()))

# print("dict_name_before", dict_name_before)
# print ("\n")
# print("dict_name_after", dict_name_after)

print("Your data has been extracted and control and experimental values are now being calculated and compared.\nPlease be patient. This may take a while.\n\n")
######################################################################################
# compare the keys in the control dictionary to they keys in the RIF dictionary
# if the keys match, print out the average expression value under each treatment
# next to each other
#####################################################################################


for key_cont, value_cont in dict_name_before.items():
    # print(key_cont, value_cont)
    key_cont_split = key_cont.split("_")
    if "CONT" in dict_key or "RIF" in dict_key:
        continue
    key_before = ("_").join(key_cont_split[0:3]) #+ "_" + (("_").join(key_cont_split[2:])).strip("L1_D2_d1_F5_f5")  #in real data the it should be [0:2]
    # print (key_before, "before")
    for key_rif, value_rif in dict_name_after.items():
        key_rif_split = key_rif.split("_")
        key_after = ("_").join(key_rif_split[0:3]) #+ "_" + ("_").join(key_rif_split[2:]).strip("L1_D2_d1_F5_f5")
        # print(key_after, "after")
        # print(int(key_before.split("_")[1].strip("M")), int(key_after.split("_")[1].strip("M"))-1 )

        if key_before.split("_")[0] == key_after.split("_")[0] \
        and int(key_before.split("_")[1].strip("M")) == (int(key_after.split("_")[1].strip("M")) -1) \
        and key_before.split("_")[2] == key_after.split("_")[2]:
            op_name = f'{key_cont.split("_")[0]}'
            strain_name_CONT = f'{"_".join(key_cont.split("_")[1:4])}'
            strain_name_RIF = f'{"_".join(key_rif.split("_")[1:4])}'
            op_ave_CONT = round(mean(value_cont),2)
            op_ave_RIF = round(mean(value_rif),2)
            FC =  str(round((op_ave_RIF - op_ave_CONT)/(op_ave_CONT),2)) + "X"  #fold change
            outline = f'{op_name} {strain_name_CONT} {strain_name_RIF} {op_ave_CONT} {op_ave_RIF} {FC}\n'
            outfile.write(outline)
            # print(outline)


# ########################################################################
# # The code below is for the files that have an easy
# # RIF or CONT in the file_name
# # such as C6_RIF etc.
# ########################################################################


for key_cont, value_cont in dict_name_before.items():
    if "CONT" in dict_key or "RIF" in dict_key:
        # print(file_name, strain, op, key_cont, value_cont)
        key_cont_split = key_cont.split("_")
        key_before = ("_").join(key_cont_split[0:2]) + "_" + (("_").join(key_cont_split[2:])).strip("CONT")  #in real data the it should be [0:2]

        for key_rif, value_rif in dict_name_after.items():
            key_rif_split = key_rif.split("_")
            key_after = ("_").join(key_rif_split[0:2]) + "_" + ("_").join(key_rif_split[2:]).strip("RIF")
            if key_before == key_after:

                op_name = f'{key_cont.split("_")[0]}'
                strain_name_CONT = f'{"_".join(key_cont.split("_")[1:4])}'
                strain_name_RIF = f'{"_".join(key_rif.split("_")[1:4])}'
                op_ave_CONT = round(mean(value_cont),2)
                op_ave_RIF = round(mean(value_rif),2)
                FC =  str(round((op_ave_RIF - op_ave_CONT)/(op_ave_CONT),2)) + "X"  #fold change
                outline = f'{op_name} {strain_name_CONT} {strain_name_RIF} {op_ave_CONT} {op_ave_RIF} {FC}\n'
                outfile.write(outline)
                # print(outline)

print("The file" , outfile_name,  "has been printed to the directory where this script is kept")
operon_file .close()
opened_files.close()
outfile.close()



#############################################################################################################################
