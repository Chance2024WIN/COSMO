
###################################################################################################################################
# Script name: ave_cov_operons.py

##### Python script to compare a wiggle file to the coordinates of UTRs operon genes.
##### If the position lies within the operon UTR, it prints out the position with its coverage to a file.


# Requirements:
# 1. Path to wiggle file eg. /home/tracey/Desktop/PhD_Project/data/Python_scripts/M2_cov_files/Unnormalized/plus/M2_Forward_L006_001_UNNORM.wig
# 2. Path to operon coordinates, eg.
# /home/tracey/Desktop/PhD_Project/data/Python_scripts/operon_list_plus_strand.txt

## NB!!!
# For Forward strand wiggle list, use plus strand operon list

###################################################################################################################################

import itertools
import collections
import numpy
from collections import defaultdict

################################################
# locate all the necessary files and open them
################################################

path_file_wiggle = input("Enter the path to the wiggle file:\n")
path_operons = input("\nEnter the path to the operon list:\n")

###################################################################################
#Format of file
#Rv0097-Rv0102	105323	119699
##################################################################################

open_wiggle_file = open(path_file_wiggle)
wiggle_list = []
for entries in open_wiggle_file:
    if entries.startswith('variableStep'):
        continue
    list_line = entries.strip().split('\t')
    #print (list_line[0:2])
    wiggle_list.append(list_line)
    # e.g. of wiggle_list  = [['1', '10.32']]

operon_list = []
open_operon_file = open(path_operons)
for operon_file in open_operon_file:
    for data in open_operon_file:
        if data.startswith("betw") or data.startswith("Gene") or data.startswith("#"):
            continue
        operon_list_lines = data.strip().split('\t')
        operon_list.append(operon_list_lines)
    #print(operon_list[0:1], out_file_name)
#e.g. of operon_list = [['Rv0097-Rv0101', '105323', '117539'], ['Rv0101-Rv0102', '117539', '119699'],....]]


count = 0
out_file_name = path_file_wiggle.split("/")[-1] + "_UTR_read_coverages.txt"
output_file = open(out_file_name, "w")

print("\n This may longer than 10 minutes. Please be patient.\n\n The wiggle file is millions of lines long.\n")

for item in wiggle_list:
    position = int(item[0])
    coverage = item[1]
    #print (position, coverage)
    for bin in operon_list:
        operon_name = bin[0]
        operon_start = int(bin[1])
        operon_end = int(bin[2])
        # print (operon_start, operon_end)

        if (position >= (operon_start-300) and position < operon_start):
            pos_name = "5UTR_" + operon_name
            count += 1
            output = operon_name + "\t" + pos_name + "\t" + str(operon_start) + "\t" + str(operon_end) + "\t" +  str(position) + "\t" + str(coverage) + "\n"

            output_file.write(output)
        elif (position > operon_end  and position <= (operon_end + 300)):
            count += 1
            pos_name = "3UTR_" + operon_name
            output = operon_name + "\t" + pos_name + "\t" + str(operon_start) + "\t" + str(operon_end) + "\t" +  str(position) + "\t" + str(coverage) + "\n"
            output_file.write(output)
#
output_file.close()

####################################################################
# Open the file with the read counts for all genomic regions
# place data into list to be used by dictionary later
####################################################################


list_coverages = []
open_operon_cover_file = open(out_file_name)
averages_file_name = path_file_wiggle.split("/")[-1].strip(".txt") + "_averages_UTRs.txt"
header = "Operon" + "\t" + "Region_name" + "\t" +  "Start" + "\t" +  "End"  + "\t" + "Position" + "\t" + "Average_coverage" + "\n"
#out_file_name + "_averages_ALL_genomic_regions.txt"
averages_file = open(averages_file_name, "w")
averages_file.write(header)

average_dict = {}
for line_list in open_operon_cover_file:
    cov_item = line_list.strip().split("\t")
    operon_name = cov_item[0]
    genomic_pos_name = cov_item[1]
    pos_start = cov_item[2]
    pos_end = cov_item[3]
    position_coord = cov_item[4]
    pos_coverage = float(cov_item[5])
    key_for_dict = f"{operon_name}:{genomic_pos_name}:{pos_start}:{pos_end}"
    list_coverages.append(genomic_pos_name)
    list_coverages.append(pos_coverage)
# print(average_dict)

####################################################################
# Use a dictionary to calculate the average of each genomic region
####################################################################
    if key_for_dict not in sorted(average_dict):
        average_dict[key_for_dict] = [pos_coverage]
    else:
        average_dict[key_for_dict].append(pos_coverage)

for k,v in average_dict.items():
    all_averages = k.split(":")[0] + "\t" + k.split(":")[1] + "\t" + k.split(":")[2] + "\t" + k.split(":")[3] + "\t" + str(round(numpy.mean(v))) + "\n"
    averages_file.write(all_averages)
    # print (all_averages)

print("\nYour file has successfully printed to the directory where this script is stored\n")
