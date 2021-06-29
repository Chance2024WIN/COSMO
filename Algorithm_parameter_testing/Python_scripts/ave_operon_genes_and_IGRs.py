###################################################################################################################################
# Script name: ave_operon_genes_and_IGRs.py

##### Python script to compare a wiggle file to the coordinates of operon genes and intergenic regions between operon genes.
##### If the position lies within the operon, it prints out the position with its coverage to a file with the operon_name.
##### It also prints out the average coverage to STDOUT

# Requirements:
# 1. Path to wiggle file eg. /home/tracey/Desktop/PhD_Project/data/Python_scripts/M2_cov_files/Unnormalized/plus/M2_Forward_L006_001_UNNORM.wig
# 2. Path to genomic region coordinates, eg.
# /home/tracey/Desktop/PhD_Project/data/Python_scripts/input_data_files/genic_IGR_coordinates_operon_plus.txt
#operon_list_plus_strand.txt

## NB!!!
# For Forward strand wiggle list, use plust strand operon list

###################################################################################################################################
import itertools
import collections
import numpy
from collections import defaultdict

path_file_wiggle = input("Enter the path to the wiggle file:\n")
path_genomic_regions = input("\nEnter the path to the operon list:\n")

###################################################################################
#Format of file
#Rv0101-Rv0101	117539	117539
##################################################################################

#dist_betw_operons_plus_strand.txt"

open_wiggle_file = open(path_file_wiggle)
open_genomic_regions_file = open(path_genomic_regions)

operon_list = []
for data in open_genomic_regions_file:
    if not data.startswith("Rv"): #or data.startswith("Gene"):
        continue

    operon_list_lines = data.strip().split('\t')
    operon_list.append(operon_list_lines)

wiggle_dict = {}
wiggle_list = []
for entries in open_wiggle_file:
    if entries.startswith('variableStep'):
        continue
    list_line = entries.strip().split('\t')
    #print (list_line[0:2])
    wiggle_list.append(list_line)


count = 0
out_file_name = path_file_wiggle.split("/")[-1] + "_operon_read_coverages.txt"
output_file = open(out_file_name, "w")

print("\n This may longer than 10 minutes. Please be patient.\n\n The wiggle file is millions of lines long.\n")

for item in wiggle_list:
    position = int(item[0])
    coverage = item[1]
    #print (position, coverage)
    for bin in operon_list:
        if (position >= int(bin[2]) and position <= int(bin[3])) or (position == int(bin[2])) or (position == int(bin[3])):
            count += 1
            output = str(bin[0]) + "\t"  + str(bin[1]) + "\t" + str(bin[2]) + "\t" +  str(bin[3]) + "\t" + str(position) + "\t" + str(coverage) + "\n"
            output_file.write(output)

output_file.close()
# print (count)

####################################################################
# Open the file with the read counts for all genomic regions
# place data into list to be used by dictionary later
####################################################################


list_coverages = []
open_operon_cover_file = open(out_file_name)
averages_file_name = out_file_name + "_averages_CDS_IGRs.txt"
header = "Operon" + "\t" + "Region_name" + "\t" +  "Start" + "\t" +  "End"  + "\t" + "Average_coverage" + "\n"
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
    pos_cov = cov_item[4]
    pos_coverage = float(cov_item[5])
    key_for_dict = f"{operon_name}:{genomic_pos_name}:{pos_end}:{pos_end}"
    list_coverages.append(genomic_pos_name)
    list_coverages.append(pos_coverage)

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
