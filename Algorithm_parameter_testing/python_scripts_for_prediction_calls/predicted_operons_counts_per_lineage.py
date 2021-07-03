## script name: predicted_operons_counts_per_lineage.py

#######################################################################################

# script that goes through all COSMO predicted files

## OUTPUT
# Output file with extention "_calls_&_cov_all_lineages.txt" : The Prediction (TP, FP, FN) **PER OPERON** for all isolates/strains
# Output file with extention "_counts_per_isolate" : Prediction **PER STRAIN ** with count of how many TPs, FPs and FNs were predicted per strain

# Requirements:
# 1. The csv files produced by COSMO containing the operon predictions
# 2. The operon list file containing all Mtb operons and their positions and genes



########################################################################################

import string
import os.path
import glob

print("Please add the path to the directory containing all the csv files produced by COSMO\n e.g. /home/user/Desktop/COSMO_files/* \n")
dir_pred_operon_files = input()

print("Please add the path to the file named: '50_combined_operon_list.txt' \n e.g. /home/user/Desktop/COSMO_files/50_combined_operon_list.txt \n")
valid_operon_list = input()

pred_operons_path = glob.glob(dir_pred_operon_files)
outfile_name = dir_pred_operon_files.split("/")[-2]  + dir_pred_operon_files.split("/")[-1].split("*")[1] + "_calls_&_cov_all_lineages.txt"
outfile_name_counts = dir_pred_operon_files.split("/")[-2] + dir_pred_operon_files.split("/")[-1].split("*")[1] + "_counts_per_isolate.txt"
outfile = open(outfile_name, "w")
outfile_counts = open(outfile_name_counts, "w")

###############################################################
# add the lines from each file with VALIDATED operons to a list
##############################################################
operon_lists, validated_operons = [],[]
list_of_list_operon_genes = []
list_operon_genes = []

header = "Call" + "\t" + "Operon" + "\t" + "Predicted_operon" + "\t" + "Operon_length" + "\t" + "Pred_operon_length" + "\t" + "Coverage" + "\t" + "Strain" + "\n"
outfile.write(header)


open_operon_files = open(valid_operon_list)
for operon_line in open_operon_files:
    if not operon_line.startswith("Operon"):
        operon_line = operon_line.strip().split("\t")
        operon_lists.append(operon_line)
        validated_operons.append(operon_line[0])

###############################################################
# add the lines from each file with PREDICTED operons to a list
##############################################################

Total_TPs, Total_FNs, Total_FN_counterpart, Total_FPs, Total_FP_counterpart = [],[],[],[],[]

for pred_operon_files in pred_operons_path:

    TP_list, TP_operon_counterpart = [],[]
    FP_list, FP_operon_counterpart = [],[]
    FN_list, FN_operon_counterpart = [],[]
    FN_single_gene_list, FN_single_gene_counterpart = [], []
    pred_operon_list = []

    pred_operon_filename  = os.path.basename(pred_operon_files)#.split(".")[0]
    strain = pred_operon_filename.rstrip(".csv")
    CDS_cov = pred_operon_filename.split("_")[-4].strip("D")
    IGR_cov = pred_operon_filename.split("_")[-3].strip("d")
    CDS_cov_diff = pred_operon_filename.split("_")[-2].strip("F")
    IGR_cov_diff = pred_operon_filename.split("_")[-1].lstrip("f").rstrip(".txt")
    # print (strain, CDS_cov, IGR_cov, CDS_cov_diff, IGR_cov_diff)
    pred_operons = open(pred_operon_files)
    for row in pred_operons:
        row = row.strip().split(",")
        pred_operon_list.append(row)

    for operon_column in operon_lists:
        operon = operon_column[0]
        operon_len = operon_column[3]
        operon_start = int(operon_column[4])
        operon_end = int(operon_column[5])
        operon_list = operon_column[6].strip().split(",")
        list_of_list_operon_genes.append(operon_list)
        for split_list in list_of_list_operon_genes:
            list_operon_genes.append(split_list)

# ###########################################################
# # Loop through the list containing the real operons and
# # check if they were predicted by algorithm
# ###########################################################

    # print (operon,pred_operon, strain)
    # break
        for elem in pred_operon_list:
            if elem[0]:
                pred_operon = elem[0]
            if elem[4]:
                strand = elem[4]
            if elem[5]:
                coverage = elem[5]
            if elem[6]:
                pred_gene_count = elem[6]

# ##################################################################
# # if the real operon matches the predicted operon, print TPs
    # ##################################################################

                if operon == pred_operon:
                    # print(operon, pred_operon, strain, operon_len, pred_gene_count, coverage)
                    TP_list.append(operon)
                    TP_operon_counterpart.append(operon) # not really necessary, but just keeping with format
                    Total_TPs.append(operon)
                    output_TP = "TP" + "\t" + operon + "\t" + pred_operon + "\t" + str(operon_len) + "\t" + str(pred_gene_count) + "\t" + coverage + "\t" + strain + "\n"
                    outfile.write (output_TP)

    # ################################################################################################
    #         # if the opredicted operon is longer than the validated operon, print FPs
    #         # the entire operon must be contained in it through
    #         # the predicted operon must contain all the genes of the real operon
    #         # AND have extra genes up/dowstream of the real operon
    # #################################################################################################
                elif elem [0] and not coverage == "NOT EXPRESSED" and not "EBG" in elem[0] and int(pred_gene_count) > 1:

                    pred_operon_prel = elem[0]
                    pred_operon_start_prel = pred_operon_prel.split(" - ")[0]
                    pred_operon_end_prel = pred_operon_prel.split(" - ")[1]
                    if pred_operon_start_prel[0].isalpha or pred_operon_start_prel[-1].isalpha:
                        pred_operon_start = int(pred_operon_start_prel.strip(string.ascii_letters))
                        str_pred_operon_start  = pred_operon_start_prel.strip(string.ascii_letters)
                    if pred_operon_end_prel[0].isalpha or pred_operon_end_prel[-1].isalpha:
                        pred_operon_end = int(pred_operon_end_prel.strip(string.ascii_letters))
                        str_pred_operon_end = pred_operon_end_prel.strip(string.ascii_letters)
                        pred_operon_len = ((pred_operon_end - pred_operon_start) +1)

                        if (pred_operon_start < operon_start and pred_operon_end > operon_end) or (pred_operon_start < operon_start and operon_end == pred_operon_end) or (pred_operon_start == operon_start and pred_operon_end > operon_end):
                            FP_list.append(pred_operon)
                            FP_operon_counterpart.append(operon)
                            Total_FPs.append(pred_operon)
                            Total_FP_counterpart.append(operon)
                            output_FP = "FP" + "\t" + operon + "\t" + pred_operon + "\t" + str(operon_len) + "\t" + str(pred_gene_count) + "\t" +  coverage + "\t" + strain + "\n"
                            outfile.write (output_FP)
                            # Total_FPs.append(operon)
    ## ##############################################################################################
    #         # if the predicted operon is shorter than it should be, print FNs
    # ##############################################################################################

                        elif (pred_operon_start > operon_start and pred_operon_end < operon_end) or (pred_operon_start > operon_start and operon_end == pred_operon_end) or (pred_operon_start == operon_start and pred_operon_end < operon_end):

                            FN_list.append(pred_operon)
                            Total_FNs.append(pred_operon)
                            FN_operon_counterpart.append(operon)
                            Total_FN_counterpart.append(operon)
                            output_FN = "FN" + "\t" + operon + "\t" + pred_operon + "\t" + str(operon_len) + "\t" + str(pred_gene_count) + "\t" +  coverage + "\t" + strain + "\n"
                            outfile.write (output_FN)


                        elif pred_operon_start_prel == operon.split(" - ")[0] or pred_operon_start_prel == operon.split(" - ")[1] or pred_operon_end_prel == operon.split(" - ")[0] or pred_operon_end_prel == operon.split(" - ")[1] or (operon.split(" - ")[0] < pred_operon_start_prel < operon.split(" - ")[1] and pred_operon_end_prel > operon.split(" - ")[1]) or ( operon.split(" - ")[0] < pred_operon_end_prel < operon.split(" - ")[1] and pred_operon_start_prel < operon.split(" - ")[0]):
                            if not pred_operon_start_prel in TP_list and not pred_operon_start_prel in FP_operon_counterpart and not pred_operon_start_prel in FN_operon_counterpart and not pred_operon_end_prel in TP_list and not pred_operon_end_prel in FP_operon_counterpart and not pred_operon_end_prel in FN_operon_counterpart :


                                FN_list.append(pred_operon)
                                Total_FNs.append(pred_operon)
                                FN_operon_counterpart.append(operon)
                                Total_FN_counterpart.append(operon)
                                output_FN_overlap = "FN" + "\t" + operon + "\t" + pred_operon + "\t" + str(operon_len) + "\t" + str(pred_gene_count) + "\t" +  coverage + "\t" + strain + "\n"
                                outfile.write (output_FN_overlap)


                elif elem [0] and not coverage == "NOT EXPRESSED" and not "EBG" in elem[0] and int(pred_gene_count) == 1:
                    pred_operon_single = elem[0]
                    if pred_operon_single[-1].isalpha or pred_operon_single[0].isalpha:
                        pred_operon_single_gene_int = int(pred_operon_single.strip(string.ascii_letters))
                    pred_operon_gene_len_single = "1"
                    if pred_operon_single in operon or (operon_start < pred_operon_single_gene_int < operon_end):
                        FN_single_gene_list.append(pred_operon_single)
                        FN_list.append(pred_operon_single)
                        Total_FNs.append(pred_operon_single)
                        FN_single_gene_counterpart.append(operon)
                        FN_operon_counterpart.append(operon)
                        Total_FN_counterpart.append(operon)
                        output_FN_single_gene = "FN_single_gene" + "\t" + operon + "\t" + pred_operon_single + "\t" + str(operon_len) + "\t" + pred_operon_gene_len_single + "\t" +  coverage + "\t" + strain + "\n"
                        outfile.write (output_FN_single_gene)


                elif elem [0] and not "EBG" in elem[0] and coverage == "NOT EXPRESSED":
                    pred_operon_single_unexp = elem[0]
                    pred_operon_gene_len_single_unexp = "0"
                    if pred_operon_single_unexp [-1].isalpha or pred_operon_single_unexp[0].isalpha:
                        pred_operon_single_gene_int_unexp  = int(pred_operon_single_unexp.strip(string.ascii_letters))
                    if pred_operon_single_unexp in operon or (operon_start < pred_operon_single_gene_int_unexp < operon_end):
                        FN_single_gene_list.append(pred_operon_single_unexp)
                        FN_list.append(pred_operon_single_unexp)
                        Total_FNs.append(pred_operon_single_unexp)
                        FN_single_gene_counterpart.append(operon)
                        FN_operon_counterpart.append(operon)
                        Total_FN_counterpart.append(operon)
                        output_FN_unexpressed = "FN_NOT_EXPRESSED" + "\t" + operon + "\t" + pred_operon_single_unexp + "\t" + str(operon_len) + "\t" + pred_operon_gene_len_single_unexp + "\t" +  coverage + "\t" + strain + "\n"
                        outfile.write (output_FN_unexpressed)
                    # outfile.write("\n\n")

######################################################################################################

    header_total_TPs = f"Predictions made by COSMO per isolate {strain} \n"
    strain_title = strain + "\n\n"
    outfile_counts.write(strain_title)

    outfile_counts.write(header_total_TPs)
    for TP in set(Total_TPs):
        out_TP = TP + "\t" + "TP" + "\t" + strain + "\n"
        outfile_counts.write(out_TP)
    outfile_counts.write("\n")


    for FP in set(Total_FPs):
        out_FP = FP + "\t" + "FP" + "\t" + strain + "\n"
        outfile_counts.write(out_FP)
    outfile_counts.write("\n")

    for FN in set(Total_FNs):
        out_FN = FN + "\t" + "FN" + "\t" + strain + "\n"
        outfile_counts.write(out_FN)
    outfile_counts.write("\n\n\n")
#
# ###############################################
    set_TPs = set(TP_list)

    set_FPs = set(FP_list)
    set_FP_counterpart = set(FP_operon_counterpart)

    set_FNs = set(FN_list)
    set_FN_counterpart = set(FN_operon_counterpart)

    TP_counts = str(len(set_TPs))  + "\t" + strain +  "\t" + "TP" + "\n"
    FP_counts = str(len(set_FP_counterpart)) + "\t" + strain + "\t" + "FP" + "\n"
    FN_counts = str(len(set_FN_counterpart)) + "\t" + strain + "\t" + "FN" + "\n"
    header_total_calls = "Summary of predictions made by COSMO per Isolate/strain\n"
    outfile_counts.write(header_total_calls)
    outfile_counts.write(TP_counts)
    outfile_counts.write(FP_counts)
    outfile_counts.write(FN_counts)
    end_line = "##################################################################################\n\n"
    outfile_counts.write("\n\n")
    outfile_counts.write(end_line)

print("The file:", outfile_name, "has been successfully printed")
print("The file:", outfile_name_counts, "has been successfully printed")

outfile.close()
outfile_counts.close()
