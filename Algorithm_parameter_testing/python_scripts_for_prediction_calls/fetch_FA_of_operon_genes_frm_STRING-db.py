## script name: fetch_FA_of_operon_genes_frm_STRING-db.py

# edited from: https://string-db.org/cgi/help.pl?subpage=api%23mapping-identifiers


#!/usr/bin/env python3

##########################################################
## For a given list of proteins the script
## finds its functional annotation
## and prints out the mapping on screen in the TSV format

###########################################################

import requests ## python -m pip install requests
import string

string_api_url = "https://string-db.org/api"
output_format = "tsv-no-header"
method = "get_string_ids"

##
## Set parameters
##

params = {

    "identifiers" : "\r".join(['Rv0046c', 'Rv0047c', 'Rv0096', 'Rv0097', 'Rv0098', 'Rv0099', 'Rv0100', 'Rv0101', 'Rv0102', 'Rv0166', 'Rv0167', 'Rv0168', 'Rv0169', 'Rv0170', 'Rv0171', 'Rv0172', 'Rv0173', 'Rv0174', 'Rv0175', 'Rv0176', 'Rv0177', 'Rv0178', 'Rv0287', 'Rv0288', 'Rv0490', 'Rv0491', 'Rv0586', 'Rv0587', 'Rv0588', 'Rv0589', 'Rv0590', 'Rv0590A', 'Rv0591', 'Rv0592', 'Rv0593', 'Rv0594', 'Rv0676c', 'Rv0677c', 'Rv0735', 'Rv0736', 'Rv0902c', 'Rv0903c', 'Rv0928', 'Rv0929', 'Rv0930', 'Rv0933', 'Rv0934', 'Rv0935', 'Rv0936', 'Rv0967', 'Rv0968', 'Rv0969', 'Rv0970', 'Rv0986', 'Rv0987', 'Rv0988', 'Rv1138c', 'Rv1139c', 'Rv1161', 'Rv1162', 'Rv1163', 'Rv1164', 'Rv1285', 'Rv1286', 'Rv1303', 'Rv1304', 'Rv1305', 'Rv1306', 'Rv1307', 'Rv1308', 'Rv1309', 'Rv1310', 'Rv1311', 'Rv1312', 'Rv1334', 'Rv1335', 'Rv1336', 'Rv1410c', 'Rv1411c', 'Rv1460', 'Rv1461', 'Rv1462', 'Rv1463', 'Rv1464', 'Rv1465', 'Rv1466', 'Rv1477', 'Rv1478', 'Rv1483', 'Rv1484', 'Rv1660', 'Rv1661', 'Rv1806', 'Rv1807', 'Rv1808', 'Rv1809', 'Rv1908c', 'Rv1909c', 'Rv1964', 'Rv1965', 'Rv1966', 'Rv1966', 'Rv1967', 'Rv1968', 'Rv1969', 'Rv1970', 'Rv1971', 'Rv2243', 'Rv2244', 'Rv2245', 'Rv2246', 'Rv2247', 'Rv2358', 'Rv2359', 'Rv2430c', 'Rv2431c', 'Rv2481c', 'Rv2482c', 'Rv2483c', 'Rv2484c', 'Rv2592c', 'Rv2593c', 'Rv2594c', 'Rv2686c', 'Rv2687c', 'Rv2688c', 'Rv2743c', 'Rv2744c', 'Rv2745c', 'Rv2871', 'Rv2872', 'Rv2873', 'Rv2874', 'Rv2875', 'Rv2877c', 'Rv2878c', 'Rv2930', 'Rv2931', 'Rv2932', 'Rv2933', 'Rv2934', 'Rv2935', 'Rv2936', 'Rv2937', 'Rv2938', 'Rv2958c', 'Rv2959c', 'Rv3083', 'Rv3084', 'Rv3085', 'Rv3086', 'Rv3087', 'Rv3088', 'Rv3089', 'Rv3132c', 'Rv3133c', 'Rv3134c', 'Rv3145', 'Rv3146', 'Rv3147', 'Rv3148', 'Rv3149', 'Rv3150', 'Rv3151', 'Rv3152', 'Rv3153', 'Rv3154', 'Rv3155', 'Rv3156', 'Rv3157', 'Rv3158','Rv3417c', 'Rv3418c', 'Rv3419c', 'Rv3420c', 'Rv3421c', 'Rv3422c', 'Rv3423c', 'Rv3493c', 'Rv3494c', 'Rv3495c', 'Rv3496c', 'Rv3497c', 'Rv3498c', 'Rv3499c', 'Rv3500c', 'Rv3501c', 'Rv3516', 'Rv3517', 'Rv3612c', 'Rv3613c', 'Rv3614c', 'Rv3615c', 'Rv3616c', 'Rv3792', 'Rv3793', 'Rv3794', 'Rv3795', 'Rv3874', 'Rv3875', 'Rv3917c', 'Rv3918c', 'Rv3919c', 'Rv3921c', 'Rv3922c', 'Rv3923c', 'Rv3924c']), # your protein list
    "species" : 83332, # species NCBI identifier
    "limit" : 1, # only one (best) identifier per input protein
    "echo_query" : 1, # see your input identifiers in the output
    "caller_identity" : "Tracey" # your app name

}

##
## Construct URL
##

request_url = "/".join([string_api_url, output_format, method])

##
## Call STRING
##

results = requests.post(request_url, data=params) # built in function of STRING-db

##
## Read and parse the results
##
outfile_name = "functions_operon_genes.txt"
outfile = open(outfile_name, "w")
heading = "Genes of 50 Operons and their annotations from STRING-db\n\n"
outfile.write(heading)

for line in results.text.strip().split("\n"):
    l = line.split("\t")
    # print(l)
    input_identifier, gene_name, FA = l[0], l[5], l[6]

    if FA.startswith("Rv"): # sometimes FA (functional annotation) was not in l[6]
        annotation = "\t".join([FA.split("aa.")[1].split(",")[0], FA.split("aa.")[1].split(",")[1]])

    elif "len:" in FA  and ";" in FA and not "kDa protein" in FA:
        description = FA.split("aa.")[0]
        annotation = description.split(",")[0].split(";")[0]


    elif "len:" in FA  and ";" in FA and "kDa protein" in FA:
        annotaion = FA.split(",")[3]

    elif not "len:" in FA and ";" in FA:
        protein_family = FA.split(";")[0]
        annotation = FA.split(";")[0] #"\t".join([FA.split(";")[0]), FA.strip().split(";")[0].split(".")[0]])


    outline = f"{input_identifier}    {annotation}\n" #{gene_name}
    outfile.write(outline)
