#! /usr/bin/env python3

import os, argparse, sys

parser = argparse.ArgumentParser(description='----- most_abundant_OTU.py v. 1.0, M. Buczek, 12th Feb 2021 -----\n'
             'This script produce four files, from COI zOTU table:\n' 
             '- barcode.txt that contains info about most abundand COI barcode, taxonomy and bacteria presence\n' 
             '- barcode.fasta, containing most abundand Eucaryotic COI barcode per each library/sample\n'
             '- euc.txt, containing information about Eucaryotic COI sequences that represents at least 5% of total Eucaryotic reads per sample\n'
             '- bac.txt, containing information about bacterial COI sequences\n'
             '     e.g.,  most_abundant_OTU.py all_zotu_table_expanded.txt\n')

# positional arguments
parser.add_argument("path_to_input_file", metavar='<zOTU_table>', help="file containing zOTU table, see eg. in the GitHub repository")

# optional arguments
parser.add_argument("-lower", metavar='<lower_abundance_treshold>', type=float, help="minimum OTU abundance treshold(default: 0.05)", default=0.05)
parser.add_argument("-upper", metavar='<upper_abundance_treshold>', type=float, help="minimum OTU abundance treshold(default: 0.00)", default=0.00)
parser.add_argument("-reads", metavar='<reads_treshold>', type=int, help="minimum number o reads per library treshold(default: 20)", default=20)

# if no arguments were given, printing the help message (args = "--help")
args = parser.parse_args(args=None if sys.argv[1:] else ['--help'])

# print(args.lower)

INPUT = open(args.path_to_input_file, "r")
OUTPUT_BAC = open("bac_" + os.path.splitext(args.path_to_input_file)[0] + ".txt", mode= "w") # bac.txt
OUTPUT_EUC = open("euc_" + os.path.splitext(args.path_to_input_file)[0] + ".txt", mode= "w") # euc.txt
OUTPUT = open("barcode_" + os.path.splitext(args.path_to_input_file)[0] + ".txt", mode= "w") # barcode.txt
OUTPUT2 = open("barcode_" + os.path.splitext(args.path_to_input_file)[0] + ".fasta", mode= "w") # barcode.fasta

# read input file as list of lists
TABLE = []
for line in INPUT:
    TABLE.append(line.strip().split())    
INPUT.close()

ROW_NO = len(TABLE) 
COL_NO = len(TABLE[0])

### Divide table to Bacteria (BAC) and Eucaryote (EUC)
HEAD = TABLE[0]

EUC = [HEAD]
BAC = [HEAD]
BAC_YES = []
NO_BAC_index = []
for row_no in range(1,ROW_NO):
    if TABLE[row_no][2].startswith("Bacteria"): # if taxonomy starts with "Bacteria" save whole row in output_bac
        BAC.append(TABLE[row_no])
        for col_no in range(5, COL_NO): # this applies only to "Bacteria" rows
            if int(TABLE[row_no][col_no]) > 0: # if number of reads is larger than 0, add sample name to the list BAC_YES
                BAC_YES.append(HEAD[col_no])

    elif not 'Spikein' in TABLE[row_no][2]:
        EUC.append(TABLE[row_no]) # if taxonomy do not starts with "Bacteria" add the row to EUC table as a list

    BAC_YES = list(set(BAC_YES)) # list of samples that have at least one bacterial zotu, takes unique values from BAC_YES

for col_no in range(5, COL_NO):
    count = 0 # number of zOTUs with bacterial reads per sample
    for row_no in range(1, len(BAC)):
        if 0 < int(BAC[row_no][col_no]): # I do not know why it has to be int() here
            count += 1
    if count == 0:
        NO_BAC_index.append(col_no)

BAC_2 = [x[:] for x in BAC] # deep copy, believe me it was needed
for row_no in range(0,len(BAC)):
    for item in sorted(NO_BAC_index, reverse = True):  
        del BAC_2[row_no][item] # removes samples that do not contain bacteria
    for item in BAC_2[row_no][:-1]:
        print(item, file=OUTPUT_BAC, end="\t")
    print(BAC_2[row_no][-1], file=OUTPUT_BAC) # it is adding last argument in the row from BAC_2 and ending with \n

### Adding the TOTAL row to Eucaryotic table --- added-up counts of all reads for the same sample
TOTALS = TABLE[0][:4]
ROW_NO = len(EUC)
for col_no in range(4,COL_NO):
    Total = 0
    for row_no in range(1,ROW_NO):
        Total += int(EUC[row_no][col_no])
    TOTALS.append(Total)
EUC.append(TOTALS)

### Removing libraries with low number of reads
reads = args.reads
for row_no in range(0,len(EUC)):
    for col_no in reversed(range(4, COL_NO)):
        if EUC[-1][col_no] < reads:
            del EUC[row_no][col_no] # removes samples that have less reads than treshold

COL_NO = len(EUC[0]) # some columns may drop out from dataframe due to above code, so updating of the variable is needed

### Translating counts into relative abundances, with keeping both tables
EUC_ABUND = [x[:] for x in EUC] # deep copy
for col_no in range(4, COL_NO):
    for row_no in range(1, ROW_NO):    # Note: does NOT include the TOTAL row!
        if EUC[-1][col_no] > 0:
            EUC_ABUND[row_no][col_no] = round(float(EUC[row_no][col_no])/EUC[-1][col_no], 4)
        else:
            EUC_ABUND[row_no][col_no] = 0

### Creating new table with most abundand barcode per library (sample) and fasta file
NEW_HEAD = ["Sample", "zOTU", "OTU", "Sequence", "Abundance", "reads", "Abundance > 5%", "Bacteria", "Taxonomy"] # new col names
for item in NEW_HEAD[:-1]:
    print(item, file=OUTPUT, end="\t")
print(NEW_HEAD[-1], file=OUTPUT)

EUC_YES = []
EUC_YES_sample = []
most_abundand = []

for col_no in range(5, COL_NO):
    first = 0 # item to store highest abundance value
    first_index = 0
    treshold = args.lower # abundance treshold
    count = 0 # number of zOTUs with abundance above treshold value
    for row_no in range(1, ROW_NO):
        if first <= EUC_ABUND[row_no][col_no]:
            first = EUC_ABUND[row_no][col_no]
            first_index = row_no # stores row number of otu with highest abundance
        if treshold < EUC_ABUND[row_no][col_no] and 0.5 > EUC_ABUND[row_no][col_no]:
            count += 1
            EUC_YES.append(row_no) # collecting rows (zOTUz) that fulfill above requirements
            EUC_YES_sample.append(col_no) # collecting amplicon library (sample) name matching OTU one line above

    if HEAD[col_no] in BAC_YES:
        bacteria = "yes"
    else:
        bacteria = "no"

    most_abundand.append(first_index) # stores row number (OTUs) with highest abundance
    most_abundand = list(set(most_abundand)) # ?

    if first >= args.upper:
        to_add = [HEAD[col_no], EUC[first_index][0], EUC[first_index][1], EUC[first_index][3], first, EUC[first_index][col_no], count, bacteria, EUC[first_index][2]]
        to_fasta = to_add[:4]
        print(">", to_fasta[0], ",", to_fasta[1], ",", to_fasta[2], "\n", to_fasta[3], file=OUTPUT2, sep="" )
        for item in to_add[:-1]:
            print(item, file=OUTPUT, end="\t")
        print(to_add[-1], file=OUTPUT)

### Creating new table with potential Parasitoids and Contamination
NEW_HEAD_EUC = ["Sample", "zOTU", "OTU", "Sequence", "Abundance", "reads", "Parasitoid/Contamination", "Taxonomy"] # new col names
for item in NEW_HEAD_EUC[:-1]:
    print(item, file=OUTPUT_EUC, end="\t")
print(NEW_HEAD_EUC[-1], file=OUTPUT_EUC)

# get list of most abundand OTUs names (not row numbers)
most_abundand_otu = []
for item in most_abundand:
    most_abundand_otu.append(EUC_ABUND[item][1])

most_abundand_otu = sorted(list(set(most_abundand_otu))) # not needed

for row_no, name in zip(EUC_YES, EUC_YES_sample): # iteration on two lists at the same time
    if EUC_ABUND[row_no][1] in most_abundand_otu: # check if OTU from EUC_YES list is among most_abundand OTUs
        to_add = [HEAD[name], EUC[row_no][0], EUC[row_no][1], EUC[row_no][3], EUC_ABUND[row_no][name], EUC[row_no][name], "Contamination", EUC[row_no][2]]
    else:
        to_add = [HEAD[name], EUC[row_no][0], EUC[row_no][1], EUC[row_no][3], EUC_ABUND[row_no][name], EUC[row_no][name], "Parasitoid", EUC[row_no][2]]
    for item in to_add[:-1]:
        print(item, file=OUTPUT_EUC, end="\t")
    print(to_add[-1], file=OUTPUT_EUC)

INPUT.close()
OUTPUT_BAC.close() # bac.txt
OUTPUT_EUC.close() # euc.txt
OUTPUT.close() # barcode.txt
OUTPUT2.close() # barcode.fasta
