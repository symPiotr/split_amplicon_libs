#! /usr/bin/env python3

import os, sys, re

nucl_ambigs = {"R": "[AG]", "Y": "[CT]", "S": "[GC]", "W": "[AT]", "K": "[GT]", "M": "[AC]", "B": "[CGT]", "D": "[AGT]", "H": "[ACT]", "V": "[ACG]", "N": "[ACGT]"}


if len(sys.argv) != 4:
	sys.exit('\nsplit_amplicon_libs_VarLenIns.py v. 1.0. Piotr Lukasik, 23rd Dec 2020\n'
	         '-------------------------------------------------------------\n'
	         'This script splits up paired fastq files corresponding to a multiplex amplicon library\n'
	         'into a collection of fastq file pairs corresponding to different targeted regions,\n'
	         'using the sequences of primers for alternative targets provided. The script assumes\n'
	         'that the primer sequences may be preceded by up to four other bases (variable-length insert).\n'
	         '    Usage: ./split_amplicon_libs.py <list_of_targets> <list_of_libraries> <output_dir> \n'
	         '    e.g., ./split_amplicon_libs.py targets.txt samples.txt ~/test_data')
	         
Script, Target_List, Sample_List, Output_Dir = sys.argv

# Formatting the Output_Dir info
if Output_Dir.endswith("/"):
    Output_Dir = Output_Dir[:len(Output_Dir)-1]
if not os.path.exists(Output_Dir):
    os.makedirs(Output_Dir)


print('\nsplit_amplicon_libs.py v. 1.0. Piotr Lukasik, 23rd Dec 2020\n'
	   '-------------------------------------------------------------\n')


###################
### Block 1. Reading and converting amplicon target list ....
### The expected format of each line:
### COI-BF3BR2      CCHGAYATRGCHTTYCCHCG    TCDGGRTGNCCRAARAAYCA
###################

print("Reading the list of targets .....", end="")

Target_List_Imported = []
TARGET_LIST = open(Target_List, 'r', encoding="utf-8")
for line in TARGET_LIST:
    if not line.startswith("#") and len(line) > 0:
        if len(line.strip().split()) != 3:
            print("There is a problem with this line, which will be skipped:\n    ", end = "")
            print(line.strip())
        else:
            Target_List_Imported.append(line.strip().split())
TARGET_LIST.close()
# Target_List_Imported = [[target_name, primer_F, primer_R], [target_name2, primer_F, primer_R]...]
# contains lists of targeted amplicon regions: target name, forward primer, reverse primer

		
Target_Dict = {}
### Now, I am going to convert primers (with ambiguous bases) to search terms for re.match(), such as 'AG[AG]CT[CTG]A'
for target in Target_List_Imported:
    if len(target) == 0:
        sys.exit("Problem with the imported Target List :|")

    primer_R1_ambig = ".{0,4}" # assumes that there can be up to four of whatever prior to the sequence!
    for base in target[1][1:]:     # skips the first base in the primer!
        if base in nucl_ambigs.keys():
            primer_R1_ambig += nucl_ambigs[base] 
        else: 
            primer_R1_ambig += base 

    primer_R2_ambig = ".{0,4}"   # assumes that there can be up to four of whatever prior to the sequence!
    for base in target[2][1:]:       # skips the first base in the primer! Helpful for COI!
        if base in nucl_ambigs.keys():
            primer_R2_ambig += nucl_ambigs[base] 
        else: 
            primer_R2_ambig += base 
                
    Target_Dict[target[0]] = [primer_R1_ambig, primer_R2_ambig, 0]
    # target_name: [primer_F, primer_R, seq_count]

print("OK! %d targets provided:\n    " % len(Target_Dict), end = "")
for target in sorted(list(Target_Dict.keys())):
    print(target, end = ", ")
print("")



### Creating a bunch of output folders --- one for each target
for target in Target_Dict.keys():
    if not os.path.exists('%s/%s' % (Output_Dir, target)):
        os.makedirs('%s/%s' % (Output_Dir, target))

### Creating the framework for storing read numbers for each target and sample
Target_Counts_per_Sample = [['']]
for target in sorted(list(Target_Dict.keys())):
    Target_Counts_per_Sample[0].append(target) # The top row will contain target labels
Target_Counts_per_Sample[0].append("Unclassified") # With "Unclassified" at the end
# The subsequent lines will contain library name and counts



###################
### Block 2. Reading sample list ....
###################

print("\nReading the list of samples .....", end ="")

Sample_List_Imported = []
SAMPLE_LIST = open(Sample_List, 'r', encoding="utf-8")
for line in SAMPLE_LIST:
    if not line.startswith("#") and len(line) > 0:
        if len(line.strip().split()) != 3:
            print("There is a problem with this line, which will be skipped:\n    ", end = "")
            print(line.strip())
        else:
            Sample_List_Imported.append(line.strip().split())
SAMPLE_LIST.close()
# Sample_List_Imported = [[libr1, libr1_R1.fastq, libr1_R2.fastq], ...]
# contains lists of targeted amplicon regions: target name, forward primer, reverse primer

Sample_count = len(Sample_List_Imported)

print("OK! %d samples will be sorted" % Sample_count)


###################
### Block 3. Processing reads
###################


Sample_index = 0

for library in Sample_List_Imported:
    Sample_index += 1
    print("Sorting data for sample %d / %d:   %s ................." % (Sample_index, Sample_count, library[0]), end = "")
    
    Unclassified_ct = 0
    
    ### Now, reading in R1 and R2 fastq files. 
    READ1 = open(library[1], "r")
    READ2 = open(library[2], "r")
    
    ### Here, we will store data for reads classified to different targets. Skipping Unclassified!
    reads_by_target = {}
    for target in Target_Dict.keys():
        reads_by_target[target] = []
    
    Unclassified_ct = 0   
    
    
    ### Now, reading input fastqs: read pair by read pair, saving details to reads_by_target dictionary
    for line in READ1:
        R1_Heading = line.strip()
        R1_Seq = READ1.readline().strip()
        READ1.readline()
        R1_Qual = READ1.readline().strip()

        R2_Heading = READ2.readline().strip()
        R2_Seq = READ2.readline().strip()
        READ2.readline()
        R2_Qual = READ2.readline().strip()
        
        Unclassified = 1
        for target in Target_Dict.keys():
            #print(target, "\n", R1_Seq[:30], Target_Dict[target][0], re.match(Target_Dict[target][0], R1_Seq))
            #print(R2_Seq[:30], Target_Dict[target][1], re.match(Target_Dict[target][1], R2_Seq))
            if bool(re.match(Target_Dict[target][0], R1_Seq)) and bool(re.match(Target_Dict[target][1], R2_Seq)): # if the primer seqs match the beginnings of target seqs
                reads_by_target[target].append([R1_Heading, R1_Seq, R1_Qual, R2_Heading, R2_Seq, R2_Qual])
                #print("HIT!")
                Unclassified = 0
        
        Unclassified_ct += Unclassified


    ### Now, exporting contents of reads_by_target dictionary
    for target in Target_Dict.keys():
        target_R1 = open("%s/%s/%s_%s_R1.fastq" % (Output_Dir, target, target, library[0]), "w")
        for seq in reads_by_target[target]:
            print(seq[0],seq[1],"+",seq[2], sep = "\n", file = target_R1)
        target_R1.close()
        
        target_R2 = open("%s/%s/%s_%s_R2.fastq" % (Output_Dir, target, target, library[0]), "w")
        for seq in reads_by_target[target]:
            print(seq[3],seq[4],"+",seq[5], sep = "\n", file = target_R2)
        target_R2.close()


    ### Preparing the line for Target_Counts_per_Sample, adding it
    Read_total = Unclassified_ct
    
    Target_Counts_per_Sample_line = [library[0]]
    for target in sorted(list(Target_Dict.keys())):
        Target_Counts_per_Sample_line.append(len(reads_by_target[target]))
        Read_total += len(reads_by_target[target])
    Target_Counts_per_Sample_line.append(Unclassified_ct)
    
    Target_Counts_per_Sample.append(Target_Counts_per_Sample_line)


    ### Looks like it's all OK for this sample!
    print("OK! %d paired-end sequences classified!" % Read_total)



###################
### Block 4. Almost done, just printing summaries :)
###################

Summary_file = open(Output_Dir + "/000_splitting_summary.txt", "w")

print("\nJob complete!\n")
print("Target counts in different samples:")
for line in Target_Counts_per_Sample:
    for item in line:
        print(item, end = "\t")
        print(item, end = "\t", file = Summary_file)
    print("")
    print("", file = Summary_file)

print("DONE!")            
