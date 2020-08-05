# split_amplicon_libs
Script split_amplicon_libs.py separates Illumina amplicon datasets (paired FASTQ files) into files corresponding to different amplification targets

### Requirements: 
Unix environment, Python 3, Biopython [details of versions it was tested with...]

### Usage:
split_amplicon_libs_VarLenIns.py <list_of_targets> <list_of_libraries> <output_dir>
e.g., 
split_amplicon_libs_VarLenIns.py standard_primers.txt sample_list.txt ~/test_data

Where sample_list.txt contains library IDs and file names ---
```
piotr.lukasik@fsm:~/data_20200730/libs$ head sample_list.txt
A_ACACON	A-ACACON_S234_L001_R1_001.fastq	A-ACACON_S234_L001_R2_001.fastq
A_CALBON1	A-CALBON1_S244_L001_R1_001.fastq	A-CALBON1_S244_L001_R2_001.fastq
A_CALBON2	A-CALBON2_S245_L001_R1_001.fastq	A-CALBON2_S245_L001_R2_001.fastqhead sample_list.txt
...
```

Where target_list.txt contains target IDs and primer sequences ---
```
piotr.lukasik@fsm:~/data_20200730/libs$ cat ~/scripts/standard_primers.txt 
### Different microbial marker genes - our standard combo
#region_name	#R1_primer_sequence	#R2_primer_sequence
COI-BF3BR2      CCHGAYATRGCHTTYCCHCG    TCDGGRTGNCCRAARAAYCA
16S-v4	GTGYCAGCMGCCGCGGTAA	GGACTACNVGGGTWTCTAAT
16S-V1V2  AGMGTTYGATYMTGGCTCAG  TGCTGCCTCCCGTAGGAGT
16S-V3V4  CCTACGGGNGGCWGCAG  GACTACHVGGGTATCTAATCC
ITS1	TGGTCATTTAGAGGAAGTAA	GCTGCGTTCTTCATCGATGC
ITS5_8  TGGTCATTTAGAGGAAGTAA  GCGTTCTTCATCGAT
ITS2  GTGARTCATCGAATCTTTG  CCTCCGCTTATTGATATGC
```

...and output_dir is wherever you want your post-splitting files created.
