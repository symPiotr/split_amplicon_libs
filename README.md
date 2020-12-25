# split_amplicon_libs.py
Script split_amplicon_libs.py splits Illumina amplicon datasets (paired FASTQ files, un-gzipped) into files corresponding to different amplification targets, checking for primer sequences within the initial portions of R1 and R2 reads. Check [https://github.com/symPiotr/amplicon_analysis_pipeline](https://github.com/symPiotr/amplicon_analysis_pipeline) for info on the organization of multi-target libraries prepared according to the Symbiosis Evolution Group protocols, that are likely to be the main source of data for this script!
&nbsp;  
  
**split_amplicon_libs.py** reads in the list of amplification targets, where 
It also reads the list of sample and fastq file names.
  
Within the provided working directory, for each of the samples, for each pairs of reads, the script searchers for primer sequences corresponding to each target within the initial portions of the R1 and R2 reads. It allows the primer sequences to be preceded by up to four additional bases (corresponding to our variable length inserts) but it does not currently allow for any mismatches within the primer sequences. The script then outputs reads corresponding to each of the targets into files "TARGET/TARGET_SAMPLE_R1.fastq" and "TARGET/TARGET_SAMPLE_R2.fastq" within the working directory. Finally, it prints the splitting summary (for each sample, the number of reads assigned to each of the targets, as well as those Unclassified) to a file "000_splitting_summary.txt" within the working directory.

### Requirements: 
Unix environment, Python 3

### Usage:
```
split_amplicon_libs.py <list_of_targets> <list_of_libraries> <output_dir>
e.g., 
split_amplicon_libs.py /mnt/matrix/symbio/db/references/standard_primers.txt sample_list.txt ~/test_data
```    
  
The required file **target_list.txt** contains target IDs and primer sequences. Any non-commented-out lines in the script should go as follows:   
**TARGET_NAME <tab> FORWARD_PRIMER_SEQ <tab> REVERSE_PRIMER_SEQ**   
Ambiguous bases within the primer sequences are allowed: they will be converted into search terms following IUPAC codes.  
Note that the script allows the primer sequences within reads to be preceded by up to four bases (any sequence), corresponding to variable length inserts used as a part of the Symbiosis Evolution Group workflow. However, the script does not currently allow for any mismatches between the provided primers and the reads.  
```
piotr.lukasik@fsm:~/scripts$ cat /mnt/matrix/symbio/db/references/standard_primers.txt
##region_name	#R1_primer_sequence	#R2_primer_sequence
COI-BF3BR2      CCHGAYATRGCHTTYCCHCG    TCDGGRTGNCCRAARAAYCA
16S-V4	GTGYCAGCMGCCGCGGTAA	GGACTACNVGGGTWTCTAAT
16S-V1V2  AGMGTTYGATYMTGGCTCAG  TGCTGCCTCCCGTAGGAGT
ITS2  GTGARTCATCGAATCTTTG  CCTCCGCTTATTGATATGC
ITS1a  TGGTCATTTAGAGGAAGTAA  GCGTTCTTCATCGAT
```  
&nbsp;  
The required file **sample_list.txt** contains library IDs and file names, where any non-commented-out lines should be organized as follows:  
**SAMPLE_NAME <tab>  R1_FILE_NAME.fastq <tab> R1_FILE_NAME.fastq**  
Note that the files need to be uncompressed --- fastq.gz would not do!  
```
piotr.lukasik@fsm:~/data_20200730/libs$ head sample_list.txt
A_ACACON	A-ACACON_S234_L001_R1_001.fastq	A-ACACON_S234_L001_R2_001.fastq
A_CALBON1	A-CALBON1_S244_L001_R1_001.fastq	A-CALBON1_S244_L001_R2_001.fastq
A_CALBON2	A-CALBON2_S245_L001_R1_001.fastq	A-CALBON2_S245_L001_R2_001.fastq
...
```  
How to create such sample list? You can use Excel or whatever you please, of course. But if your samples have the names as above, a convenient way could be a Unix pipe:  
```
for file in *_R1_001.fastq; do
    SampleName=`basename $file _L001_R1_001.fastq `
    SampleNameMod=$(echo "$SampleName" | sed 's/-/_/g' | sed 's/_S[0-9]\+$//g')
    echo $SampleNameMod "$SampleName"_L001_R1_001.fastq "$SampleName"_L001_R2_001.fastq >> sample_list.txt
done
```  
  
&nbsp;  

The output_dir is wherever you want to move your post-splitting files!
