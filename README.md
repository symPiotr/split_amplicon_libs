# split_amplicon_libs.py
Script [split_amplicon_libs.py](split_amplicon_libs.py) splits Illumina amplicon datasets (paired FASTQ files, **un-gzipped**) into files corresponding to different amplification targets, checking for primer sequences within the initial portions of R1 and R2 reads.   
A version of the script, [split_gzipped_amplicon_libs.py](split_gzipped_amplicon_libs.py) does pretty much the same, except that it inputs and outputs gzipped files.  
  
Check [https://github.com/symPiotr/amplicon_analysis_pipeline](https://github.com/symPiotr/amplicon_analysis_pipeline) for info on the organization of multi-target libraries prepared according to the Symbiosis Evolution Group protocols, that will often be the main source of data for this script!
&nbsp;    
  
The script reads in the list of amplification targets and the list of sample and fastq file names. In each of the samples, in each pairs of reads, the script then searchers for primer sequences corresponding to each target within the initial portions of the R1 and R2 reads. It assign the read pair to a given target in the case of a match. It allows the primer sequences to be preceded by up to four additional bases (corresponding to our variable length inserts) but it does not currently allow for any mismatches within the primer sequences.  
&nbsp;  
  
The script then outputs reads corresponding to each of the targets into files "TARGET/TARGET_SAMPLE_R1.fastq" and "TARGET/TARGET_SAMPLE_R2.fastq" within the working directory. Finally, it prints the splitting summary (for each sample, the number of reads assigned to each of the targets, as well as those Unclassified) to a file **000_splitting_summary.txt** within the working directory.  
&nbsp  
  
### Requirements: 
Unix environment, Python 3

### Usage:
```
split_amplicon_libs.py <list_of_targets> <list_of_libraries> <output_dir>
e.g., 
split_amplicon_libs.py /mnt/matrix/symbio/db/references/standard_primers.txt sample_list.txt ~/test_data
```
OR
```
split_gzipped_amplicon_libs.py <list_of_targets> <list_of_libraries> <output_dir>
e.g., 
split_gzipped_amplicon_libs.py ~/references/standard_primers.txt sample_list.txt ~/test_data
```
  
The required file **target_list.txt** contains target IDs and primer sequences. Any non-commented-out lines in the script should go as follows:   
**TARGET_NAME <tab> FORWARD_PRIMER_SEQ <tab> REVERSE_PRIMER_SEQ**   
Ambiguous bases within the primer sequences are allowed: they will be converted into search terms following IUPAC codes.  
Note that the script allows the primer sequences within reads to be preceded by up to six bases (any sequence), corresponding to variable length inserts used as a part of the Symbiosis Evolution Group workflow. However, the script does not currently allow for any mismatches between the provided primers and the reads.  
```
piotr.lukasik@azor:~$ cat /mnt/qnap/users/symbio/references/primer_list_standard.txt
#### Standard primers
COI_BF3HYMBR2      CCHGAYATRKCHTWYCCHCG    TCDGGRTGNCCRAARAAYCA
16S_V4  GTGYCAGCMGCCGCGGTAA     GGACTACNVGGGTWTCTAAT
16S_V1V2  AGMGTTYGATYMTGGCTCAG  TGCTGCCTCCCGTAGGAGT
ITS2  GTGARTCATCGAATCTTTG  CCTCCGCTTATTGATATGC
ITS1a  TGGTCATTTAGAGGAAGTAA  GCGTTCTTCATCGAT
nu_SSU	CGATAACGAACGAGACCT	ANCCATTCAATCGGTANT
```  
&nbsp;  
The required file **sample_list.txt** contains library IDs and file names, where any non-commented-out lines should be organized as follows:  
**SAMPLE_NAME <tab>  R1_FILE_NAME.fastq <tab> R1_FILE_NAME.fastq**  
Note that the files need to be uncompressed --- fastq.gz would not do!  
```
(base) piotr.lukasik@azor:~/split_gunzipped_tests$ head sample_list.txt
IPA0680 IPA0680_R1.fq.gz IPA0680_R2.fq.gz
IPA0681 IPA0681_R1.fq.gz IPA0681_R2.fq.gz
IPA0682 IPA0682_R1.fq.gz IPA0682_R2.fq.gz
...
```  
How to create such sample list? You can use some combination of Excel and REGEX, or whatever you please, of course. But if your samples have the names as above, a convenient way could be a Unix pipe:  
```
### assuming that your input files are in the directory where you execute the script, and have name format SampleName_R1.fq.gz, SampleName_R2.fq.gz:
for file in *_R1.fq.gz; do
    SampleName=`basename $file _R1.fq.gz `
    SampleNameMod=$(echo "$SampleName" | sed 's/-/_/g' | sed 's/_S[0-9]\+$//g')
    echo $SampleNameMod "$SampleName"_R1.fq.gz "$SampleName"_R2.fq.gz >> sample_list.txt
done
```  
  
&nbsp;  

The output_dir is wherever you want to move your post-splitting files!
