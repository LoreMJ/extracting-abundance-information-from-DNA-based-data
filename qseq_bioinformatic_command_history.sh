#!/bin/bash
set -e
set -u
set -o pipefail
##########################################################################################################
##########################################################################################################
# shell script for qseq data
##########################################################################################################
##########################################################################################################

# Usage: bash qseq_bioinformatic_command_history.sh

# PIPESTART=$(date)

# run script from ~/PATH_TO/qseq/scripts
# fq files in qseq/data/
# info files about primers, tags and PCR_sets in qseq/info/
# analysis outputs in qseq/analysis/

# set variables
HOMEFOLDER="/Users/luomingjie/qseq/" # do not have a ~ in the path
Begum=/Users/luomingjie/src/Begum/src/
DAME=/Users/luomingjie/src/DAMe/bin/
PAIRFQ=/Users/luomingjie/src/Pairfq/bin/

cd "${HOMEFOLDER}"

if [ ! -d analysis ] # if directory analysis does not exist
then
	mkdir analysis
fi

if [ ! -d PandaSeq_outputs ] # if directory PandaSeq_outputs does not exist
then
	mkdir PandaSeq_outputs
fi
mv PandaSeq_outputs analysis/

if [ ! -d Sort_outputs ] # if directory Sort_outputs does not exist
then
	mkdir Sort_outputs
fi
mv Sort_outputs analysis/

# read in library folder list and make a bash array
find data/* -maxdepth 0 -type d | sed 's/^.*\///g' > librarylist.txt # find all libraries in the data folder
library_info=librarylist.txt # put librarylist.txt into variable
library_names=($(cut -f 1 "${library_info}" | uniq)) # convert variable to array this way
echo "there are" "${#library_names[@]}" "libraries will be processed."
 # echo number of elements in the array
INDEX=1

for library in ${library_names[@]} # ${library_names[@]} is the full bash array
do
  echo "Now on library" ${INDEX} of ${#library_names[@]}
  INDEX=$((INDEX+1))

  mkdir ${library}_working
  cp data/${library}/*.fq.gz ${library}_working
  cd ${library}_working
  
  # Use AdapterRemoval to trim Illumina sequencing adapters
  # brew install adapterremoval
  # or download from https://github.com/MikkelSchubert/adapterremoval
  AdapterRemoval --file1 MJ${library}_*_1.fq.gz --file2 MJ${library}_*_2.fq.gz --output1 ${library}_R1_ADtrimed.fq.gz --output2 ${library}_R2_ADtrimed.fq.gz --gzip --trimns
  rm MJ${library}_*.fq.gz

  # Sickle error correction strategies and identified quality trimming
  # brew install sickle
  # or download from https://github.com/najoshi/sickle
  sickle pe -f ${library}_R1_ADtrimed.fq.gz -r ${library}_R2_ADtrimed.fq.gz -t sanger -o ${library}_R1_ADStrimed.fastq.gz -p ${library}_R2_ADStrimed.fastq.gz -s ${library}_trimmed_singles_file.fastq.gz -g >> sickle.out
  rm -f *_ADtrimed.fq.gz

  # run bfc to do error correction
  # brew install bfc
  # or download from https://github.com/lh3/bfc
  bfc -s 3g -t16 ${library}_R1_ADStrimed.fastq.gz | gzip -9 > ${library}_R1_ADStrimed_bfc.fastq.gz
  bfc -s 3g -t16 ${library}_R2_ADStrimed.fastq.gz | gzip -9 > ${library}_R2_ADStrimed_bfc.fastq.gz
  rm -f *_ADStrimed.fastq.gz

  # use pandaseq to merge pair-end reads
  # brew install homebrew/science/pandaseq
  # or download from https://github.com/neufeld/pandaseq/releases
  	# because bfc_outputs can't be used directly by pandaseq, pairfq has to be used to pair reads before running pandaseq
  	# download pairfq from https://github.com/sestaton/Pairfq
  gunzip -d ${library}_R*_ADStrimed_bfc.fastq.gz
  cat ${library}_R1_ADStrimed_bfc.fastq | tr "\t" " " > temp.fastq; rm -f ${library}_R1_ADStrimed_bfc.fastq; sed 's/ .*/ 1:N:0:0/g' temp.fastq > ${library}_R1_ADStrimed_bfc.fastq
  cat ${library}_R2_ADStrimed_bfc.fastq | tr "\t" " " > temp.fastq; rm -f ${library}_R2_ADStrimed_bfc.fastq; sed 's/ .*/ 2:N:0:0/g' temp.fastq > ${library}_R2_ADStrimed_bfc.fastq
  rm -f temp.fastq

  
  pairfq makepairs -f ${library}_R1_ADStrimed_bfc.fastq -r ${library}_R2_ADStrimed_bfc.fastq -fp R1_paired.fastq -rp R2_paired.fastq -fs r1_sing.fastq -rs r2_sing.fastq -stats

  rm -f r*.fastq #delete all the unpaired reads data files

  #run pandaseq
  pandaseq -f R1_paired.fastq -r R2_paired.fastq -A simple_bayesian -B -F -d bfsrk -N -g ${library}_pandaseq_log.txt -w ${library}_ADStrimed_bfc_merged.fastq

  rm -f R*_paired.fastq
  gzip -9 ${library}_ADStrimed_bfc_merged.fastq
  
  cd ${HOMEFOLDER}
  mv ${library}_working/${library}_ADStrimed_bfc_merged.fastq.gz analysis/PandaSeq_outputs
  rm -rf ${library}_working
  
  cd analysis/PandaSeq_outputs/
  seqkit fq2fa ${library}_ADStrimed_bfc_merged.fastq.gz -w 0 -j 4 -o ${library}_merged.fasta # this is for the latter steps in generating OTU-tables of mock soup
  
  ######################################### Begum sort #################################
  # prepare primer.txt, tags.txt and ${library}PCRset_Info.txt, Pool${library}_Info.txt; put them in folder $HOMEFOLDER/info/
  # beginning of sort
  # Download Begum from https://github.com/shyamsg/Begum. Begum needs python 2.7.16
  cd "${HOMEFOLDER}"
  mkdir pool${library}
  cp analysis/PandaSeq_outputs/${library}_ADStrimed_bfc_merged.fastq.gz pool${library}
  cd pool${library}
  #sorting is quite slow. better run on HPC. It took more than 24 hours on my imac for each library.
  python2 ${Begum}Begum.py sort -p "${HOMEFOLDER}"/info/primers.txt -t "${HOMEFOLDER}"/info/tags.txt -s "${HOMEFOLDER}"/info/${library}PCRset_Info.txt -l "${HOMEFOLDER}"/info/pool${library}_Info.txt -pm 2 -d . -o 'sorted'
  mv sorted_${library}.tagInfo "${HOMEFOLDER}"/analysis/Sort_outputs
  cd "${HOMEFOLDER}"
  rm -rf pool${library}
  
done
 
##################################### Begum filter#####################################
#prepare PCRset_Info.txt
cd "${HOMEFOLDER}"/analysis/
mkdir Filter_outputs
cd Sort_outputs
 
# python2 ${Begum}Begum.py filter -h
# usage: Begum filter [-h] -i InputPrefix -s SampleInformationFile [-p propPCRs]
#                     [-m minTimes] [-l minLength] [-d OutDirectory]
#                     [-o OutPrefix]
# 
# optional arguments:
#   -h, --help            show this help message and exit
#   -i InputPrefix, --inputPrefix InputPrefix
#                         Information file with prefix information for the sort
#                         tagInfo files.
#   -s SampleInformationFile, --sampleInfo SampleInformationFile
#                         File with tag combo and pool for each sample (Format:
#                         Sample FwdTagName RevTagName PoolName)
#   -p propPCRs, --propPCRs propPCRs
#                         Minimum proportion of PCR replicates a sequence should
#                         be present in.
#   -m minTimes, --minOccurence minTimes
#                         Minimum number of times a sequence should be present,
#                         in a PCR replicate to be consider a true sequence.
#   -l minLength, --minLength minLength
#                         Minimum length of the amplicon sequence - in case of
#                         single end or merged sequences, it is the length of
#                         the sequence, and in case of paired end reads, it is
#                         the sum of the length of the 2 reads.
#   -d OutDirectory, --output_directory OutDirectory
#                         Output directory
#   -o OutPrefix, --output_prefix OutPrefix
#                         Prefix for output files

# we only accept the reads that appeared in at least 2 of 3 PCRs (-p 0.6) at a min-copy number of 4 reads per PCR (-m 4), with min-length of 300 bp (-l 300)
python2 ${Begum}Begum.py filter -i 'sorted' -s "${HOMEFOLDER}"/info/PCRset_Info.txt -p 0.6 -m 4 -l 300 -d "${HOMEFOLDER}"/analysis/Filter_outputs -o BegumFilter_0.6_copies4
 
#################################### mock soup OTU table ###################################
# for mock soup data, we need to compare the results with and without UMI info.
# But begum pipeline can't extract the information of UMI, so we have to write our own scripts to deal with mock soup data
# we will generate two OTU tables for mock soup data, seperately with UMIs' numbers and reads' numbers

# this step is for transfering UMIs from the reads to their IDs
cd "${HOMEFOLDER}"
bash UMI_to_header.sh

# to obtain the UMI's info, we need to do picking OTUs in each PCRset's subdata
# pick OTUs by using QIIME's pick_otus.py
# macqiime needs to be installed following the instruction. http://www.wernerlab.org/software/macqiime/macqiime-installation
macqiime #open macqiime
bash pick_otus.sh
exit #back to normal shell
 
# generate the otu-table with UMIs' numbers
bash for_OTUtable_withUMI.sh
 
# generate the otu-table with reads' numbers
bash for_OTUtable_withREAD.sh
 
############################################### Malaise-trap samples OTU table ################################
cd "${HOMEFOLDER}"
# ${DAME}convertToUSearch.py -h
# usage: convertToUSearch.py [-h] -i InputFasta [-lmin minLength]
#                            [-lmax maxLength] [-u] [-s]
# 
# Take fasta file and fix it for input to USEARCH/Sumaclust
# 
# optional arguments:
#   -h, --help            show this help message and exit
#   -i InputFasta, --inFasta InputFasta
#                         Input fasta file.
#   -lmin minLength, --minLength minLength
#                         Minimum length of read
#   -lmax maxLength, --maxLength maxLength
#                         Maximum length of read
#   -u, --usearch         Generate Usearch input (default: sumaclust input)
#   -s, --sampleFastas    Make sample level fastas in SampleFastas directory

cd analysis
mkdir MalaiseTrap_outputs
cd MalaiseTrap_outputs

${DAME}convertToUSearch.py -i "${HOMEFOLDER}"/analysis/Filter_outputs/BegumFilter_0.6_copies4.fna -lmin 300 -lmax 330
gsed 's/ count/;size/' FilteredReads.forsumaclust.fna > FilteredReads.forvsearch.fna
vsearch --sortbysize FilteredReads.forvsearch.fna --output FilteredReads.forvsearch_sorted.fna
vsearch --uchime_denovo FilteredReads.forvsearch_sorted.fna --nonchimeras FilteredReads.forvsearch_sorted_nochimeras.fna #this took more than 2 hours on my imac
gsed 's/;size/ count/' FilteredReads.forvsearch_sorted_nochimeras.fna > Filtered.forsumaclust.nochimeras.fna
sumaclust -t 0.97 -e Filtered.forsumaclust.nochimeras.fna > OTUs_97_sumaclust.fna
${DAME}tabulateSumaclust.py -i OTUs_97_sumaclust.fna -o sumaclust_97_otutalbe.txt -blast #OTU table

################################ make matchlist for lulu #######################################################
# brew install blast
makeblastdb -in sumaclust_97_otutalbe.txt.blast.txt -parse_seqids -dbtype nucl
blastn -db sumaclust_97_otutalbe.txt.blast.txt -num_threads 4 -outfmt '6 qseqid sseqid pident' -out sumaclust_matchlist.txt -qcov_hsp_perc .90 -perc_identity .84 -query sumaclust_97_otutalbe.txt.blast.txt


######################################### check OTU representative seqs ##########################################
#after lulu
cd ~/analysis/MalaiseTrap_outputs

seqkit grep sumaclust_97_otutalbe.txt.blast.txt -f dilution_lulu_otuID.txt -w 0 > dilution_lulu_otu_seqs.fasta #dilution_lulu_otuID.txt is generated from R 
seqkit translate dilution_lulu_otu_seqs.fasta -T 5 -f 2 -w 0 -o dilution_lulu_otu_translated.fasta
seqkit grep dilution_lulu_otu_translated.fasta -s -p "*" -i -w 0 -o lulu_stopcondons_frame2.txt #check stopcodons

######################################## find spike-in ############################################
vsearch --search_exact sumaclust_97_otutalbe.txt.blast.txt --db spike-in.fasta --blast6out spike_matches.txt --strand plus 
