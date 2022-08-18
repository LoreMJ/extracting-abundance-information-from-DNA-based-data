#!/bin/bash
set -e
set -u
set -o pipefail
###########################################################################################
###########################################################################################
# a shell script to pick OTUs in QIIME
###########################################################################################
###########################################################################################

# Usage: bash pick_otus.sh
# run on my mac pro

# need to prepare the Sanger sequences of those individuals used here as ref (Sanger_ref.fas), put it in Folder info/

PIPESTART=$(date)

HOMEFOLDER=$(pwd)

echo "Home folder is "${HOMEFOLDER}""

# set variables
HOMEFOLDER="/Users/apple/src/qseq/" # do not have a ~ in the path
INDEX=1

cd "${HOMEFOLDER}"/analysis/

if [ ! -d mocksoup_outputs ]
then
	mkdir mocksoup_outputs
fi
cd mocksoup_outputs

if [ ! -d uclust_ref_99 ]
then
	mkdir uclust_ref_99
fi

# read a list and make a bash array
cd "${HOMEFOLDER}"/analysis/split_filter_outputs/
find * -type f -name "*_filtered_IDN.fasta.gz" | sed 's/_filtered_IDN.fasta.gz//g'> "${HOMEFOLDER}"/list.txt
cd "${HOMEFOLDER}"
PCRset_info=list.txt
PCRset_names=($(cut -f 1 "$PCRset_info" | uniq))
for PCRset in ${PCRset_names[@]}
do
	echo  "Now on PCRset" $INDEX
	INDEX=$((INDEX+1))

	mkdir ${PCRset}_temp
	cd ${PCRset}_temp
	cp "${HOMEFOLDER}"/analysis/split_filter_outputs/${PCRset}_filtered_IDN.fasta.gz .
	gunzip -d ${PCRset}_filtered_IDN.fasta.gz
	pick_otus.py -i ${PCRset}_filtered_IDN.fasta -r "${HOMEFOLDER}"/info/Sanger_ref.fas -m uclust_ref -s 0.99 -C -o ${PCRset}_uclust_ref_99
	mv ${PCRset}_uclust_ref_99/${PCRset}_filtered_IDN_otus.txt "${HOMEFOLDER}"/analysis/mocksoup_outputs/uclust_ref_99
	cd .. #now in countN folder
	rm -rf ${PCRset}_temp
done
