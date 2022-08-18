#!/bin/bash
set -e
set -u
set -o pipefail
###########################################################################################
###########################################################################################
# a shell script to put UMI on header
###########################################################################################
###########################################################################################

# Usage: bash UMI_to_header.sh
# run on my mac pro

# need to prepare tagF.txt, tagR.txt, mock_PCRset_Info.txt. Put them in Folder info/

PIPESTART=$(date)

HOMEFOLDER=$(pwd)

echo "Home folder is "${HOMEFOLDER}""

# set variables
HOMEFOLDER="/Users/apple/src/qseq/" # do not have a ~ in the path
INDEX=1

. ~/.linuxify
cd "${HOMEFOLDER}"/analysis/

if [ ! -d split_filter_outputs ] # if directory split_filter_outputs does not exist
then
	mkdir split_filter_outputs
fi

cd Filter_outputs
cat BegumFilter_0.6_copies4.fna | tr '\t' '#' | tr '\n' ' ' | sed 's/ >/\n>/g' > FilteredReadsOneline.txt
cd "${HOMEFOLDER}"/analysis/

# read in folder list and make a bash array
cat "${HOMEFOLDER}"/info/mock_PCRset_Info.txt | tr "\t" "_" > PCRs_list.txt
PCRset_info=PCRs_list.txt # put folderlist.txt into variable
PCRset_names=($(cut -f 1 "$PCRset_info" | uniq)) # convert variable to array this way

for PCRset in ${PCRset_names[@]}  # ${PCRset_names[@]} is the full bash array
do
	echo "Now on PCRset" ${INDEX} of ${#PCRset_names[@]}". Moved back to starting directory:"
	sample=$(sed -n "${INDEX}, ${INDEX}p" "${HOMEFOLDER}"/info/mock_PCRset_Info.txt | cut -f 1)
	FtagName=$(sed -n "${INDEX}, ${INDEX}p" "${HOMEFOLDER}"/info/mock_PCRset_Info.txt | cut -f 2)
	RtagName=$(sed -n "${INDEX}, ${INDEX}p" "${HOMEFOLDER}"/info/mock_PCRset_Info.txt | cut -f 3)
	pool=$(sed -n "${INDEX}, ${INDEX}p" "${HOMEFOLDER}"/info/mock_PCRset_Info.txt | cut -f 4)
	INDEX=$((INDEX+1))
	pwd
	
	mkdir ${PCRset}_temp
	cd ${PCRset}_temp
	Ftagseq=$(grep -w "${FtagName}" "${HOMEFOLDER}"/info/tagF.txt | cut -f 2) #Tagseq of this tag/PCRset
	Rtagseq=$(grep -w "${RtagName}" "${HOMEFOLDER}"/info/tagR.txt | cut -f 2)

	if grep -q -m 1 ">${sample}#" "${HOMEFOLDER}"/analysis/Filter_outputs/FilteredReadsOneline.txt
	then
    	# make a ${sample}_ref.fna file from Begum filter's results
    	grep ">${sample}#" "${HOMEFOLDER}"/analysis/Filter_outputs/FilteredReadsOneline.txt | tr '\n' ' ' | sed 's/ /\n/g' > ${sample}_ref.fna
    	# obtain the reads belonging to this PCRset of this sample
    	cat "${HOMEFOLDER}"/analysis/PandaSeq_outputs/${pool}_merged.fasta | seqkit grep -s -i -p ${Ftagseq} | seqkit grep -s -i -p ${Rtagseq} -w 0 > ${PCRset}_temp.fasta
    	# because the function of search_exact in VSEARCH requests that the input reads and the ref sequences should have the same length,  so we have to delete the primers and tags before --search_exact.
    	sed 's/^.*GG.AC.GG.TGAAC.GT.TA.CC.CC//g' ${PCRset}_temp.fasta | sed 's/TG.TT.TT.GG.CA.CC.GA.GT.TA.*//g' > ${PCRset}_temp_noprimertag.fasta
    	# filter the reads based on the ref file from Begum filter's results
    	vsearch --search_exact ${PCRset}_temp_noprimertag.fasta -db ${sample}_ref.fna -matched hits.fasta -strand plus
    	# obtain the whole length for these filtered reads
    	grep -E "^>" hits.fasta | gsed 's/>//g' > hits_IDs.txt
    	seqtk subseq ${PCRset}_temp.fasta hits_IDs.txt > ${PCRset}_filtered.fasta
    	# because the filtered reads' IDs are from the raw data and don't have the sample information, here we make a new ID for each filtered read which include the sample's name and a code
		seq $(grep -c ">" ${PCRset}_filtered.fasta) | sed -e "s/^/>${sample}_/g" > ID.txt
		# obtain the UMI info
		grep -Ev "^>" ${PCRset}_filtered.fasta > seqs.txt
		sed 's/^.*GG.AC.GG.TGAAC.GT.TA.CC.CC//g' seqs.txt | sed 's/TG.TT.TT.GG.CA.CC.GA.GT.TA.*//g' > seqs_no_primertag.txt
		sed "s/.${Ftagseq}.*${Rtagseq}.//g" seqs.txt | sed "s/^/N/g" | sed "s/$/N/g" > Ncontents.txt
		# put the UMI info into the new IDs and make a filtered reads' file with new IDs including UMIs
		paste ID.txt Ncontents.txt | tr '\t' '_' > ID_N.txt
		paste ID_N.txt seqs_no_primertag.txt | grep -E "_N.........N" | tr '\t' '\n' > ${PCRset}_filtered_IDN.fasta
		gzip -9 ${PCRset}_filtered_IDN.fasta
    	mv ${PCRset}_filtered_IDN.fasta.gz  "${HOMEFOLDER}"/analysis/split_filter_outputs/
		cd ..
    	rm -rf ${PCRset}_temp
	else
		cd ../
		rm -rf ${PCRset}_temp
		continue
	fi
#
done

echo "UMI to header done"
