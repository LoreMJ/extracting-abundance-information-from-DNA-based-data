#!/bin/bash
set -e
set -u
set -o pipefail
###########################################################################################
###########################################################################################
# a shell script to generate the OTU table with reads' numbers
###########################################################################################
###########################################################################################

# Usage: bash for_OTUtable_withREAD.sh
# run on my mac pro

PIPESTART=$(date)

HOMEFOLDER=$(pwd)

echo "Home folder is "${HOMEFOLDER}""

# set variables
HOMEFOLDER="/Users/apple/src/qseq/" # do not have a ~ in the path
INDEX=1

# read a list and make a bash array
cd "${HOMEFOLDER}"/
PCRset_info=list.txt
PCRset_names=($(cut -f 1 "$PCRset_info" | uniq))

cd  analysis/mocksoup_outputs/

if [ ! -d for_OTUtable_withREAD ]
then
	mkdir for_OTUtable_withREAD
fi

. ~/.linuxify
echo "OTUID" > OTUIDs.txt
grep ">" "${HOMEFOLDER}"/info/Sanger_ref.fas | sed 's/>//g' >> OTUIDs.txt

for PCRset in ${PCRset_names[@]}
do
	echo "Now on PCRset" ${INDEX} of ${#PCRset_names[@]}
	INDEX=$((INDEX+1))

	mkdir ${PCRset}_temp
	cp "${HOMEFOLDER}"/analysis/mocksoup_outputs/uclust_ref_99/${PCRset}_filtered_IDN_otus.txt ${PCRset}_temp/
	cd ${PCRset}_temp
	echo "${PCRset}" > ${PCRset}_countReads.txt
	sed '1d' ../OTUIDs.txt > OTUIDs.txt
	otuid_info=OTUIDs.txt
	otuid_names=($(cut -f 1 "$otuid_info" | uniq))
	NUMBER=1
	for otuid in ${otuid_names[@]}
	do
  		NUMBER=$((NUMBER+1))
  		if grep -q -m1 "^${otuid}" ${PCRset}_filtered_IDN_otus.txt
  		then
  			grep "^${otuid}" ${PCRset}_filtered_IDN_otus.txt > ${otuid}_0.txt
  			cat ${otuid}_0.txt | tr '\t' ' ' | sed 's/ /\n/g' > ${otuid}_1.txt
  			sed '1d' ${otuid}_1.txt | sed "s/^.*_//g" > ${otuid}_2.txt
			grep -c "^N" ${otuid}_2.txt >> ${PCRset}_countReads.txt # count reads for each OTU in each PCRset
			rm -f ${otuid}*.txt
  		else
			echo "0" >> ${PCRset}_countReads.txt
			continue
		fi
	done

	mv ${PCRset}_countReads.txt ../for_OTUtable_withREAD/
	cd ../
	rm -rf ${PCRset}_temp


done

paste OTUIDs.txt for_OTUtable_withREAD/* > mock_countReads_OTUtable.txt
rm -f OTUIDs.txt
rm -rf for_OTUtable_withREAD
	

echo "Pipeline started at $PIPESTART"
NOW=$(date)
echo "Pipeline ended at   $NOW"

