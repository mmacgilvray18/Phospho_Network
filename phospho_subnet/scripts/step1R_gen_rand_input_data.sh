#!/bin/bash
# Processes files from Matt for use in ILP.
# This script works on the submodules derived from random subnets.
# First, generates the submodules file
# Generates simple IDs for submodules:
# Example: Induced_.....s..S..._cdc14_Induced_Amplified
# Pattern (_-delimited)
# 	Induced/repressed
#	Motif string			
#	Source					
#	Source knockout phenotype: Induced/Repressed_Amplified/Defective
# New pattern:
#	<source><I/R><A/D><number>	
# Example: cdc14IA0
# No phenotype module: np instead of source
set -eux

mapper=map_module_names.py

# map source name to orf
naclmap='s/cdc14/YFR028C/; s/hog1/YLR113W/; s/pde2/YOR360C/;'

nrand=4 # index of highest random nw
cond=NaCl 
LOC=../rand_ilp/input_data/${cond}_RAND

# loop over the random submodule connections
for i in `seq 0 $nrand`
do
	# locate the connections file
	siffileOrig=${LOC}/${cond}_SIF_Randomized_${i}_orig.txt
	siffile=${LOC}/${cond}_SIF_Randomized_${i}.txt

	# remove empty lines
	awk 'NF>1' ${siffileOrig} > $siffile

	# first, create the submodules file
	# sif format:
	#	Interactor_A	Edge_Type	Interactor_B	Annotation	SI	FDR_Score	Match_FDR
	#	Induced_.....s..S..._cdc14_Induced_Amplified	Constituent	YMR031C	
	# Submodule format:
	#	subModule	Source
	#	Induced_.....s..S..._cdc14_Induced_Amplified	cdc14

	modfile=${LOC}/${cond}_Submodules_Randomized_${i}.txt
	umod=${LOC}/${cond}_Submodules_Randomized_${i}_unique.txt
	grep Constituent ${siffile} | cut -f1 | sort -u > $umod

	printf "subModule\tSource\n" > ${modfile}
	paste <(cut -f 1 $umod) <(cut -f 3 -d"_" $umod) | sed 's/\tNo//' | sort -k2,2 >> ${modfile}


	modoutfn=${LOC}/${cond}_module_names_Randomized_${i}.txt
	pairfn=${LOC}/${cond}_source_module_pairs_Randomized_${i}.txt
	sifoutfn=${LOC}/${cond}_mappedsif_Randomized_${i}.txt
	nwoutfn=${LOC}/${cond}_module_edges_Randomized_${i}.txt
	ntypefn=${LOC}/${cond}_node_types_Randomized_${i}.txt
	
	n=1 #module counter

	printf "source\tmodule_id\tmodule_pattern\n" > ${modoutfn}
	while read mod src
	do		
		p1=i	# induced	
		p2=a	# amplified
			
		if [[ $mod == *"No_Phenotype"* ]];
		then
			src="np"	
			p2=""	
		fi
	
		if [[ $mod == *"Repressed"* ]];
		then
			p1=r
		fi
	
		if [[ $mod == *"Defective"* ]];
		then
			p2=d
		fi
	
		name=$src$p1$p2$n
		
		if [[ $cond == "DTT" ]]; then
			name=${name/_2}
		fi
		echo -e "$src\t$name\t$mod" | sed $'s/\r//g' >> ${modoutfn}
		
		#name=${src}${p1}${p2}${n}
		n=$((n+1))
	done < <(tail -n +2 ${modfile})

	echo "Created $modoutfn"
	wc -l $modfile
	wc -l $modoutfn
	
	# map source name to ORF
	smap=$naclmap
	
	tail -n +2 $modoutfn | sed "$smap" | grep -v "No_Phenotype" > ${pairfn}
	
	# map names in Matt's SIF file to these IDs 
	python $mapper $modoutfn $siffile 0 2 > ${sifoutfn}
	echo "Created ${sifoutfn}"
	wc -l $siffile
	wc -l $sifoutfn
	
	# Convert to my format
	# get possible match types
	alltypes=$(cut -f 2 $siffile | tail -n +2 | sort -u)
	etypes=$(echo $alltypes | sed 's/ /|/g')
	
	# header
	head -n 1 $sifoutfn | sed "s/Edge_Type/matchType=CatSet(${etypes})/" | awk '{print $2"\t"$1"\t"$3"\tDir\tSign\tAnnotation=CatSet(Kinase|Phosphatase|Yes)"}' > $nwoutfn
	# remove dashes here
	tail -n +2 $sifoutfn | awk '{print $2"\t"$1"\t"$3"\t1\t0\t"$4"\t"$5}' | sed 's/Kinase\tYes/Kinase|Yes/; s/Phosphatase\tYes/Phosphatase|Yes/; s/-//g' >> $nwoutfn
	
	wc -l $nwoutfn
	
	# make file with source and module node types
	#printf "name\tntype=Discrete(module|source)\n" > $ntypefn
	cut -f 1 $pairfn | sort -u | awk '{print $1"\tsource"}' > $ntypefn
	# get module names from original map; want to include the NPs
	cut -f 2 $modoutfn | tail -n +2 | sort -u | awk '{print $1"\tmodule"}' >> $ntypefn

	wc -l $pairfn
	wc -l $modoutfn
	wc -l $ntypefn

done # loop over random nets

