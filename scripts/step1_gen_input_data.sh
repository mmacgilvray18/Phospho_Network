#!/bin/bash
# Generates short module identifiers for use in subnetwork inference programs.
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


LOC=../input_data/
mapper=map_module_names.py

# map source name to orf
naclmap='s/cdc14/YFR028C/; s/hog1/YLR113W/; s/pde2/YOR360C/;'
dttmap='s/ire1/YHR079C/; s/mkk1_2/MKK1MKK2/;'

for cond in DTT #NaCl DTT
do 
	# locate the module filename
	modfile=$(ls ${LOC}/${cond}/*sub*odules.txt)
	siffile=$(ls ${LOC}/${cond}/*SIF*.txt)
	
	modoutfn=${LOC}/${cond}/${cond}_module_names.txt
	pairfn=${LOC}/${cond}/${cond}_source_module_pairs.txt
	sifoutfn=${LOC}/${cond}/${cond}_mappedsif.txt
	nwoutfn=${LOC}/${cond}/${cond}_module_edges.txt
	ntypefn=${LOC}/${cond}/${cond}_node_types.txt
	
	n=1
	
	# remove if already exists
	if [ -e $modoutfn ]; then
		rm $modoutfn
	fi
	
	echo -e "source\tmodule_id\tmodule_pattern" > ${modoutfn}
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
		
		if [[ $cond == DTT ]]; then
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
	if [[ $cond == "NaCl" ]]; then
		echo "NaCl"
		smap=$naclmap
	elif [[ $cond == "DTT" ]]; then
		echo "dtt"
		smap=$dttmap
		
		# also make the name map for DTT
		printf "YOR231W\tMKK1MKK2\nYPL140C\tMKK1MKK2\n" > ${LOC}/${cond}/DTT_merge_map.txt
		echo "Created name map ${LOC}/${cond}/DTT_merge_map.txt"
	else
		echo "fail"
		exit
	fi
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
	
done
