#!/bin/bash
# Runs pathfinder on randomized input data.
# Will create directories named NaCl_rand${i}
set -eux

cond=NaCl
template=../rand_ilp/config/${cond}_rand_template.config
nethome=../rand_ilp/input_data/network/
datahome=../rand_ilp/input_data/NaCl_RAND
outhome=../rand_ilp/pathfinder

# we need to fill in these variables:
# {BGNET} -- random bg net
# {MODEDGES} -- module edge SIF file
# {NTYPES} -- source/submodule node info
# {PAIRFILE} -- source-target pair file
# {OUT}  -- directory to put output in
# {OUTPREF} -- file prefix for output

nrand=4
for i in `seq 1 $nrand`;
do

	echo "Running pathfinding for $cond"
	date

	OUT=${outhome}/${cond}_rand${i}
	if [[ ! -d $OUT ]]; then
		mkdir -p $OUT
	fi
	OUTPREF=${OUT}

	# generate config file
	CONFIG=${OUT}/${OUTPREF}_config.txt
	echo "Generating $CONFIG"
	
	BGNET=$(ls ${nethome}/rewired-BG-network-*_${i}_bgnet.txt | head -n1)
	MODEDGES=${datahome}/${cond}_module_edges_Randomized_${i}.txt
	PAIRFILE=${datahome}/${cond}_source_module_pairs_Randomized_${i}.txt
	NTYPES=${datahome}/${cond}_node_types_Randomized_${i}.txt

	sed "s:{BGNET}:${BGNET}:; s:{MODEDGES}:${MODEDGES}:; s:{NTYPES}:${NTYPES}:; s:{PAIRFILE}:${PAIRFILE}:; s:{OUT}:${OUT}:; s:{OUTPREF}:${OUTPREF}:;" ${template} > ${CONFIG}
	
	java -jar phosphonet_v4_pathfinder.jar ${CONFIG} > ${OUT}/run_${cond}_rand${i}_pathfinder.log
	echo "Check output in ${OUT}"
	echo "See ${OUT}/run_${cond}_rand${i}_pathfinder.log"
	date

done # end loop over bgnets

