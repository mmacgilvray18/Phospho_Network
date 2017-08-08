#!/bin/bash
# Runs pathfinder for DTT and NaCl.
# Will put output into a test subdirectory directly in this same directory (DTT_test, NaCl_test)

for cond in DTT #NaCl
do
	echo "Running pathfinding for $cond"
	date
	if [ ! -d ${cond}_test ]; then
		mkdir ${cond}_test	
	fi
	
	java -jar phosphonet_v4_pathfinder.jar ../pathfinder/${cond}/${cond}_phosphonet_v4.config > run_${cond}_pathfinder.log
	echo "Check output in ${cond}_test"
	date
done
