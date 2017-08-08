#!/bin/bash
# Example script for extracting paths for specific nodes and confidence values.
# Note that genes must be referred to by YORF without dashes.


# USAGE: python get_paths.py pathfile_with_confs.tab [min_conf ORF1 ... ] > yourfilename.sif
python ../scripts/get_paths.py NaCl_subnet_results/NaCl_3_paths_with_confs.tab 0.75 > NaCl_3_0.75.sif
python ../scripts/get_paths.py NaCl_subnet_results/NaCl_3_paths_with_confs.tab 0.75 YLR248W > NaCl_3_0.75_RCK2.sif

# Example for augmenting with background network.
for cond in NaCl DTT; do
	for f in ${cond}_subnet_results/${cond}_3_*0.75.sif;
	do
		python ../scripts/augment.py $f ${cond}_cytoscape_files/${cond}_3_v4_background.sif > ${f/.sif/_plusbg.tab}
	echo "Augmented $f with background: ${f/.sif/_plusbg.tab}"
	done	
done
