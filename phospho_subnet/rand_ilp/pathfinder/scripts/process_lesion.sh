#!/bin/bash
# Processes gams output files into confidence files and cytoscape-friendly files
# Args: top dir of analysis, name of path contents file (eg, data/nacl_v2_all_paths.tab),
#	prefix for output files
# Source nodes hard-coded for NaCl and DTT 
# You should call this script from a parent directory that contains the script directory AND gams output directory

if [[ $# != 3 ]]; then
	echo "USAGE: bash scripts/process.sh analysis_directory path_file.tab output_prefix"
	echo "analysis_directory is where you ran GAMS"
	echo "path_file.tab contains the node/edge contents of paths"
	echo "output_prefix is a short string that identifies this analysis (no directory names needed)"
	exit 2
fi
topdir=$1
pathfile=$2
prefix=$3


if [[ $pathfile == *"NaCl"* ]]; then
	sources="YFR028C YLR113W YOR360C"
elif [[ $pathfile == *"DTT"* ]]; then
	sources="YHR079C MKK1MKK2"
else
	echo "I don't recognize the condition for $pathfile"
	exit
fi

# look for hidden nodes

#prefix=${main/.gms/}
#prefix=${prefix/run_/}
#prefix=${prefix/_sampler/}
#echo $prefix, $main

wd=$(pwd)
cd $topdir

if [[ ! -e ${wd}/scripts/add_path_confs.py ]]; then
	echo "Cannot find scripts for processing -- you should call this program from above the scripts directory (or update the paths)"
	exit
fi	
if [[ ! -e hidden_nodes.tab ]]; then
	echo "Can't find required hidden_nodes.tab file in $topdir"
	exit
fi

for d in gdx result cyto; do
	if [[ ! -d $d ]]; then
		mkdir $d
	fi
done

mv *.gdx gdx
# CHECK SCRIPT PATHS HERE
ls
bash ${wd}/scripts/dump_information_lesion.sh gdx result path hidden_nodes.tab 

python ${wd}/scripts/add_path_confs.py ${wd}/${pathfile} result/path_sigma.tab > cyto/${prefix}_paths_with_confs.tab 


for f in 0.50 0.75; do
	python ${wd}/scripts/make_sif_update.py cyto/${prefix}_paths_with_confs.tab $f > cyto/${prefix}_${f}.sif
	# for each source...
	for s in ${sources}; do
		python ${wd}/scripts/make_sif_update.py cyto/${prefix}_paths_with_confs.tab $f $s > cyto/${prefix}_${s}_${f}.sif
	done	

	echo "Made ${topdir}/cyto/${prefix}_${f}.sif"
done

cd $wd


