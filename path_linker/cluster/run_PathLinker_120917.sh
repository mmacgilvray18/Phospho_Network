#!/bin/bash
# Run PathLinker using the parameters specified below

echo _CONDOR_JOB_IWD $_CONDOR_JOB_IWD
echo Cluster $cluster
echo Process $process
echo RunningOn $runningon

# PathLinker doesn't support a named output directory
cd ../results/cdc14_hog1_pde2_sources_pl_120917

plpath=../../../path-linker/PathLinker
# Number of source-target paths to find
k=100000
inpath=../../data/PathLinker
netfile=${inpath}/PathLinker_Background_Network.txt

# Must have already created the path-linker environment from environment.yml
source activate path-linker

# Iterate over the sources
for source in Cdc14 Hog1 Pde2
do
	# Create the PathLinker command
	CMD="python $plpath/run.py \
		$netfile \
		${inpath}/${source}_Source_Target.txt \
		--PageRank \
		--output ${source}_ \
		-k $k"

	echo $CMD
	$CMD
done

# Track the contents of the conda environment
echo "path-linker conda environment"
conda list
echo

# Track the status of the gasch-nacl-phospho-net repostiory, which can be used
# to recover the input data files
echo "gasch-nacl-phospho-net version:"
git log -1
echo

# Track which version was run with the git commit log
echo "Path Linker version:"
cd $plpath
git log -1
