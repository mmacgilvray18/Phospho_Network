#!/bin/bash
# Run PCSF using the parameters set in the wrapper file

echo _CONDOR_JOB_IWD $_CONDOR_JOB_IWD
echo Cluster $cluster
echo Process $process
echo RunningOn $runningon

prizefile=${prizepath}/${prizetype}.txt
outlabel=${prizetype}_beta${beta}_omega${omega}_mu${mu}_seed${seed}

# Create the Forest command
# Do not run with noisy edges because the r parameter in the configuration
# file is already being used to add edge noisy with msgsteiner if noise
# is needed
CMD="python $oipath/scripts/forest.py \
	-p $prizefile \
	-e $edgefile \
	-c $conf \
	-d $sources \
	--msgpath=$msgsteinerpath \
	--outpath=$outpath \
	--outlabel=$outlabel \
	--cyto30 \
	--noisyEdges=0 \
	-s $seed"

echo $CMD
$CMD

# Track the status of the gasch-nacl-phospho-net repostiory, which can be used
# to recover the input data files
echo "gasch-nacl-phospho-net version:"
git log -1
echo

# Track which version was run with the git commit log
echo "Omics Integrator version:"
cd $oipath
git log -1
