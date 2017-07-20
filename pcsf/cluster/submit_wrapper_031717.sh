#!/bin/bash
# Submit PCSF runs with one of three deleted genes as the source and phospho
# submodules of that deletion as as targets.  Run PCSF separatley for each
# source and its respective targets.  Run with different seeds without random
# edge noise because all edges have the same confidence.  Use multiple values
# of beta and mu.  Sets up a new configuration file and uses environment
# variables to pass other arguments to forest.py.  Creates a summarization
# script to be run manually.

# Set the output directory for PCSF and HTCondor output
export outpath=../results/cdc14_hog1_pde2_sources_031717
mkdir -p $outpath

# Set the code paths for the Omics Integrator and msgsteiner dependencies
## This is the directory that contains scripts/forest.py
export oipath=/mnt/ws/virology/home/agitter/projects/OmicsIntegrator
## This is the path to the msgsteiner executable, including the executable name
export msgsteinerpath=/mnt/ws/virology/shared/lab_folders/GitterLab/progs/msgsteiner-1.3/msgsteiner

# Fixed parameters, sweep over multiple mu and beta
# Depth from root of tree
D=10
# Convergence parameter
g=1e-3
# Edge noise to compute a family of solutions
# Adding edge noise may make the solutions too dependent on the random
# edges that have lower cost in a particular run.  Instead rely on different
# beta and mu to obtain diverse networks.
r=0
# Cost of adding a new tree
# Use a small value because this is irrelevant when there is only one source node
w=0.01

# The path to the prize file above
export prizepath=../data
# The background network, including the path
export edgefile=../data/background_network.txt
# The source and prizes vary so set them below

# Try different values of mu, which penalizes high-degree nodes
for m in $(seq 0 0.01 0.1)
do

	# Try different prize scaling factors, which reward including prize nodes
	for b in 0.1 $(seq 0.5 0.5 5)
	do

		# Create the configuration file, removing an older copy of the file if it exists
		mkdir -p ${outpath}/conf
		filename=${outpath}/conf/conf_w${w}_b${b}_D${D}_m${m}_r${r}_g${g}.txt
		rm -f $filename
		touch $filename
		printf "w = ${w}\n" >> $filename
		printf "b = ${b}\n" >> $filename
		printf "D = ${D}\n" >> $filename
		printf "mu = ${m}\n" >> $filename
		printf "r = ${r}\n" >> $filename
		printf "g = ${g}\n" >> $filename

		# Set the remaining environment variables
		export conf=$filename
		export beta=$b
		export mu=$m
		export omega=$w
		
		# Iterate over the sources
		for source in cdc14 hog1 pde2
		do

			# Prize filename prefix (assume a .txt extension follows, e.g. prize.txt)
			# It will be used to create an output prefix for the Steiner forest networks
			export prizetype=${source}_prizes
			# A file listing the protein names that should be treated as source nodes,
			# including the path
			export sources=../data/${source}_sources.txt

			# Use different seeds for each run
			# Sample 10 forests per parameter combination
			for s in $(seq 1 10)
			do
				# Set the seed
				export seed=$s

				# Submit the job to HTCondor with the configuration file, params, and seed
				condor_submit submit_PCSF_031717.sub
			done
		done
	done
done

# Generate a wrapper script to summarize the family of forests
# This must be run after all PCSF runs terminate
# HTCondor can manage these dependencies with the Directed Acyclic Graph Manager
# but this strategy generalizes to other setups
# Set the name of the summarization script, which overwrites an existing
# file with the same name
# The script assumes that the summarization Python code resides in the same
# directory
sumscript=summarize_forests_031717.sh
rm -f $sumscript
touch $sumscript
printf "#!/bin/bash\n" >> $sumscript
printf "#Summarize a collection of Steiner forest samples\n" >> $sumscript
# Include all .sif files in the results directory, merging networks from the different sources
# Don't specify the --prizefile argument because different prizes were used for each source
printf "python summarize_sif.py --indir ${outpath} --pattern *optimalForest.sif --outfile ${outpath}/cdc14_hog1_pde2_summary\n" >> $sumscript
chmod u+x $sumscript
