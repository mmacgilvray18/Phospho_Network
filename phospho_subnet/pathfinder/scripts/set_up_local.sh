#!/bin/bash
# Sets up a run with choice of run-type (pool/lesion) and choice of condition (DTT/NaCl).
# Add a directory name of your choice.
# You will then descend into that directory to run GAMS on the model.
# USAGE: ./set_up.sh condition runtype outdir

if [[ $# != 4 ]]; then
	echo "USAGE: ./set_up.sh condition pathlen runtype outdir"
	echo "Where condition = {DTT, NaCl}"
	echo "pathlen = { 3,4}"
	echo "runtype = {lesion, pool}"
	echo "outdir = your choice of output directory (must be new)"
	exit
fi

cond=$1
pathlen=$2
runtype=$3
outdir=$4

wd=$(pwd)

datafile=${wd}/${cond}/${cond}_${pathlen}_v4.gms
if [[ ! -e $datafile ]]; then
	echo "Cannot find data file $datafile"
	exit 2
fi

template=${wd}/model/template_run_${runtype}.gms
if [[ ! -e $template ]]; then
	echo "Cannot find template gams file $template"
	exit 2
fi

if [[ -d $outdir ]]; then
	echo "$outdir already exists."
	exit 2
fi

#########

# make directory and copy in relevant stuff -- template and CPLEX options
mkdir -p $outdir
sed "s:{DATA}:${datafile}:" $template > $outdir/run_${runtype}.gms
cp model/cplex* ${outdir}

echo "Directory $outdir is ready for you to run GAMS on ${outdir}/run_${runtype}.gams"

