# Uses gdxdump to dumps node, edge, path relevance variables
# from gdx files in directory $1
# into a directory specified by $2
gdxdir=$1
outdir=$2
pref=$3
hidden=$4

if [[ $# != 4 ]]; then
	echo "dump_info gdxdir outdir prefix hidden_node_file"
	exit
fi

for f in $hidden $outdir $gdxdir; do
	if [[ ! -e $f ]]; then
		echo "Can't find $f"
		exit
	fi
done

echo "Deleting existing dump files"
echo $(pwd)

rm -f ${outdir}/*_dump

inpref=""
for fn in ${gdxdir}/${pref}*gdx
do
	f=${fn##*/}  # remove path
	f=${f%.*} # remove extension
	inpref=$f
	if [[ "$inpref" == "solnpool" ]]; then
		continue
	fi
	for i in sigma x y 
	do
		gdxdump ${gdxdir}/${f}.gdx symb=${i} Format=csv | grep -v ",0" | grep -v "Val" | awk -v i="${i}" -F"," '{print 	i"\t"$1"\t"$2}' >> ${outdir}/${f}_dump
		
	done	
done

# last sol ID is number of solutions
totsols=$(tail -n 1 $hidden | awk '{printf "%d", $2}')

python ../scripts/gather_path_solution_info_lesion.py "${outdir}/${pref}*_dump" ${outdir}/${pref} $hidden $totsols

