# Uses gdxdump to dumps node, edge, path relevance variables
# from gdx files in directory $1
# into a directory specified by $2
gdxdir=$1
outdir=$2
pref=$3

if [[ $# != 3 ]]; then
	echo "dump_info gdxdir outdir prefix"
	exit
fi

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

python /mnt/ws/home/dchasman/opt-cluster/phosphonet_v4/scripts/gather_path_solution_info.py "${outdir}/${pref}*_dump" ${outdir}/${pref}

