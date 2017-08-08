"""
Collects the variable info in a directory of dumped gdx files and 
produces three files: path confidence, node confidence, and edge confidence.

Special lesion version: look at all solutions, but take into account
which solutions have a particular node held aside.

gather_path_solution_info_lesion.py gdx_dump_pattern output_prefix hidden_nodes.tab
"""
import sys
import glob

# argv 1: pattern for matching dump files; eg: path_sol*dump
ls = glob.glob(sys.argv[1])
outpref = sys.argv[2]
print "Writing confs to files with prefix %s" % outpref

# third argument: hidden_nodes file (node sol "hidden" -- weird spacing due to GAMS)
hiddenFN=sys.argv[3]
totsols=int(sys.argv[4])

hidden={}	# {node : [solIDs hidden]}
solsfound=set()
with open(hiddenFN) as f:
	for line in f:
		line=line.strip()
		if len(line)==0:
			continue
		sp=line.split()	# any whitespace
		node=sp[0]
		sol=int(float(sp[1]))
		#solsfound.add(sol)
		if sol>=totsols:
			print >> sys.stderr, "Solution %f is bogus (> totsols %d); skipping b/c probably last in file and artifact" % (sol,totsols)
			continue
			
		if node not in hidden:
			hidden[node]=[]
		hidden[node].append(sol)
		
solIDs=range(0, len(ls))	# all IDs to consider

# number of times we see a path
paths={}
nodes={}
edges={}
#dirs={}
gomap = {"sigma":paths, "x":edges, "y":nodes} #,"d":dirs}

tot=0
for fn in ls:
	#print "#", fn
	bad=False
	
	# assume filename in format blah_blah_ID_dump
	sol=int(fn.split("_")[-2])
	
	if sol not in solIDs:
		#print "skipping %d" % sol
		continue
	#print "using %d" % sol
	
	seen=set()
	with open(fn) as f:
		for line in f:
			# didn't work
			if "Symbol not found" in line:
				print >> sys.stderr, "Bad file: %s" % fn
				bad=True
				break
			
			# var\t"name"\tsetting
			sp=line.strip().split("\t")
			symbol = sp[1][1:-1]

			if symbol in seen:
				continue
			try:
				setting = int(sp[2])		
			except IndexError, ValueError:
				print >> sys.stderr, "Bad format", fn, sp
				bad=True
				break

			themap = gomap[sp[0]]
			themap[symbol] = themap.get(symbol,0)+1
		if not bad:
			tot+=1

#print "# %d total solutions" % tot

total=float(tot)

for var in gomap.keys():
	with open("%s_%s.tab" % (outpref, var), "w") as f:
		f.write("#item\tconf(%s) (out of %d total sols)\n" % (var, tot))
		for (s, count) in gomap[var].items():
			# if nodes, divide by total minus number in which the node is hidden
			if var=="y":
				hideTot=len(hidden.get(s,[]))
				
				if count>total>1:
					print >> sys.stderr, "%s" % s
				if count/(total-hideTot) > 1:
					print >> sys.stderr, "%s\t%f\t%d\t%d\n" % (s, count/(total-hideTot), count, hideTot)
				f.write("%s\t%f\t%d\t%d\n" % (s, count/(total-hideTot), count, hideTot))
			else:
				f.write("%s\t%f\n" % (s, count/total))
			
			#print s, count, total
			if (count/total) > 1:
				print >> sys.stderr, "Too many entries per symbol in your dump files. (%s, %d, %d)" % (s, count, total)





