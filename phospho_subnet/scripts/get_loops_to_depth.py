"""
Given a sif file, enumerates loops (cyclic paths) up to a given depth.
Optionally the user can provide a list of start nodes. Otherwise, we will
search starting from all nodes.

NOTE:
This code will stop searching once a start node is reached, so if we 
have a loop like this:
YDR293C -> NPR74 -> YDR293C
We will not be able to identify loops of this format:
YDR293C -> NPR74 -> ... lots of other nodes ... -> YDR293C


USAGE: python get_loops_to_depth.py sif_file maxDepth [start_nodes]

sif_file is a network file in sif format (nodeA etype nodeB),

maxDepth is an integer marking the number of nodes allowed in the loop (counting start node twice),

and start_nodes is an optional filename for a tab-delim file with a list of
start nodes in the first column (eg, you could use it to indicate all kinases,
or any node of interest).

The loop length will really be one plus unique node count.
So, the loop YDR293C -> NPR74 -> YDR293C is considered length 3.

examples:

python get_loops_to_depth.py test_v3_maxall_0.75.sif 3 
finds all loops of format A -> B -> A

python get_loops_to_depth.py test_v3_maxall_0.75.sif 3 start_nodes.tab

chasman@cs.wisc.edu
"""
import sys

USAGE="USAGE: python get_loops_to_depth.py network.sif int(maxDepth) [start_nodes.tab]"

def main(argv):
	""" 
	Main function: 
		reads network
		parses max depth
		optionally reads start node file
	"""
	if len(argv)<3 or len(argv) > 4:
		print USAGE
		return 2

	(nw, allnodes)=read_nw(argv[1])
	maxDepth=int(argv[2])
	starts=allnodes
	if len(argv)==4:
		starts=read_nodes(argv[3])
	
	alloops=[]
	for s in starts:
		# start from node s, search for all loops up to max depth.
		loops=search([s], nw, maxDepth)
		alloops.extend(loops)
		
	# print the nodes in the found loops, separated by spaces.
	for loop in alloops:
		print " ".join(loop)
		
def search(loop, nw, maxDepth):
	"""
	Searches for loops: paths that start and end with the same node.
	Loops are otherwise acyclic (can't contain duplicates other than
	start/end).
	"""
	# return condition 1: same first and last node. we found a loop!
	# stop searching.
	if len(loop)>1 and loop[0]==loop[-1]:
		#print "\tFOUND LOOP:", loop
		return [loop]
	
	# return condition 2: hit max depth without finding loop.
	# return empty list as we are empty-handed.
	if len(loop)==maxDepth:
		#print "\tDepth exceeded: ABORT", loop
		return []
	
	
	# otherwise, start extending out from last node in loop
	curNode=loop[-1]
	nexts=nw.get(curNode,[])
	
	furtherLoops=[]
	
	#print loop, "to:"
	
	# for each outgoing neighbor, search it to completion.
	for n in nexts:
		# hit a cycle involving any non-start node? don't bother searching
		if n in loop[1:]:
			continue
			
		# make a new candidate loop by adding the current node to the 
		# current loop.
		newLoop=loop+[n]
		
		# recursive call to same function to continue searching.
		# this will return:
		#	a loop if we hit the same node again,
		# 	an empty list if we
		#		- exceed depth limit
		#		- can't find anymore places to search
		nextLoops=search(newLoop, nw, maxDepth)
		furtherLoops.extend(nextLoops)
	
	return furtherLoops
	
	
		
def read_nw(fn):
	"""
	Reads network from sif file.
	Assumes the network is directed:
	nodeA d/u nodeB
	
	Returns network as adjacency map { a : [b1, b2, ...] } 
	and set of all nodes.
	"""
	nw={}
	allnodes=[]
	with open(fn) as f:
		for row in f:
			# split on spaces
			sp=row.strip().split()
			nodeA=sp[0]
			nodeB=sp[2]
			
			if nodeA not in nw:
				nw[nodeA]=[]
			nw[nodeA].append(nodeB)		
			
			# store all the nodes
			allnodes.append(nodeA)
			allnodes.append(nodeB)
	
	return nw, set(allnodes)
		
def read_nodes(fn):
	"""
	Reads nodes from tab-delim file. Takes first column.
	"""
	nodes=[]
	with open(fn) as f:
		for row in f:
			sp=row.strip().split("\t")
			nodes.append(sp[0])
	
	return set(nodes)


if __name__=="__main__":
	sys.exit(main(sys.argv))
