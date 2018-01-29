** Main program **
** This version doesn't check for subsets of SIs
** Maximize source-target connections, max matched SIs, minimize nodes, maximize paths.
** Lesion version for generating multiple sols: 
**		Run whole optimization process.
**		Hide a random selection of non-module/non-source nodes CHOSEN in the last solution.
**		Repeat x 1000
** Template: Replace /mnt/ws/home/dchasman/nacl/phosphonet/phospho_v4/pathfinder/NaCl_rand2/NaCl_rand2.gms with your data file path

$phantom null

$offlisting
$offsymxref
option
        limrow=0,
        limcol=0,
	solprint = off;


* Load data file
$Include /mnt/ws/home/dchasman/nacl/phosphonet/phospho_v4/pathfinder/NaCl_rand2/NaCl_rand2.gms 

alias (node,node1);



* make pair paths	
Set pairPath(node,node,path)	"paths that satisfy receptor-source pairs";
pairPath(node,node1,path)=no;
loop( pair(node,node1),
	loop(pstart(node,path)$pnode(node1,path), pairPath(node,node1,path)=yes;););

* Load model file - updated for new GAMS
$Include ../model/model_2018.gms


* First, find max number of pairs that we can connect
Free variable
	conn	"Count number of connected pairs";
Binary variable 	
	sat(node,node1) "pairs connected";

Equation	
	satPair(node,node)	"Assess pair connection - at least one path active"
	count	"count up pairs";
	
satPair(pair(node,node1)) .. sat(node,node1) =l= sum( pairPath(node,node1,path), sigma(path) );
count .. conn =e= sum( pair(node,node1), sat(node,node1));

* Make sure that at least this fraction of pairs are connected
*Scalar satFrac /0.0/;
*Equation mustConnect	"require that only a few pairs are not connected";
*mustConnect .. sum( pair(node,node1), sat(node,node1) ) =g= satFrac*card(pair);

Model maxConn /all/;
maxConn.optcr=0.0;
maxConn.optca=0.0;
maxConn.reslim=100000;
maxConn.optfile=1;

Scalar maxConnSol	"maximum connections possible" /0/;

display maxConnSol;
 
Equation setMaxConn	"Set maximum connections (no tolerance)";
setMaxConn .. conn =g= 1.0*maxConnSol;

* Maximize SI edges
Set si(edge)    "SI-type edges";
si(edge)=no;
* si(edge)$(motifmatchValueE(edge) or motifunmatchedValueE(edge) or unknownrecognitionmotifValueE(edge) or SharedInteractionValueE(edge))=yes;
si(edge)$(SharedInteractionValueE(edge))=yes;
Free variable matchCount        "total SI edges ";
Equation countMatch     "count up active SI edges";
countMatch .. matchCount =l= sum(si(edge), x(edge));

Model maxMatchModel /all/;
maxMatchModel.optcr=0.0;
maxMatchModel.optca=0.0;
maxMatchModel.reslim=100000;
maxMatchModel.optfile=1;

* Fix the value
Scalar maxMatches "max matches calculated"	/0/;

Equation fixMatchCount "Set matches active within tolerance";
fixMatchCount .. matchCount =g= 1.0*maxMatches;

* Minimize nodes - find one solution
Model minNodeModel /all/ ;
minNodeModel.optca=0;
minNodeModel.reslim=100000;
minNodeModel.optfile=1;

* Maximize paths to improve interpretation.
Model maxPathModel /all/ ;
maxPathModel.optcr=0;
maxPathModel.optca=0;
maxPathModel.reslim=100000;
maxPathModel.optfile=1;

Set solNode(node)	"nodes on in a solution";

Set canHide(node)	"non-module, non-source, nodes in paths";
canHide(node)=no;
loop(pnode(node, path)$(not moduleValueN(node) and not sourceValueN(node)),
	canHide(node)=yes;);
Scalar numhide	/0/;
numhide=card(canHide);
display numhide;
display canHide;
	

Set 	soln 	possible solutions in the solution pool 
	/file1*file1000/;
file fsol / hidden_nodes.tab / ;

* don't hide anything the first time
Scalar s	"count solns"		/0/;
Scalar nFrac	"fraction of last sol's nodes to hide next sol"	/ 0.05 /;
loop(soln,	
	solve maxConn using mip max conn;

	maxConnSol=conn.l;

	solve maxMatchModel using mip max matchCount;

	maxMatches=matchCount.l;
	solve minNodeModel using mip min nodeCount;

	solNode(node)=no;
	solNode(node)$(y.l(node)>0)=yes;
	
	y.fx(node)=0;
	y.fx(solNode)=1;	
	
	solve maxPathModel using mip max pathCount;
	put fsol;
	put_utility 'gdxout' / 'path_sol_' s:0:0 '.gdx';
	execute_unload sigma,x,y,d,satPair;	
	s=s+1;
	
* reset variables that might affect later solutions
	y.lo(node)=0; y.up(node)=1;
	x.lo(edge)=0; x.up(edge)=1;
	sigma.lo(path)=0; sigma.up(path)=1;
	d.lo(ppi)=0; d.up(ppi)=1;
	satPair.lo(pair)=0; satPair.up(pair)=1;	
	

* randomly hides some % of the canHide nodes that WERE chosen last time.
	put fsol;	
	loop(canHide(node)$solNode(node),
		if (uniform(0,1) < nFrac, y.fx(node)=0; put node.tl '':8 s '':8 'hidden'/;););
	
);
