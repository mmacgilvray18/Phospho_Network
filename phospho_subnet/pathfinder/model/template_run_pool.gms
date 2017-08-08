** Main program **
** Maximize source-target connections, max matched SIs, minimize nodes, maximize paths.
** Solution pool version: Get 1000 solutions.
** Template: Replace {DATA} with location of your data file.

$phantom null

$offlisting
$offsymxref offsumlist
option 
	limrow=0,
	limcol=0;


* Load data file
$Include {DATA}

alias (node,node1);

* make pair paths	
Set pairPath(node,node,path)	"paths that satisfy receptor-source pairs";
pairPath(node,node1,path)=no;
loop( pair(node,node1),
	loop(pstart(node,path)$pnode(node1,path), pairPath(node,node1,path)=yes;););

* Load model file
$Include ../model/model.gms


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

* Commented out -- Make sure that at least this fraction of pairs are connected
*Scalar satFrac /0.0/;
*Equation mustConnect	"require that only a few pairs are not connected";
*mustConnect .. sum( pair(node,node1), sat(node,node1) ) =g= satFrac*card(pair);

Model maxConn /all/;
maxConn.optcr=0.0;
maxConn.optca=0.0;
maxConn.reslim=100000;
maxConn.optfile=1;
solve maxConn using mip max conn;

Scalar maxConnSol	"maximum connections possible" /0/;
maxConnSol=conn.l;

display maxConnSol;
 
Equation setMaxConn	"Set maximum connections within tolerance";
setMaxConn .. conn =g= 0.95*maxConnSol;
** This is it!!!

* Maximize SI edges
Set si(edge)	"SI-type edges";
si(edge)=no;
si(edge)$(motifmatchValueE(edge) or motifunmatchedValueE(edge) or unknownrecognitionmotifValueE(edge) or SharedInteractionValueE(edge))=yes;
Free variable matchCount	"total match edges ";
Equation countMatch	"count up active SI edges";
countMatch .. matchCount =l= sum(si(edge), x(edge));

Model maxMatchModel /all/;
maxMatchModel.optcr=0.0;
maxMatchModel.optca=0.0;
maxMatchModel.reslim=100000;
maxMatchModel.optfile=1;
solve maxMatchModel using mip max matchCount;

* Fix the value
Scalar maxMatches "max matches calculated"	/0/;
maxMatches=matchCount.l;

Equation fixMatchCount "Set matches active within tolerance";
fixMatchCount .. matchCount =g= 0.95*maxMatches;

* Minimize nodes - find multiple solutions.
Model minNodeModel /all/ ;
minNodeModel.optca=5;
minNodeModel.reslim=100000;
minNodeModel.optfile=2;

solve minNodeModel using mip min nodeCount;

* Maximize paths to improve interpretation.
Model maxPathModel /all/ ;
maxPathModel.optcr=0;
maxPathModel.optca=0;
maxPathModel.reslim=100000;
maxPathModel.optfile=1;

Set solNode(node)	"nodes on in a solution";

Set 	soln 	possible solutions in the solution pool 
	/file1*file1000/;
set	solnpool(soln) actual solutions;
file fsol;

execute_load 'solnpool.gdx', solnpool=Index;

Scalar s	"count solns"		/0/;
loop(solnpool(soln),
	put_utility fsol 'gdxin' / solnpool.te(soln):0:0;
	putclose;
	execute_loadpoint;

	solNode(node)=no;
	solNode(node)$(y.l(node)>0)=yes;
	
	y.fx(node)=0;
	y.fx(solNode)=1;	
	
	solve maxPathModel using mip max pathCount;
	
	put_utility 'gdxout' / 'path_sol_' s:0:0 '.gdx';
	execute_unload sigma,x,y;	
	s=s+1;
);
