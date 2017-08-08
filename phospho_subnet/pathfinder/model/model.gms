** Simple model ***

Binary Variables
	x(edge)	"edge relevance"
	y(node)	"node relevance"
	d(ppi)	"edge direction"
	sigma(path)	"path relevance";
	
Free Variables
	nodeCount	"total relevant nodes"
	pathCount	"total relevant paths";
		
Equations
	usePath1(edge,path)	"path must contain all relevant edge"
	usePath2(edge)	"edge must be in a relevant path";
	
usePath1(pedge(edge,path)) .. sigma(path) =l= x(edge);
usePath2(edge) .. x(edge) =l= sum( pedge(edge,path), sigma(path));

Equations
	pathDir1(ppi,path)	"fwd path-edge pairs"
	pathDir2(ppi,path) "back path-edge pairs";
	
pathDir1(fwd(ppi,path)) .. sigma(path) =l= d(ppi);
pathDir2(back(ppi,path)) .. sigma(path) =l= 1 - d(ppi);


Equations
	useNode1(edge,node,node)	"edge must have two relevant nodes"	
	useNode2(edge,node,node)	"ditto"
	useNode3(node)	"node must have a relevant edge";
	
useNode1(enode(edge,node,node1)) .. x(edge) =l= y(node);
useNode2(enode(edge,node,node1)) .. x(edge) =l= y(node1);
useNode3(node) .. y(node) =l= sum( (edge,node1)$(enode(edge,node,node1) or enode(edge,node1,node)), x(edge));


Equation countActiveNodes	"count active nodes";
countActiveNodes .. nodeCount =e= sum(node, y(node));

Equation countActivePaths	"count active paths";
countActivePaths .. pathCount =e= sum(path, sigma(path));
