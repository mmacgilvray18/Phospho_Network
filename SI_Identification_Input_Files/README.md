# Shared interactor input files

### `Annotation.csv`
This file contains two columns "systematic_name_dash_removed" and "systematic
name". The only difference between these columns is that dashes for the gene
systematic names have been removed. This file is used to change the systematic
name for all proteins in submodules to the non-dash version. This was necessary
based on how we identify shared interactors. Shared interactors are proteins
that are enriched for interactions with submodule constituent proteins that
exist in the background network, compared to chance. Since the background
network does NOT contain dashes; we would miss interactions for proteins that
contain a dash in their systematic name. After shared interactors are identified
using the non-dash systematic name, we change the name back to the appropriate
yeast systematic name (with dashes), using this file.  

### `Background_Network.csv`
This file contains 25,682 protein-protein interactions in yeast that were
curated in Chasman and Ho *et al* (2014). These interactions are in SIF file
format, thus interactions are directed from the protein in the "Protein 1"
column towards the protein in the "Protein 2" column, unless otherwise noted.
This file contains 51,364 rows, representing all interactions being included
twice, but in their reverse orientation (ie, A:B is also listed as B:A). This
was done because when we search for a given submodule protein's interactors, we
only search in the "Protein 1 column" and don't want to miss any interactions.
Interaction types (ie, kinase-substrate, etc) are listed in the "Interaction"
column. "ppi" stands for "Protein Protein Interaction", which are interactions
with no known directionality. Within the "Interaction" column, "both" indicates
that an interaction is bidirectional. "Reverse" on the end of an interaction
indicates a reverse interaction pair, as previously described. All directed
interactions, in a single direction, are denoted by a "1" in the "Directed"
column, whereas bidirectional and interactions lacking directionality are
denoted by "0".

### `Number_Interactions_Each_Protein.csv`
This file contains the number of unique proteins each protein in the background
network interacts with in the background network. This number is displayed in
the "Total" column. The number of proteins each protein interacts is critical
for identifying enriched shared interactors according to the hypergeometric
test.

## References
Chasman D, Ho YH, Berry DB, et al. (2014) Pathway connectivity and signaling
coordination in the yeast stress-activated signaling network. Mol Syst Biol 10: 759.
