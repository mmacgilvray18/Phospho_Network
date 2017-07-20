# Prize-collecting Steiner forest
Network analysis of the yeast salt stress phosphorylation response with the
[Omics Integrator](https://github.com/fraenkel-lab/OmicsIntegrator) software.
The Omics Integrator site provides instructions for installing Omics Integrator
and the msgsteiner dependency.

## cluster
The scripts in this subdirectory are provided for transparency and are intended
to show in detail how Omics Integrator was run.  They are not runnable on a new
system without minor modifications.  The scripts assume that Omics Integrator
and msgsteiner have been installed.  They use the
[HTCondor](https://research.cs.wisc.edu/htcondor/) queuing system, but this can
be substituted by replacing the `condor_submit` calls in
`submit_wrapper_031717.sh` with the appropriate submit commands for a different
queuing system.

The overall workflow is to submit a series of HTCondor jobs with
`submit_wrapper_031717.sh`. These use the job specifications in
`submit_PCSF_031717.sub` and execute instances of `run_PCSF_031717.sh`.  After
all the jobs complete, the  scripts `summarize_forests_031717.sh` and
`summarize_single_source_forests_031717.sh` are used to call `summarize_sif.py`
and summarize the Steiner forest ensembles.

## data
The data subdirectory contains the input files for Omics Integrator.
`background_network.txt` contains the same network used for the ILP network
analysis with different formatting.  The Prize-collecting Steiner forest was
run independently for each source so each of the three sources has its own
sources file and prizes files.

## results
The `cdc14_hog1_pde2_sources_031717` subdirectory contains the output from the
Omics Integrator runs and summarization scripts.  There are summaries of all
Steiner forests in the ensemble as well as the subsets of forests that are
specific to a single target.  The [Cytoscape](http://www.cytoscape.org/) file
`cdc14_hog1_pde2_summary_union.cys` visualizes these summaries and was created
with Cytoscape version 3.2.  `*_edgeAnnotation.txt` files contain the frequency
of edges in the ensemble.  `*nodeAnnotation.txt` contains the node ensemble
frequency.  `_size.txt` shows the sizes (number of nodes) of the individual
Steiner forests in the ensemble.  `_union.sif` and `_union.tsv` contain the
ensemble edges in different formats.
`cdc14_hog1_pde2_summary_union_NP_sumod_added.txt` extends
`cdc14_hog1_pde2_summary_union.sif` by adding the no phenotype submodules into
the ensemble network, and this is the version of the Steiner forest network used
to compare with the ILP network.

## References
Tuncbag N, Gosline SJ, Kedaigle A, Soltis AR, Gitter A, Fraenkel E. (2016)
Network-based interpretation of diverse high-throughput datasets through the
Omics Integrator software package.
PLoS Comput Biol. 12(4):e1004879.
[doi:10.1371/journal.pcbi.1004879](https://doi.org/10.1371/journal.pcbi.1004879).
