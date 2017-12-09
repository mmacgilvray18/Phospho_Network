# PathLinker
This subdirectory contains network analysis of the yeast salt stress phosphorylation response with the [PathLinker](https://github.com/Murali-group/PathLinker) software.
PathLinker was run in the conda environment specified by `environment.yml`.

## cluster
The script in this subdirectory is provided for transparency and is intended to show in detail how PathLinker was run.
It is not runnable on a new system without minor modifications to the paths.
The scripts assume that PathLinker has been installed and the `path-linker` conda environment is available.

## data
The data subdirectory contains the input files for PathLinker.
`PathLinker_Background_Network.txt` contains the same network used for the ILP  and PCSF network analysis but converts undirected edges to a pair of directed edges.
PathLinker was run independently for each source so each of the three sources has its own input file.

## results
The `cdc14_hog1_pde2_sources_pl_120917` subdirectory contains the PathLinker output.
See the [PathLinker repository](https://github.com/Murali-group/PathLinker) for a description of the output file formats.

## References
Ritz A, Poirel CL, Tegge AN, Sharp N, Simmons K, Powell A, Kale SD, Murali TM.
(2016)
Pathways on demand: automated reconstruction of human signaling networks.
NPJ Syst Biol Appl. 2:16002.
[doi:10.1038/npjsba.2016.2](https://doi.org/10.1038/npjsba.2016.2).
