#!/bin/bash
# Summarize a collection of Steiner forest samples
# Summarize the subset for forest for each source separately
echo Summarize cdc14
python summarize_sif.py --indir ../results/cdc14_hog1_pde2_sources_031717 --pattern cdc14*optimalForest.sif --prizefile ../data/cdc14_prizes.txt --outfile ../results/cdc14_hog1_pde2_sources_031717/cdc14_summary

echo Summarize hog1
python summarize_sif.py --indir ../results/cdc14_hog1_pde2_sources_031717 --pattern hog1*optimalForest.sif --prizefile ../data/hog1_prizes.txt --outfile ../results/cdc14_hog1_pde2_sources_031717/hog1_summary

echo Summarize pde2
python summarize_sif.py --indir ../results/cdc14_hog1_pde2_sources_031717 --pattern pde2*optimalForest.sif --prizefile ../data/pde2_prizes.txt --outfile ../results/cdc14_hog1_pde2_sources_031717/pde2_summary
