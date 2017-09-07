# pathophy
Scripts for aiding in pathogen phylogenetics analysis

Phylogenetic Inclusion python script:<br/>
This script is used to take a transmission cluster (e.g. defined by SNP cut-offs) and extend it to include all other related strains that share the common ancestor in a phylogenetic tree.<br/>
Inputs:<br/>

A cluster file: one cluster per line, separated by tabs<br/>
A tree file in newick format<br/>
NOTE: requires dendropy to be installed<br/>
usage<br/>
python phylogeneticInclusion.py -t treeFile -c clusterFile<br/>
