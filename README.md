# pathophy
Scripts for aiding in pathogen phylogenetics analysis

phylogeneticInclusion python script:
This script takes a tree file and a cluster file and extends clusters to include all tips that are from the common ancestor of all members of a cluster
The cluster file is one cluster per line, separated by tabs
The tree file is in newick format
NOTE: requires dendropy to be installed
usage
python phylogeneticInclusion.py -t treeFile -c clusterFile
