# pathophy
Scripts for aiding in pathogen phylogenetics analysis

phylogeneticInclusion python script:<br/>

This script takes a tree file and a cluster file and extends clusters to include all tips that are from the common ancestor of all members of a cluster<br/>
The cluster file is one cluster per line, separated by tabs<br/>
The tree file is in newick format<br/>
NOTE: requires dendropy to be installed<br/>
usage<br/>
python phylogeneticInclusion.py -t treeFile -c clusterFile<br/>
