#! /usr/bin/env python
import sys
import re
import dendropy
import copy
from dendropy.calculate import treemeasure
import argparse

"""
Author: Conor Meehan 07/09/17

This script takes a tree file and a cluster file and extends clusters to include all tips that are from the common ancestor of all members of a cluster
The cluster file is one cluster per line, separated by tabs
The tree file is in newick format
NOTE: requires dendropy to be installed
usage
python phylogeneticInclusion.py -t treeFile -c clusterFile

"""

#parse the inputs
parser = argparse.ArgumentParser()
parser.add_argument('-t', required=True, help='Newick formatted tree file')
parser.add_argument('-c', required=True, help='Cluster file with each cluster as a tab separated list of taxon names. 1 cluster per line.')
args = parser.parse_args()


#get the tree name and cluster file
try:
	tree=args.t
except IOError:
	print "\n tree file not found in directory."
	sys.exit()
try:
	clusterF=open(args.c, 'r')
except IOError:
	print "\n cluster file not found in directory."
	sys.exit()

#create save file 
try:
	save=open("ClustersExtendedOnTree.txt",'w')
except IOError:
	print 'no room for save file'
	sys.exit()	


#read in the tree, keeping underscores in the names
tree = dendropy.Tree.get(path=tree,schema='newick',preserve_underscores=True)
pdm = treemeasure.PatristicDistanceMatrix(tree)
taxa=tree.taxon_namespace

#go cluster by cluster, get the mrca of them all, then output all the tip child nodes of that MRCA as the new cluster
clusters=[]
while 1:
	s=clusterF.readline()
	if not s:
		break
	s=s.rstrip()
	sections=s.split("\t")
	mrca = tree.mrca(taxon_labels=sections)
	newClusterNodes=mrca.leaf_nodes()
	newCluster=[]
	for leaf in newClusterNodes:
		newCluster.append(re.sub("\'","",str(leaf.taxon)))
	newCluster.sort()	
	clusters.append(set(newCluster))

#there is a possibility that some overlaps were missed due to extension encompassing an existing cluster so do one extra round of additions
looseClusters_final=[]
skip=[]
for lcNum in range(len(clusters)):
	 if lcNum in skip:
	 	continue
	 else:
	 	lc1=clusters[lcNum]
	 	for lcNum2 in range(lcNum+1,len(clusters)):
	 		lc2=clusters[lcNum2]
	 		if lc1.intersection(lc2):
	 			lc1=lc1.union(lc2)
				skip.append(lcNum2)
		if lc1 not in looseClusters_final:	#if case identical cluster already made it in, skip it
			looseClusters_final.append(lc1)	
#create a unique list for each loose cluster and save as strings
looseClustersStrings=[]
for slc in looseClusters_final:
	lc = list(set(slc))
	lc.sort()
	looseClustersStrings.append("\t".join(lc))
looseClustersStrings.sort()	
			

#save to file the clusters
for lcs in looseClustersStrings:
	save.write(lcs+"\n")
save.close()

sys.exit()
