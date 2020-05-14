#! /usr/bin/env python
import sys
import re
from collections import defaultdict

#script to cluster sequences based on a pairwise distance file
#Script takes the output from pairwiseDist.py and a cut-off
#Output is 
#A file listing all sequences that are within the cut-off from another sequence
#A file where each line is a list of sequences within cut-off distance of each other (tight cluster)
#A file where each line is a list of sequences where at each sequence is within the cut-off from at least 1 other (e.g. if cutoff is 5 then seq1 and seq2 could be 5 apart and seq3 could be 10 from seq2 and 3 from seq1 but all on same line) (loose cluster)
#NOTE: cut-off is used as <= so if cutoff is 5 then those a distance of 5 will be counted as clustered
#NOTE2: NA can be an entry (e.g. for a tip distance that is not supported by bootstrap or unequal lengths of patterns) so if it is, count as not within cutoff
#usage
#python pairwiseCluster.py pairwiseDistFile cut-off 

#subroutine to remove duplicates from an array
def unique(array):
	uniq = list(set(array))

	seen = {}
	for item in array:
   	 seen[item] = seen.get(item, 0) + 1

	uniq = seen.keys()

	seen = {}
	uniq = []
	for item in array:
		count = seen.get(item, 0)
		if count == 0:
			uniq.append(item)
		seen[item] = count + 1
	return uniq	

#read in file
try:
	pwMxF=open(sys.argv[1], 'r')
except IndexError:
	print("\n pairwise distance matrix not supplied.")
	sys.exit()
except IOError:
	print("\n pairwise distance matrix not found in directory.")
	sys.exit()
#get the cut-off
try:
	cutoff = float(sys.argv[2]) 
except IndexError:
	print("\ncluster distance cut-off not supplied.")
	sys.exit()	
except ValueError:
	print("\ncluster distance cut-off is not a number.")
	sys.exit()

#create save files
try:
	saveList=open("samplesInCluster_"+sys.argv[2]+"cutoff.txt",'w')
except IOError:
	print('no room for save file')
	sys.exit()
try:
	saveTightCluster=open("tightClusters_"+sys.argv[2]+"cutoff.txt",'w')
except IOError:
	print('no room for save file')
	sys.exit()
try:
	saveLooseCluster=open("looseClusters_"+sys.argv[2]+"cutoff.txt",'w')
except IOError:
	print('no room for save file')
	sys.exit()

#read in the pairwise matrix and for each sequence keep a dictionary of what others they are in the cluster with
#keep a list of all that are in a cluster
first=0
names=0
clusters=defaultdict(list)
inCluster=[]
while 1:
	s=pwMxF.readline()
	if not s:
		break
	s=s.rstrip()
	sections=s.split("\t")
	if first==0:#header line
		names=sections
		names.pop(0)
		first=1
		continue
	dists=list(sections)
	seqName=dists.pop(0)
	seqPos=names.index(seqName)
	for i in range(len(dists)):
		if i==seqPos: #skip the self hits
			continue
		elif dists[i]=="NA": #is a distance that is not applicable so skip
			continue
		elif float(dists[i])<=cutoff:
			clusters[seqName].append(names[i])
			inCluster.append(names[i])
pwMxF.close()

#create a unique list of those in clusters and print to file
inCluster=unique(inCluster)		
inCluster.sort()
saveList.write("\n".join(inCluster)+"\n")
saveList.close()

#turn the dictionary into a unique sorted list of sorted strings (tight clusters) and sorted lists(loose clusters)
#the strings are made unique and these are the tight clusters
stringClust=[]
for seq in clusters:
	temp=list(clusters[seq])
	temp.append(seq)
	temp.sort()
	stringClust.append("\t".join(temp))
stringClust=unique(stringClust)
stringClust.sort()

#save to file the clusters
for sc in stringClust:
	saveTightCluster.write(sc+"\n")
saveTightCluster.close()

#loose groups are formed by going sample by sample and adding all the other samples they are in a cluster with to a loose cluster. If none are in an existing cluster, a new one is started
looseClusters=[]
for sample in inCluster:
	newClust=list(clusters[sample])
	newClust.append(sample)
	sNC=frozenset(newClust) #create a set version of the list
	if looseClusters==[]: #no existing clusters
		looseClusters.append(sNC)
	else:
		#check all samples in this cluster against samples already in the loose clusters list. If any are there, add all to that loose cluster
		found=0
		for lcNum in range(len(looseClusters)):
			slc=looseClusters[lcNum]
			if sNC.intersection(slc):
				newLC=sNC.union(slc)
				looseClusters[lcNum]=newLC
				found=1
				break
		if found==0:
			looseClusters.append(sNC)
#there is a possibility that some overlaps were missed due to order of adding them together so do one extra round of additions
looseClusters_final=[]
skip=[]
for lcNum in range(len(looseClusters)):
	 if lcNum in skip:
	 	continue
	 else:
	 	found=0
	 	lc1=looseClusters[lcNum]
	 	for lcNum2 in range(lcNum+1,len(looseClusters)):
	 		lc2=looseClusters[lcNum2]
	 		if lc1.intersection(lc2):
	 			newLC=lc1.union(lc2)
				looseClusters_final.append(newLC)	
				found=1
				skip.append(lcNum2)
				break
		if found==0:
			looseClusters_final.append(lc1)			

#create a unique list for each loose cluster and save as strings
looseClustersStrings=[]
for slc in looseClusters_final:
	lc=unique(list(slc))
	lc.sort()
	looseClustersStrings.append("\t".join(lc))
looseClustersStrings.sort()	
			

#save to file the clusters
for lcs in looseClustersStrings:
	saveLooseCluster.write(lcs+"\n")
saveLooseCluster.close()
				