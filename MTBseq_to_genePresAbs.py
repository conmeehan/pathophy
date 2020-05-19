#! /usr/bin/env python
import sys
import re
import argparse
import glob
from collections import defaultdict

"""
Author: Conor Meehan 18/02/2019

Script to find genes that are not present in samples from an MTBseq run output (Position Table files)


Input:
Folder containing Position Tables from MTBseq (default is Position_Tables)
The genes file from the reference used for MTBseq (can be found in the var/ref/ folder of the MTBseq installation)
The minimum proportion of positions absent to be counted as an absent gene (default is 0.95)
The minimum number of reads to consider a position to be covered (default is 8, same as MTBseq)

Output:
A presence/absence map of the genes in a matrix. Names are the position table names with everythign after the first _ removed (similar to MTBseq)
The samples are listed in the first column; the genes are listed in the first row. There is a 1 for a present gene and a 0 for an absent gene

Usage:
python MTBseq_to_genePresAbs.py --folder Position_TablesFolder --genes genesFile --cutoff proportionAbsent --minreads coverageMinumum

"""

#parse the inputs
parser = argparse.ArgumentParser(description= 'Script to find genes that are not present in samples from an MTBseq run output (Position Table files)')
parser.add_argument('--folder', required=False, default="Position_Tables", help='Folder containing Position Tables from MTBseq (default is Position_Tables)')
parser.add_argument('--genes', required=True, help='The genes file from the reference used for MTBseq (can be found in the var/ref/ folder of the MTBseq installation)')
parser.add_argument('--cutoff', required=False, default="0.95", help='The minimum proportion of positions absent to be counted as an absent gene (default is 0.95)')
parser.add_argument('--minreads', required=False, default="8", help='The minimum number of reads to consider a position to be covered (default is 8, same as MTBseq)')

args = parser.parse_args()

co=float(args.cutoff)
minReads=int(args.minreads)

#read in the genes list and create a dictionary of name to start/end positions (stored as a list)
genes=defaultdict(list)
geneNames=[]
#open the genes map file
try:
	mapF=open(args.genes, 'r')
except IOError:
	print("\n gene mapping file not found in folder.")
	sys.exit()

while 1:
	line=mapF.readline()
	if not line:
		break
	line=line.rstrip()
	sections=line.split("\t")
	if sections[0]=="# ID":	#header line
		continue
	
	#some genes are in reverse so if the end position is before the start position, reverse them
	if int(sections[3])>int(sections[2]):
		start=int(sections[2])
		end=int(sections[3])
	else:
		start=int(sections[3])
		end=int(sections[2])
		
	genes[sections[0]]=[start,end]
	geneNames.append(sections[0])
mapF.close()

#go through each position table and create a list fo the genes it contains
#this is done by going from the first to last position of the gene, counting the number of positions that have no coverage and then dividing that by the total length
#if this number is greater than the cut-off, count as absent
tables=glob.glob(args.folder+'/*.tab')
presAbs={}
for table in tables:
	print("Processing "+table)
	presAbs[table]={}
	#open the table file
	try:
		tableF=open(table, 'r')
	except IOError:
		print(table+" file not found in folder.")
		sys.exit()
	#read the table in as all positions with a 1 if present and a 0 if coverage is absent
	#this is done by adding up all the reads in columns 3-14 and if this is below the minimum read count, this is marked as absent
	positions={}
	while 1:
		line=tableF.readline()
		if not line:
			break
		line=line.rstrip()
		sections=line.split("\t")
		if sections[0]=="#Pos":	#header line
			continue
		coverage=sections[3:15]
		totalReads=0
		for c in coverage:
			totalReads+=int(c)
		if totalReads < minReads:
			positions[int(sections[0])]=0
		else:
			positions[int(sections[0])]=1
	tableF.close()				
	
	
	#go through the genes list and count the positions that have a 0
	for gene in genes:
		noCov=0
		length=0
		for pos in range(genes[gene][0],genes[gene][1]+1):
			if positions[pos]==0:
				noCov+=1
			length+=1	
		#get the proportion that is not covered

		prop=float(noCov)/length
		if prop>=co:
			presAbs[table][gene]=0
		else:
			presAbs[table][gene]=1		



#create the output
print("Processing done, creating output")
#create the save file
try:
	save=open("positionTables_genepresAbs.txt",'w')
except IOError:
	print('no room for save file')
	sys.exit()


#save the gene names in the header
header="\t"+"\t".join(geneNames)+"\n"
save.write(header)
for entry in presAbs:
	#take the folder name and everything after the _ form the name
	name=re.sub(args.folder+"/","",entry)
	name=re.sub("_.*","",name)
	save.write(name)
	for gene in geneNames:
		save.write("\t"+str(presAbs[entry][gene]))
	save.write("\n")
save.close()
sys.exit()