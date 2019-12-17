#! /usr/bin/env python
import sys
import re
import argparse
import glob
import os
from collections import defaultdict
from Bio import SeqIO

"""
Author: Conor Meehan 17/08/2018

This script is used to combine files of gene alignments into a concatenated alignment (assumes all species/isolates in all gene alignments)

Input:
A folder of fasta gene alignmentds where the sequence names is the species/isolate name (identical between files)
NOTE: assumes all end in .fasta

Output:
A concatenated alignment of all genes, using the species name as the sequence header (pulled from the first input file)
A file that outlines the gene name and the start/stop positions of the gene in the concatenated alignment
e.g. STENO_1, 1-1345

NOTE: uses biopython
NOTE: Assumes all genes have all species/isolates. If not, will mess up the alignment

Usage:
python genesToConcatAlign.py --folder folderName

"""
#parse the inputs
parser = argparse.ArgumentParser(description= 'This script is used to combine files of gene alignments into a concatenated alignment (assumes all species/isolates in all gene alignments)')
parser.add_argument('--folder', required=True, help='Folder containing the gene alignment fasta formatted files')
args = parser.parse_args()

#create the save file
try:
	saveAlign=open("concatAlign.fasta",'w')
except IOError:
	print 'no room for save file'
	sys.exit()
try:
	saveGene=open("concatAlignGenePartition.txt",'w')
except IOError:
	print 'no room for save file'
	sys.exit()
#Go through the folder file by file, extracting the gene sequences into the relevant portions of the dictionary
files=glob.glob(args.folder+'/*.fasta')

align=defaultdict(str)
isolates=[]
genes={}
geneOrder=[]
position=0
firstFile=0
for file in files:
	#get the gene name
	geneName=re.sub(".*\/","",file)
	geneName=re.sub(".fasta","",geneName)
	length=0
	first=0
	
	#parse the fasta file, isolate by isolate and add entries to the dictionary
	for record in SeqIO.parse(file, "fasta"):
		if firstFile==0:
			isolates.append(record.id)
		if first==0: #first entry, get the length of the gene
			first=1
			length=len(record.seq)
		align[record.id]+=(str(record.seq))
	genes[geneName]=str(position+1)+"-"+str(position+length)
	geneOrder.append(geneName)
	position=position+length
	firstFile=1	
	print geneName+" parsed"
#print the alignment and gene partition information to files
print "Writing out alignment and gene partition file"
for isolate in isolates:
	saveAlign.write(">"+isolate+"\n"+align[isolate]+"\n")

for gene in geneOrder:
	saveGene.write(gene+" = "+genes[gene]+"\n")

saveAlign.close()
saveGene.close()
sys.exit()