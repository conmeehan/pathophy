#! /usr/bin/env python
import sys
import re
import argparse
import glob
import os
from collections import defaultdict
from Bio import SeqIO

"""
Author: Conor Meehan 16/08/2018

This script is used to combine files of individual gene sequences, per isolate, into gene-specific alignment files

Input:
A folder of isolate files where the isolate name is everything before the first _.
NOTE: assumes all end in .fasta
Inside each file is a set of fasta formatted sequences where the gene name is xxx_xxx (i.e. based around 1st _ only)
e.g. >STENO_1_6, the gene name is STENO_1

Output:
In a folder named geneAlignments:
Per gene files where the file name is the gene name (e.g. STENO_1.fasta)
Inside each file the sequences are named as isolate_geneName (e.g. 677_STENO_1)

NOTE: uses biopython
NOTE: assumes only 1 entry per gene in each isolate. Will warn if double and only keep last one

Usage:
python isolateFilesToGeneFiles.py --folder folderName

"""

#parse the inputs
parser = argparse.ArgumentParser(description= 'This script is used to combine files of individual gene sequences, per isolate, into gene-specific alignment files')
parser.add_argument('--folder', required=True, help='Folder containing the isolate fasta formatted files')
args = parser.parse_args()

#Go through the folder file by file, extracting the gene sequences into the relevant portions of the dictionary
files=glob.glob(args.folder+'/*.fasta')
if not os.path.exists("geneAlignments"):
    os.makedirs("geneAlignments")

genes=defaultdict(dict) #for each gene name there is a dictionary where the key is the sample name and the value is the sequence
isolates=[] #list of isolate names so that always printed out in the same order

for file in files:
	try:
		F=open(file, 'rU')
	except IOError:
		print "\n "+file+" not found in folder."
		sys.exit()
	#get the sample name
	sampleName=re.sub(".*\/","",file)
	sampleName=re.sub(".fasta","",sampleName)
	sampleName=re.sub("_.*","",sampleName)
	isolates.append(sampleName)
	
	#parse the fasta file, gene by gene and add entries to the dictionary
	for record in SeqIO.parse(file, "fasta"):
		#keep only data to the left and right of 1st _
		sections=record.id.split("_")
		name=sections[0]+"_"+sections[1]
		if sampleName in genes[name]:
			print "double entry for gene "+name+" in isolate "+sampleName+". Only retaining the final entry in isolate file."
		genes[name][sampleName]=record.seq

#print a per gene file of all the associated sequences
for gene in genes:
	try:
		save=open("geneAlignments/"+gene+".fasta",'w')
	except IOError:
		print 'no room for save file'
		sys.exit()
	
	for isolate in genes[gene]:
		save.write(">"+isolate+"\n"+str(genes[gene][isolate])+"\n")
	
	save.close()