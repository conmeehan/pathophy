#! /usr/bin/env python
import sys
import re
import argparse

"""
Author: Conor Meehan 21/07/2020

This script takes a BLAST output in tabular format and the input query fasta file and filters out all entries under a certain similarity percentage and percentage of length cut-off

Input:
BLAST table (outfmt 6 from blast+)
The original query fasta file
A percentage similarity/identity above which to retain hits (is >=, not just >)
A percentage length of the query sequence above which to retain a hit (is >=, not just >)

#NOTE: assumes proper formatting of query names, so there is an exact match between query name and what is in the BLAST table
#NOTE:fasta file should be a single line per sequence, not multi line
Output:
A subset of results based on the filtering above


Usage:
python blastHitsFilterIdLength.py --blast <blast file> --query <query fasta file> --id <% id cut-off> --length <% query length cut-off>

"""
#parse the inputs
parser = argparse.ArgumentParser(description= 'This script takes a BLAST output in tabular format and the input query fasta file and filters out all entries under a certain similarity percentage and percentage of length cut-off')
parser.add_argument('--blast', required=True, help='BLAST table file (outfmt 6 from blast+)')
parser.add_argument('--query', required=True, help='The original query fasta file')
parser.add_argument('--id', required=True, help='A percentage similarity/identity above which to retain hits (is >=, not just >)')
parser.add_argument('--length', required=True, help='A percentage length of the query sequence above which to retain a hit (is >=, not just >)')

args = parser.parse_args()


#read in the files and get the %id and %length
try:
	blastF = open(args.blast, 'r') 
except IOError:
	print("\nblast file not found.")
	sys.exit()

try:
	queryF = open(args.query, 'r') 
except IOError:
	print("\nquery fatsa file not found.")
	sys.exit()
try:
	idco = float(args.id) 
except IOError:
	print("\n% id cut-off not found.")
	sys.exit()
	
try:
	lenco = float(args.length) 
except IOError:
	print("\n% length cut-off not found.")
	sys.exit()	
	
#create a save file
try:
	save = open("blastIdLengthFiltered.txt",'w') 
except IndexError:
	print("\n no room for save file")
	sys.exit()	

#Get the length of each query sequence
lengths={}
while 1:
	line=queryF.readline()
	if not line:
		break
	
	if line.startswith(">"):
		name=re.sub(">","",line)
		name=name.rstrip()
	else:
		seq=line.rstrip()
		querylen=len(seq)
		lengths[name]=querylen
queryF.close()			

#read the blast file in and for each line do the processing
while 1:
	line=blastF.readline()
	if not line:
		break
	line=line.rstrip()
	sections =line.split("\t")
	#first filter on the ID percentage
	if float(sections[2])<idco:
		continue
	
	#get the length of the query sequence for this hit
	qlen=lengths[sections[0]]
	#get the length covered by the hit (end query - start query). Add one so it includes that start position
	hitlen=abs(float(int(sections[7])-int(sections[6])+1))

	#get the length as a percentage of the total query length
	pclen=(hitlen/qlen)*100
	if(pclen<lenco):
		continue
	
	#any hits that are at this stage have passed so save to file
	save.write(line+"\n")
	
				
blastF.close()
save.close()
sys.exit()