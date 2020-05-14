#! /usr/bin/env python
import sys
import re

#script to take an alignment and get pairwise distance between all sequences as a matrix, with allowed exclusions
#Script also takes a list of characters not to count as divergent if encountered
#e.g. if the excluded characters are N- then if a comparison has a N or - this will not be counted as different
#NOTE: is not case sensitive so N and n are the same thing
#usage
#python pairwiseDist.py fastaAlignment excludedCharacters

def pairwiseDist(seq1,seq2):
	total=0
	for pos in range(len(seq1)):
		char1=seq1[pos].upper()
		char2=seq2[pos].upper()
		
		skip=0
		for skipChar in skipChars:
			if char1==skipChar or char2==skipChar:
				skip=1
		if skip==1:
			continue
		
		if char1!=char2:
			total+=1
	return total		


#get file
try:
	fastaF=open(sys.argv[1], 'r')
except IndexError:
	print "\n fasta file not supplied."
	sys.exit()
except IOError:
	print "\n fasta file not found in directory."
	sys.exit()
#get the characters to not be counted as variants
try:
	skipCharStr = sys.argv[2] 
	skipChars=list(skipCharStr)
except IndexError:
	print "\ncharacters to be skipped as variants missing so assuming none are to be counted"
	skipChars=['']


#create save file 
try:
	save=open("pairwise_distances.txt",'w')
except IOError:
	print 'no room for save file'
	sys.exit()


#read all sequences into an array
align=[]
names=[]
seq=''
pos=0
name=''
while 1:
	s=fastaF.readline()
	
	if not s:
		break
	s=s.rstrip()
	if re.match(">",s):
		if pos==1:
			align.append(seq)
			names.append(re.sub(">","",name))
			seq=''
		pos=1
		name=s
	else:
		seq+=s

#do for last entry
align.append(seq)
names.append(re.sub(">","",name))

#get the pairwise distances as a matrix
dists=[[0 for col in range(len(names))] for row in range(len(names))]
for seq1Num in range(len(align)):
	print "processing "+names[seq1Num]
	for seq2Num in range(seq1Num+1,len(align)):
		dist=pairwiseDist(align[seq1Num],align[seq2Num])
		dists[seq1Num][seq2Num]=dist
		dists[seq2Num][seq1Num]=dist

#print the matrix to file
save.write("\t"+"\t".join(names)+"\n")
for rowNum in range(len(dists)):
	save.write(names[rowNum])
	for entry in dists[rowNum]:
		save.write("\t"+str(entry))
	save.write("\n")
save.close()
sys.exit()		