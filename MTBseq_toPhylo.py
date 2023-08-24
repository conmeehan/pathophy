#! /usr/bin/env python
import sys
import re
import collections
import argparse
import glob
import os.path

"""
Author: Conor Meehan 04/04/2018

This script takes the output of MTBseq and prepares the SNP alignment and invariant site count for use in RAxML-NG
Input:
Amend folder with *_amended_u95_phylo_w12.plainIDs.fasta and *_amended_u95_phylo_w12.tab files inside (can be any starting text but must end with these suffixes)
The reference genome used for the assembly in MTBseq
Any characters to be ignored in terms of invariance. E.g. if the skipping characters are - and N a site such as --A-N will be removed as only 1 useable character is left

The script goes through the tab file, finds all positions that are SNPs and counts all other (invariant) positions in the genome
The fasta file output from MTBseq may have some invariant positions (e.g. the same SNP in all samples compared to reference will be a SNP in the dataset, but invariant in the alignment)
These sites are removed from the SNP fasta file and the invariant counts are modified accordingly.

All these files are then placed into a suggested RAxML-NG command.
This command will run 20 starting trees, select the best one, then run 100 bootstraps and place on the best tree. Site-repate is turned on for optimal SNP processing and threads are set to 3, the often optimal for MTB alignments.
If run, the best tree with BS support on the branches will be in <projectName>.raxml.support

Output:
A SNP file for use in phylogenetics (e.g. RAxML-NG). This file is named <projectName>_SNPalignment.fasta 
A count of the invariant sites in RAxML-NG format (i.e. countA/countC/countG/countT). This file is named <projectName>_invariantSiteCount.txt
A list of the site sin the SNP file that were removed in the 2nd step of the script. This file is named <projectName>_invariantSitesRemovedFromSNPfile.txt
The suggested RAxML-NG command. This is placed in file <project>_suggestedRAxML-NG.sh
(Project name is derived from the name before the first _ in the tab file above)

Usage:
python  MTBseq_toPhylo.py --amend <apathToAmendFolder> --ref <pathToReferenceSequence> --skip <charactersNotToBeCountedAsVariants>
If given with no options, assumes that the Amend folder is in the current folder, the reference sequence is in the current folder and named M._tuberculosis_H37Rv_2015-11-13.fasta (default for MTBseq) and the skipped characters are N- (case insensitive)

"""
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



#parse the inputs
parser = argparse.ArgumentParser()
parser.add_argument('--amend', required=False, default="Amend", help='folder containing fasta and tab file (default is Amend in current folder)')
parser.add_argument('--ref', required=False, default="M._tuberculosis_H37Rv_2015-11-13.fasta", help='Reference genome fasta file (default is M._tuberculosis_H37Rv_2015-11-13.fasta in current folder)')
parser.add_argument('--skip', required=False, default="N-", help='Characters not to be counted as variants in SNP file (default is N-)')

args = parser.parse_args()


#read in files and characters to be skipped
try:
	tabFile=glob.glob(os.path.join(args.amend,"*_amended_u95_phylo_w12.tab"))[0]
	tabF=open(tabFile, 'rU')
except IOError:
	print "\n tab variant file not found in directory."
	sys.exit()
try:
	refF=open(args.ref, 'rU')
except IOError:
	print "\n reference fasta file not found in directory."
	sys.exit()
try:
	fastaFile=glob.glob(os.path.join(args.amend,"*_amended_u95_phylo_w12.plainIDs.fasta"))[0]
	fastaF=open(fastaFile, 'rU')
except IOError:
	print "\n fasta file not found in directory."
	sys.exit()
try:
	remCharsStr=args.skip
	remChars=list(remCharsStr)
except IOError:
	print "\n Characters to be skipped, not supplied."
	sys.exit()

#get the project name from the tab file
#Do this by removing everything before the last / and everything after the first _
project=re.sub(".*\/","",tabFile)
project=re.sub("_.*","",project)

#create save files
try:
	saveFasta = open(project+"_SNPalignment.fasta",'w') 
except IOError:
	print "\n no room for save file"
	sys.exit()	
try:
	saveRemPos = open(project+"_invariantSitesRemovedFromSNPfile.txt",'w') 
except IOError:
	print "\n no room for save file"
	sys.exit()	
try:
	saveInvCount = open(project+"_invariantSiteCount.txt",'w') 
except IOError:
	print "\n no room for save file"
	sys.exit()	
try:
	saveScript = open(project+"_suggestedRAxML-NG.sh",'w') 
except IOError:
	print "\n no room for save file"
	sys.exit()

print "Reading in reference"
#place the reference sequence into a string and then a list
refStr=''
while 1:
	line=refF.readline()
	if not line:
		break
	line=line.rstrip()	
	if re.match(">",line):
		continue
	else:
		refStr+=line
refF.close()
refStr=refStr.upper()
ref=list(refStr)

#go through the tab file and for each position, change that position to a - (position is tab position-1 for mapping to list)
header=0
print "Removing variants based on tab file"
while 1:
	line=tabF.readline()
	if not line:
		break
	if re.match("#Position",line):
		header=1
		continue
	if header==1:
		line=line.rstrip()
		sections=line.split("\t")
		#remove that position from the reference file
		position=int(sections[0])
		ref[position-1]="-"
tabF.close()
#count the A, C, G, T in the reference (all upper from command above) and save to a dictionary
countsRef=collections.Counter(ref)

#read in the SNP alignment file. Allows for sequences to be multiple line per sample (needed?)
print "Modifying SNP file (if needed)"
align=[]
seq=''
pos=0
name=''
while 1:
	line=fastaF.readline()
	if not line:
		break
	line=line.rstrip()	
	
	if re.match(">",line):
		if pos==1:
			align.append([name,seq])
			seq=''
		pos=1
		name=line
	else:
		seq+=line
#do for last entry
align.append([name,seq])

fastaF.close()

#assumed all sequences are the same length so get the length of the first entry
al=len(align[0][1])
#go through alignment column by column and record if a site changes away from that the character in sequence 1. 
#If so, skip, if none do, record as invariant
remSites=[]
remSitesChars=[]
for i in range(al):
	#go through all samples and get the character in this column
	charsInColumn=[]
	for entry in align:
		character=entry[1][i].upper()
		if character.upper() in remChars: #skip it as a variant site if it in the list of sites to ignore
			continue
		charsInColumn.append(character) #if not skipped, add to character list
		
	uniqueChars=unique(charsInColumn) #once finished all sites, get a unique list of the characters
	if len(uniqueChars)<2: #only 1 character found at this position, so remove it
		remSitesChars.append(uniqueChars[0])
		remSites.append(i)
#convert to set for quicker lookup
srs=set(remSites)	

#go sequence by sequence and remove all the sites listed in remSites
nalign=[]
for entry in align:
	newseq=''
	oldseq=entry[1]
	for pos in range(len(oldseq)):
		if pos not in srs:
			newseq+=oldseq[pos]
	nalign.append([entry[0],newseq])

#output the new sequences to file
for entry in nalign:
	saveFasta.write(entry[0]+"\n"+entry[1]+"\n")		
saveFasta.close()

#save all the positions removed
#increase the number by 1 for alignment positions (ie start at 1, not 0)
for i in range(len(remSites)):
	pos=remSites[i]+1
	saveRemPos.write(str(pos)+"\t"+",".join(remSitesChars[i])+"\n")
saveRemPos.close()

#Count up those removed positions and add to the A/C/G/C counts from reference sequence, ave as strings
print "Counting invariant sites"
countsSNP=collections.Counter("".join(remSitesChars))
counts=[str(countsRef['A']+countsSNP['A']), str(countsRef['C']+countsSNP['C']), str(countsRef['G']+countsSNP['G']), str(countsRef['T']+countsSNP['T'])]
saveInvCount.write("/".join(counts))
saveInvCount.close()

#save the suggested RAxML-NG command to file
print "Writing suggested RAxML-NG command to file"
saveScript.write("raxml-ng --model GTR+G+ASC_STAM{"+"/".join(counts)+"} --site-repeat on -all --msa "+project+"_SNPalignment.fasta --prefix "+project+" --threads 3")
saveScript.close()

print 'Done'
sys.exit()	




