#! /usr/bin/env python
import sys
import re
import os
import glob
from collections import defaultdict
import argparse

"""
Author: Conor Meehan 31/05/17

This pipeline converts outputs from snippy into a whole genome alignment, a SNP alignment and an invariant count of positions for use in ascertainment bias correction in phylogenetics

Input:
A folder in which there are snippy output folders. In particular the following files must be in each
	snps.vcf	This is used to get all snps and indels
	snps.aligned.fa	This is used for the Ns and large areas of either deletions or unaligned regions (e.g. insertion sites; multiple mapping sites). These two latter cannot be distinguished within this script.
	these files can be named other things (e.g. not snps but another prefix, as done by the --prefix in snippy) but must be PREFIX.snps and PREFIX.aligned.fa with no other aligned.fa file
The fasta reference sequence used to create the snippy outputs (assumed same for all). E.g. this is the ref.fa file in the reference folder of a snippy output.
An interval masking file (optional; see below)
A set of characters not to count as variable (optional; see below)
A flag for whether to include the refrence sequence (optional)

Output:
Whole genome alignment: includes all SNPs, indels, Ns and large deletions/unaligned regions
SNP alignment: N's, SNPs and indels are included
Invariant site counts: A count of A, G, C and T that are invariant in this whole geome alignment. This allows for input to RAxML or BEAST (or others) along with SNP alignment for whole geome phylogenetics and ascertainment bias correction.
	There are two invariant site count files: the definitive counts of invariant sites (i.e. only A, C, G, T) and the count where sites that are invariant also if excluded characters are present (e.g. N characters in an invariant site).
A map of genome position to SNP alignment position. Insertions will be listed as all at the same position
All sequences are named per the snippy output folder name

NOTES:
VCF file parsing:
	4 of 5 types of changes (snp, mnp, ins, and del) are handled
	the vcf file is first passed to vcfallelicprimitives to create the simple input. 
	Thus, vcfallelicprimitives must be installed and in the path

INDEL handling:
	If an indel has a variable numbe rof positions in different sequences (e.g. 2 extra in one sequence and 3 in another), alignment is not undertaken, the two will be placed at the start

Interval file
	Intervals can be defined to replace sections of all sequences with N's (i.e. hard masking)
	Interval file should be start of interval<tab>end of interval per line. Alternatively, if only one position on a line, that single site will be masked.
	Assumed interval numbers are 'regular' numbers. i.e. start at 1, not 0
	Intervals are inclusive (start and end positions will be masked)
	Reference (1st sequence in output) is also masked (if included with flag)

Characters to skip for SNP alignment
	Certain characters such as ambiguous or missing data can be excluded as variants when creating the SNP alignment.
	For example: If a column is --A-N, RAxML and other programs will this is as invariant because N is not counted as a SNP. Thus, such sites should be listed as invariant and not included in the SNP alignment.
	The characters to be ignored in terms of invariance are listed at prompt.
	Thus if sites of characters are variable only with gaps but the user wishes them to be removed they shall place -
	Often for WGS invariant removal -nN is needed
	For full removal of potential problem characters, use nNoOxX?KMRYSWBVHD-
	NOTE: don't put the - first, the program will be confused. It is ok if it is by itself
	Removal characters are case sensitive so n and N are not the same
	If a column is all excluded characters this will also be removed. Thus if you have N- as the exclusion set and the site is --NN- then this site will be removed

Inclusion of reference sequence
	Flag is added to include (1) or exclude (0) the reference sequence form the alignment
	It is assumed that the reference is only 1 sequence.

usage
python snippyVcfFilesWithAmbigToAlignmentsWithPosMap.py -f folderContainingSnippyFolders -r referenceGenome -i intervalFile (optional) -c charactersToRemove (optional) -a addReferenceFlag(0/1) (optional)
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
parser.add_argument('--f', required=True, help='Folder containing snippy output folders')
parser.add_argument('--r', required=True, help='Reference genome file')
parser.add_argument('--i', required=False, default="NONE", help='Interval file (optional)')
parser.add_argument('--c', required=False, default="NONE", help='Characters to be excluded from counting as variable (optional)')
parser.add_argument('--a', required=False, default=0, help='1 to add reference genome (optional)')
args = parser.parse_args()


#read in files
try:
	snippyFolder = args.f 
except IndexError:
	print "\nsnippy vcf files folder not found."
	sys.exit()
try:
	refF=open(args.r, 'rU')
except IOError:
	print "\n reference fasta file not found in directory."
	sys.exit()

intervalF=args.i
if intervalF!="NONE":
	try:
		intF=open(intervalF, 'rU')
	except IOError:
		print "\n interval file not found in directory."
		sys.exit()	
	
#get the characters to not be counted as variants
remCharStr = args.c
if remCharStr.upper()=="NONE":
	remChars=['']
else:
	remChars=list(remCharStr)
#get the flag for whether the reference is to be included or not
refInc = args.a 
	
print "Reading in reference"
#place the reference sequence into a string and then array
refStr=''
while 1:
	s=refF.readline()
	if not s:
		break
	s=s.rstrip()	
	if re.match(">",s):
		continue
	else:
		refStr+=s	
refF.close()
ref=list(refStr)

print "Creating SNP modified genomes"
#go through the SNP vcf file for each genome and first simplify
#modify the reference as needed per genome to create a whole genome sequence for each entry
seqs=[]
names=[]

insertions=defaultdict(dict) #keep track of insertions to do all of them at the end
insPos=[] #keep track of all places there is an insertion
folders=os.listdir(snippyFolder)
os.chdir(snippyFolder)
#remove the mac finder file
if ".DS_Store" in folders:
	folders.pop(folders.index(".DS_Store"))

listPos=-1
prefixDict={} #keep a list of the prefix used for each folder
for Folder in folders:
	names.append(Folder)
	os.chdir(Folder)
	#get the prefix of the .vcf and .aligned.fa files
	prefix=re.sub(".aligned.fa","",glob.glob("*.aligned.fa")[0])
	prefixDict[Folder]=prefix
	os.system("vcfallelicprimitives "+prefix+".vcf >snps_simple.vcf")
	try:
		snippyFile = open("snps_simple.vcf") 
	except IOError:
		print "\nError reading "+file
		sys.exit()
	
	
	seqs.append(list(ref))
	listPos+=1
	insertions[listPos]={} #create a blank section for the list positions
	while 1:
		s=snippyFile.readline()
		if not s:
			break
		s=s.rstrip()
		if re.match("#",s):#on an info line so skip
			continue
		#have reached a change line so get the position and type
		sections=s.split("\t")
		pos=int(sections[1])-1
		type=sections[7] #the type is listed at the end of section 7
		type=re.sub(".*TYPE=","",type)
		
		#if it is a snp
		if type=="snp":
			seqs[listPos][pos]=sections[4].upper()
		elif type=="del": #is always simple deletion so retains the first position and deletes all the ones after
			refEntry=list(sections[3])
			for extraNum in range(len(refEntry)):
				seqs[listPos][pos+extraNum]="-"
			seqs[listPos][pos]=	sections[4].upper()
		elif type=="ins":
			insertions[listPos][pos]=sections[4].upper()
			insPos.append(pos)
		elif type=="complex":#currently cannot handle. Code is here as legacy to maybe handle in the future
			#refEntry=sections[3]
			#sampleEntry=sections[4]
			#if len(refEntry)<len(sampleEntry):#its an insertion
			#	insertions[listPos][pos]=sections[4].upper()
			#	insPos.append(pos)
			#elif len(refEntry)>len(sampleEntry):#its a deletion so change snp and put in -
			#	for extraNum in range(len(refEntry)):
			#		seqs[listPos][pos+extraNum]="-"
			#	seqs[listPos][pos]=	sections[4][0].upper()	
			#else: #unknown complex type
			#	print "complex type found but unknown how to handle. Entry is:\n"+s
			continue
		else: #not a snp, insertion, deletion or complex so unknown how to handle
			print "entry type is not a snp, insertion, deletion or complex so unknown how to handle. Entry is:\n"+s			
	snippyFile.close()
	os.chdir("..")
#place the large gaps and Ns onto the alignment
print "Adding ambiguous characters and large deletions/unaligned portions"
ambigAlign=[]
for Folder in folders:
	prefix=prefixDict[Folder]
	path=os.path.join(Folder,prefix+".aligned.fa")
	try:
		nFile = open(path) 
	except IOError:
		print "\nError reading "+prefix+".aligned.fa file from folder "+Folder
		sys.exit()
	#get the sequence
	seq=""
	while 1:
		s=nFile.readline()
		if not s:
			break
		s=s.rstrip()	
		if	not re.match(">",s):
			seq+=s
	ambigAlign.append(seq)
	nFile.close()		
#compare the ambigAlign and core align and replace the align nucleotides with N and - as needed
for i in range(len(seqs)):
	main=list(seqs[i])
	ambig=list(ambigAlign[i])
	
	for pos in range(len(main)):
		if ambig[pos]=="N":
			main[pos]="N"
		elif ambig[pos]=="-":
			main[pos]="-"
	seqs[i]=main
			
#go through insertions and modify all sequences appropriately
insPos=unique(insPos)
posIns=defaultdict(list)
insKeys=insertions.keys()
insKeys.sort()

#create an array of positions where each position has at least 1 insertion and each sequence is accounted for
for seqNum in insKeys:
	seqIns=insertions[seqNum]
	for position in insPos:
		if position in seqIns:
			posIns[position].append(seqIns[position])#keep the insertion
		else:
			posIns[position].append(seqs[seqNum][position])#keep the original snp if an insertion isnt made in that sequence

#insert the insertions into the alignment and keep an original position to new position reference map

posInsKeys=posIns.keys()
posInsKeys.sort()
offset=0
for pos in posInsKeys:
	finalPos=pos+offset
	lenIns=0	
	for ins in posIns[pos]:
		if len(ins)>lenIns:
			lenIns=len(ins)		
	offset+=lenIns-1 #-1 as it contains the unchanged snp beforehand
	for seqNum in insKeys:
		seqIns= posIns[pos][seqNum]
		gapInc=lenIns-len(seqIns)
		if gapInc>0 and gapInc<lenIns-1:
			print "Gap with multiple length insertion found at reference position "+str(pos)+" ("+str(finalPos)+" in final alignment). Shorter insertion placed at start of whole insertion."
		gap=""
		for g in range(gapInc):
			gap+="-"		
		insertionNew=seqIns+gap
		seqs[seqNum][pos]=insertionNew
	#modify the reference
	refGap=""
	for g in range(lenIns-1):	
		refGap+="-"
	refIns=ref[pos]+refGap
	ref[pos]=refIns


if intervalF!="NONE":
	print 'Masking interval'
	#mask the intervals as required. Any insertions created inside an interval will also be masked
	while 1:
		s=intF.readline()
		s=s.rstrip()
		if not s:
			break
		sections=s.split("\t")
		start=int(sections[0])-1
		if len(sections)==1:#single position
			end=start
		else:	
			end=int(sections[1])-1
		for pos in range(start,end+1):
			ref[pos]=re.sub(".","N",ref[pos])
			for seqNum in range(len(seqs)):
				seqs[seqNum][pos]=re.sub(".","N",seqs[seqNum][pos])
	intF.close()

#create save file for whole genome 
try:
	path=os.path.join("..","SNP_wholeGenome.fasta")
	saveGenome=open(path,'w')
except IOError:
	print 'no room for save file'
	sys.exit()

print 'Writing whole genome alignment'
#first write the reference, if wanted
if refInc=="1":
	saveGenome.write(">reference\n")
	seqStr="".join(ref)
	saveGenome.write(seqStr+"\n")
#turn each sequence into fasta format and output
for seqNum in range(len(names)):
	saveGenome.write(">"+names[seqNum]+"\n")
	seqStr="".join(seqs[seqNum])
	saveGenome.write(seqStr+"\n")
saveGenome.close()

print 'Finding and counting invariant sites'
#go through the 2D list above and if the elements in the 2nd dimension are all the same, skip it (but count if its an ACGT), checking for if they would be the same without the skipped characters
#Elements longer than 1 character (i.e. insertions), need to be expanded and checked character by character
#If the reference is to be included, also check the reference sequence list 
#keep a list of the genome positions where variants are retained in the SNP alignment
baseCount=[0,0,0,0]
baseCountWithSkipped=[0,0,0,0]
basePattern=["A","C","G","T"]
snpAlign=[]
genPosMap=[]
for i in range(len(seqs[0])):
	column=[]
	if refInc=="1": #include reference, if asked for
		column.append(ref[i])	
	for seqNum in range(len(names)):
		column.append(seqs[seqNum][i])
	if len(column[0])==1: #non-insertion
		#check for skipped sites and remove them 
		skippedSites=0
		patterns=set(column)
		for skip in remChars:
			if skip in patterns:
				patterns.remove(skip)
				skippedSites=1
		if len(patterns)==0: #Contained only skipped sites so ignore this site
			continue		
		if len(patterns)==1: #invariant site so count the base. If there was no skipped sites, count in 'pure' base count too
			base=patterns.pop()
			if base in basePattern:
				baseCountWithSkipped[basePattern.index(base)]+=1
				if skippedSites==0:
					baseCount[basePattern.index(base)]+=1
				continue
		#any site that remains is a variable site so take note of its genome position and add to the SNP alignment	
		genPosMap.append([len(snpAlign)+1, i+1])
		snpAlign.append(column)	
		continue
	
	#length longer than 1 means an insertion site. Split the insertions up and do as above
	for x in range(len(column[0])):
		miniCol=[]
		for y in range(len(column)):
			miniCol.append(column[y][x])
		
		#check for skipped sites and remove them 
		skippedSites=0
		patterns=set(miniCol)
		for skip in remChars:
			if skip in patterns:
				patterns.remove(skip)
				skippedSites=1
		if len(patterns)==0: #Contained only skipped sites so ignore this site
			continue		
		if len(patterns)==1: #invariant site so count the base. If there was no skipped sites, count in 'pure' base count too
			base=patterns.pop()
			if base in basePattern:
				baseCountWithSkipped[basePattern.index(base)]+=1
				if skippedSites==0:
					baseCount[basePattern.index(base)]+=1
				continue
		#any site that remains is a variable site so take note of its genome position and add to the SNP alignment	
		genPosMap.append([len(snpAlign)+1, i+1])
		snpAlign.append(miniCol)	
		continue

#transpose the SNP alignment 2d list to get it per sequence, not per column
transposeSnpAlign=[list(i) for i in zip(*snpAlign)]

#create save file for SNP alignment
try:
	path=os.path.join("..","SNPonly.fasta")
	saveSNP=open(path,'w')
except IOError:
	print 'no room for save file'
	sys.exit()

print "Writing variant site alignment"
#output the new sequences to file, plus reference if asked for
if refInc=="1":
	saveSNP.write(">reference\n")
	saveSNP.write("".join(transposeSnpAlign[0])+"\n")
	transposeSnpAlign.pop(0)	
for num in range(len(names)):
	saveSNP.write(">"+names[num]+"\n"+"".join(transposeSnpAlign[num])+"\n")		
saveSNP.close()

print "Writing base counts for Stamatakis correction method. Only sites that are pure invariant are in baseCount.txt. Sites that are invariant if skipped characters are not counted are in baseCountWithSkipped.txt"
#create save files for counts
try:
	path=os.path.join("..","baseCount.txt")
	saveCount=open(path,'w')
except IOError:
	print 'no room for save file'
	sys.exit()
try:
	path=os.path.join("..","baseCountWithSkipped.txt")
	saveCountWS=open(path,'w')
except IOError:
	print 'no room for save file'
	sys.exit()

#write out the A, C, G, T counts for the stamatakis method
for base in baseCount:
	saveCount.write(str(base)+" ")
saveCount.write("\n")
saveCount.close()	
for base in baseCountWithSkipped:
	saveCountWS.write(str(base)+" ")	
saveCountWS.write("\n")	
saveCountWS.close()

print "Writing SNP alignment positions to genome position map"
#create save file for map
try:
	path=os.path.join("..","snpAlignToGenomeMap.txt")
	saveMap=open(path,'w')
except IOError:
	print 'no room for save file'
	sys.exit()

for pos in genPosMap:
	saveMap.write(str(pos[0])+"\t"+str(pos[1])+"\n")
saveMap.close()	

print "Done"
sys.exit()


