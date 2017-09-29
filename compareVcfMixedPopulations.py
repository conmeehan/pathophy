#! /usr/bin/env python
import sys
import re
import argparse


"""
Author: Conor Meehan 29/09/17

This script is used to try estimate mixed populations from whole genome sequencing VCF files.

Input is 2 vcf files output from a SNP calling pipeline such as snippy but at different minumum proportion for variant evidence (fractions)
E.g. run snippy with --minfrac 0.3 and compare to default (0.9)

Output are the SNPs (genome position, ref base, alt base) which are present in the lower fraction cut-off vcf but not the upper

User may supply a specific set of positions to output into a separate filtered list. These may be drug resistant sites, lineage defining SNPs etc. original output file will still contain these, these positions willm just appear in both files.
Each line should be a different should be start of interval<tab character>end of interval per line (inclusive). Alternatively, if only one position on a line, that single site will be filtered.


NOTE: in theory, this should work with all vcf files but only has been tested with those output from snippy
NOTE: only SNPs are output. Indels, mnp and complex are skipped (for now)
Usage:
python compareVcfMixedPopulations.py --lower vcfLowerFraction --upper vcfUpperFraction --filter positionFile (optional)

"""

#parse the inputs
parser = argparse.ArgumentParser()
parser.add_argument('--lower', required=True, help='VCF file output from pipeline at lower fraction')
parser.add_argument('--upper', required=True, help='VCF file output from pipeline at higher fraction')
parser.add_argument('--filter', required=False, default="NONE", help='List of genome positions to output separately (optional)')
args = parser.parse_args()

#read in files
try:
	lowF=open(args.lower, 'rU')
except IOError:
	print "\n Lower fraction VCF file not found in directory."
	sys.exit()
try:
	upF=open(args.upper, 'rU')
except IOError:
	print "\n Upper fraction VCF file not found in directory."
	sys.exit()

#create save file for all positions
try:
	saveList=open("mixedPositions.txt",'w')
except IOError:
	print 'no room for save file'
	sys.exit()

#get the filter file, if needed
if args.filter!="NONE":
	try:
		filtF=open(args.filter, 'rU')
	except IOError:
		print "\n Upper fraction VCF file not found in directory."
		sys.exit()
	#create save file for filtered positions
	try:
		saveFilt=open("mixedPositionsFiltered.txt",'w')
	except IOError:
		print 'no room for save file'
		sys.exit()



#get the genome position, ref base and alt base from the lower VCF file
lowerDict={}
while 1:
	line=lowF.readline()
	if not line:
		break
	if re.match("#",line):#on an info line so skip
		continue
	sections=line.split("\t")
	type=sections[7] #the type is listed at the end of section 7
	type=re.sub(".*TYPE=","",type)
	type=re.sub(";.*","",type)
	if type.lower()=='snp':
		lowerDict[int(sections[1])]=[sections[3], sections[4]]
lowF.close()
lowerSet=set(lowerDict.keys())
		
#get the genome position from the upper VCF file
upperSet=set()
while 1:
	line=upF.readline()
	if not line:
		break
	if re.match("#",line):#on an info line so skip
		continue
	sections=line.split("\t")
	type=sections[7] #the type is listed at the end of section 7
	type=re.sub(".*TYPE=","",type)
	type=re.sub(";.*","",type)
	if type.lower()=='snp':
		upperSet.add(int(sections[1]))
upF.close()

	
#get the positions that are in the lower set but not the upper set
inLow=lowerSet.difference(upperSet)

#print that set to file
for pos in inLow:
	saveList.write(str(pos)+"\t"+"\t".join(lowerDict[pos])+"\n")
saveList.close()	

#save filtered list, if requested
if args.filter!="NONE":
	while 1:
		line=filtF.readline()
		if not line:
			break
		line=line.rstrip()
		sections=line.split("\t")
		start=int(sections[0])
		if len(sections)==1:#single position
			end=start
		else:	
			end=int(sections[1])-1	
		for pos in range(start,end+1):
			if pos in inLow:
				saveFilt.write(str(pos)+"\t"+"\t".join(lowerDict[pos])+"\n")
	saveFilt.close()		
sys.exit()