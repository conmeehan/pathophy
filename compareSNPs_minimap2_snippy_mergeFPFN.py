#! /usr/bin/env python
import sys
import re
import argparse
import copy
from collections import defaultdict
import pysam
from types import MappingProxyType
import gzip

"""
Author: Conor Meehan 11/04/2025

Script takes two VCF files, one from minimap2/paftools.js with all data (snps, indels, large gaps etc) and one from snippy
The reference genome for both of these should be the same
This is used to compare SNP calling (no other types) between two pipelines on the same reference genome.
Output are metrics (FP, TP, FN) and also distance to a indel/large gap and other descriptives for each SNP
Outputs are modified so that FP and FN within a certain distance of each other are merged into a TP, listed at the FN position

NOTE: The contig name in both VCF files should be the same
NOTE2: Only SNPs will be compared between the VCFs (i.e. single character (not -) in the REF and ALT sections of the VCF)
NOTE3: Requires pysam to process VCF
NOTE4: The minimap2 output stores indels of any size as the letter *before* the indel in the ref or alt and then the nucleotides in the indel.
	This means that indels are stored as the VCF position +1 and the length is the insertion (either in the ref or alt) length -1

Input:
--minimap2VCF The VCF file from minimap2
--snippyVCF The VCF file that is to be compared to the reference from snippy
--dist The maximum distance to consider FP and FN for merging
--prefix A string to put at the start of each output file

Output:
A semi-confusion matrix in the format:
True positives	False positives	False negatives

A set of metrics:
Precision (also known as Positive predictive value (PPV)) : (\frac{TP}{TP + FP})
Recall (Sensitivity) (also known as True positive rate (TPR)): (\frac{TP}{TP + FN})
F1 Score: (2 \times \frac{\text{Precision} \times \text{Recall}}{\text{Precision} + \text{Recall}})
False negative rate (FNR): (\frac{FN}{TP + FN})



A tab-delimited file of TP, FP and FN (in that order) in this format:
Chromosome Name	Site	Type	Closest Gap	Distance to closest gap	Gap size	Closest SNP	Distance to closest SNP	Closest SNP Type
Where type is TP, FP or FN, dist is the distance to the start, end of the nearest indel or large gap and the size of that indel
The file name has _siteMetrics.txt


To run:
python compareSNPs_minimap2_snippy.py --minimap2VCF  minimap2VCFfile --snippyVCF snippyVCFfile --dist maxDistForMerging --prefix StringBeforeOutputs

"""	

#parse the inputs
parser = argparse.ArgumentParser(description= 'Script takes two VCF files, one from minimap2 (as groudn truth) and snippy, and creates a semi-confusion matrix (true positive, false positive, true negative) and other metrics of the sample VCF SNPs vs the reference SNPs')
parser.add_argument('--minimap2VCF', required=True, help='The VCF file from minimap2')
parser.add_argument('--snippyVCF', required=True, help='The VCF file that is to be compared to the reference from snippy')
parser.add_argument('--prefix', required=True, help='A string to put at the start of each output file')
parser.add_argument('--dist', required=True, help='A maximum distance to consider merging FP and FNs')

args = parser.parse_args()

#read in the files
try:
	rVCFF=pysam.VariantFile(args.minimap2VCF)
except IOError:
	print("\n Minimap2 VCF file not found.")
	sys.exit()
try:
	sVCFF=pysam.VariantFile(args.snippyVCF)
except IOError:
	print("\n Snippy VCF file not found.")
	sys.exit()

#Get the prefix
prefix = str(args.prefix)

#Get the distance
maxDist= int(args.dist)

#Get all the SNPs from the minimap2 ref VCF file. Store them as a contig to SNP position to SNP contents map
#Alongside this get all the indels (a - in the ref or alt) and store their position and length
refVCF_Map=defaultdict(lambda: defaultdict(str))
refVCF_indel=defaultdict(lambda: defaultdict(int))


for vcfRecord in rVCFF:
	contig=vcfRecord.chrom
	contigPos=int(vcfRecord.pos)
	ref=vcfRecord.ref
	alt=vcfRecord.alts[0]
	#check for indels and if found save position and length otherwise check for SNP and save
	if len(ref)==1 and len(alt)==1: #is a true SNP and not MNP etc
		refVCF_Map[contig][str(contigPos)]=alt.upper()
	elif len(ref)>1 and len(alt)==1: #Is a deletion relative to reference
		refVCF_indel[contig][str(contigPos+1)]=len(ref)-1
	elif len(ref)==1 and len(alt)>1: #Is an insertion relative to reference
		refVCF_indel[contig][str(contigPos+1)]=len(alt)-1	

#vcfMap is complete so stop more elements being added
refVCF_Map_set = MappingProxyType(refVCF_Map)

rVCFF.close()

#Get all the SNPs from the snippy VCF file. Store them as a contig to SNP position to SNP contents map
sampleVCF_Map=defaultdict(lambda: defaultdict(str))

for vcfRecord in sVCFF:
	contig=vcfRecord.chrom
	contigPos=str(vcfRecord.pos)
	ref=vcfRecord.ref
	SNP=vcfRecord.alts[0]
	if len(ref)==1 and len(SNP)==1: #is a true SNP and not MNP etc
		sampleVCF_Map[contig][contigPos]=SNP.upper()

#vcfMap is complete so stop more elements being added
sampleVCF_Map_set = MappingProxyType(sampleVCF_Map)

sVCFF.close()


#Compare the two maps in the following ways:
#SNP in both: true positive
#SNP in sample VCF but not ref VCF: false positive
#SNP in ref VCF but not sample VCF: false negative
#Save them as counts and also as nested dictionaries of type to then chromosome name to then position to then distance to nearest gap and the length of that gap
TP=0
FP=0
FN=0
typeMap=defaultdict(lambda: defaultdict(lambda: defaultdict(str)))

#First loop through all the SNPs in the minimap2 file and populate the true positives and false negatives

for rVCF_Chrom in refVCF_Map_set:
	for rVCF_Pos in refVCF_Map_set[rVCF_Chrom]:
		if rVCF_Chrom in sampleVCF_Map_set:
			if rVCF_Pos in sampleVCF_Map_set[rVCF_Chrom]:
				if refVCF_Map_set[rVCF_Chrom][rVCF_Pos] == sampleVCF_Map_set[rVCF_Chrom][rVCF_Pos]:
					TP+=1 #SNP is found in both so is a TP

					#Go through the indel map and find the one closest to this SNP
					rVCF_Pos_int=int(rVCF_Pos)
					closest=[0,0,0]
					for indel_Pos in refVCF_indel[rVCF_Chrom]:
						indel_Pos_int=int(indel_Pos)
						distIndel=abs(indel_Pos_int-rVCF_Pos_int)
						if distIndel<abs(indel_Pos_int-closest[0]):
							closest=[indel_Pos_int,distIndel,refVCF_indel[rVCF_Chrom][indel_Pos]]
					typeMap["TP"][rVCF_Chrom][rVCF_Pos]=closest
				else:
					FN+=1 #SNP position in both but not the same, so its a FN

					#Go through the indel map and find the one closest to this SNP
					rVCF_Pos_int=int(rVCF_Pos)
					closest=[0,0,0]
					for indel_Pos in refVCF_indel[rVCF_Chrom]:
						indel_Pos_int=int(indel_Pos)
						distIndel=abs(indel_Pos_int-rVCF_Pos_int)
						if distIndel<abs(indel_Pos_int-closest[0]):
							closest=[indel_Pos_int,distIndel,refVCF_indel[rVCF_Chrom][indel_Pos]]
					typeMap["FN"][rVCF_Chrom][rVCF_Pos]=closest


			else:
				FN+=1 #SNP position not found in VCF so FN
				
				
				#Go through the indel map and find the one closest to this SNP
				rVCF_Pos_int=int(rVCF_Pos)
				closest=[0,0,0]
				for indel_Pos in refVCF_indel[rVCF_Chrom]:
					indel_Pos_int=int(indel_Pos)
					distIndel=abs(indel_Pos_int-rVCF_Pos_int)
					if distIndel<abs(indel_Pos_int-closest[0]):
						closest=[indel_Pos_int,distIndel,refVCF_indel[rVCF_Chrom][indel_Pos]]
				typeMap["FN"][rVCF_Chrom][rVCF_Pos]=closest

		else:
			print("Chromosome "+rVCF_Chrom+" found in reference VCF but not sample VCF so this is being skipped and wont be in any comparisons")		

#Repeat with the VCF file but only look for those that dont match to catch the false positives
for sVCF_Chrom in sampleVCF_Map_set:
	for sVCF_Pos in sampleVCF_Map[sVCF_Chrom]:
		if sVCF_Chrom in refVCF_Map_set:
			if sVCF_Pos in refVCF_Map_set[sVCF_Chrom]:
				if refVCF_Map_set[sVCF_Chrom][sVCF_Pos] == sampleVCF_Map_set[sVCF_Chrom][sVCF_Pos]:
					continue #SNP is found in both so is a TP but already counted so continue
				else:
					FP+=1 #SNP position in both but not the same, so its a FP
					
					#Go through the indel map and find the one closest to this SNP
					sVCF_Pos_int=int(sVCF_Pos)
					closest=[0,0,0]
					for indel_Pos in refVCF_indel[sVCF_Chrom]:
						indel_Pos_int=int(indel_Pos)
						distIndel=abs(indel_Pos_int-sVCF_Pos_int)
						if distIndel<abs(indel_Pos_int-closest[0]):
							closest=[indel_Pos_int,distIndel,refVCF_indel[sVCF_Chrom][indel_Pos]]
					typeMap["FP"][sVCF_Chrom][sVCF_Pos]=closest
			else:
				FP+=1 #SNP position not found in mauve so FP
				
				#Go through the indel map and find the one closest to this SNP
				sVCF_Pos_int=int(sVCF_Pos)
				closest=[0,0,0]
				for indel_Pos in refVCF_indel[sVCF_Chrom]:
					indel_Pos_int=int(indel_Pos)
					distIndel=abs(indel_Pos_int-sVCF_Pos_int)
					if distIndel<abs(indel_Pos_int-closest[0]):
						closest=[indel_Pos_int,distIndel,refVCF_indel[sVCF_Chrom][indel_Pos]]
				typeMap["FP"][sVCF_Chrom][sVCF_Pos]=closest

		else:
			print("Chromosome "+sVCF_Chrom+" found in sample VCF output but not reference VCF so this is being skipped and wont be in any comparisons")		

#Go through the FP and FN data to merge based on distance
#Create a loop from 1 to the given maximum distance and for each FP, if the closest type is an FN and less than that distance, merge them and change the FN to a TP and remove the FP
#this ensures that the FP will always be merged with its closest FN

#Keep lists of the positions that have been merged for later skipping
mergedFP=[]
mergedFN=[]
for testDist in range(1,maxDist+1):
	for typeChrom in typeMap["FP"]:
		for typePos in typeMap["FP"][typeChrom]:
			typePos_int=int(typePos)
			
			#Skip if this FP has already been merged with an FN
			if [typeChrom,typePos] in mergedFP:
				continue
				
			for typeChrom2 in typeMap["FN"]:
				for typePos2 in typeMap["FN"][typeChrom2]:
					typePos2_int=int(typePos2)
					
					#Skip if this FN has already been merged with an FP
					if [typeChrom2,typePos2] in mergedFN:
						continue
					
					distSNP=abs(typePos_int-typePos2_int)
					if distSNP<=testDist:
						typeMap["TP"][typeChrom2][typePos2]=typeMap["FN"][typeChrom2][typePos2]
						mergedFP.append([typeChrom,typePos])
						mergedFN.append([typeChrom2,typePos2])

#Remove the FN and FP entries
for FPremove in mergedFP:
	del typeMap["FP"][FPremove[0]][FPremove[1]]
for FNremove in mergedFN:
	del typeMap["FN"][FNremove[0]][FNremove[1]]


#Go through each SNP and append to the value in the dictionary the closest SNP and its type
for typeName in typeMap:
	for typeChrom in typeMap[typeName]:
		for typePos in typeMap[typeName][typeChrom]:		
			typePos_int=int(typePos)
			closest=[0,0,""]
			for typeName2 in typeMap:
				for typeChrom2 in typeMap[typeName2]:
					for typePos2 in typeMap[typeName2][typeChrom2]:
						typePos2_int=int(typePos2)
						distSNP=abs(typePos_int-typePos2_int)
						if distSNP==0: #self hit so skip
							continue 
						if closest[0]==0: #its the first one so just set this first SNP as 'closest' in case our case SNP is closer to 0 than any other SNP
							closest=[typePos2_int,distSNP,typeName2]
						if distSNP<abs(typePos_int-closest[0]):
							closest=[typePos2_int,distSNP,typeName2]
			typeMap[typeName][typeChrom][typePos]=typeMap[typeName][typeChrom][typePos]+closest


#Modify the TP, FN and FP counts
TP=TP+len(mergedFN)
FP=FP-len(mergedFP)
FN=FN-len(mergedFN)

#Create the three save files
try:
	saveMx = open(prefix+"_VCFcomparisons_semi-confusion-matrix.txt", 'w')   # open file
except IndexError:
	print("\nno room for save file.")
	sys.exit()
try:
	saveMetrics = open(prefix+"_VCFcomparisons_SNP-comparison-metrics.txt", 'w')   # open file
except IndexError:
	print("\nno room for save file.")
	sys.exit()	

try:
	saveSites = open(prefix+"_siteMetrics.txt", 'w')   # open file
except IndexError:
	print("\nno room for save file.")
	sys.exit()	

#Output the confusion matrix
saveMx.write("True Positive\tFalse Positive\tFalse Negative\n")
saveMx.write(str(TP)+"\t"+str(FP)+"\t"+str(FN)+"\n")
saveMx.close()

#Calculate and output the metrics
precision= TP/(TP+FP)
recall= TP/(TP+FN)
F1= 2*((precision*recall)/(precision+recall))
FNR = FN/(TP+FN)

saveMetrics.write("Precision\tRecall\tF1 Score\tFalse negative rate\n")
saveMetrics.write(str(precision)+"\t"+str(recall)+"\t"+str(F1)+"\t"+str(FNR)+"\n")
saveMetrics.close()

#Output the site statistics
saveSites.write("Chromosome Name\tSite\tType\tClosest Gap\tDistance to closest gap\tGap size\tClosest SNP\tDistance to closest SNP\tClosest SNP Type\n")
#Do the TP first, then FP, then FN
for typeChrom in typeMap["TP"]:	
	typePosList=list(typeMap["TP"][typeChrom].keys())
	typePosList.sort(key=int)
	for typePos in typePosList:
		saveSites.write(typeChrom+"\t"+typePos+"\tTP\t"+"\t".join(str(x) for x in typeMap["TP"][typeChrom][typePos])+"\n")
for typeChrom in typeMap["FP"]:	
	typePosList=list(typeMap["FP"][typeChrom].keys())
	typePosList.sort(key=int)
	for typePos in typePosList:
		saveSites.write(typeChrom+"\t"+typePos+"\tFP\t"+"\t".join(str(x) for x in typeMap["FP"][typeChrom][typePos])+"\n")
for typeChrom in typeMap["FN"]:	
	typePosList=list(typeMap["FN"][typeChrom].keys())
	typePosList.sort(key=int)
	for typePos in typePosList:
		saveSites.write(typeChrom+"\t"+typePos+"\tFN\t"+"\t".join(str(x) for x in typeMap["FN"][typeChrom][typePos])+"\n")
saveSites.close()

sys.exit()	
	