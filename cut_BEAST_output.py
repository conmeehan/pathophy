#! /usr/bin/env python
import sys
import re
import argparse

"""
Author: Conor Meehan 18/05/2018
Script takes a MCMC log or trees file from BEAST and will cut them at a specific step
E.g. the file may have more than 40,000 steps and you wish to cut at that step (step requested will be included)

Inputs
File to cut
File type (log or trees)
Last step to include

Output
New log or trees file. Will be named same as input with _cut<number>
e.g. if input it Dataset.log and the last step is 40000, the output is Dataset_cut40000.log

Usage
python cut_BEAST_output.py --file FileToCut --type <log/trees> --step <finalStepToInclude>
"""

parser = argparse.ArgumentParser()
parser.add_argument('--file', required=True, help='MCMC log or trees file to cut')
parser.add_argument('--type', required=True, help='"log" for log file input or "trees" for trees file input (Case insensitive)')
parser.add_argument('--step', type=int, required=True, help='Last step to include in output')

args = parser.parse_args()

#get the file and ensure the type is lower case
try:
	file=open(args.file, 'rU')
except IOError:
	print "\n file not found in directory."
	sys.exit()
type=args.type
type=type.lower()

#If its a log file, cut at the position which will be from the start of the line until the first tab
if type=="log":
	#create the save file
	outName=re.sub(".log","",args.file)
	try:
		save=open(outName+"_cut"+str(args.step)+".log",'w')
	except IOError:
		print 'no room for save file'
		sys.exit()
	
	#parse the file
	while 1:
		line=file.readline()
		if not line:
			break
		if re.match("#", line) or re.match("Sample", line):
			save.write(line)
			continue
				
		sections=line.split("\t")
		if int(sections[0])>args.step:
			break
		save.write(line)	
	save.close()
	sys.exit()	

#if its a trees file, its tree STATE_<step> =
if type=="trees":
	#create the save file
	outName=re.sub(".trees","",args.file)
	try:
		save=open(outName+"_cut"+str(args.step)+".trees",'w')
	except IOError:
		print 'no room for save file'
		sys.exit()
	
	#parse the file
	while 1:
		line=file.readline()
		if not line:
			break
		if not re.match("tree STATE_", line):
			save.write(line)
			continue
				
		currentStep=re.sub("tree STATE_","",line)
		currentStep=int(re.sub(" .*","",currentStep))
		if currentStep>args.step:
			break
		save.write(line)
	save.close()
	sys.exit()

print "Type was neither log nor trees so nothing was done"
sys.exit()		
