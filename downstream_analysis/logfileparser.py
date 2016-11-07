#!/usr/bin/env python
import os
import sys
import argparse




parser = argparse.ArgumentParser(description='Aim of this code is to read Structural Variant Files and find TMPRSS2:ERG fusion. To run this command logfileparser.txt -i cotainerfile.txt -o complited.txt')
parser.add_argument('-i',type=str, help='path of the top folder which contains all of the batch directories')
parser.add_argument('-o',type=str,help='Path to Reference dictionaries.')
args=parser.parse_args()

### checks arguments ###
print "Welcome to the mutational database creator."
for item in vars(args): ###
	if getattr(args,item) == None:
		raise ValueError , "All arguments must be set"



if args.i == args.o:
	raise ValueError , 'Warning! inplut and output file path is same !!!!!' ## to avaid erasing input file	

args.i= os.getcwd() + '/' + args.i
args.o= os.getcwd() + '/' + args.o



filein=args.i
with open(filein,'rb') as myfile:
	content = [i.rstrip().split('\t')  for i in myfile if '##' not in i]

memory=[]
storage=[]
for pivot in content:
	temp=[]
	
	if pivot[0] in memory:
		continue
	for i in content:
		if pivot[0] == i[0]:
			temp.append(i)
			
	storage.append(temp)
	memory.append(pivot[0])


fileout=args.o
with open(fileout,'w') as log:
	for patient in storage:
		for element in patient[0][0:2]:
			log.write(element + '\t')
		
		if len(patient) > 1:
			for specimen in patient:
				for element in specimen[2::]:
					log.write(element + '\t')
		else:
			for specimen in patient:
				for element in specimen[2::]:
					log.write(element +'\t')
			
		log.write('\n')
print "Done!"
			
		
		