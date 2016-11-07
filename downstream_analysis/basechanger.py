#!/usr/bin/env python

import os
import re
from guppy import hpy
h=hpy()

### this code creates txt files that can be parsed by basechangestats.py ### without this code oth ohter wont work. ###
### CREATE A "basechange" DIRECTORY WITH CONTAINING "input" AND "output" DIRECTORY###
### direct to maindirectory path to the input directry which is in the basechange directory ###
maindirectory= '/Volumes/blackhole/ICGC_other_than_PCA_consensus/ICGC_other_than_PCA_AR250/ICGC_other_than_PCA_AR250_final_output_from_all_batches'
os.chdir(maindirectory)
topdir=maindirectory

def manyakreader(filepath):
	df=[]
	with open(filepath,'rb') as myfile:
		for line in myfile:
			if '#' in line:
				pass
			else:
				temp=line.rstrip().split("\t")[0:2] + line.rstrip().split("\t")[3:5] + [line.rstrip().split("\t")[7]]
				df.append(temp)
				
		return df

				
#print manyakreader('/Volumes/blackhole/ICGC_other_than_PCA_consensus/ICGC_other_than_PCA_AR250/ICGC_other_than_PCA_AR250_annotated_output_from_all_batches/OV-US/DO27978/annotated_full.0ae2193f-0d68-485a-b8c2-7568cbcce33e.indel_consensus')
				
listdir=os.listdir(maindirectory)

for project in listdir:
	filepath='../basechange/input/' + project +'_all_snv_txt'
	snvlog= open(filepath,'w')
	print 'Single Project Before MEMORY ALLOCATION CHECK'
	print h.heap()
	if project == '.DS_Store': ### to get around with a stupid error 
		continue
	
	print project
	for patient in os.listdir("/".join([os.getcwd(),project])):
		if patient == '.DS_Store': ### to get around with a stupid error 
				continue
		print patient
		for files in os.listdir("/".join([os.getcwd(),project,patient])):
			if files == '.DS_Store': ### to get around with a stupid error 
				continue
			print files
			snvlog.flush()
			if "random_vcf" in files:
				pass
			if ".txt" in files:
				pass
			if 'snv_mnv' in files:
					for line in manyakreader("/".join([topdir,project,patient,files])):
						for i in line:
							snvlog.write(i + '\t')
						snvlog.write('\n')
					break
			print "Single Patient After MEMORY ALLOCATION CHECK"
			print h.heap()
			
				
	snvlog.close()


			
		
