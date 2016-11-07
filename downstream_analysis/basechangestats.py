#!/usr/bin/env python

import os
import matplotlib
import collections
from guppy import hpy




### HOW TO RUN IT ?###
### CREATE A BASECHANGE DIRECTORY CONTAINING INPUT AND OUTPUT SUBDIRECTORIES ###
### GIVE MAINDIRECTORY TO INPUT DIRECTORY. CODE WILL PUSH OUTPUTS TO THE OUTPUT DIRECTORY ###

maindirectory='/Users/morova/Google Drive/TuncProjectStep2-3/Yunus_HPC/PCA_ICGC/AR_250widen/PCa_concensus_AR250widen/basechange_filtered/input'
os.chdir(maindirectory)
listdir=os.listdir(maindirectory)

def basechangereader(filepath):
	df=[]
	with open(filepath,'rb') as myfile:
		for line in myfile:
			if '#' in line:
				pass
			else:
				temp=line.rstrip().split("\t")
				df.append(temp)
				
		return df
	

def	variantclassificationcounter(basechangetype):
	#basechangenetype denen list, her elemani 2 elemandan olusan bir list of list. biriinci elemani chromosome info digeri ise mutation tipi.
	outputchrom = [ 'chr'+i[0] for i in basechangetype]
	counterchrom = collections.Counter(outputchrom)
	
	outputSV = [i[1] for i in basechangetype]
	counterSV = collections.Counter(outputSV)
	
	
	return [counterchrom,counterSV]
	
for item in listdir:
	h=hpy()
	if '.DS' in item:
		continue

	#df= basechangereader('/Volumes/blackhole/ICGC_other_than_PCA_consensus/ICGC_other_than_PCA_AR250/basechange_output/input/BLCA-US_all_snv.txt')
	
	df = basechangereader(item)
	A_T=0 ; A_G=0 ; A_C=0;
	A_T_vc=[] ; A_G_vc=[] ; A_C_vc=[];
	
	T_A=0 ; T_G=0 ; T_C=0;
	T_A_vc=[] ; T_G_vc=[] ; T_C_vc=[];
	
	C_A=0 ; C_T=0 ; C_G=0;
	C_A_vc=[]; C_T_vc=[]; C_G_vc=[];
	
	G_A=0 ; G_T=0 ; G_C=0;
	G_A_vc=[]; G_T_vc=[] ; G_C_vc=[];
	
	def variantclassification(filterline):
		a=[i.split('=')[1] for i in filterline.split(';') if 'Variant_Classification' in i]
		return a[0]

	for line in df:
		if line[2] == 'A':
			
			if line[3] == 'T':
				A_T_vc.append([line[0],variantclassification(line[4])])
			elif line[3] == 'C':
				A_C_vc.append([line[0],variantclassification(line[4])])
			elif line[3] =='G':
				A_G_vc.append([line[0],variantclassification(line[4])])
			else:
				print line[3]
		
		elif line[2] == 'T':
				
			if line[3] == 'A':
				T_A_vc.append([line[0],variantclassification(line[4])])
			elif line[3] == 'C':
				T_C_vc.append([line[0],variantclassification(line[4])])
			elif line[3] =='G':
				T_G_vc.append([line[0],variantclassification(line[4])])
			
			else:
				print line[3]
			
		elif line[2] == "G":
			
			if line[3] == 'T':
				G_T_vc.append([line[0],variantclassification(line[4])])
			elif line[3] == 'C':
				G_C_vc.append([line[0],variantclassification(line[4])])
			elif line[3] =='A':
				G_A_vc.append([line[0],variantclassification(line[4])])
			else:
				print line[3]
			
		elif line[2] == 'C':
			
			if line[3] == 'T':
				C_T_vc.append([line[0],variantclassification(line[4])])
			elif line[3] == 'A':
				C_A_vc.append([line[0],variantclassification(line[4])])
			elif line[3] =='G':
				C_G_vc.append([line[0],variantclassification(line[4])])
			else:
				print line[3]
			
		else:
			print line[2]
	print "Matrix is finished"
	print h.heap() #heap is for performance estimation
	
	print item.split('_')[0]
	print ' ', 'A', 'T', 'C','G'
	print 'A', 0, len(A_T_vc), len(A_C_vc), len(A_G_vc);
	print 'T', len(T_A_vc),0,len(T_C_vc),len(T_G_vc);
	print 'C',len(C_A_vc),len(C_T_vc),0, len(C_G_vc);
	print 'G',len(G_A_vc),len(G_T_vc),len(G_C_vc),0;
	
	print "Writing Matrix..."
	os.mkdir('../output/'+item.split('_')[0])
	with open('../output/'+ item.split('_')[0]+ '/'+ item.split('_')[0]  + '_matrix.txt' , 'w') as output:
		output.write('\t'.join([' ', 'A', 'T', 'C','G']))
		output.write('\n')
		output.write('\t'.join(['A', '0', str(len(A_T_vc)), str(len(A_C_vc)), str(len(A_G_vc))]))
		output.write('\n')
		output.write('\t'.join(['T',str(len(T_A_vc)),'0',str(len(T_C_vc)),str(len(T_G_vc))]))
		output.write('\n')
		output.write('\t'.join(['C',str(len(C_A_vc)),str(len(C_T_vc)),'0', str(len(C_G_vc))]))
		output.write('\n')
		output.write('\t'.join(['G',str(len(G_A_vc)),str(len(G_T_vc)),str(len(G_C_vc)),'0']))
		output.write('\n')
	print "Writing Variant Classification vs Basechange..."	
	print h.heap()
	filename='../output/' + item.split('_')[0] +'/'+ item.split('_')[0] +  '_SV.txt' 
	with open(filename , 'w') as output:
		variables=['A_T_vc','A_G_vc','A_C_vc','T_A_vc', 'T_G_vc', 'T_C_vc','C_A_vc','C_T_vc','C_G_vc','G_A_vc', 'G_T_vc','G_C_vc']
		for i in variables:
			output.write(i + '\t')
			for key,value in zip(variantclassificationcounter(eval(i))[1].keys(),variantclassificationcounter(eval(i))[1].values()): output.write(str(key)+'|'+str(value)+';')
			output.write('\n')
		output.write('\n')
	print "Writing Chromosome vs Basechange..."	
	print h.heap()
	filename='../output/' + item.split('_')[0] +'/'+ item.split('_')[0] +'_CHROM.txt'		
	with open( filename , 'w') as output:
		variables=['A_T_vc','A_G_vc','A_C_vc','T_A_vc', 'T_G_vc', 'T_C_vc','C_A_vc','C_T_vc','C_G_vc','G_A_vc', 'G_T_vc','G_C_vc']
		for i in variables:
			output.write(i + '\t')
			for key,value in zip(variantclassificationcounter(eval(i))[0].keys(),variantclassificationcounter(eval(i))[0].values()): output.write(str(key)+'|'+str(value)+';')
			output.write('\n')
		output.write('\n')
	
	print h.heap()






	
	





