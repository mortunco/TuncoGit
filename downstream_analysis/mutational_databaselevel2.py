#!/usr/bin/env python
import os
import re
import subprocess
import cPickle as pickle
import time
import sys
import argparse
import gzip
import tarfile
import pandas

parser = argparse.ArgumentParser(description='Aim of this code is to read Structural Variant Files and find TMPRSS2:ERG fusion.')
parser.add_argument('-r',type=str,help='Path to Reference dictionaries.')
parser.add_argument('-snv',type=str,help='Path to SNV compiledfiles.')
parser.add_argument('-indel',type=str,help='Path to indel compiledfiles.')
parser.add_argument('-slist',type=str,help='Path to selected specimen list')
parser.add_argument('-output',type=str,help='Path to output dataframe location and name.')

args=parser.parse_args()


### checks arguments ###
print "Welcome to the mutational database creator."
for item in vars(args): ###
	if getattr(args,item) == None:
		raise ValueError , "All arguments must be set"

def whattimeisit():
	whattimeistit=time.localtime()
	days=['Mon','Tue','Wed','Thu','Fri','Sat','Sun']
	return '{0}/{1}/{2} - {6} - {3}:{4}:{5}'.format(time.localtime()[0],time.localtime()[1],time.localtime()[2],time.localtime()[3],time.localtime()[4],time.localtime()[5],days[time.localtime()[6]])


def mutation_numberreader(pathtofile):
	adict={}
	with open(pathtofile,'rb') as sourcefile:
		for line in sourcefile:
			adict[line.rstrip().split('\t')[2]] = line.rstrip().split('\t')[0::]
	return adict

with open('{0}/reference_dictionaries/reference_dict.p'.format(args.r),'rb') as fp:
	fileprefix_reference= pickle.load(fp)
fp.close()

with open('{0}/reference_dictionaries/FILEPREFIXbridgeCLINICAL.p'.format(args.r),'rb') as fp:
	bridge= pickle.load(fp)
fp.close()

with open('{0}/reference_dictionaries/clinical_reference_dict.p'.format(args.r),'rb') as fp:
	clinical_reference= pickle.load(fp)
fp.close()

with open('{0}/reference_dictionaries/folderprefix_reference_dict.p'.format(args.r),'rb') as fp:
	folder_prefix=pickle.load(fp)
fp.close()

with open('{0}/reference_dictionaries/consensus_dictionary.pc'.format(args.r),'rb') as fp:
	consensusdict=pickle.load(fp)
fp.close()

with open('/Users/morova/Google Drive/TuncProjectStep2-3/reference_dictionaries/consensusclinical_dictionary.pc','rb') as fp:
 	consensus_clinical = pickle.load(fp)

selected_specimen_list=[]
with open(args.slist,'rb') as sourcefile:
	for line in sourcefile:
		selected_specimen_list.append(line.rstrip())
		
		
### TMPRSS2:ERG postive filter###
sv_list=[]
with open('/Users/morova/Google Drive/TuncProjectStep2-3/Yunus_HPC/PCA_ICGC/TMPRSS2info/SV_number_nospace.txt','rb') as sv_doc:
	for line in sv_doc:
		sv_list.append(line.rstrip().split('\t'))

for i in xrange(len(sv_list)):


	if len(sv_list[i]) == 4 or 'TMPRSS2' not in sv_list[i][4]:
		sv_list[i].append('Negative')

	if len(sv_list[i]) > 5:
		while len(sv_list[i]) > 5:
			sv_list[i].pop()
	if 'TMPRSS2' not in sv_list[i][4]:
		sv_list[i][4] = 'Negative'

	sv_dict= {}

for specimen in sv_list:
	sv_dict[specimen[3]] = [specimen[4],specimen[0],specimen[1]]
		
		
#print selected_specimen_list

snv_list=mutation_numberreader(args.snv)
indel_list=mutation_numberreader(args.indel)


### clinical_reference[bridge['f9c52187-2e82-d58a-e040-11ac0d484fc4'][0]] = ['DO50427', 'PRAD-UK', 'Primary tumour - solid tissue', 'Gleason', '3+4', 'NA']
df=[]
ct=0
for specimen in selected_specimen_list:
	#clinical_reference[consensusdict['DO38671'][0:1] + ['Unknown'] + [consensusdict['DO38671'][0:1]] , snv_list[patient][2::], indel_list[patient][2::]
	if specimen not in snv_list.keys():
		snv_list[specimen]=[specimen,'consensus','-1|-1','-1']

		
	if specimen not in indel_list.keys():
		indel_list[specimen]=[specimen,'consensus','-1|-1','-1']

	# print clinical_reference[bridge[specimen][0]][::-1][3:4]
	# #print [sv_dict[specimen][0]]
	# print clinical_reference[bridge[specimen][0]][3::]
	# print snv_list[specimen][2:4]
	# print indel_list[specimen][2:4]
		
		
	if specimen in sv_dict.keys():
		temp= [specimen] + clinical_reference[bridge[specimen][0]][::-1][4::] + [sv_dict[specimen][0]] + clinical_reference[bridge[specimen][0]][2::] + snv_list[specimen][3::] + indel_list[specimen][3::] 
	else:
		try:
			if bridge[specimen][0] not in clinical_reference.keys():
				temp= [specimen] + clinical_reference[bridge[specimen][0]][::-1][4::] + ['Unknown'] + clinical_reference[bridge[specimen][0]][2::] + snv_list[specimen][3::] + indel_list[specimen][3::]
			else:
				temp= [specimen] + consensus_clinical[specimen][::-1][2::] + ['Unknown'] +  [ consensus_clinical[specimen][3],'NA','NA','NA' ] + snv_list[specimen][3::] + indel_list[specimen][3::]
		except:
			ct += 1
			print ct
			
	df.append(temp)
	

print ct

df2=pandas.DataFrame(df)
 
df2.columns=['Aliquot-ID','ProjectName','ICGCPatientID','TMPRSS2:ERG_Fusion_Status','Specimen_Type','Staging_System','Cancer_Degree','Cancer_Degree_ALT','AR250_SNV','FDR','AR250_INDEL','FDR']

df2.to_csv(args.output)


	
