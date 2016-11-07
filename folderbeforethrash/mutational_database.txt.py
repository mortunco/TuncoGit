##! /usr/bin/env python


#### eger bunu goruyorsam ####
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
parser.add_argument('-p',type=str, help='path of the top folder which contains all of the batch directories')
parser.add_argument('-r',type=str,help='Path to Reference dictionaries.')
args=parser.parse_args()

def whattimeisit():
	whattimeistit=time.localtime()
	days=['Mon','Tue','Wed','Thu','Fri','Sat','Sun']
	return '{0}/{1}/{2} - {6} - {3}:{4}:{5}'.format(time.localtime()[0],time.localtime()[1],time.localtime()[2],time.localtime()[3],time.localtime()[4],time.localtime()[5],days[time.localtime()[6]])

print "Welcome to the mutational database creator."
for item in vars(args): ###
	if getattr(args,item) == None:
		raise ValueError , "All arguments must be set"


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

### TMPRSS2:ERG postive filter###
sv_list=[]
with open('/Users/morova/Google Drive/TuncProjectStep2-3/Yunus_HPC/PCA_ICGC/SV_number_sanger.txt','rb') as sv_doc:
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

	print sv_list[i]


### CONTROL INDEL ###
control_dict={}
with open('/Users/morova/Google Drive/TuncProjectStep2-3/Yunus_HPC/PCA_ICGC/Control_Experiment/control_final_indellog_from_all_batches/Control_sanger_indel.txt','rb') as control:
	for line in control:
		control_dict[line.rstrip().split('\t')[0]] = line.rstrip().split('\t')[0::]


### AR INDEL ###
ar_dict={}
with open('/Users/morova/Google Drive/TuncProjectStep2-3/Yunus_HPC/PCA_ICGC/AR_250widen/Prostate_Cancer_AR_250widen_new_final_indellog_from_all_batches/Prostate_Cancer_AR_250widen_indel_sanger.txt','rb') as AR:
	for line in AR:
		ar_dict[line.rstrip().split('\t')[0]] = line.rstrip().split('\t')[0::]


### ['DO51068', 'PRAD-CA', 'Primary tumour - solid tissue', 'Gleason', '4+3', 'NA'] ### code beloew returns this.###
#print clinical_reference[bridge[sv_list[1][3]][0]]

###['PRAD-CA', 'DO51068', 'sanger', '45d0ccb2-641f-4348-b3a8-61f4113cd85b', '1|ERG:TMPRSS2|820'] ### sv_list structure .
df=[]
for patientindex in range(len(sv_list)):

	patient_structure=[sv_list[patientindex][0:2]+[sv_list[patientindex][4]] + clinical_reference[bridge[sv_list[patientindex][3]][0]][2::] + control_dict[sv_list[patientindex][1]][2:4] + ar_dict[sv_list[patientindex][1]][2:4] ]
	patient_structure=[i for a in patient_structure for i in a]
	df.append(patient_structure)


df2=pandas.DataFrame(df)

df2.columns=['ProjectName','ICGCPatientID','TMPRSS2:ERG_Fusion_Status','Specimen_Type','Staging_System','Cancer_Degree','Cancer_Degree_ALT','Control_BED_INDEL_count','FDR','AR_250_BED','FDR']


print 'Data frame was saved under /Users/morova/Google Drive/TuncProjectStep2-3/Yunus_HPC/PCA_ICGC/ as ICGC_PCa_dataframe.csv'
df2.to_csv('/Users/morova/Google Drive/TuncProjectStep2-3/Yunus_HPC/PCA_ICGC/ICGC_PCa_dataframe.csv')












