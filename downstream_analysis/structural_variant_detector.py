##! /usr/bin/env python
### this hpc version ###
### STABLE VERSION ###

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
parser = argparse.ArgumentParser(description='Aim of this code is to read Structural Variant Files and find TMPRSS2:ERG fusion.')
parser.add_argument('-p',type=str, help='path of the top folder which contains all of the batch directories')
args=parser.parse_args()

def whattimeisit():
	whattimeistit=time.localtime()
	days=['Mon','Tue','Wed','Thu','Fri','Sat','Sun']
	return '{0}/{1}/{2} - {6} - {3}:{4}:{5}'.format(time.localtime()[0],time.localtime()[1],time.localtime()[2],time.localtime()[3],time.localtime()[4],time.localtime()[5],days[time.localtime()[6]])

print "Welcome to the Structural Variant Detector"
for item in vars(args): ###
	if getattr(args,item) == None:
		raise ValueError , "All arguments must be set"


with open('/mnt/kufs/scratch/tmorova15/references/reference_dictionaries/reference_dict.p','rb') as fp:
	fileprefix_reference= pickle.load(fp)
fp.close()

with open('/mnt/kufs/scratch/tmorova15/references/reference_dictionaries/FILEPREFIXbridgeCLINICAL.p','rb') as fp:
	bridge= pickle.load(fp)
fp.close()

with open('/mnt/kufs/scratch/tmorova15/references/reference_dictionaries/clinical_reference_dict.p','rb') as fp:
	clinical_reference= pickle.load(fp)
fp.close()

with open('/mnt/kufs/scratch/tmorova15/references/reference_dictionaries/folderprefix_reference_dict.p','rb') as fp:
	folder_prefix=pickle.load(fp)
fp.close()


print "------------------------------------------------------------"
print 'Initializing Reference Check...'
print time.time()

print whattimeisit(),'Checking File Prefix Reference'

if fileprefix_reference['a9ec7d9e-b179-4782-a589-43c7d1642be9'][0] != 'DO51159':
	raise ValueError , 'Your File Prefix Dictionary has problems, Check dictionary or Recreate it!!!'

if bridge['a9ec7d9e-b179-4782-a589-43c7d1642be9'][1] != 'DO51159':
	raise ValueError , 'Your Bridge Dictionary has problems, Check dictionary or Recreate it !!!'

if clinical_reference['SP113005'][0] != 'DO51159':
	print 'clinical reference dictionary says:%s  but it is supposed to be DO51159' % clinical_reference['SP113005']
	raise ValueError, 'Your clinical reference dictionary has problems, Check dictionary or Recreate it !!!'

if folder_prefix['d4fa1f86-ceb6-11e5-bcb8-d4714d100685'][1] != 'PRAD-CA':
	raise ValueError, 'Your folder reference dictionary has problems, Check dictionary or Recreate it !!!'

topdirectory=args.p

print 'Checking Batch directories under {0}'.format(topdirectory)
os.chdir(topdirectory)

toplist=os.listdir(os.getcwd())
batchlist=[i for i in toplist if 'Batch' in i]
print batchlist

### This is for eleminating problematic id problem ###
project_names =[]
for batch in batchlist:
	for item in os.listdir('{0}/input/'.format(batch)):
		if item in folder_prefix.keys():
			project_names.append(folder_prefix[item][1])
		else:
			pass

project_names=list(set(project_names))

structural_log= open('{0}/SV_number.txt'.format(topdirectory),'w')
structural_log.write('#### patientid Possible TMPRSS2:ERG|FilePathinYUNUS\n')
structural_log2=open('{0}/SV_location.txt'.format(topdirectory),'w')

structural_log.flush()
structural_log2.flush()
ETSfamily=['ERG','ETV1','ETV4','ETV5']
for batch in batchlist:
	for folder in os.listdir('{0}/input/'.format(batch)):

		if len(folder.split('-')) == 5: ### There are other folders/files located in the same input folder such as log.txt of the download or manifests. This wi

			if folder not in folder_prefix.keys():

				if os.path.exists('{1}/final/Problematic_id/{0}'.format(folder,batch)):
					pass
				else:
					os.mkdir('{1}/final/Problematic_id/{0}'.format(folder,batch))

			else:

				### analysis method : sanger, dkfz, broad
				analysis_method = folder_prefix[folder][2]
				
				#indel_log.write(folder_prefix[folder][0] + '\t' + analysis_method +'\t')## write to the file which contains all of the numbers



				for name in os.listdir('{0}/input/{1}'.format(batch,folder)):
					
					if bool(re.search('somatic.sv.tar.gz', name)): ### checks if snv_mnv.vcf.gz in the iterated variable
						print 'Analysing-->', '{1}/input/{0}/'.format(folder,batch) + name
						analysis_id=name.split('.')[0]
						start=time.time()
						#print name
						#print analysis_method
						if analysis_method != 'sanger':
							continue

						### Searches through the vcf file and checks if these gene names are in the file. If any of them are not in the
						tar= tarfile.open('{1}/input/{0}/'.format(folder,batch) + name)

						for member in tar.getnames():
							if 'annot.bedpe' in member:
								SVvcf=tar.extractfile(member)
								structural_log.write(folder_prefix[folder][1]+'\t'+ folder_prefix[folder][0] + '\t' + analysis_method +'\t') ## write to the file which contains all of the numbers
								structural_log.write(analysis_id + '\t')
								
								structural_log2.write(folder_prefix[folder][1]+'\t'+ folder_prefix[folder][0] + '\t' + analysis_method +'\t') ## write to the file which contains all of the numbers
								structural_log2.write(analysis_id + '\t')									
								
								### SOME PATIENTS HAVE MULTIPLE TMPRSS ERG FUSIONS IN THEIR GENOME THEREFORE. AFTER A  SPECIMEN THERE MIGHT BE MORE THAN SINGLE TMPRSS2:ERG FUSION ENTRIES.###
								for line in SVvcf:
									TMPRSS2count=0
									ETScount=0
									line=line.rstrip().split('\t')

									if '#' in line[0]:
										continue

									if line[27] == line[37]:
										continue
									
									if ('TMPRSS2' in line[27]) or ('TMPRSS2' in line[37]):
										TMPRSS2count += 1
									
									if any(line[27] == i for i in ETSfamily) or any(line[37] == i for i in ETSfamily):
										ETScount += 1

									if TMPRSS2count * ETScount > 0:

										structural_log.write(str(TMPRSS2count * ETScount) + '|' + line[27]+':'+line[36] +'|' + line[45] + '\t')
										structural_log2.write(line[45] +'|'+ line[27] +':'+ line[0] + ':' + line[1] + ':' + line[2] + '|' + line[36] +':'+ line[3] + ':' + line[4] + ':' +  line[5]+ '\t')

								structural_log.write('\n')
								structural_log2.write('\n')





				structural_log.write('\n')
				structural_log2.write('\n')


structural_log.close()
structural_log2.close()


print whattimeisit(),'Process has finalized without any error ! Well Done !!'