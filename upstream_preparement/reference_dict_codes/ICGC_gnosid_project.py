### This is the reference that gives information about which file prefix is corresponded to which donor ####

import cPickle as pickle
import json


f3= open('/Users/morova/Google Drive/TuncProjectStep2-3/ICGC/release_may2016.v1.tsv.txt','r')

reference={}
for line in f3:
	temp = line.rstrip().split('\t')
	project_name, donor_id, sanger, dkfz, broad, muse = temp[1] , temp[3], temp[38], temp[43], temp[48], temp[53]

	### This process is applied to all of the entries because, \ ###
	### some patients may have multiple samples which will have multiple mutation calls. ###

	sanger = sanger.split(',')
	dkfz = dkfz.split(',')
	broad = broad.split(',')
	muse = muse.split(',')


	### This commanded part was written for the test the algorithm ###
	for file_prefix in sanger:
		reference[file_prefix] = [donor_id,project_name,'sanger']
	for file_prefix in dkfz:
		reference[file_prefix] = [donor_id,project_name,'dkfz']
	for file_prefix in broad:
		reference[file_prefix] = [donor_id,project_name,'broad']
	for file_prefix in muse:
		reference[file_prefix] = [donor_id,project_name,'muse']
	#if len(sanger) > 1:
	#	print sanger
	#	for item in sanger:
	#		print '%s ----> %s' % (item, reference[item])
f3.close()



data=[]
with open('/Users/morova/Google Drive/TuncProjectStep2-3/ICGC/release_may2016.v1.4.with_consensus_calls.jsonl.txt') as jsonldata:

	for line in jsonldata:
		data.append(json.loads(line)) ### every line is a one json. 



for i in range(len(data)):
	try:
		for specimen in xrange(len(data[i]["consensus_somatic_variant_calls"]['snv_mnv'])):
			reference[data[i]["consensus_somatic_variant_calls"]['snv_mnv'][specimen]['gnos_id']]= [data[i]["icgc_donor_id"], data[i]['dcc_project_code'],'consensus' ]
		
		for specimen in xrange(len(data[i]["consensus_somatic_variant_calls"]['indel'])):
			reference[data[i]["consensus_somatic_variant_calls"]['indel'][specimen]['gnos_id']]= [data[i]["icgc_donor_id"] , data[i]['dcc_project_code'],'consensus' ]
	
	except KeyError:
		continue



with open('/Users/morova/Google Drive/TuncProjectStep2-3/reference_dictionaries/folderprefix_reference_dict.p','wb') as fp:
	data = pickle.dump(reference,fp)
	


### TEST ###

if reference['d4fa1f86-ceb6-11e5-bcb8-d4714d100685'][1] != 'PRAD-CA':
	raise ValueError, 'Your folder reference dictionary has problems, Check dictionary or Recreate it !!!'

if reference['ff464a6e-d7a3-4b32-a247-6fd984e162bb'][1] != 'PRAD-UK':
	raise ValueError, 'Your folder reference dictionary has problems, Check dictionary or Recreate it !!!'

print 'All Tests are finished, reference dictionary is working fine!. Good Luck !!!'


#with open('/Volumes/blackhole/reference_dict.p','rb') as fp:
#	data= pickle.load(fp)