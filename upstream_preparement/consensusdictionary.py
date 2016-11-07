#! /usr/bin/env python
import json
import cPickle as pickle

data=[]
with open('/Users/morova/Google Drive/TuncProjectStep2-3/ICGC/release_may2016.v1.4.with_consensus_calls.jsonl.txt') as jsonldata:

	for line in jsonldata:
		data.append(json.loads(line)) ### every line is a one json.
		

PATHHH='/Volumes/blackhole'
with open('{0}/reference_dictionaries/reference_dict.p'.format(PATHHH),'rb') as fp:
	fileprefix_reference= pickle.load(fp)
fp.close()

with open('{0}/reference_dictionaries/FILEPREFIXbridgeCLINICAL.p'.format(PATHHH),'rb') as fp:
	bridge= pickle.load(fp)
fp.close()

with open('{0}/reference_dictionaries/clinical_reference_dict.p'.format(PATHHH),'rb') as fp:
	clinical_reference= pickle.load(fp)
fp.close()

with open('{0}/reference_dictionaries/folderprefix_reference_dict.p'.format(PATHHH),'rb') as fp:
	folder_prefix=pickle.load(fp)
fp.close()


print data[5]["consensus_somatic_variant_calls"]['snv_mnv'][0]['specimen_type']
consensusdict={}
for i in range(len(data)):
	try:
		for specimen in xrange(len(data[i]["consensus_somatic_variant_calls"]['snv_mnv'])):
			consensusdict[data[i]["consensus_somatic_variant_calls"]['snv_mnv'][specimen]['aliquot_id']]= [data[i]["icgc_donor_id"], data[i]['dcc_project_code'],'consensus' ,data[5]["consensus_somatic_variant_calls"]['snv_mnv'][0]['specimen_type']]
		
		for specimen in xrange(len(data[i]["consensus_somatic_variant_calls"]['indel'])):
			consensusdict[data[i]["consensus_somatic_variant_calls"]['indel'][specimen]['aliquot_id']]= [data[i]["icgc_donor_id"] , data[i]['dcc_project_code'],'consensus' ,data[5]["consensus_somatic_variant_calls"]['snv_mnv'][0]['specimen_type']]
	
	except KeyError:
		continue


print consensusdict['3e7ccab5-5b1d-4147-b907-77cab8f0837e']

with open('/Users/morova/Google Drive/TuncProjectStep2-3/reference_dictionaries/consensusclinical_dictionary.pc','wb') as fp:
 	data = pickle.dump(consensusdict,fp)





# consensusdict={}
# for i in range(len(data)):
# 	try:
# 		
# 		consensusdict[data[i]["icgc_donor_id"] ]= [ data[i]["consensus_somatic_variant_calls"]['snv_mnv'][0]['icgc_specimen_id'] , data[i]['dcc_project_code'],'consensus' ]
# 		
# 	except KeyError:
# 		continue
# print consensusdict['DO51159']
# 
# with open('/Users/morova/Google Drive/TuncProjectStep2-3/reference_dictionaries/consensus_dictionary.pc','wb') as fp:
# 	data = pickle.dump(consensusdict,fp)
# 	
# 	
# 
# # o=0
# # for i in range(len(data)):
# # 	try:
# # 		a = data[i]["consensus_somatic_variant_calls"]['snv_mnv'][0]['icgc_specimen_id']
# # 		print clinical_reference[a]
# # 		o += 1
# # 	except:
# # 		KeyError
# # 		print 'siktigim'
# # 		
# # print o, len(data) 
# 		
# 
# 
# 
