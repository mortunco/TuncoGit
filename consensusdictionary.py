import json
import cPickle as pickle

data=[]
with open('/Users/morova/Downloads/release_may2016.v1.4.with_consensus_calls.jsonl.txt') as jsonldata:

	for line in jsonldata:
		data.append(json.loads(line)) ### every line is a one json. 


consensusdict={}
for i in range(len(data)):
	try:
		for specimen in xrange(len(data[i]["consensus_somatic_variant_calls"]['snv_mnv'])):
			consensusdict[data[i]["consensus_somatic_variant_calls"]['snv_mnv'][specimen]['gnos_id']]= [data[i]["icgc_donor_id"], data[i]['dcc_project_code'],'consensus' ]
		
		for specimen in xrange(len(data[i]["consensus_somatic_variant_calls"]['indel'])):
			consensusdict[data[i]["consensus_somatic_variant_calls"]['indel'][specimen]['gnos_id']]= [data[i]["icgc_donor_id"] , data[i]['dcc_project_code'],'consensus' ]
	
	except KeyError:
		continue


print 