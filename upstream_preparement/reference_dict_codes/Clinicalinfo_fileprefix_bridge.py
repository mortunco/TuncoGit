#Is built for identfying which fileprefix belongs to which specimen. Some of the patients may have multiple specimen variant called thats why
# this dictionary is required #
import cPickle as pickle


metafilepath='/Users/morova/Google Drive/ICGC/pcawg_sample_sheet.tsv.txt'

ICGC_metadict={}

with open(metafilepath,'r') as f3:
	for line in f3:
		temp = line.rstrip('\n').split('\t')

		fileprefix, donor_id,project_name, icgc_sample_id= temp[4], temp[2], temp[3], temp[6]

		ICGC_metadict[fileprefix] = [icgc_sample_id,donor_id,project_name]



f3.close()

with open('/Volumes/blackhole/FILEPREFIXbridgeCLINICAL.p','wb') as fp:
	print 'Clinical Reference is written...'
	data = pickle.dump(ICGC_metadict,fp)

fp.close()