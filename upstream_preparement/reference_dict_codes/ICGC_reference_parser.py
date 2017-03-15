### This is the reference that gives information about which file prefix is corresponded to which donor ####

import cPickle as pickle


f3= open('/Users/morova/Google Drive/ICGC/release_may2016.v1.tsv.txt','r')

reference={}
for line in f3:
	temp = line.rstrip().split('\t')
	project_name, donor_id, sanger, dkfz, broad, muse = temp[1] , temp[3], temp[39], temp[44], temp[49], temp[54]

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


with open('/Volumes/blackhole/reference_dict.p','wb') as fp:
	data = pickle.dump(reference,fp)


### This
#with open('/Volumes/blackhole/reference_dict.p','rb') as fp:
#	data= pickle.load(fp)

