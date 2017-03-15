#!/usr/bin/env python

import os
import re


os.mkdir('{0}_final_snvlog_from_all_batches'.format(args.name))
os.mkdir('{0}_variant_classification_snv'.format(args.name))
os.mkdir('{0}_final_indellog_from_all_batches'.format(args.name))
os.mkdir('{0}_variant_classification_indel'.format(args.name))


 
 os.chwd('somedirectory')
 
 for directory in os.listdir(os.getcwd()):
 
	if bool(re.search('*final_snvlog_from_all_bathces', directory)):
		
		for filename in os.listdir(os.getcwd() + '/' +directory ):
			if not bool(re.search('snv_numberBatch*.txt',filename)):
				with open(os.getcwd()+'/'+directory+'/'+filename,'rb') as snvlogfile:
					 snvlog=[i.rstrip('\n').split('\t') for i in snvlogfile]
						
			
		
 