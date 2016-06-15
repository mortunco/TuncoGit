##!/usr/bin/env python
### this hpc version ###
### v.01 ###

#### eger bunu goruyorsam ####
import os
import re
import subprocess
import cPickle as pickle


with open('/mnt/home/tmorova15/RS/references/reference_dictionaries/reference_dict.p','rb') as fp:
	fileprefix_reference= pickle.load(fp)
fp.close()

with open('/mnt/home/tmorova15/RS/references/reference_dictionaries/FILEPREFIXbridgeCLINICAL.p','rb') as fp:
	bridge= pickle.load(fp)
fp.close()

with open('/mnt/home/tmorova15/RS/references/reference_dictionaries/clinical_reference_dict.p','rb') as fp:
	clinical_reference= pickle.load(fp)
fp.close()

with open('/mnt/home/tmorova15/RS/references/reference_dictionaries/folderprefix_reference_dict.p','rb') as fp:
	folder_prefix=pickle.load(fp)
fp.close()

print 'Initializing Reference Check...'


print "------------------------------------------------------------"
print 'Initializing Reference Check...'

print 'Checking File Prefix Reference'

if fileprefix_reference['a9ec7d9e-b179-4782-a589-43c7d1642be9'][0] != 'DO51159':
	raise ValueError , 'Your File Prefix Dictionary has problems, Check dictionary or Recreate it!!!'

if bridge['a9ec7d9e-b179-4782-a589-43c7d1642be9'][1] != 'DO51159':
	raise ValueError , 'Your Bridge Dictionary has problems, Check dictionary or Recreate it !!!'

if clinical_reference['SP113005'][0] != 'DO51159':
	print 'clinical reference dictionary says:%s  but it is supposed to be DO51159' % clinical_reference['SP113005']
	raise ValueError, 'Your clinical reference dictionary has problems, Check dictionary or Recreate it !!!'

if folder_prefix['d4fa1f86-ceb6-11e5-bcb8-d4714d100685'][1] != 'PRAD-CA':
	raise ValueError, 'Your folder reference dictionary has problems, Check dictionary or Recreate it !!!'

topdirectory='/mnt/home/tmorova15/RS/vcf_analysis'
os.chdir(topdirectory)

project_names = list(set([folder_prefix[item][1] for item in os.listdir('input/')]))

print '%s  projects detected, single folder for each project is being created' % len(project_names)
print project_names

for project_id in project_names:
	os.mkdir('final/%s' % project_id)


bed_file_input="/mnt/kufs/scratch/tmorova15/references/highconf_250wider.txt"
SnpSiftjarpath='/mnt/kufs/scratch/tmorova15/softwares/snpEff/SnpSift.jar'
SnpEffjarpath='/mnt/kufs/scratch/tmorova15/softwares/snpEff/snpEff.jar'
annotation_file_path='/mnt/kufs/scratch/tmorova15/references/common_all_20160601.vcf'
snp_filter_option = "\"(! ID =~ 'rs' ) & (isHet(GEN[TUMOUR]))\""
indel_filter_option = "\"(! ID =~ 'rs' )\""







os.mkdir('./final/Problematic_id')

for folder in os.listdir('input/'):

	if folder not in folder_prefix.keys():
		os.mkdir('./final/Problematic_id/%s' % folder)

	else:
		analysis_method = folder_prefix[folder][2]
		os.mkdir('./final/%s/%s' % (folder_prefix[folder][1],folder_prefix[folder][0]))

		output_directory = './final/%s/%s/' % (folder_prefix[folder][1],folder_prefix[folder][0])
		### Bootstrappingi burada yapacagiz ###

		for name in os.listdir('input/%s' % folder):

			if bool(re.search('snv_mnv.vcf.gz$', name)): ### checks if snv_mnv.vcf.gz in the iterated variable
				analysis_id=name.split('.')[0]

				vcf_isolation_bash_script=["vcftools","--gzvcf", './input/%s/' % folder  + name , "--bed" ,  bed_file_input , "--recode" ,"--recode-INFO-all", "--out", output_directory + "unannotated." + analysis_id + ".snv_mnv_" + analysis_method ]
				r1=subprocess.call(vcf_isolation_bash_script)

				intermediate_file= output_directory + "unannotated." + analysis_id + ".snv_mnv_" + analysis_method + '.recode.vcf' ### This file keeps the record of spesific mutations filtered by the bed file

				annotation_script=' '.join(['java', '-jar' , SnpSiftjarpath , 'annotate' , '-id' ,annotation_file_path, intermediate_file, '|', 'java', '-Xmx4g','-jar',SnpEffjarpath,'eff', 'GRCh37.75' , '-','|', \
											'java' ,'-jar', SnpSiftjarpath,'filter',snp_filter_option,'>' ,output_directory+'final.' + analysis_id+'.snv_mnv_'+analysis_method+ '.vcf'])

				subprocess.call(annotation_script,shell=True)



				#p1=subprocess.Popen(["java", "-jar", SnpSiftjarpath, 'annotate','-id', annotation_file_path, intermediate_file], stdout= subprocess.PIPE)
				#p2=subprocess.Popen(["java" , "-Xmx4g", "-jar", SnpEffjarpath, "eff",'GRCh37.75','-'], stdin=p1.stdout,stdout= subprocess.PIPE)
				#final_output = open(output_directory +"final." + analysis_id + ".snv_mnv_" + analysis_method+".vcf", 'w')
				#p3=subprocess.call(["java" ,"-jar",SnpSiftjarpath, "filter", snp_filter_option], stdin=p2.stdout,stdout = final_output)




			elif bool(re.search('indel.vcf.gz$', name)):


				analysis_id=name.split('.')[0]

				vcf_isolation_bash_script=["vcftools","--gzvcf", './input/%s/' % folder  + name , "--bed" ,  bed_file_input , "--recode" ,"--recode-INFO-all", "--out", output_directory + "unannotated." + analysis_id + ".indel_" + analysis_method ]
				r1=subprocess.call(vcf_isolation_bash_script)

				intermediate_file= output_directory+"unannotated." + analysis_id + ".indel_" + analysis_method + ".recode.vcf"

				annotation_script=' '.join(['java', '-jar' , SnpSiftjarpath , 'annotate' , '-id' ,annotation_file_path, intermediate_file, '|', 'java', '-Xmx4g','-jar',SnpEffjarpath,'eff', 'GRCh37.75' , '-','|', \
											'java' ,'-jar', SnpSiftjarpath,'filter',indel_filter_option,'>' ,output_directory+'final.' + analysis_id+'.indel_'+analysis_method+ '.vcf'])

				subprocess.call(annotation_script,shell=True)

