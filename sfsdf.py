#!/usr/bin/env python
### this hpc version ###
### deneme deneme ###

#### eger bunu goruyorsam ####
import os
import re
import subprocess
import cPickle as pickle


with open('/Volumes/blackhole/reference_dictionaries/reference_dict.p','rb') as fp:
	fileprefix_reference= pickle.load(fp)
fp.close()

with open('/Volumes/blackhole/reference_dictionaries/FILEPREFIXbridgeCLINICAL.p','rb') as fp:
	bridge= pickle.load(fp)
fp.close()

with open('/Volumes/blackhole/reference_dictionaries/clinical_reference_dict.p','rb') as fp:
	clinical_reference= pickle.load(fp)
fp.close()

with open('/Volumes/blackhole/reference_dictionaries/folderprefix_reference_dict.p','rb') as fp:
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

topdirectory='/Volumes/blackhole/icgc'
os.chdir(topdirectory)

project_names = list(set([folder_prefix[item][1] for item in os.listdir('input/')]))

print '%s  projects detected, single folder for each project is being created' % len(project_names)
print project_names

for project_id in project_names:
	os.mkdir('final/%s' % project_id)


bed_file_input="/Volumes/blackhole/highconf_250wider.txt"
SnpSiftjarpath='/Users/morova/snpEff/SnpSift.jar'
SnpEffjarpath='/Users/morova/snpEff/snpEff.jar'
annotation_file_path='/Volumes/blackhole/common_all_20160601.vcf'
snp_filter_option = "\"(! ID =~ 'rs' ) & (isHet(GEN[TUMOUR]))\""
indel_filter_option = "\"(! ID =~ 'rs' )\""



#directory_path='/Volumes/blackhole/0b140a8c-d2e0-11e5-9065-4a655c6878de'
#filenames=os.listdir(directory_path)
#os.chdir(directory_path)
#snv_file_count=0
#indel_file_count=0

snv_log= open('./final/snv_number.txt','w')
snv_log.write('#### patientid found_mutations|number_of_total_mutations\n')
indel_log =open('./final/indel_number.txt','w')
indel_log.write('#### patientid found_mutations|number_of_total_mutations\n')


for folder in os.listdir('input/'):
	analysis_method = folder_prefix[folder][2]
	os.mkdir('./final/%s/%s' % (folder_prefix[folder][1],folder_prefix[folder][0]))

	output_directory = './final/%s/%s/' % (folder_prefix[folder][1],folder_prefix[folder][0])
	### Bootstrappingi burada yapacagiz ###

	print 'new patient'
	for name in os.listdir('input/%s' % folder):

		if bool(re.search('snv_mnv.vcf.gz$', name)): ### checks if snv_mnv.vcf.gz in the iterated variable
			analysis_id=name.split('.')[0]

			print 'Annotating SNV'
			annotation_script=' '.join(['java', '-jar' , SnpSiftjarpath , 'annotate' , '-id' ,annotation_file_path,  './input/%s/' % folder  + name , '-v' ,'|', 'java', '-Xmx4g','-jar',SnpEffjarpath,'eff', 'GRCh37.75' , '-','-v','|', \
										'java' ,'-jar', SnpSiftjarpath,'filter',snp_filter_option,'>' ,output_directory+'annotated_full_' + analysis_id+'.snv_mnv_'+analysis_method])

			subprocess.check_call(annotation_script,shell=True)
			vcf_isolation_bash_script=["vcftools","--vcf", output_directory+ 'annotated_full_' + analysis_id+'.snv_mnv_'+analysis_method, "--bed" ,  bed_file_input , "--recode" ,"--recode-INFO-all", "--out", output_directory + "final." + analysis_id + ".snv_mnv_" + analysis_method ]

			print 'Cleaving SNV VCF'
			r1=subprocess.check_call(' '.join(vcf_isolation_bash_script),shell=True)

			#intermediate_file= output_directory + "unannotated." + analysis_id + ".snv_mnv_" + analysis_method + '.recode.vcf' ### This file keeps the record of spesific mutations filtered by the bed file

					### This is for reading log files  and creating a single file which contains all of the number of the mutations.###
			with open(output_directory + "final." + analysis_id + ".snv_mnv_" + analysis_method +'.log' , 'r') as logfile:
				for line in logfile:
					if 'possible' and 'After filtering, kept' in line:
						foundmuts = re.search('After filtering, kept (.+?) out of a possible (.+?) Sites', line).group(1)
						total = re.search('After filtering, kept (.+?) out of a possible (.+?) Sites', line).group(2)
						snv_log.write(foundmuts + '|' + total+ '\t')

			logfile.close()



			#p1=subprocess.Popen(["java", "-jar", SnpSiftjarpath, 'annotate','-id', annotation_file_path, intermediate_file], stdout= subprocess.PIPE)
			#p2=subprocess.Popen(["java" , "-Xmx4g", "-jar", SnpEffjarpath, "eff",'GRCh37.75','-'], stdin=p1.stdout,stdout= subprocess.PIPE)
			#final_output = open(output_directory +"final." + analysis_id + ".snv_mnv_" + analysis_method+".vcf", 'w')
			#p3=subprocess.call(["java" ,"-jar",SnpSiftjarpath, "filter", snp_filter_option], stdin=p2.stdout,stdout = final_output)




		elif bool(re.search('indel.vcf.gz$', name)):


			analysis_id=name.split('.')[0]
			print 'Annotating Indel'
			annotation_script=' '.join(['java', '-jar' , SnpSiftjarpath , 'annotate' , '-id' ,annotation_file_path,  './input/%s/' % folder  + name ,'-v', '|', 'java', '-Xmx4g','-jar',SnpEffjarpath,'eff', 'GRCh37.75' , '-','-v','|', \
										'java' ,'-jar', SnpSiftjarpath,'filter',indel_filter_option,'>' ,output_directory+'annotated_full_' + analysis_id+'.indel_'+analysis_method])

			subprocess.check_call(annotation_script,shell=True)
			vcf_isolation_bash_script=["vcftools","--gzvcf", output_directory+ 'annotated_full_' + analysis_id+'.snv_mnv_'+analysis_method, "--bed" ,  bed_file_input , "--recode" ,"--recode-INFO-all", "--out", output_directory + "final." + analysis_id + ".indel_" + analysis_method ]

			print 'Cleaving Indel VCF'
			r1=subprocess.check_call(' '.join(vcf_isolation_bash_script),shell=True)
			#subprocess.call(annotation_script,shell=True)

			with open(output_directory + "final." + analysis_id + ".indel_" + analysis_method +'.log' , 'r') as logfile:
				for line in logfile:
					if 'possible' and 'After filtering, kept' in line:
						foundmuts = re.search('After filtering, kept (.+?) out of a possible (.+?) Sites', line).group(1)
						total = re.search('After filtering, kept (.+?) out of a possible (.+?) Sites', line).group(2)
						indel_log.write(foundmuts + '|' + total+ '\t')

			logfile.close()

	snv_log('\n')
	indel_log('\n')


snv_log.close()
indel_log.close()