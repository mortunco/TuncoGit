##!/usr/bin/env python
### this hpc version ###
### v.01 ###

#### eger bunu goruyorsam ####
import os
import re
import subprocess
import cPickle as pickle


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

topdirectory='/mnt/kufs/scratch/tmorova15/vcf_analysis'
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

broad_snp_filter_option = "\"(! ID =~ 'rs' )\""
sanger_snp_filter_option= "\"(! exists SNP ) & (HE='1')\""
dkfz_snp_filter_option ="\"(! ID =~ 'rs') & (HE='1') \""

sanger_indel_filter_option = "\"(! ID =~ 'rs' )\""
dkfz_indel_filter_option = "\"(! ID =~ 'rs' ) & (HE='1')\""
broad_indel_filter_option=  "\"(! ID =~ 'rs' )\""






snv_log= open('./final/snv_number.txt','w')
snv_log.write('#### patientid found_mutations|number_of_total_mutations\n')
indel_log =open('./final/indel_number.txt','w')
indel_log.write('#### patientid found_mutations|number_of_total_mutations\n')

os.mkdir('./final/Problematic_id')

for folder in os.listdir('input/'):

	if folder not in folder_prefix.keys():
		os.mkdir('./final/Problematic_id/%s' % folder)



	else:

		### Python gives a file already existed error, so this line is for to prevent that error.
		#### If a patient folder is created for another caller, code does not create a new dir and write to the currently existed dir
		if os.path.exists('./final/%s/%s' % (folder_prefix[folder][1],folder_prefix[folder][0])):
			pass
		else:
			os.mkdir('./final/%s/%s' % (folder_prefix[folder][1],folder_prefix[folder][0]))

		### analysis method : sanger, dkfz, broad
		analysis_method = folder_prefix[folder][2]
		snv_log.write(folder_prefix[folder][0] + '\t') ## write to the file which contains all of the numbers

		output_directory = './final/%s/%s/' % (folder_prefix[folder][1],folder_prefix[folder][0])
		### Bootstrappingi burada yapacagiz ###

		for name in os.listdir('input/%s' % folder):

			if bool(re.search('somatic.snv_mnv.vcf.gz$', name)): ### checks if snv_mnv.vcf.gz in the iterated variable
				analysis_id=name.split('.')[0]


				print "Annotating %s SNV" % folder_prefix[folder][0]

				if analysis_method == 'sanger':

					annotation_script =' '.join(['java','-jar',SnpSiftjarpath ,'gt', './input/%s/' % folder  + name , '|' ,'java','-jar',SnpSiftjarpath,'filter',sanger_snp_filter_option, '>', output_directory+'annotated_full_' + analysis_id+'.snv_mnv_'+analysis_method])

				elif analysis_method == 'dkfz':

					annotation_script = ' '.join(['java','-jar',SnpSiftjarpath ,'gt', './input/%s/' % folder  + name , '|' ,'java','-jar',SnpSiftjarpath,'filter',dkfz_snp_filter_option, '>' ,output_directory+'annotated_full_' + analysis_id+'.snv_mnv_'+analysis_method])

				#elif belki muse eklenir buraya:

				else:
					### broad icin ###
					annotation_script=' '.join(['java', '-jar' , SnpSiftjarpath , 'annotate' , '-id' ,annotation_file_path,  './input/%s/' % folder  + name , '|', 'java', '-Xmx4g','-jar',SnpEffjarpath,'eff', 'GRCh37.75' , '-','|', \
										'java' ,'-jar', SnpSiftjarpath,'filter',broad_snp_filter_option,'>' ,output_directory+'annotated_full_' + analysis_id+'.snv_mnv_'+analysis_method])


				subprocess.check_call(annotation_script, shell=True)


				print "Cleaving spesific regions given by %s" % bed_file_input
				vcf_isolation_bash_script=["/mnt/kufs/scratch/tmorova15/softwares/vcftools_0.1.13/bin/vcftools","--gzvcf", output_directory+ 'annotated_full_' + analysis_id+'.snv_mnv_'+analysis_method, "--bed" ,  bed_file_input , "--recode" ,"--recode-INFO-all", "--out", output_directory + "final." + analysis_id + ".snv_mnv_" + analysis_method ]
				r1=subprocess.check_call(' '.join(vcf_isolation_bash_script),shell=True)


				### This is for reading log files  and creating a single file which contains all of the number of the mutations.###
				with open(output_directory + "final." + analysis_id + ".snv_mnv_" + analysis_method +'.log' , 'r') as logfile:
					for line in logfile:
						line = line.rstrip('\n')
						if ('After filtering' in line) and ('possible' in line):

							foundmuts = re.search('After filtering, kept (.+?) out of a possible (.+?) Sites', line).group(1)
							total = re.search('After filtering, kept (.+?) out of a possible (.+?) Sites', line).group(2)
							snv_log.write(foundmuts + '|' + total+ '\t')

				logfile.close()





				#p1=subprocess.Popen(["java", "-jar", SnpSiftjarpath, 'annotate','-id', annotation_file_path, intermediate_file], stdout= subprocess.PIPE)
				#p2=subprocess.Popen(["java" , "-Xmx4g", "-jar", SnpEffjarpath, "eff",'GRCh37.75','-'], stdin=p1.stdout,stdout= subprocess.PIPE)
				#final_output = open(output_directory +"final." + analysis_id + ".snv_mnv_" + analysis_method+".vcf", 'w')
				#p3=subprocess.call(["java" ,"-jar",SnpSiftjarpath, "filter", snp_filter_option], stdin=p2.stdout,stdout = final_output)




			elif bool(re.search('somatic.indel.vcf.gz$', name)):
				analysis_id=name.split('.')[0]


				print "Annotating %s Indel" % folder_prefix[folder][0]

				if analysis_method == 'sanger':

					annotation_script=' '.join(['java', '-jar' , SnpSiftjarpath , 'annotate' , '-id' ,annotation_file_path,  './input/%s/' % folder  + name , '|', 'java', '-Xmx4g','-jar',SnpEffjarpath,'eff', 'GRCh37.75' , '-','|', \
										'java' ,'-jar', SnpSiftjarpath,'filter',sanger_indel_filter_option,'>' ,output_directory+'annotated_full_' + analysis_id+'.indel_'+analysis_method])


				elif analysis_method == 'dkfz':

					annotation_script = ' '.join(['java','-jar',SnpSiftjarpath ,'gt', './input/%s/' % folder  + name , '|' ,'java','-jar',SnpSiftjarpath,'filter',dkfz_indel_filter_option, '>' ,output_directory+'annotated_full_' + analysis_id+'.indel_'+analysis_method])


				#elif belki muse eklenir buraya:

				else:

					annotation_script=' '.join(['java', '-jar' , SnpSiftjarpath , 'annotate' , '-id' ,annotation_file_path,  './input/%s/' % folder  + name , '|', 'java', '-Xmx4g','-jar',SnpEffjarpath,'eff', 'GRCh37.75' , '-','|', \
										'java' ,'-jar', SnpSiftjarpath,'filter',broad_indel_filter_option,'>' ,output_directory+'annotated_full_' + analysis_id+'.indel_'+analysis_method])




				subprocess.check_call(annotation_script,shell=True)


				print "Cleaving spesific regions given by %s" % bed_file_input
				vcf_isolation_bash_script=["/mnt/kufs/scratch/tmorova15/softwares/vcftools_0.1.13/bin/vcftools","--gzvcf", output_directory+ 'annotated_full_' + analysis_id+'.indel_'+analysis_method, "--bed" ,  bed_file_input , "--recode" ,"--recode-INFO-all", "--out", output_directory + "final." + analysis_id + ".indel_" + analysis_method ]
				r1=subprocess.check_call(' '.join(vcf_isolation_bash_script),shell=True)

				with open(output_directory + "final." + analysis_id + ".indel_" + analysis_method +'.log' , 'r') as logfile:
					for line in logfile:
						line = line.rstrip('\n')
						if ('After filtering' in line) and ('possible' in line):

							foundmuts = re.search('After filtering, kept (.+?) out of a possible (.+?) Sites', line).group(1)
							total = re.search('After filtering, kept (.+?) out of a possible (.+?) Sites', line).group(2)
							indel_log.write(foundmuts + '|' + total+ '\t')

				logfile.close()


			




	snv_log.write(analysis_method+'-SNV'+'\n')
	indel_log.write(analysis_method + '-Indel' + '\n')


snv_log.close()
indel_log.close()