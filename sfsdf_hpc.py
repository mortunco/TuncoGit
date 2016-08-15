##! /usr/bin/env python
### this hpc version ###
### v.03 ### after paralelization ###

#### eger bunu goruyorsam ####
import os
import re
import subprocess
import cPickle as pickle
import time
import sys
import argparse
parser = argparse.ArgumentParser(description='This code runs annotation of VCFs and filtration of the annotated file based on a given location file')
parser.add_argument('-p',type=str, help='path of the top folder which contains input and final directories.')
args=parser.parse_args()

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

topdirectory=args.p

print topdirectory



#topdirectory='/mnt/kufs/scratch/tmorova15/vcf_analysis'
os.chdir(topdirectory)
print os.getcwd()
#project_names = list(set([folder_prefix[item][1] for item in os.listdir('input/')]))

project_names =[]

for item in os.listdir('input/'):
	if item in folder_prefix.keys():
		project_names.append(folder_prefix[item][1])
	else:
		pass

project_names=list(set(project_names))


print '%s  projects detected, single folder for each project is being created' % len(project_names)
print project_names

for project_id in project_names:
	os.mkdir('final/%s' % project_id)


bed_file_input="/mnt/kufs/scratch/tmorova15/references/highconf_250wider.txt"
chromsizes='/mnt/kufs/scratch/tmorova15/references/GRCh37.genome.txt'
bedtoolspath='/mnt/kufs/scratch/tmorova15/bcbio/bin/bedtools'
randombedcount=1000
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

snv_log.flush()
indel_log.flush()

os.mkdir('./final/Problematic_id')
repeat_cheacker=[]
for folder in os.listdir('input/'):


	if len(folder.split('-')) == 5: ### There are other folders/files located in the same input folder such as log.txt of the download or manifests. This wi


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
			snv_log.write(folder_prefix[folder][0] + '\t' + analysis_method +'\t') ## write to the file which contains all of the numbers
			indel_log.write(folder_prefix[folder][0] + '\t' + analysis_method +'\t')## write to the file which contains all of the numbers

			output_directory = './final/%s/%s/' % (folder_prefix[folder][1],folder_prefix[folder][0])

			start=time.time()
			if not os.path.exists('./final/randombeds'): #creates random bed directory
				os.mkdir('./final/randombeds')
				### Random bed generation if not exists ###
				print 'Genarating Random Bed Files...'
				for i in range(randombedcount):
					random_commandline=' '.join([bedtoolspath,'shuffle', '-chrom','-i',bed_file_input,'-g',chromsizes,'>', './final/randombeds/random_' + '%s' % str(i)])
					subprocess.check_call(random_commandline,shell= True)

			print 'Took', time.time() - start,'to run....'

			### Random Vcf Genaration ###
			if not os.path.exists(output_directory + 'random_vcfs'):
				os.mkdir(output_directory + 'random_vcfs')


			for name in os.listdir('input/%s' % folder):

				if bool(re.search('somatic.snv_mnv.vcf.gz$', name)): ### checks if snv_mnv.vcf.gz in the iterated variable
					analysis_id=name.split('.')[0]

					start=time.time()
					print "**Annotating %s SNV" % folder_prefix[folder][0]

					if analysis_method == 'sanger':

						annotation_script =' '.join(['java','-jar',SnpSiftjarpath ,'gt', './input/%s/' % folder  + name , '|' ,'java','-jar',SnpSiftjarpath,'filter',sanger_snp_filter_option, '>', output_directory+'annotated_full.' + analysis_id+'.snv_mnv_'+analysis_method])

					elif analysis_method == 'dkfz':

						annotation_script = ' '.join(['java','-jar',SnpSiftjarpath ,'gt', './input/%s/' % folder  + name , '|' ,'java','-jar',SnpSiftjarpath,'filter',dkfz_snp_filter_option, '>' ,output_directory+'annotated_full.' + analysis_id+'.snv_mnv_'+analysis_method])

					#elif belki muse eklenir buraya:

					else:
						### broad icin ###
						annotation_script=' '.join(['java', '-jar' , SnpSiftjarpath , 'annotate' , '-id' ,annotation_file_path,  './input/%s/' % folder  + name , '|', 'java', '-Xmx6g','-jar',SnpEffjarpath,'eff', 'GRCh37.75' , '-','|', \
											'java' ,'-jar', SnpSiftjarpath,'filter',broad_snp_filter_option,'>' ,output_directory+'annotated_full.' + analysis_id+'.snv_mnv_'+analysis_method])


					subprocess.check_call(annotation_script, shell=True)

					print 'Annotation of SNV Took', time.time() - start,'to run....'

					print "***Cleaving spesific regions given by %s" % bed_file_input
					start=time.time()
					#vcf_isolation_bash_script=["/mnt/kufs/scratch/tmorova15/softwares/vcftools_0.1.13/bin/vcftools","--gzvcf", output_directory+ 'annotated_full.' + analysis_id+'.snv_mnv_'+analysis_method, "--bed" ,  bed_file_input , \
					#						   "--recode" ,"--recode-INFO-all", "--out", output_directory + "final." + analysis_id + ".snv_mnv_" + analysis_method]

					### takes only the mutations within given bed region ###
					vcf_isolation_bash_script =['java','-jar',SnpSiftjarpath,'interval', '-i' , output_directory+ 'annotated_full.' + analysis_id+'.snv_mnv_'+analysis_method, bed_file_input, '>',  output_directory + "final." + analysis_id + ".snv_mnv_" + analysis_method + '.vcf']
					subprocess.check_call(' '.join(vcf_isolation_bash_script),shell=True)

					### counts the mutations
					mutation_count_line=[bedtoolspath,'intersect', '-a' ,  output_directory + 'annotated_full.' + analysis_id + ".snv_mnv_" + analysis_method ,'-b',bed_file_input, '|' ,'wc']

					### count is the count of mutation counts around given bed regions ###
					count = subprocess.check_output(' '.join(mutation_count_line),shell = True)

					### Gives the total mutations from the vcf ###
					total_mutation_count_line=['egrep','-v','^#',  output_directory+ 'annotated_full.' + analysis_id+'.snv_mnv_'+analysis_method , '|', 'wc']
					total = subprocess.check_output(' '.join(total_mutation_count_line),shell=True)

					count=count.split()
					total=total.split()
					print count,total,'snv'
					print 'Cleavage of SNV Took', time.time() - start,'to run....'


					print count[0], 'SNV mutations found in', total[0]

					### This is for reading log files  and creating a single file which contains all of the number of the mutations.###

					snv_log.write(str(count[0]) + '|' + str(total[0]) + '\t')



					### bootstrap ###
					start=time.time()
					random_counter=0
					for n,randombed in enumerate(os.listdir('./final/randombeds')):
						sys.stdout.write("\r{0}".format(randombed))
						sys.stdout.flush()
						random_isolation_bash_script = ['java','-jar',SnpSiftjarpath,'interval', '-i' , output_directory+ 'annotated_full.' + analysis_id+'.snv_mnv_'+analysis_method, './final/randombeds/'+ randombed, '|',\
						bedtoolspath,'intersect', '-a' ,"-",'-b',bed_file_input, '|' ,'wc']

						random_mutation_count=subprocess.check_output(' '.join(random_isolation_bash_script),shell=True)
						if count < random_mutation_count:
							random_counter += 1
					snv_log.write(str(random_counter/float(randombedcount)) + '\t')

					snv_log.flush()
					print 'Bootstraping of SNV Took', time.time() - start,'to run....'
				elif bool(re.search('somatic.indel.vcf.gz$', name)):
					analysis_id=name.split('.')[0]
					start=time.time()

					print "**Annotating %s Indel" % folder_prefix[folder][0]

					if analysis_method == 'sanger':

						annotation_script=' '.join(['java', '-jar' , SnpSiftjarpath , 'annotate' , '-id' ,annotation_file_path,  './input/%s/' % folder  + name , '|', 'java', '-Xmx6g','-jar',SnpEffjarpath,'eff', 'GRCh37.75' , '-','|', \
											'java' ,'-jar', SnpSiftjarpath,'filter',sanger_indel_filter_option,'>' ,output_directory+'annotated_full.' + analysis_id+'.indel_'+analysis_method])


					elif analysis_method == 'dkfz':

						annotation_script = ' '.join(['java','-jar',SnpSiftjarpath ,'gt', './input/%s/' % folder  + name , '|' ,'java','-jar',SnpSiftjarpath,'filter',dkfz_indel_filter_option, '>' ,output_directory+'annotated_full.' + analysis_id+'.indel_'+analysis_method])


					#elif belki muse eklenir buraya:

					else:

						annotation_script=' '.join(['java', '-jar' , SnpSiftjarpath , 'annotate' , '-id' ,annotation_file_path,  './input/%s/' % folder  + name , '|', 'java', '-Xmx6g','-jar',SnpEffjarpath,'eff', 'GRCh37.75' , '-','|', \
											'java' ,'-jar', SnpSiftjarpath,'filter',broad_indel_filter_option,'>' ,output_directory+ 'annotated_full.' + analysis_id+'.indel_'+analysis_method])


					print 'Tunc you are so good man. So good....'

					subprocess.check_call(annotation_script,shell=True)
					print 'Annotation of Indel Took', time.time() - start,'to run....'

					print "***Cleaving spesific regions given by %s" % bed_file_input
					# vcf_isolation_bash_script=["/mnt/kufs/scratch/tmorova15/softwares/vcftools_0.1.13/bin/vcftools","--gzvcf", output_directory+ 'annotated_full.' + analysis_id+'.indel_'+analysis_method, "--bed" ,  bed_file_input , "--recode" ,"--recode-INFO-all", "--out", output_directory + "final." + analysis_id + ".indel_" + analysis_method ]
					# r1=subprocess.check_call(' '.join(vcf_isolation_bash_script),shell=True)

					vcf_isolation_bash_script =['java','-jar',SnpSiftjarpath,'interval', '-i' , output_directory+ 'annotated_full.' + analysis_id+'.indel_'+analysis_method, bed_file_input, '>',  output_directory + "final." + analysis_id + ".indel_" + analysis_method + '.vcf']
					subprocess.check_call(' '.join(vcf_isolation_bash_script),shell=True)

					### counts the mutations
					mutation_count_line=[bedtoolspath,'intersect', '-a' ,  output_directory + 'annotated_full.' + analysis_id+'.indel_'+analysis_method,'-b',bed_file_input, '|' ,'wc']

					### count is the count of mutation counts around given bed regions ###
					count = subprocess.check_output(' '.join(mutation_count_line),shell = True)

					### Gives the total mutations from the vcf ###
					total_mutation_count_line=['egrep','-v','^#', output_directory+ 'annotated_full.' + analysis_id+'.indel_'+analysis_method, '|', 'wc']
					total = subprocess.check_output(' '.join(total_mutation_count_line),shell=True)

					count=count.split()
					total=total.split()
					print count[0], 'Indel mutations found in', total[0]
					indel_log.write(str(count[0]) + '|' + str(total[0]) + '\t')
					print 'Cleavage of Indel Took', time.time() - start,'to run....'

					start=time.time()
					random_counter=0
					for n,randombed in enumerate(os.listdir('./final/randombeds')):
						sys.stdout.write("\r{0}".format(randombed))
						sys.stdout.flush()
						random_isolation_bash_script = ['java','-jar',SnpSiftjarpath,'interval', '-i' , output_directory+ 'annotated_full.' + analysis_id+'.indel_'+analysis_method, './final/randombeds/'+ randombed, '|',\
						bedtoolspath,'intersect', '-a' ,"-",'-b',bed_file_input, '|' ,'wc']

						random_mutation_count=subprocess.check_output(' '.join(random_isolation_bash_script),shell=True)
						if count < random_mutation_count:
							random_counter += 1

					indel_log.write(str(random_counter/float(randombedcount)) + '\t')
					indel_log.flush()

					print 'Bootstraping of Indel Took', time.time() - start,'to run....'







		snv_log.write('\n')
		indel_log.write('\n')


snv_log.close()
indel_log.close()

print 'Process has finalized without any error ! Well Done !!'