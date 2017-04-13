##! /usr/bin/env python
### this hpc version ###
### STABLE VERSION ###


import os
import re
import subprocess
import cPickle as pickle
import time
import sys
import argparse
import collections
parser = argparse.ArgumentParser(description='This code runs annotation of VCFs and filtration of the annotated file based on a given location file')
parser.add_argument('-p',type=str, help='path of the top folder which contains input and final directories.')
parser.add_argument('-l',type=str, help='bed formatted full file path that tell in which locations we should survey.')
parser.add_argument('-it',type=int,help='Number of iterations will be run by bootstrapping step. Default = 1000', default=1000)
args=parser.parse_args()

def whattimeisit():
	whattimeistit=time.localtime()
	days=['Mon','Tue','Wed','Thu','Fri','Sat','Sun']
	return '{0}/{1}/{2} - {6} - {3}:{4}:{5}'.format(time.localtime()[0],time.localtime()[1],time.localtime()[2],time.localtime()[3],time.localtime()[4],time.localtime()[5],days[time.localtime()[6]])

def variantclassificationcounter(filename,textfilename):
	script=' '.join(['cat',filename,'|', 'grep','-o','Variant_Classification=.*','|' 'cut','-f','2','-d','='])
	output=subprocess.check_output(script,shell=True).split('\n')
	counter=collections.Counter(output)
	for key,value in zip(counter.keys(),counter.values()): textfilename.write(str(key)+'|'+str(value)+';')
	textfilename.write('\n')
	textfilename.flush()
	return

for item in vars(args): ###
	if getattr(args,item) == None:
		raise ValueError , "All arguments must be set"


def BEDbootstrap(logfilename,SNV,verbose):		
	'''This code handles bootstrapping of the spesific mutation. It returns False Discovery Rate probability by comparing with random locations in the genome
	It requires relative location of the randombed files from the top folder. logfilename is to determine in which file it will write the output. Type of mutation that will determine which annotated vcf file.(.snv_mnv_ or .indel_)
	and finally verbose option to for the debugging. Verbose option write the number of the mutation found in each random bed file.'''
	random_counter=0
	
	if SNV == True:
		mutation_type='.snv_mnv_'
	else:
		mutation_type='.indel_'
	
	for n,randombed in enumerate(os.listdir('./final/randombeds')):
							
		sys.stdout.write("\r{0}".format(randombed))
		sys.stdout.flush()
		random_isolation_bash_script = ['java','-jar',SnpSiftjarpath,'interval', '-i' , output_directory+ 'annotated_full.' + analysis_id+ mutation_type +analysis_method, './final/randombeds/'+ randombed, '|',\
		bedtoolspath,'intersect', '-a' ,"-",'-b','./final/randombeds/'+ randombed, '|' ,'wc']
		#print random_isolation_bash_script

		random_mutation_count=subprocess.check_output(' '.join(random_isolation_bash_script),shell=True)
		random_mutation_count=int(random_mutation_count.rstrip().split()[0])
		if int(count[0]) <= int(random_mutation_count):
			#print 'snv'
			#print 'count',count[0]
			#print 'random_no', random_mutation_count[0]
			random_counter += 1
		if verbose == True: logfilename.write(str(random_mutation_count) + '\t')
	
	return random_counter



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
print time.time()

print whattimeisit(),'Checking File Prefix Reference'

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
print whattimeisit(),os.getcwd()
#project_names = list(set([folder_prefix[item][1] for item in os.listdir('input/')]))

project_names =[]

for item in os.listdir('input/'):
	if item in folder_prefix.keys():
		project_names.append(folder_prefix[item][1])
	else:
		pass

project_names=list(set(project_names))


print whattimeisit(),'%s  projects detected, single folder for each project is being created' % len(project_names)
print whattimeisit(),project_names

for project_id in project_names:
	if not os.path.exists('final/%s' % project_id):
		os.mkdir('final/%s' % project_id)


bed_file_input=args.l
chromsizes='/mnt/kufs/scratch/tmorova15/references/GRCh37.genome.txt'
bedtoolspath='/mnt/kufs/scratch/tmorova15/bcbio/bin/bedtools'
randombedcount=args.it
SnpSiftjarpath='/mnt/kufs/scratch/tmorova15/softwares/snpEff/SnpSift.jar'
SnpEffjarpath='/mnt/kufs/scratch/tmorova15/softwares/snpEff/snpEff.jar'
annotation_file_path='/mnt/kufs/scratch/tmorova15/references/common_all_20160601.vcf'


### SNP ###
sanger_snp_filter_option = "\"((FILTER = 'PASS') && (! exists SNP)) && ((GEN[1].GT == '1|0') || (GEN[1].GT == '0|1'))\"" ### filtering NON PASS option is added.
#sanger_snp_filter_option = "\"(! exists SNP ) && ((GEN[1].GT == '1|0') || (GEN[1].GT == '0|1'))\""
##sanger_snp_filter_option= "\"(! exists SNP ) & (HE='1')\"" #This was old because snpsifts gt option removes format column from the file.
dkfz_snp_filter_option = "\"(FILTER == 'PASS') && ((! ID =~ 'rs' ) && ((GEN[1].GT == '1/0') || (GEN[1].GT == '0/1')))\""
#dkfz_snp_filter_option = "\"(! ID =~ 'rs' ) && ((GEN[1].GT == '1/0') || (GEN[1].GT == '0/1'))\""
##dkfz_snp_filter_option ="\"(! ID =~ 'rs') & (HE='1') \"" #This was old because snpsifts gt option removes format column from the file.
broad_snp_filter_option = "\"(! ID =~ 'rs' )\""
consensus_snp_filter_option="\"(FILTER == '') && (! exists dbsnp)\""

### INDEL ###
sanger_indel_filter_option = "\"(! ID =~ 'rs' ) && (FILTER == 'PASS')\"" ### filtering NON PASS option is added.
#sanger_indel_filter_option = "\"(! ID =~ 'rs' )\""
dkfz_indel_filter_option= "\"(FILTER == 'PASS') && ((! ID =~ 'rs' ) && ((GEN[1].GT == '1/0') || (GEN[1].GT == '0/1')))\""
##dkfz_indel_filter_option = "\"(! ID =~ 'rs' ) && ((GEN[1].GT == '1/0') || (GEN[1].GT == '0/1'))\""
###dkfz_indel_filter_option = "\"(! ID =~ 'rs' ) & (HE='1')\"" This was old because snpsifts gt option removes format column from the file.
broad_indel_filter_option=  "\"(! ID =~ 'rs' )\""
consensus_indel_filter_option="\"(FILTER == '')&& (! exists dbsnp)\""





snv_log= open('./final/snv_number.txt','w')
snv_log_annotated = open('./final/snv_variant_classification_annotated.txt','w')
snv_log_final= open('./final/snv_variant_classification_final.txt','w')
snv_log.write('#### patientid found_mutations|number_of_total_mutations\n')
snv_randomized_log=open('./final/snv_randomized_mutation_counts.txt','w')

indel_log =open('./final/indel_number.txt','w')
indel_log_annotated  = open('./final/indel_variant_classification_annotated.txt','w')
indel_log_final  = open('./final/indel_variant_classification_final.txt','w')
indel_log.write('#### patientid found_mutations|number_of_total_mutations\n')
indel_randomized_log=open('./final/indel_randomized_mutation_counts.txt','w')

snv_log.flush()
indel_log.flush()
for project_id in project_names:
	if not os.path.exists('final/Problematic_id'):
		os.mkdir('./final/Problematic_id')

repeat_cheacker=[]
for folder in os.listdir('input/'):


	if len(folder.split('-')) == 5: ### There are other folders/files located in the same input folder such as log.txt of the download or manifests. This wi


		if folder not in folder_prefix.keys():
			if not os.path.exists('final/Problematic_id/{0}'.format(folder)):
				os.mkdir('./final/Problematic_id/%s' % folder)
		
			else:
				pass



		else:

			### Python gives a file already existed error, so this line is for to prevent that error.
			#### If a patient folder is created for another caller, code does not create a new dir and write to the currently existed dir
			if os.path.exists('./final/%s/%s' % (folder_prefix[folder][1],folder_prefix[folder][0])):
				pass
			else:
				os.mkdir('./final/%s/%s' % (folder_prefix[folder][1],folder_prefix[folder][0]))

			### analysis method : sanger, dkfz, broad or consensus.
			analysis_method = folder_prefix[folder][2]

			


			output_directory = './final/%s/%s/' % (folder_prefix[folder][1],folder_prefix[folder][0])

			start=time.time()
			if not os.path.exists('./final/randombeds'): #creates random bed directory
				os.mkdir('./final/randombeds')
				### Random bed generation if not exists ###
				print whattimeisit(), 'Genarating Random Bed Files...'
				for i in range(randombedcount):
					random_commandline=' '.join([bedtoolspath,'shuffle', '-chrom','-i',bed_file_input,'-g',chromsizes,'-excl',gapped_genome_file,'>', './final/randombeds/random_' + '%s' % str(i)])
					subprocess.check_call(random_commandline,shell= True)

			print 'Took', time.time() - start,'to run....'

			### Random Vcf Genaration ###
			# if not os.path.exists(output_directory + 'random_vcfs'):
			# 	os.mkdir(output_directory + 'random_vcfs')


			for name in os.listdir('input/%s' % folder):

				if bool(re.search('somatic.snv_mnv.vcf.gz$', name)): ### checks if snv_mnv.vcf.gz in the iterated variable
					analysis_id=name.split('.')[0]
					
					snv_log.write(folder_prefix[folder][0] + '\t' + analysis_method +'\t' +  analysis_id + '\t') ## write to the file which contains all of the numbers
					snv_log_annotated.write(folder_prefix[folder][0] + '\t' + analysis_method +'\t' + analysis_id + '\t') ## write to the file which contains all of the numbers
					snv_log_final.write(folder_prefix[folder][0] + '\t' + analysis_method +'\t' +  analysis_id + '\t') ## write to the file which contains all of the numbers
					snv_randomized_log.write(folder_prefix[folder][0] + '\t' + analysis_method +'\t' +  analysis_id + '\t') ## write to the file which contains all of the randomized mutations.
					start=time.time()




					if analysis_method == 'sanger':

						annotation_script =' '.join(['java','-jar',SnpSiftjarpath ,'gt','-u','./input/%s/' % folder  + name , '|' ,'java','-jar',SnpSiftjarpath,'filter',sanger_snp_filter_option, '>', output_directory+'annotated_full.' + analysis_id+'.snv_mnv_'+analysis_method])

					elif analysis_method == 'dkfz':

						annotation_script = ' '.join(['java','-jar',SnpSiftjarpath ,'gt', '-u','./input/%s/' % folder  + name , '|' ,'java','-jar',SnpSiftjarpath,'filter',dkfz_snp_filter_option, '>' ,output_directory+'annotated_full.' + analysis_id+'.snv_mnv_'+analysis_method])

					elif analysis_method == 'consensus':
						annotation_script=' '.join(['java','-jar',SnpSiftjarpath,'filter', consensus_snp_filter_option, './input/%s/' % folder  + name, '>' ,output_directory+'annotated_full.' + analysis_id+'.snv_mnv_'+analysis_method])


					else:
						### broad icin ###
						annotation_script=' '.join(['java', '-jar' , SnpSiftjarpath , 'annotate' , '-id' ,annotation_file_path,  './input/%s/' % folder  + name , '|', 'java', '-Xmx6g','-jar',SnpEffjarpath,'eff', 'GRCh37.75' , '-','|', \
											'java' ,'-jar', SnpSiftjarpath,'filter',broad_snp_filter_option,'>' ,output_directory+'annotated_full.' + analysis_id+'.snv_mnv_'+analysis_method])


					if not os.path.exists(output_directory+'annotated_full.' + analysis_id+'.snv_mnv_'+analysis_method): ### This prevents useless annotation.
						print whattimeisit(), "**Annotating %s SNV" % folder_prefix[folder][0]
						subprocess.check_call(annotation_script, shell=True)
						print whattimeisit(), 'Annotation of SNV Took', time.time() - start,'to run....'
						variantclassificationcounter(output_directory+'annotated_full.' + analysis_id+'.snv_mnv_'+analysis_method,snv_log_annotated) ### we would like to see how variants separated through different mutation types 
					
				
					
					print whattimeisit(),"**Patient %s SNV has already annotated, cleaving step will be initiated" % folder_prefix[folder][0]
					print whattimeisit(), "***Cleaving spesific regions given by %s" % bed_file_input
					start=time.time()
					#vcf_isolation_bash_script=["/mnt/kufs/scratch/tmorova15/softwares/vcftools_0.1.13/bin/vcftools","--gzvcf", output_directory+ 'annotated_full.' + analysis_id+'.snv_mnv_'+analysis_method, "--bed" ,  bed_file_input , \
					#						   "--recode" ,"--recode-INFO-all", "--out", output_directory + "final." + analysis_id + ".snv_mnv_" + analysis_method]

					### takes only the mutations within given bed region ###
					vcf_isolation_bash_script =['java','-jar',SnpSiftjarpath,'interval', '-i' , output_directory+ 'annotated_full.' + analysis_id+'.snv_mnv_'+analysis_method, bed_file_input, '>',  output_directory + "final." + analysis_id + ".snv_mnv_" + analysis_method + '.vcf']
					subprocess.check_call(' '.join(vcf_isolation_bash_script),shell=True)
					
					variantclassificationcounter(output_directory + "final." + analysis_id + ".snv_mnv_" + analysis_method + '.vcf',snv_log_final) ### variant classification is written.

					### counts the mutations
					mutation_count_line=[bedtoolspath,'intersect', '-a' , output_directory + 'annotated_full.' + analysis_id + ".snv_mnv_" + analysis_method ,'-b',bed_file_input, '|' ,'wc']

					### count is the count of mutation counts around given bed regions ###
					count = subprocess.check_output(' '.join(mutation_count_line),shell = True)

					### Gives the total mutations from the vcf ###
					total_mutation_count_line=['egrep','-v','^#',  output_directory+ 'annotated_full.' + analysis_id+'.snv_mnv_'+analysis_method , '|', 'wc']
					total = subprocess.check_output(' '.join(total_mutation_count_line),shell=True)

					count=count.split()
					total=total.split()

					print whattimeisit(), 'Cleavage of SNV Took', time.time() - start,'to run....'


					print whattimeisit(), count[0], 'SNV mutations found in', total[0]

					### This is for reading log files  and creating a single file which contains all of the number of the mutations.###

					snv_log.write(str(count[0]) + '|' + str(total[0]) + '\t')



					### bootstrap ###
					start=time.time()
					# random_counter=0
					# for n,randombed in enumerate(os.listdir('./final/randombeds')):
					# 	sys.stdout.write("\r{0}".format(randombed))
					# 	sys.stdout.flush()
					# 	random_isolation_bash_script = ['java','-jar',SnpSiftjarpath,'interval', '-i' , output_directory+ 'annotated_full.' + analysis_id+'.snv_mnv_'+analysis_method, './final/randombeds/'+ randombed, '|',\
					# 	bedtoolspath,'intersect', '-a' ,"-",'-b',bed_file_input, '|' ,'wc']
					# 
					# 	random_mutation_count=subprocess.check_output(' '.join(random_isolation_bash_script),shell=True)
					# 	if count < random_mutation_count:
					# 		random_counter += 1
					# snv_log.write(str(random_counter/float(randombedcount)) + '\t')
					
					random_counter=BEDbootstrap(logfilename = snv_randomized_log, SNV=True, verbose=True)
					snv_randomized_log.write('\n') ### to go to the next patient
					snv_randomized_log.flush()
					snv_log.write(str(float(random_counter) / float(randombedcount)) + '\t')
					snv_log.flush()
					print whattimeisit(), 'Bootstraping of SNV Took', time.time() - start,'to run....'

				elif bool(re.search('somatic.indel.vcf.gz$', name)): #### INDEL PART ####
					analysis_id=name.split('.')[0]
					start=time.time()
					
					indel_log.write(folder_prefix[folder][0] + '\t' + analysis_method +'\t' + analysis_id + '\t')## write to the file which contains all of the numbers
					indel_log_annotated.write(folder_prefix[folder][0] + '\t' + analysis_method +'\t' + analysis_id + '\t') ## write to the file which contains all of the numbers
					indel_log_final.write(folder_prefix[folder][0] + '\t' + analysis_method +'\t' + analysis_id + '\t') ## write to the file which contains all of the numbers
					indel_randomized_log.write(folder_prefix[folder][0] + '\t' + analysis_method +'\t' +  analysis_id + '\t')



					if analysis_method == 'sanger':

						annotation_script=' '.join(['java', '-jar' , SnpSiftjarpath , 'annotate' , '-id' ,annotation_file_path,  './input/%s/' % folder  + name , '|', 'java', '-Xmx6g','-jar',SnpEffjarpath,'eff', 'GRCh37.75' , '-','|', \
											'java' ,'-jar', SnpSiftjarpath,'filter',sanger_indel_filter_option,'>' ,output_directory+'annotated_full.' + analysis_id+'.indel_'+analysis_method])


					elif analysis_method == 'dkfz':

						annotation_script = ' '.join(['java','-jar',SnpSiftjarpath ,'gt', '-u','./input/%s/' % folder  + name , '|' ,'java','-jar',SnpSiftjarpath,'filter',dkfz_indel_filter_option, '>' ,output_directory+'annotated_full.' + analysis_id+'.indel_'+analysis_method])


					elif analysis_method == 'consensus':
						annotation_script=' '.join(['java','-jar',SnpSiftjarpath,'filter', consensus_indel_filter_option, './input/%s/' % folder  + name, '>' ,output_directory+'annotated_full.' + analysis_id+'.indel_'+analysis_method])


					else:

						annotation_script=' '.join(['java', '-jar' , SnpSiftjarpath , 'annotate' , '-id' ,annotation_file_path,  './input/%s/' % folder  + name , '|', 'java', '-Xmx6g','-jar',SnpEffjarpath,'eff', 'GRCh37.75' , '-','|', \
											'java' ,'-jar', SnpSiftjarpath,'filter',broad_indel_filter_option,'>' ,output_directory+ 'annotated_full.' + analysis_id+'.indel_'+analysis_method])


					if not os.path.exists(output_directory+ 'annotated_full.' + analysis_id+'.indel_'+analysis_method):  ### This prevents useless annotation.
						print whattimeisit(),"**Annotating %s Indel" % folder_prefix[folder][0]
						subprocess.check_call(annotation_script,shell=True)
						print whattimeisit(),'Annotation of Indel Took', time.time() - start,'to run....'
						variantclassificationcounter(output_directory+'annotated_full.' + analysis_id+'.indel_'+analysis_method,indel_log_annotated) ### we would like to see how variants separated through different mutation types 
						
					
					
					print whattimeisit(),"**Patient %s Indel has already annotated, cleaving step will be initiated"  % folder_prefix[folder][0]
					print whattimeisit(),"***Cleaving spesific regions given by %s" % bed_file_input
					# vcf_isolation_bash_script=["/mnt/kufs/scratch/tmorova15/softwares/vcftools_0.1.13/bin/vcftools","--gzvcf", output_directory+ 'annotated_full.' + analysis_id+'.indel_'+analysis_method, "--bed" ,  bed_file_input , "--recode" ,"--recode-INFO-all", "--out", output_directory + "final." + analysis_id + ".indel_" + analysis_method ]
					# r1=subprocess.check_call(' '.join(vcf_isolation_bash_script),shell=True)

					vcf_isolation_bash_script =['java','-jar',SnpSiftjarpath,'interval', '-i' , output_directory+ 'annotated_full.' + analysis_id+'.indel_'+analysis_method, bed_file_input, '>',  output_directory + "final." + analysis_id + ".indel_" + analysis_method + '.vcf']
					subprocess.check_call(' '.join(vcf_isolation_bash_script),shell=True)
					
					variantclassificationcounter(output_directory + "final." + analysis_id + ".indel_" + analysis_method + '.vcf',indel_log_final) ### variant classification is written.

					### counts the mutations
					mutation_count_line=[bedtoolspath,'intersect', '-a' ,  output_directory + 'annotated_full.' + analysis_id+'.indel_'+analysis_method,'-b',bed_file_input, '|' ,'wc']

					### count is the count of mutation counts around given bed regions ###
					count = subprocess.check_output(' '.join(mutation_count_line),shell = True)

					### Gives the total mutations from the vcf ###
					total_mutation_count_line=['egrep','-v','^#', output_directory+ 'annotated_full.' + analysis_id+'.indel_'+analysis_method, '|', 'wc']
					total = subprocess.check_output(' '.join(total_mutation_count_line),shell=True)

					count=count.split()
					total=total.split()
					print whattimeisit(), count[0], 'Indel mutations found in', total[0]
					indel_log.write(str(count[0]) + '|' + str(total[0]) + '\t')
					print whattimeisit(), 'Cleavage of Indel Took', time.time() - start,'to run....'

					start=time.time()
					# random_counter=0
					# for n,randombed in enumerate(os.listdir('./final/randombeds')):
					# 	sys.stdout.write("\r{0}".format(randombed))
					# 	sys.stdout.flush()
					# 	random_isolation_bash_script = ['java','-jar',SnpSiftjarpath,'interval', '-i' , output_directory+ 'annotated_full.' + analysis_id+'.indel_'+analysis_method, './final/randombeds/'+ randombed, '|',\
					# 	bedtoolspath,'intersect', '-a' ,"-",'-b',bed_file_input, '|' ,'wc']
					# 
					# 	random_mutation_count=subprocess.check_output(' '.join(random_isolation_bash_script),shell=True)
					# 	if count < random_mutation_count:
					# 		random_counter += 1

					# indel_log.write(str(random_counter/float(randombedcount)) + '\t')
					# indel_log.flush()
					random_counter=BEDbootstrap(logfilename = indel_randomized_log, SNV=False, verbose=True)
					indel_randomized_log.write('\n') ### to go to the next patient 
					indel_randomized_log.flush()
					indel_log.write(str( float(random_counter) / float(randombedcount) ) + '\t')
					indel_log.flush()
					


					print whattimeisit(), 'Bootstraping of Indel Took', time.time() - start,'to run....'







		snv_log.write('\n')
		indel_log.write('\n')


snv_log.close()
snv_log_annotated.close()
snv_randomized_log.close()
indel_log.close()
indel_log_annotated.close()
indel_randomized_log.close()




print whattimeisit(),'Process has finalized without any error ! Well Done !!'