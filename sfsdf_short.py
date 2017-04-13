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

print '### Input files must be in zipped format. Otherwise, code will not be reliable ###'
print 'Bootstrap Iteration is set up to {0}'.format(args.it)

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
		### checks the number of mutations found in given bed proximity for each randomized Bed.### 
		#random_isolation_bash_script = ['java','-jar',SnpSiftjarpath,'interval', '-i' ,output_directory+'annotated_full.' + os.path.splitext(name)[0] , './final/randombeds/'+ randombed, '|',\
		#bedtoolspath,'intersect', '-a' ,"-",'-b','./final/randombeds/'+ randombed, '|' ,'wc']
		random_isolation_bash_script = [bedtoolspath, 'intersect','-a' ,output_directory+'annotated_full.' + os.path.splitext(name)[0] , '-b','./final/randombeds/'+ randombed, '|','wc']
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

print 'Warning you are running the SFSDF script with any given VCF input. In this mode, there will be no annotation done unless it was set so.'
topdirectory=args.p

print topdirectory

os.chdir(topdirectory)
print whattimeisit(),os.getcwd()

bed_file_input=args.l
chromsizes='/mnt/kufs/scratch/tmorova15/references/GRCh37.genome.txt'
bedtoolspath='/mnt/kufs/scratch/tmorova15/bcbio/bin/bedtools'
randombedcount=args.it
SnpSiftjarpath='/mnt/kufs/scratch/tmorova15/softwares/snpEff/SnpSift.jar'
SnpEffjarpath='/mnt/kufs/scratch/tmorova15/softwares/snpEff/snpEff.jar'
annotation_file_path='/mnt/kufs/scratch/tmorova15/references/common_all_20160601.vcf'
gapped_genome_file='/mnt/kufs/scratch/tmorova15/references/gap_regions_grch37.bed'


### SNP ###
sanger_snp_filter_option = "\"((FILTER = 'PASS') && (! exists SNP)) && ((GEN[1].GT == '1|0') || (GEN[1].GT == '0|1'))\"" ### filtering NON PASS option is added.
#sanger_snp_filter_option = "\"(! exists SNP ) && ((GEN[1].GT == '1|0') || (GEN[1].GT == '0|1'))\""
##sanger_snp_filter_option= "\"(! exists SNP ) & (HE='1')\"" #This was old because snpsifts gt option removes format column from the file.
dkfz_snp_filter_option = "\"(FILTER == 'PASS') && ((! ID =~ 'rs' ) && ((GEN[1].GT == '1/0') || (GEN[1].GT == '0/1')))\""
#dkfz_snp_filter_option = "\"(! ID =~ 'rs' ) && ((GEN[1].GT == '1/0') || (GEN[1].GT == '0/1'))\""
##dkfz_snp_filter_option ="\"(! ID =~ 'rs') & (HE='1') \"" #This was old because snpsifts gt option removes format column from the file.
broad_snp_filter_option = "\"(! ID =~ 'rs' )\""
#consensus_snp_filter_option="\"(FILTER == '') && (! exists dbsnp)\"" ##Deprecated after getting an error states if dbsnp is now existed in a VCF header, it will stop the process ##
consensus_snp_filter_option="\"(! ID =~ 'rs')\""

### INDEL ###
sanger_indel_filter_option = "\"(! ID =~ 'rs' ) && (FILTER == 'PASS')\"" ### filtering NON PASS option is added.
#sanger_indel_filter_option = "\"(! ID =~ 'rs' )\""
dkfz_indel_filter_option= "\"(FILTER == 'PASS') && ((! ID =~ 'rs' ) && ((GEN[1].GT == '1/0') || (GEN[1].GT == '0/1')))\""
##dkfz_indel_filter_option = "\"(! ID =~ 'rs' ) && ((GEN[1].GT == '1/0') || (GEN[1].GT == '0/1'))\""
###dkfz_indel_filter_option = "\"(! ID =~ 'rs' ) & (HE='1')\"" This was old because snpsifts gt option removes format column from the file.
broad_indel_filter_option=  "\"(! ID =~ 'rs' )\""
consensus_indel_filter_option="\"(FILTER == '')&& (! exists dbsnp)\""


snv_log=open('./final/snv_number.txt','w')
snv_log_annotated = open('./final/snv_variant_classification_annotated.txt','w')
snv_log_final= open('./final/snv_variant_classification_final.txt','w')
snv_log.write('#### patientid found_mutations|number_of_total_mutations\n')
snv_randomized_log=open('./final/snv_randomized_mutation_counts.txt','w')
snv_log.flush()

for folder in os.listdir('input/'):


	

		### Python gives a file already existed error, so this line is for to prevent that error.
		#### If a patient folder is created for another caller, code does not create a new dir and write to the currently existed dir
		if os.path.exists('./final/%s/' % (folder)):
			pass
		else:
			os.mkdir('./final/%s' % (folder))

		output_directory = './final/%s/' % (folder)

		start=time.time()
		if not os.path.exists('./final/randombeds'): #creates random bed directory
			os.mkdir('./final/randombeds')
			### Random bed generation if not exists ###
			print whattimeisit(), 'Genarating Random Bed Files...'
			for i in range(randombedcount):
				random_commandline=' '.join([bedtoolspath,'shuffle', '-chrom','-i',bed_file_input,'-g',chromsizes,'-excl',gapped_genome_file,'>', './final/randombeds/random_' + '%s' % str(i)])
				subprocess.check_call(random_commandline,shell= True)

		print 'Took', time.time() - start,'to run....'
		
		### Bundan altina bakmadim ###

		for name in os.listdir('input/%s' % folder):

			if bool(re.search('vcf.gz$', name)): ### checks if snv_mnv.vcf.gz in the iterated variable
							
				snv_log.write(folder + '\t' + name +'\t' +  'AnalysisName' + '\t') ## write to the file which contains all of the numbers
				snv_log_annotated.write(folder + '\t' + name +'\t' +  'AnalysisName' + '\t') ## write to the file which contains all of the numbers
				snv_log_final.write(folder + '\t' + name +'\t' +  'AnalysisName' + '\t') ## write to the file which contains all of the numbers
				snv_randomized_log.write(folder + '\t' + name +'\t' +  'AnalysisName' + '\t') ## write to the file which contains all of the randomized mutations.
				start=time.time()



				### This summer edit this line for new upcoming vcf types. 
				#if analysis_method == 'sanger':
				
				### Annotation scriptin sonuna os.path.splitext(name)[0] ibaresini koymamin sebebi, snpEff gz gorunce gunzipliyor ama dosya daha ocneden ziplenmis durumdaysa error veriyor. Gz den kurtulmak icin koyuldu. ###
				annotation_script =' '.join(['java','-jar',SnpSiftjarpath ,'gt','-u','./input/%s/' % folder  + name , '|' ,'java','-jar',SnpSiftjarpath,'filter',consensus_snp_filter_option, '>', output_directory+'annotated_full.' + os.path.splitext(name)[0] ])
			

				#elif analysis_method == 'consensus':
				#	annotation_script=' '.join(['java','-jar',SnpSiftjarpath,'filter', consensus_snp_filter_option, './input/%s/' % folder  + name, '>' ,output_directory+'annotated_full.' + analysis_id+'.snv_mnv_'+analysis_method])



				if not os.path.exists(output_directory+'annotated_full.' + os.path.splitext(name)[0]): ### This prevents useless annotation.
					print whattimeisit(), "**Annotating %s " % name 
					subprocess.check_call(annotation_script, shell=True)
					print whattimeisit(), 'Annotation of vcf Took', time.time() - start,'to run....'
					variantclassificationcounter(output_directory+'annotated_full.' + os.path.splitext(name)[0],snv_log_annotated) ### we would like to see how variants separated through different mutation types 
				
			
				
				print whattimeisit(),"**%s vcf file has already annotated, cleaving step will be initiated" % name
				print whattimeisit(), "***Cleaving spesific regions given by %s" % bed_file_input
				start=time.time()
				#vcf_isolation_bash_script=["/mnt/kufs/scratch/tmorova15/softwares/vcftools_0.1.13/bin/vcftools","--gzvcf", output_directory+ 'annotated_full.' + analysis_id+'.snv_mnv_'+analysis_method, "--bed" ,  bed_file_input , \
				#						   "--recode" ,"--recode-INFO-all", "--out", output_directory + "final." + analysis_id + ".snv_mnv_" + analysis_method]

				### takes only the mutations within given bed region ###
				vcf_isolation_bash_script =['java','-jar',SnpSiftjarpath,'interval', '-i' , output_directory+ 'annotated_full.' + os.path.splitext(name)[0], bed_file_input, '>',  output_directory + "final." + os.path.splitext(name)[0] ]
				subprocess.check_call(' '.join(vcf_isolation_bash_script),shell=True)
				
				variantclassificationcounter(output_directory + "final." + os.path.splitext(name)[0]  ,snv_log_final) ### variant classification is written.

				### counts the mutations
				mutation_count_line=[bedtoolspath,'intersect', '-a' , output_directory+ 'annotated_full.' + os.path.splitext(name)[0]  ,'-b',bed_file_input, '|' ,'wc']

				### count is the count of mutation counts around given bed regions ###
				count = subprocess.check_output(' '.join(mutation_count_line),shell = True)

				### Gives the total mutations from the vcf ###
				total_mutation_count_line=['egrep','-v','^#', output_directory+ 'annotated_full.' + os.path.splitext(name)[0] , '|', 'wc']
				total = subprocess.check_output(' '.join(total_mutation_count_line),shell=True)

				count=count.split()
				total=total.split()

				print whattimeisit(), 'Cleavage  Took', time.time() - start,'to run....'


				print whattimeisit(), count[0], 'mutations found in', total[0]

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
				snv_log.write('\n')



snv_log.close()
snv_log_annotated.close()
snv_randomized_log.close()




print whattimeisit(),'Process has finalized without any error ! Well Done !!'