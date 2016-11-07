import argparse
import os
import subprocess
import shutil
parser = argparse.ArgumentParser(description='This code gathers the information in Batch directories and merge in single finalized folder.')
parser.add_argument('-p',type=str, help='path of the project directory which contains Batch directories of spesific run.')
parser.add_argument('-name',type=str,help='name of the project you want to put')
args=parser.parse_args()


for item in vars(args): ###
	if getattr(args,item) == None:
		raise ValueError , "All arguments must be set"


topdirectory=args.p

os.chdir(topdirectory)

toplist=os.listdir(os.getcwd())

batchlist= [i for i in toplist if 'Batch' in i]

os.mkdir('{0}_final_output_from_all_batches'.format(args.name))
os.mkdir('{0}_annotated_from_all_batches'.format(args.name))
os.mkdir('{0}_final_snvlog_from_all_batches'.format(args.name))
os.mkdir('{0}_variant_classification_snv'.format(args.name))
os.mkdir('{0}_final_indellog_from_all_batches'.format(args.name))
os.mkdir('{0}_variant_classification_indel'.format(args.name))


print batchlist


for folder in batchlist:
	print folder
	
	### SNV Mutation No#
	os.rename('{0}/final/snv_number.txt'.format(folder),folder+'/final/'+'snv_number'+folder+'.txt')
	shutil.move(folder+'/final/'+'snv_number'+folder+'.txt','{0}_final_snvlog_from_all_batches/'.format(args.name)+'snv_number'+folder+'.txt' )
	
	### SNV all Variant Distribution ###
	os.rename('{0}/final/snv_variant_classification_annotated.txt'.format(folder),folder+'/final/'+'snv_variant_classification_annotated'+folder+'.txt')
	shutil.move(folder+'/final/'+'snv_variant_classification_annotated'+folder+'.txt','{0}_variant_classification_snv/'.format(args.name)+'snv_variant_classification_annotated_'+folder+'.txt')
	
	### SNV final Variant Distribution ###
	os.rename('{0}/final/snv_variant_classification_final.txt'.format(folder),folder+'/final/'+'snv_variant_classification_final'+folder+'.txt')
	shutil.move(folder+'/final/'+'snv_variant_classification_final'+folder+'.txt','{0}_variant_classification_snv/'.format(args.name)+'snv_variant_classification_final_'+folder+'.txt')
	
	### Indel Mutation No#
	os.rename('{0}/final/indel_number.txt'.format(folder),folder+'/final/'+'indel_number'+folder+'.txt')
	shutil.move(folder+'/final/'+'indel_number'+folder+'.txt','{0}_final_indellog_from_all_batches/'.format(args.name)+'indel_number'+folder+'.txt')
	
	### Indel all Variant Distribution ###
	os.rename('{0}/final/indel_variant_classification_annotated.txt'.format(folder),folder+'/final/'+'indel_variant_classification_annotated'+folder+'.txt')
	shutil.move(folder+'/final/'+'indel_variant_classification_annotated'+folder+'.txt','{0}_variant_classification_indel/'.format(args.name)+'indel_variant_classification_annotated_'+folder+'.txt')
	
	### Indel final Variant Distributio ###
	os.rename('{0}/final/indel_variant_classification_final.txt'.format(folder),folder+'/final/'+'indel_variant_classification_final'+folder+'.txt')
	shutil.move(folder+'/final/'+'indel_variant_classification_final'+folder+'.txt','{0}_variant_classification_snv/'.format(args.name)+'indel_variant_classification_final_'+folder+'.txt')
	
	


	### For processesed indels
	subprocess.check_call(" ".join(['rsync','-av','--exclude="randombeds/"','--exclude="*/*/annotated_full.*"',folder+'/final/','{0}_final_output_from_all_batches/'.format(args.name)]) ,shell=True) ### keeps final directory
	##subprocess.check_call(" ".join(['rm','-rf','--',folder+'/final']) ,shell=True )
	##os.mkdir(folder + '/final')

	#### For input files
	subprocess.check_call(" ".join(['rsync','-av','--exclude="randombeds/"','--exclude="*/*/final.*"',folder+'/final/','{0}_annotated_output_from_all_batches/'.format(args.name)]) ,shell=True) ### keeps annotated files
