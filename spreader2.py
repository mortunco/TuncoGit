import os
import subprocess
import shutil
import sys
import argparse

parser = argparse.ArgumentParser(description='This code splits a sigle folder to the sub directories.')
parser.add_argument("-p",type=str,help='Path to the directory which contains the directory that contains all patients files')
parser.add_argument("--dirname",help='name of the directory that contains patient info.')
parser.add_argument('-n',type=int,help='Number of cores allocated in total')
parser.add_argument('-q',type=str,help='Full path of the qsub file')


args=parser.parse_args()

for item in vars(args): ###
	if getattr(args,item) == None:
		raise ValueError , "All arguments must be set"




n=args.n
os.chdir(args.p)
all_file_list=os.listdir("./{0}".format(args.dirname))

qsub_file_name=args.q

print len(all_file_list) , "Sub folders found"


for item in range(1,n+1):

	os.mkdir("./Batch{0}".format(item))
	os.mkdir("./Batch{0}/input".format(item))
	os.mkdir("./Batch{0}/final" .format(item))
	shutil.copy2("./{0}".format(qsub_file_name),"./Batch{0}/example_final".format(item) )

while len(all_file_list) > 0:

	for item in range(1,n+1):
		if len(all_file_list) <1:
			break
		shutil.move("./{0}/{1}" .format(args.dirname,all_file_list[-1]), "./Batch{0}/input/{1}".format(item,all_file_list[-1]))
		all_file_list.pop()
		all_file_list=os.listdir('./{0}'.format(args.dirname))






