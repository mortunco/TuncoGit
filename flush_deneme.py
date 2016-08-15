import argparse

parser = argparse.ArgumentParser(description='This code splits a sigle folder to the sub directories.')
parser.add_argument("-p",type=str,help='The path of the main folder which contains all the patients folders')
parser.add_argument('-n',type=int,help='Number of cores allocated in total')

args=parser.parse_args()

# print args.p
# print args.n

for item in vars(args):
	if getattr(args,item) == None:
		raise ValueError , 'All arguments must be setted'