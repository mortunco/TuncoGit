
f3=open('librarygercek.fa','wb')

with open('/Users/morova/libarygercek.csv') as sourcefile:
	for line in sourcefile:
		#print line
		f3.write('>')
		f3.write(line.rstrip().split(',')[1] + '\n')
		f3.write(line.rstrip().split(',')[2] + '\n')
		
		
f3.close()