import subprocess
import collections


filen='/Users/morova/Desktop/vcfclinicalanlysis/concensuscalls/ff464a6e-d7a3-4b32-a247-6fd984e162bb/0bfd1068-3fd8-a95b-e050-11ac0c4860c3.consensus.20160830.somatic.snv_mnv.vcf'
script=' '.join(['cat',filen,'|', 'grep','-o','Variant_Classification=.*','|' 'cut','-f','2','-d','='])

a=subprocess.check_output(script,shell=True)

a=a.split('\n')

#counter = collections.Counter(a)
#f3.write(str(counter))


f3=open('manyakleleve.txt','w')

def variantclassificationcounter(filename,textfilename):
	script=' '.join(['cat',filename,'|', 'grep','-o','Variant_Classification=.*','|' 'cut','-f','2','-d','='])
	output=subprocess.check_output(script,shell=True).split('\n')
	counter=collections.Counter(output)
	for key,value in zip(counter.keys(),counter.values()): textfilename.write(str(key)+'|'+str(value)+';')
	textfilename.write('\n')
	textfilename.flush()
	return

variantclassificationcounter(filen,f3)


f3.close()