

df=[]
with open('/Users/morova/Desktop/vcfclinicalanlysis/structural_variation_intervalcheck/SV_location.txt','r') as data:
	for item in data:
		temp=item.rstrip().split('\t')
		df.append(temp)
		

print df[6]
peakno_1139 = (39613382 + 39612569) / 2
peakno_1140 = (39675741 + 39674943) / 2
peakno_1141 = (40279595 + 40278811) / 2
peakno_1395 = (42718670 + 42717859) / 2
peakno_1396 = (42806084 + 42805320) / 2
peakno_1398 = (42893039 + 42892275) / 2
peakno_1399 = (42898848 + 42898084) / 2
peakno_1403 = (42955717 + 42954947) / 2

peaklist=[peakno_1139,peakno_1140,peakno_1141,peakno_1395,peakno_1396,peakno_1398,peakno_1399,peakno_1403]


f3=open('/Users/morova/Desktop/vcfclinicalanlysis/structural_variation_intervalcheck/location_summary.txt','w')
for aliquot in df:
	f3.write('\t'.join(aliquot[0:3]) + '\t')
	if len(aliquot) > 5:
		for fusion in aliquot[4:]:
			print fusion
			quality, firstone, secondone = fusion.split('|')[0] , fusion.split('|')[1] , fusion.split('|')[2]
			middpoint_first=[firstone.split(':')[0],(int(firstone.split(':')[3]) + int(firstone.split(':')[2])) / 2]
			midpoint_second=[secondone.split(':')[0], (int(secondone.split(':')[3]) + int(secondone.split(':')[2]) ) / 2 ]
			
			midpoint_first_distance=min([ abs(i - middpoint_first[1] + 1) for i in peaklist])
			midpoint_second_distance=min([ abs(i - midpoint_second[1] + 1) for i in peaklist])
			
			f3.write(str(quality) + '|' + str(middpoint_first[0]) + ':' + str(midpoint_first_distance) + '|' + str(midpoint_second[0]) + ':' + str(midpoint_second_distance) + '\t' )
	
	elif len(aliquot) == 5:
	
		quality, firstone, secondone = fusion.split('|')[0] , fusion.split('|')[1] , fusion.split('|')[2]
		middpoint_first=[firstone.split(':')[0],(int(firstone.split(':')[3]) + int(firstone.split(':')[2])) / 2]
		midpoint_second=[secondone.split(':')[0], (int(secondone.split(':')[3]) + int(secondone.split(':')[2]) ) / 2 ]
			
		midpoint_first_distance=min([ abs(i - middpoint_first[1] + 1) for i in peaklist])
		midpoint_second_distance=min([ abs(i - midpoint_second[1] + 1) for i in peaklist])
			
		f3.write(str(quality) + '|' + str(middpoint_first[0]) + ':' + str(midpoint_first_distance) + '|' + str(midpoint_second[0]) + ':' + str(midpoint_second_distance) + '\t' )
		
	else:
		pass
	
	f3.write('\n')
	



f3.close()
			
			
			
			
			
			
			
			
			
			
			
			
			
			
			
			
			
			
			




