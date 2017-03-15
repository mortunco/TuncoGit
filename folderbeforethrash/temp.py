#!/usr/bin/env python

### bu kodun amaci, yaptigim mutational database'de bir seferlik hazirlik yapilmis listeyi eslestirmektir. Bu sekilde direk patient id ile klinik bilgiye ulasacagiz. ###
import pandas

with open('/Users/morova/Google Drive/TuncProjectStep2-3/Yunus_HPC/PCA_ICGC/ICGC_PCa_dataframe_1.csv') as myfile:
	content = [i.rstrip().split(',')[1::] for i in myfile]
	


mydict={}
for i in range(len(content)):
	
	mydict[content[i][1]] =  [content[i][0]] +content[i][2::]
	

### TMPRSS2:ERG postive filter###
sv_list=[]
with open('/Users/morova/Google Drive/TuncProjectStep2-3/Yunus_HPC/PCA_ICGC/SV_number_sanger.txt','rb') as sv_doc:
	for line in sv_doc:
		sv_list.append(line.rstrip().split('\t'))

for i in xrange(len(sv_list)):


	if len(sv_list[i]) == 4 or 'TMPRSS2' not in sv_list[i][4]:
		sv_list[i].append('Negative')

	if len(sv_list[i]) > 5:
		while len(sv_list[i]) > 5:
			sv_list[i].pop()
	if 'TMPRSS2' not in sv_list[i][4]:
		sv_list[i][4] = 'Negative'

	
sv_dict={}

for i in sv_list:
	sv_dict[i[1]] = i[4]
	

	
	

# ### CONTROL INDEL ###
# control_dict={}
# with open('/Users/morova/Google Drive/TuncProjectStep2-3/Yunus_HPC/PCA_ICGC/Control_Experiment/control_final_indellog_from_all_batches/Control_sanger_indel.txt','rb') as control:
# 	for line in control:
# 		control_dict[line.rstrip().split('\t')[0]] = line.rstrip().split('\t')[0::]


### AR SNV ###
ar_dict_snv={}
with open('/Users/morova/Google Drive/TuncProjectStep2-3/Yunus_HPC/PCA_ICGC/AR_250widen/PCa_concensus_AR250widen/PCa_consensus_AR250Wide_new_final_snvlog_from_all_batches/snv_number_compiled.txt','rb') as AR:
	for line in AR:
		ar_dict_snv[line.rstrip().split('\t')[0]] = line.rstrip().split('\t')[0::]
### AR Indel ###
ar_dict_indel={}
with open('/Users/morova/Google Drive/TuncProjectStep2-3/Yunus_HPC/PCA_ICGC/AR_250widen/PCa_concensus_AR250widen/PCa_consensus_AR250Wide_new_final_indellog_from_all_batches/indel_compiled.txt','rb') as AR:
	for line in AR:
		ar_dict_indel[line.rstrip().split('\t')[0]] = line.rstrip().split('\t')[0::]
		
		
df=[]
for i in ar_dict_indel.keys():
	try:
		temp = [ [i] + mydict[i][0:6] + ar_dict_snv[i][2:4] + ar_dict_indel[i][2:4] ]
		df.append(temp[0])
	except:
		pass
for i in df:
	if len(i) > 11:
		print i
	

df2=pandas.DataFrame(df)

print len(['ProjectName','ICGCPatientID','TMPRSS2:ERG_Fusion_Status','Specimen_Type','Staging_System','Cancer_Degree','Cancer_Degree_ALT','SNV','FDR','INDEL','FDR'])
df2.columns=['ProjectName','ICGCPatientID','TMPRSS2:ERG_Fusion_Status','Specimen_Type','Staging_System','Cancer_Degree','Cancer_Degree_ALT','SNV','FDR','INDEL','FDR']


#print 'Data frame was saved under /Users/morova/Google Drive/TuncProjectStep2-3/Yunus_HPC/PCA_ICGC/ as ICGC_PCa_dataframe.csv'
df2.to_csv('/Users/morova/Google Drive/TuncProjectStep2-3/Yunus_HPC/PCA_ICGC/ICGC_PCa_denem2.csv')

		
		
	
		
		


