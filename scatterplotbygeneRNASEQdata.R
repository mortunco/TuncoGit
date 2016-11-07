
###Hierchical

setwd('/Volumes/MyPassport/MyProjects/Machine_Learning/PRAD/')
load('PRAD_data_matrix.RData')

### primary cancer list ####
b=0
only_primary=c()
for (patient_index in 1:nrow(total_mat)){
  
  temp=unlist(strsplit(rownames(total_mat)[patient_index],split = '-',fixed = T))[4]
  if (!(xor(grepl('11',temp),grepl('6',temp)))){
    #print(temp)
    b=b+1 
    only_primary=c(only_primary,patient_index)}
}

cat(paste('Number of primary cancer patients:',b))
total_mat= total_mat[only_primary,]

### entrez thing ###
a=colnames(total_mat)
a=toupper(unlist(lapply(a, function(x)  unlist((strsplit(x,'|',fixed = TRUE)))[1]))) #ask mehmet hoca if this is OK ? 
colnames(total_mat)=a


###Scattter + Abline ###
setwd('~/Desktop/')
x="SKP2"
y="KLK3"
pdf(file = paste(x,y,'scatter_plots.pdf',sep ='|'))
plot(log2(total_mat[,x]+1),log2(total_mat[,y]+1),pch=16,col='red',
     xlab=paste(x,'log2(Expression Count)',sep = ' '),ylab=paste(y,'log2(Expression Count)',sep = ' '))
reg=lm(log2(total_mat[,y]+1) ~ log2(total_mat[,x]+1))
summary(reg)$r.squared
abline(lm(log2(total_mat[,y]+1) ~ log2(total_mat[,x]+1)),col='blue',cex=3)
dev.off()


plot(1:10,1:10)
abline(0,1)
