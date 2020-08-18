`+` <- function(e1, e2) {
  if (is.character(e1) | is.character(e2)) {
    paste0(e1, e2)
  } else {
    base::`+`(e1, e2)
  }
}

library(survival)
library(survminer)
library(dplyr)

cancer='UCEC'             #cancer name
val=90                    #%size of sampling
iter=100                   # num of iterations

out_file='/Users/macbook/Desktop/UCEC/'+cancer+'_valid_90.csv'
beta_file='/Users/macbook/Desktop/APOP_PAN_CANCER/HR_APOP/'+cancer+'_HRmedian_raw.csv'
tcga_datafile="/Users/macbook/Desktop/webserver-pan-can/data_basic/RDS_files/"+cancer+".RDS"

beta<-read.csv(file=beta_file,header =TRUE, sep = ",", dec = ".")
data1<-readRDS(tcga_datafile)
genes <- scan(file = '/Users/macbook/Desktop/UCEC/'+cancer+'_genes.csv', what = 'character', sep = ',')


#select(LUSC_C, gep[i])[,1]
PI=0
for(i in seq(from=1, to=length(genes), by=1))#no of cancers
{
  k=beta[beta$Gene==genes[i],][2][,1]
  if (k>0)
  {
    b=1*(select(data1, genes[i])[,1]>median(select(data1, genes[i])[,1]))}
  
  if (k<0)
  {b=1*(select(data1, genes[i])[,1]<median(select(data1, genes[i])[,1]))}
  
  PI=PI+b
}
data1$PI=PI
#data1$PI=PI+1*(data1$age_at_diagnosis>60)
#g=length(genes)/2

HR=c()
C_ind=c()
j=1
while(j!=iter) {
  
  data2=sample_n(data1, ceiling((val/100)*nrow(data1)))
  
  surv_object <- Surv(time = data2$OS.time/30, event = data2$vital_status)
  
  #survival analysis: fits cox ph model to find HR for mean cut
  fit1 <- survfit(surv_object~data2$PI>7); 
  #summary(fit1);
  #ggsurvplot(fit1, data=data1)
  
  
  fit1.coxph <- coxph(surv_object~data2$PI>7)
  # summary(fit1.coxph);
  first <- coef(summary(fit1.coxph))
  #first
  a=summary(fit1.coxph)
  #a
  #max U and Q, min p value
  
  #fit1.coxph$concordance[6]
  #gg<-ggsurvplot(fit1, data=data1,legend=c(0.8,0.8),xlab = "Time (Months)",axes.offset=FALSE,legend.labs = c("Low Risk", "High Risk"))
  
  #ggsave('/Users/macbook/Desktop/APOP_PAN_CANCER/top_ten/KMplots/'+cancer+'voting.jpg', plot = print(gg),height = 5 ,width = 7, dpi = 1000)
  if((!is.na(first[5]))&&(first[5]<0.05))
  {
    HR[j]=first[2]
    C_ind[j]=a$concordance[1]
    j=j+1
  }
  #j=j+1
}
print(mean(HR))
print(sqrt(var(HR)))

print(mean(C_ind))
print(sqrt(var(C_ind)))


{write.table(cbind('HR','C_index'),
             file=out_file,row.names=F,col.names=F,sep = ',',append = T);#output file
}

{write.table(cbind(HR,C_ind),
             file=out_file,row.names=F,col.names=F,sep = ',',append = T);#output file
}


boxplot(HR1,HR2,HR3,ylab="Hazard Ratio (HR)", xlab="Sampling size",names = c('50%','70%','90%'),col="grey")
#boxplot(C1,C2,C3,ylab="C-Index", xlab="Sampling size",names = c('50%','70%','90%'),col="grey")
