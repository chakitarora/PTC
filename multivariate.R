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


cancer='THCA'
beta_file='/Users/macbook/Desktop/APOP_PAN_CANCER/HR_APOP/THCA_HRmedian_raw.csv'
tcga_datafile="/Users/macbook/Desktop/THCA/THCA.RDS"

beta<-read.csv(file=beta_file,header =TRUE, sep = ",", dec = ".")
data1<-readRDS(tcga_datafile)
genes <- scan(file = '/Users/macbook/Desktop/THCA/THCA_genes.csv', what = 'character', sep = ',')

PI=0
v=0
for(i in seq(from=1, to=length(genes), by=1))#no of cancers
{
  if (beta[beta$Gene==genes[i],][4][,1]<0.05)
  {
    v=v+1
    k=beta[beta$Gene==genes[i],][2][,1]
    if (k>0)
    {
      b=1*(select(data1, genes[i])[,1]>median(select(data1, genes[i])[,1]))}
    
    if (k<0)
    {b=1*(select(data1, genes[i])[,1]<median(select(data1, genes[i])[,1]))}
    
    PI=PI+b
  }
}
data1$PI=PI
#data1$PI=PI+1*(data1$age_at_diagnosis>60)

#############################################
data2<-subset(data1,ajcc_pathologic_tumor_stage!='[NA]')
#data2<-subset(data2,age_at_diagnosis>45)
#data2<-subset(data2,(ajcc_metastasis_pathologic_pm!='[NA]')&(ajcc_metastasis_pathologic_pm!='MX'))
data2<-subset(data2,(residual_tumor!='[NA]')&(residual_tumor!='[Unknown]')&(residual_tumor!='RX'))
#stage. m stage, residual tumor, age>45


data2 <- mutate(data2, PI = ifelse((PI > 5), "2", "1"))
data2 <- mutate(data2, age_at_diagnosis = ifelse((age_at_diagnosis > 60), "2", "1"))
data2 <- mutate(data2, ajcc_pathologic_tumor_stage = ifelse(((ajcc_pathologic_tumor_stage=='Stage III')|(ajcc_pathologic_tumor_stage=='Stage IV')|(ajcc_pathologic_tumor_stage=='Stage IVA')|(ajcc_pathologic_tumor_stage=='Stage IVC')), "2", "1"))
data2 <- mutate(data2, ajcc_metastasis_pathologic_pm = ifelse((ajcc_metastasis_pathologic_pm=='M1'), "2", "1"))
data2 <- mutate(data2, residual_tumor = ifelse(((residual_tumor=='R1')|(residual_tumor=='R2')), "2", "1"))

data2$age_at_diagnosis<-factor(data2$age_at_diagnosis)
data2$ajcc_pathologic_tumor_stage<-factor(data2$ajcc_pathologic_tumor_stage)
data2$ajcc_metastasis_pathologic_pm<-factor(data2$ajcc_metastasis_pathologic_pm)
data2$residual_tumor<-factor(data2$residual_tumor)
data2$PI<-factor(data2$PI)



surv_object <- Surv(time = data2$OS.time/30, event = data2$vital_status)

#survival analysis: fits cox ph model to find HR for mean cut
fit1 <- survfit(surv_object~data2$PI>6); 
#summary(fit1);
ggsurvplot(fit1, data=data2,xlab = "Time (Months)")



fit1.coxph <- coxph(surv_object~data2$PI>6)

# summary(fit1.coxph);
first <- coef(summary(fit1.coxph))
#first
a=summary(fit1.coxph)
a$concordance[1]

print(c(nrow(data2),first[2],first[5],a$concordance[1],a$conf.int[3],a$conf.int[4],a$logtest[3]))

gg<-ggsurvplot(fit1, data=data2,xlab = "Time (Months)",pval = TRUE,legend.labs = c("Low Risk", "High Risk"),axes.offset=FALSE)
gg
ggsave('/Users/macbook/Desktop/THCA/riskgroups_KMplots/M1.jpg', plot = print(gg),height = 5 ,width = 7, dpi = 1000)

############################################################

# Fit a Cox proportional hazards model
fit.coxph <- coxph(surv_object ~ PI+ajcc_pathologic_tumor_stage+residual_tumor+age_at_diagnosis, 
                   data = data2)
ggforest(fit.coxph, data = data2)
