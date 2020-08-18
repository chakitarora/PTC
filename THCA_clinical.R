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
beta_file='/Users/macbook/Desktop/APOP_PAN_CANCER/HR_APOP/'+cancer+'_HRmedian_raw.csv'
tcga_datafile="/Users/macbook/Desktop/THCA/"+cancer+".RDS"

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

######################################## general ############################
#use | for what u want and & for what you don't want
#data2<- subset(data1, ((race != '[NA]')&(race != '[Unknown]')&(race != '[Not Evaluated]')))
#data2<- subset(data1, race=='WHITE')
#data2<- subset(data1, (extrathyroidal_extension=='Minimal (T3)')|(extrathyroidal_extension=='Very Advanced (T4b)')|(extrathyroidal_extension=='Moderate/Advanced (T4a)'))
data2<- subset(data1, (laterality=='Right lobe')|(laterality=='Left lobe'))
#data2<- subset(data1, ((ajcc_nodes_pathologic_pn == 'N1')|(ajcc_nodes_pathologic_pn == 'N1b')|(ajcc_nodes_pathologic_pn == 'N1a')))


library("survMisc")
fit1<-Surv(time = data2$OS.time/30, event = data2$vital_status)
fit2 <- coxph(fit1~age_at_diagnosis, data = data2)
cut <- cutp(fit2)$age_at_diagnosis

g=cut[1,1][,1]

rm(fit1)

# x=(0.06066)*data1$PLP2-(0.35719)*data1$TGFBR3-(0.19667)*data1$LYVE1-(0.10089)*data1$FXYD6-(0.01634)*data1$FABP4 
# x>median(x)
surv_object <- Surv(time = data2$OS.time/30, event = data2$vital_status)

#survival analysis: fits cox ph model to find HR for mean cut
fit1 <- survfit(surv_object~data2$age_at_diagnosis>60); 
#summary(fit1);
ggsurvplot(fit1, data=data2,xlab = "Time (Months)")



fit1.coxph <- coxph(surv_object~data2$age_at_diagnosis>60)

# summary(fit1.coxph);
first <- coef(summary(fit1.coxph))
#first
a=summary(fit1.coxph)
a$concordance[1]

print(c(nrow(data2),first[2],first[5],a$concordance[1],a$conf.int[3],a$conf.int[4],a$logtest[3]))

gg<-ggsurvplot(fit1, data=data2,xlab = "Time (Months)",pval = TRUE,legend.labs = c("Low Risk", "High Risk"),axes.offset=FALSE)
gg
ggsave('/Users/macbook/Desktop/THCA/riskgroups_KMplots/bilateral.jpg', plot = print(gg),height = 5 ,width = 7, dpi = 1000)

######################################## voting based #######################################
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

if (sum(PI)!=0)
{  g=v/2
surv_object <- Surv(time = data1$OS.time/30, event = data1$vital_status)

#survival analysis: fits cox ph model to find HR for mean cut
fit1 <- survfit(surv_object~PI>g); 
#summary(fit1);
ggsurvplot(fit1, data=data1)


fit1.coxph <- coxph(surv_object~PI>g)
# summary(fit1.coxph);
first <- coef(summary(fit1.coxph))
#first
a=summary(fit1.coxph)
a$concordance[1]
#max U and Q, min p value
}

###################################################### prognostic index ######################  
PI=0
for(i in seq(from=1, to=length(genes), by=1))#no of cancers
{
  b=beta[beta$Gene==genes[i],][2][,1]
  PI=PI+b*select(data1, genes[i])[,1]
}

surv_object <- Surv(time = data1$OS.time/30, event = data1$vital_status)

#survival analysis: fits cox ph model to find HR for mean cut
fit1 <- survfit(surv_object~PI>median(PI)); 
#summary(fit1);
#ggsurvplot(fit1, data=test1)


fit1.coxph <- coxph(surv_object~PI>median(PI))
# summary(fit1.coxph);
first <- coef(summary(fit1.coxph))
#first
a=summary(fit1.coxph)

############################## cutp : can be used for both ########################## 
library("survMisc")
fit1<-Surv(time = data1$OS.time/30, event = data1$vital_status)
fit2 <- coxph(fit1~PI, data = data1)
cut <- cutp(fit2)$PI

g=cut[1,1][,1]

rm(fit1)

surv_object <- Surv(time = data1$OS.time/30, event = data1$vital_status)

#survival analysis: fits cox ph model to find HR for mean cut
fit1 <- survfit(surv_object~PI>=g$PI); 
#summary(fit1);
ggsurvplot(fit1, data=data1)


fit1.coxph <- coxph(surv_object~PI>=g$PI)
# summary(fit1.coxph);
first <- coef(summary(fit1.coxph))
#first
a=summary(fit1.coxph)



ggsurv<-ggsurvplot(
  fit1,                            # survfit object with calculated statistics.
  data = data1,                   # data used to fit survival curves.
  risk.table = TRUE,               # show risk table.
  #title=" \t \t BCL2+BCLXL-BAX-BAK ",
  conf.int = FALSE,                 # show confidence intervals for
  # palette = c("#E7B800", "#2E9FDF"),
  xlim = c(0,150),                  # present narrower X axis, but not affect
  #ylim=c(0.1,1),
  xlab = "Time (Months)",         # customize X axis label.
  # break.time.by = 20,             # break X axis in time intervals by 500.
  ggtheme = theme_light(),         # customize plot and risk table with a theme.
  risk.table.height = 0.25,         # the height of the risk table
  risk.table.y.text = FALSE,        # show bars instead of names in text annotations
  # # in legend of risk table.
  conf.int.style = "step",         # customize style of confidence intervals
  legend=c(0.8,0.8),
  legend.labs = c("High Risk", "Low Risk"),  # change legend labels.
  # ncensor.plot = TRUE,               # plot the number of censored subjects at time t
  # ncensor.plot.height = 0.25,
  # pval = TRUE,                       # show p-value of log-rank test.
  axes.offset=FALSE,
)
ggsurv$plot <- ggsurv$plot+ 
  ggplot2::annotate("label", x = 30, y = 0.25, # x and y coordinates of the text
                    label = "HR = 2.94, p-val= 6.54e-10 \n 95%CI(2.09 - 4.15)
Wald p= 7e-10 & logrank, p= 9e-11", size = 4)
ggsurv


#fit1.coxph$concordance[6]