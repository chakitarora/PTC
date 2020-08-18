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



beta_file='/Users/macbook/Desktop/UCEC/UCEC-final.csv'

data2<-read.csv(file=beta_file,header =TRUE, sep = ",", dec = ".")
#data2$OS.time=data2$OS.time_ACTUAL

# library("survMisc")
# fit1<-Surv(time = data2$OS.time/30, event = data2$vital_status)
# fit2 <- coxph(fit1~OS.time_Linear, data = data2)
# cut <- cutp(fit2)$OS.time_Linear
# 
# g=cut[1,1][,1]
# 
# rm(fit1)

val=data2$PKD1L3
#print(median(val))

gene='PKD1L3'

surv_object <- Surv(time = data2$OS.time/30, event = data2$vital_status)

#survival analysis: fits cox ph model to find HR for mean cut
fit1 <- survfit(surv_object~val>median(val)); 
#summary(fit1);
ggsurvplot(fit1, data=data2,xlab = "Time (Months)")



fit1.coxph <- coxph(surv_object~val>median(val))
# summary(fit1.coxph);
first <- coef(summary(fit1.coxph))
#first
a=summary(fit1.coxph)
a$concordance[1]

print(c(first[2],first[5],a$concordance[1],a$conf.int[3],a$conf.int[4],a$logtest[3],median(val)))

gg<-ggsurvplot(fit1, data=data2,xlab = "Time (Months)",pval = TRUE,legend.labs = c("<=median", ">median"),axes.offset=FALSE)
gg
ggsave('/Users/macbook/Desktop/UCEC/GENE_KM'+gene+'.jpg', plot = print(gg),height = 5 ,width = 7, dpi = 1000)
