library(hgu133a.db)
library(annotate)
x <- hgu133aSYMBOL
# Get the probe identifiers - gene symbol mappings
mapped_probes <- mappedkeys(x)
# Convert to a dataframe
genesym.probeid <- as.data.frame(x[mapped_probes])
#head(genesym.probeid)

#genes <- c("ANXA1", "APCS", "C6", "COL14A1", "CYB5R3", "FAM82A1", "GBE1", "GNAI1", "HLA-DQA1", "IL16", "LTBP2", "MMS19", "MXRA5", "SEC13", "TUBAL3", "WDR1")
#genes <- scan(file = '/Users/macbook/Desktop/UCEC/UCEC_genes.csv', what = 'character', sep = ',')

genes<-c('TYRO3')
results <- genesym.probeid[which(genesym.probeid$symbol %in% genes),]
results
