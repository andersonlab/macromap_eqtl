library("DESeq2")
library("sva")

################################################################################################################################################################################################
# load data 
dds_expression<-readRDS("/lustre/scratch117/cellgen/team170/np12/MacroMap/Analysis/QC/RDS_expression_covariates/Macromap_final_expression/protein_coding_lincRNA/dds_expression_macromap_fds.rds")
covariates<-readRDS("/lustre/scratch117/cellgen/team170/np12/MacroMap/Analysis/QC/RDS_expression_covariates/Macromap_final_expression/protein_coding_lincRNA/covariates_macromap_fds.rds")


dds_expression$Stimulus_Hours<-as.factor(dds_expression$Stimulus_Hours)
dds_expression$Donor<-as.factor(dds_expression$Donor)
dds_expression$Differentiation_time_No_Days <-as.factor(dds_expression$Differentiation_time_No_Days)


dds_expression <- estimateSizeFactors(dds_expression)

dat  <- counts(dds_expression, normalized = TRUE)
idx  <- rowMeans(dat) > 1
dat  <- dat[idx, ]
mod  <- model.matrix(~ Stimulus_Hours , colData(dds_expression))
mod0 <- model.matrix(~   1, colData(dds_expression))
svseq <- svaseq(dat, mod, mod0, n.sv = 10)
saveRDS(svseq,"svaseq_10.RDS")



