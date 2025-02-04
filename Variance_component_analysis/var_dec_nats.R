library("RColorBrewer")
library("ggplot2")
library("RColorBrewer")
source("~/bin/Rscripts/lm.pval.R")
library("qvalue")
library("lme4")
library("MuMIn")
library("DESeq2")
library("GenABEL")


#dds_expression_norm_cqn_quantile<-readRDS("/lustre/scratch117/cellgen/team170/np12/MacroMap/Analysis/QC/RDS_expression_covariates/Freeze_Aug2019/dds_expression_norm_cqn_quantile_freezeAug2019.rds")
covariates<-readRDS("/lustre/scratch117/cellgen/team170/np12/MacroMap/Analysis/QC/RDS_expression_covariates/Macromap_final_expression/protein_coding_lincRNA/covariates_macromap_fds.rds")
tpms.log10<-readRDS("/lustre/scratch117/cellgen/team170/np12/MacroMap/Analysis/QC/RDS_expression_covariates/Macromap_final_expression/protein_coding_lincRNA/RDS_log10_TPMs_filtd_macromap_fds.rds")

source("/nfs/users/nfs_n/nk5/Project/R/Util/LFLMM.R")

tpms.log10<-tpms.log10[,-c(1:6)]

cov_model<-c("RunID","Donor","Stimulus_Hours","Sex","Library_prep","Date_thawed","Passage_number_at_thawing","Passage_EB_formation","IPSCs_culture_time","Total_Harvests","Differentiation_time_No_Days",
             "Purity_result_per","Estimated_cell_diameter","Cell_diameter_SD","Differentiation_media","SeV_Result")

cov_to_test<-covariates[cov_model]

cov_to_test <- data.frame(lapply(cov_to_test, as.character), stringsAsFactors=FALSE)

LFLMM.results<-LFLMM(as.matrix(tpms.log10),cov_to_test,ITRMAX=100)

saveRDS(LFLMM.results,"LFLMM.results.RDS")