
args = commandArgs(trailingOnly=TRUE) # input run

expression<-readRDS("/lustre/scratch117/cellgen/team170/np12/MacroMap/Analysis/QC/RDS_expression_covariates/Macromap_final_expression/protein_coding_lincRNA/RDS_TPMs_filtd_macromap_fds.rds")
covariates<-readRDS("/lustre/scratch117/cellgen/team170/np12/MacroMap/Analysis/QC/RDS_expression_covariates/Macromap_final_expression/protein_coding_lincRNA/covariates_macromap_fds.rds")

covariates <- data.frame(lapply(covariates, as.character), stringsAsFactors=FALSE)

expression_ctrl6_to_map<-data.frame(expression[,c(1:6)],expression[-c(1:6)][covariates$Stimulus_Hours==args[1]])

covariates_ctrl6_to_map<-covariates[covariates$Stimulus_Hours==args[1],]

# when you map ids be sure to check for matches between RNA-seq and genotypes

colnames(expression_ctrl6_to_map)[-c(1:6)]<-covariates_ctrl6_to_map$HipsciID


# This is a crucial step to match vcf file
write.table(expression_ctrl6_to_map,file=paste0("/nfs/users/nfs_n/np12/myscratch/MacroMap/Analysis/eQTLs/Macromap_fds/",args[1],"/expression/exp_TPMs_fild_",args[1],".txt"),
            quote=F,col.names=T,row.names=F,sep="\t")
write.table(covariates_ctrl6_to_map,file=paste0("/nfs/users/nfs_n/np12/myscratch/MacroMap/Analysis/eQTLs/Macromap_fds/",args[1],"/covariates/covariates_",args[1],".txt"),
            quote=F,col.names=T,row.names=F,sep="\t")
setwd(paste0("/nfs/users/nfs_n/np12/myscratch/MacroMap/Analysis/eQTLs/Macromap_fds/",args[1],"/expression/"))

comm <- paste0("Rscript /lustre/scratch117/cellgen/team170/np12/MacroMap/Analysis/eQTLs/Macromap_fds/scripts/create_input_QTL.R exp_TPMs_fild_",args[1],".txt autosomesX QTLtools")
system (sprintf(comm))
