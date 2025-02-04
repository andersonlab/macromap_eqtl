
packages = c("ashr", "mashr","RColorBrewer", "gplots","gridExtra","spatstat","dplyr","corrplot","lattice","data.table","doMC"
             ,"directlabels","rmeta","pheatmap","fst")

## Now load or install&load all
package.check <- lapply(
  packages,
  FUN = function(x) {
    if (!require(x, character.only = TRUE)) {
      install.packages(x, dependencies = TRUE)
      library(x, character.only = TRUE)
    }
  }
)

# enable parallel processing for computationally intensive operations.
registerDoMC(cores = 5)

#source("~/myscratch/MacroMap/Analysis/eQTLs/Macromap_fds/analysis/scripts/flashR_cov.R")

# Create this function in order to calculate the covariance matrices with z-score 
# cov_pca_Z<-function (data, npc, subset = NULL) 
# {
#   assert_that(npc > 1)
#   assert_that(npc <= n_conditions(data))
#   if (is.null(subset)) {
#     subset = 1:n_effects(data)
#   }
#   data$Z<-data$Bhat/data$Shat 
#   
#   res.svd = svd(data$Z[subset, ], nv = npc, nu = npc)
#   f = res.svd$v
#   Ulist = cov_from_factors(t(f), "PCA")
#   d = diag(res.svd$d[1:npc])
#   Ulist = c(Ulist, list(tPCA = f %*% d^2 %*% t(f)/length(subset)))
#   return(Ulist)
# }
# environment(cov_pca_Z) <- asNamespace('mashr')

#m.ed<-readRDS(file="/nfs/users/nfs_n/np12/myscratch/MacroMap/Analysis/eQTLs/Freeze_Feb2019/analysis/mashR/rds_final/m.ed.RDS")
#m.ed.t<-readRDS(file="/nfs/users/nfs_n/np12/myscratch/MacroMap/Analysis/eQTLs/Freeze_Feb2019/analysis/mashR/rds_final/m.ed.t.RDS")
#U.ed<-readRDS("/nfs/users/nfs_n/np12/myscratch/MacroMap/Analysis/eQTLs/Freeze_Feb2019/analysis/mashR/rds_final/U.ed.RDS")

path<-"/lustre/scratch119/realdata/mdt3/teams/gaffney/np12/MacroMap/Analysis/eQTLs/Macromap_fds/analysis/mashR/analysis/analysisV2/8_flashR_best_eQTLs_random200K_commonbaseline/"

#########################################################################################################################################################################################
# Load the nominal data 
MergedDT<-read.fst("~/myscratch/MacroMap/Analysis/eQTLs/Macromap_fds/analysis/mashR/files/merged_nominal_final.fst",as.data.table=TRUE)
#MergedDT<-data.table(MergedDT)

cond_pcs<-read.table("~/myscratch/MacroMap/Analysis/eQTLs/Macromap_fds/analysis/eQTLs_per_TSS/1MB/cond.eqtls.txt")

conditions<-sapply(strsplit(as.character(cond_pcs$V1),"_1"),"[[",1)
pcs<-cond_pcs$V2


MergedDT_std.err<-MergedDT[,.SD,.SDcols=colnames(MergedDT)[grep ("std.err",colnames(MergedDT))]]
MergedDT_beta<-MergedDT[,.SD,.SDcols=colnames(MergedDT)[grep ("Beta",colnames(MergedDT))]]
MergedDT_nominalpval<-MergedDT[,.SD,.SDcols=colnames(MergedDT)[grep ("Nominal_pvalue",colnames(MergedDT))]]

rownames(MergedDT_std.err)<-MergedDT$gene_snp
rownames(MergedDT_beta)<-MergedDT$gene_snp
rownames(MergedDT_nominalpval)<-MergedDT$gene_snp

#########################################################################################################################################################################################
# load significant results per conditions 

for (i in 1:length(conditions)) {
  cat("\nRead Input data",conditions[i],"\n");
  cond<-paste0("~/myscratch/MacroMap/Analysis/eQTLs/Macromap_fds/",conditions[i],"/analysis/1MB/PC",pcs[i],"/",conditions[i],"_0.05_1MB_PC",pcs[i],".significant.stderr.txt")
  cond_name<-paste0("condition_",conditions[i],"_significant")
  assign(cond_name,fread(paste("cat",cond)))
}

conditions_significant.list<-NULL
for (i in 1:length(conditions)) {
  conditions_significant.list[[i]]<-cond_name<-paste0("condition_",conditions[i],"_significant")
}

# Assign gene_snp_chr_pos to the data tables, set the same key to all data tables, change colnames for the binding step 

for (i in 1:length(conditions_significant.list)) {
  dat<-get(conditions_significant.list[[i]]) # 
  cat("Creating key (gene_snp) for",conditions_significant.list[[i]],"\n")
  dat<-dat[,'gene_snp':=paste(Phenotype_ID,Best_variant_in_cis,Chromosome_var,Pos_variant,sep="_")]
  cat("Setting key (gene_snp) for",conditions_significant.list[[i]],"\n")
  setkey(dat,gene_snp)
  cat("changing colnames for",conditions_significant.list[[i]],"\n")
  colnames(dat)[1:length(colnames(dat))-1]<-paste(colnames(dat)[1:length(colnames(dat))-1],conditions[i],sep="_")
  cat("assigning data table to",conditions_significant.list[[i]],"\n")
  assign(conditions_significant.list[[i]],dat)
}

list.df = mget(ls(pattern = "*significant$"))

# gene snp as rownames 
list.df<-lapply(list.df,function (x)  { rownames(x) <-paste(x$Phenotype_ID,x$Best_variant_in_cis,x$Chromosome_phe,x$Pos_variant,sep="_"); x})

# all significant eQTLs for every conditions 
eQTLs_all_conditions<-unique(unlist(lapply(list.df,function (x) x$gene_snp )))
length(unique(substr(eQTLs_all_conditions,1,15)))  # number of unique genes that have an eQTL This is a result 

sign_eQTLs_betas<-MergedDT_beta[rownames(MergedDT_beta) %in% eQTLs_all_conditions]
sign_eQTLs_std.err<-MergedDT_std.err[rownames(MergedDT_std.err) %in% eQTLs_all_conditions]
sign_eQTLs_nominalpval<-MergedDT_nominalpval[rownames(MergedDT_nominalpval) %in% eQTLs_all_conditions]

conditions.order<-substr(colnames(MergedDT_std.err),9,length(colnames(MergedDT_std.err)))

#########################################################################################################################################################################################
# load best snp per gene for all conditions 
for (i in 1:length(conditions)) {
  cond<-paste0("~/myscratch/MacroMap/Analysis/eQTLs/Macromap_fds/",conditions[i],"/analysis/1MB/PC",pcs[i],"/1MB_",pcs[i],".GENE.stderr.txt")
  cond_name<-paste0("condition_",conditions[i],"_all_significant")
  assign(cond_name,fread(paste("cat",cond)))
}

conditions_all_significant.list<-NULL
for (i in 1:length(conditions)) {
  conditions_all_significant.list[[i]]<-cond_name<-paste0("condition_",conditions[i],"_all_significant")
}

# Assign gene_snp_chr_pos to the data tables, set the same key to all data tables, change colnames for the binding step 
# Be careful here the key is the gene name in order to merge the data 

for (i in 1:length(conditions_all_significant.list)) {
  dat<-get(conditions_all_significant.list[[i]]) # 
  cat("Creating key (gene_snp) for",conditions_all_significant.list[[i]],"\n")
  dat<-dat[,'gene_snp':=paste(Phenotype_ID,Best_variant_in_cis,Chromosome_var,Pos_variant,sep="_")]
  cat("Setting key (gene_snp) for",conditions_all_significant.list[[i]],"\n")
  setkey(dat,Phenotype_ID)
  dat$z_score<- dat$Beta_regression/dat$std.err 
  dat<-dat[!is.na(dat$Corrected_pvalue)]
  cat("changing colnames for",conditions_all_significant.list[[i]],"\n")
  colnames(dat)[2:length(colnames(dat))]<-paste(colnames(dat)[2:length(colnames(dat))],conditions[i],sep="_")
  
  cat("assigning data table to",conditions_all_significant.list[[i]],"\n")
  assign(conditions_all_significant.list[[i]],dat)
}


list.df = mget(ls(pattern = "*all_significant$"))

MergedDT_all_significant = Reduce(function(...) merge(..., all = FALSE), list.df) 

MergedDT_all_significant_zscore<-MergedDT_all_significant[,.SD,.SDcols=colnames(MergedDT_all_significant)[grep ("z_score",colnames(MergedDT_all_significant))]]

# Find a better solution for this 

best_snp_gene<-c(condition_CIL_24_all_significant$gene_snp[apply(abs(MergedDT_all_significant_zscore),1,which.max)==1],
                 condition_CIL_6_all_significant$gene_snp[apply(abs(MergedDT_all_significant_zscore),1,which.max)==2],
                 condition_Ctrl_24_all_significant$gene_snp[apply(abs(MergedDT_all_significant_zscore),1,which.max)==3],
                 condition_Ctrl_6_all_significant$gene_snp[apply(abs(MergedDT_all_significant_zscore),1,which.max)==4],
                 condition_IFNB_24_all_significant$gene_snp[apply(abs(MergedDT_all_significant_zscore),1,which.max)==5],
                 condition_IFNB_6_all_significant$gene_snp[apply(abs(MergedDT_all_significant_zscore),1,which.max)==6],
                 condition_IFNG_24_all_significant$gene_snp[apply(abs(MergedDT_all_significant_zscore),1,which.max)==7],
                 condition_IFNG_6_all_significant$gene_snp[apply(abs(MergedDT_all_significant_zscore),1,which.max)==8],
                 condition_IL4_24_all_significant$gene_snp[apply(abs(MergedDT_all_significant_zscore),1,which.max)==9],
                 condition_IL4_6_all_significant$gene_snp[apply(abs(MergedDT_all_significant_zscore),1,which.max)==10],
                 condition_LIL10_24_all_significant$gene_snp[apply(abs(MergedDT_all_significant_zscore),1,which.max)==11],
                 condition_LIL10_6_all_significant$gene_snp[apply(abs(MergedDT_all_significant_zscore),1,which.max)==12],
                 condition_MBP_24_all_significant$gene_snp[apply(abs(MergedDT_all_significant_zscore),1,which.max)==13],
                 condition_MBP_6_all_significant$gene_snp[apply(abs(MergedDT_all_significant_zscore),1,which.max)==14],
                 condition_P3C_24_all_significant$gene_snp[apply(abs(MergedDT_all_significant_zscore),1,which.max)==15],
                 condition_P3C_6_all_significant$gene_snp[apply(abs(MergedDT_all_significant_zscore),1,which.max)==16],
                 condition_PIC_24_all_significant$gene_snp[apply(abs(MergedDT_all_significant_zscore),1,which.max)==17],
                 condition_PIC_6_all_significant$gene_snp[apply(abs(MergedDT_all_significant_zscore),1,which.max)==18],
                 condition_Prec_D0_all_significant$gene_snp[apply(abs(MergedDT_all_significant_zscore),1,which.max)==19],
                 condition_Prec_D2_all_significant$gene_snp[apply(abs(MergedDT_all_significant_zscore),1,which.max)==20],
                 condition_R848_24_all_significant$gene_snp[apply(abs(MergedDT_all_significant_zscore),1,which.max)==21],
                 condition_R848_6_all_significant$gene_snp[apply(abs(MergedDT_all_significant_zscore),1,which.max)==22],
                 condition_sLPS_24_all_significant$gene_snp[apply(abs(MergedDT_all_significant_zscore),1,which.max)==23],
                 condition_sLPS_6_all_significant$gene_snp[apply(abs(MergedDT_all_significant_zscore),1,which.max)==24])


bhat_matrix_best_snp_gene<-as.data.frame(MergedDT_beta[(rownames(MergedDT_beta) %in% best_snp_gene)])
rownames(bhat_matrix_best_snp_gene)<-rownames(MergedDT_beta)[(rownames(MergedDT_beta) %in% best_snp_gene)]

shat_matrix_best_snp_gene<-as.data.frame(MergedDT_std.err[(rownames(MergedDT_std.err) %in% best_snp_gene)])
rownames(shat_matrix_best_snp_gene)<-rownames(MergedDT_std.err)[(rownames(MergedDT_std.err) %in% best_snp_gene)]


data_to_test = mash_set_data(as.matrix(bhat_matrix_best_snp_gene),as.matrix(shat_matrix_best_snp_gene))
# For each condition I pick the best SNP-GENE pair that shows the highest z score for all the genes we need all the nominals for this 


best_snp_gene.all.conditions<-as.data.frame(MergedDT[(MergedDT$gene_snp %in% best_snp_gene)])

# calclulate z score and pick the condition that the gene snp pair shows higher z score 
#genes.significant_per_condition<-get_n_significant_conditions(m.ed.t,sig_fn=get_lfsr,thresh=0.05)
#best_snp_gene.all.conditions_all.info<-best_snp_gene.all.conditions[best_snp_gene.all.conditions$gene_snp %in% names(genes.significant_per_condition[genes.significant_per_condition %in% c(1:24)]),]
#rownames(best_snp_gene.all.conditions_all.info)<-best_snp_gene.all.conditions_all.info$gene_snp

#z.scores_picks<-apply(abs(best_snp_gene.all.conditions_all.info[grep("Beta",names(best_snp_gene.all.conditions_all.info))]/best_snp_gene.all.conditions_all.info[grep("std",names(best_snp_gene.all.conditions_all.info))]),1,which.max)
#########################################################################################################################################################################################
list_signif<-mget(grep("(^condition_)(?!.*all.*)",ls(),value=TRUE,perl=TRUE))

eQTLs_10FDR_250K<-t(as.data.frame(lapply(list_signif, function (x) length(x$gene_snp))))
rownames(eQTLs_10FDR_250K)<-conditions.order
colnames(eQTLs_10FDR_250K)<-"eQTLs"

eGenes_1<-eQTLs_10FDR_250K[order(eQTLs_10FDR_250K[,1],decreasing=T)]

donors_1<-read.table("~/myscratch/MacroMap/Analysis/eQTLs/Macromap_fds/analysis/eQTLs_per_TSS/power.conditions.txt",h=F)
donors_1<-donors_1[order(eQTLs_10FDR_250K[,1],decreasing=T),]
labels_eQTLs<-rownames(eQTLs_10FDR_250K)[order(eQTLs_10FDR_250K[,1],decreasing=T)]


#########################################################################################################################################################################################
### This is the actual analysis summarised 
set.seed(1234)  
rnd<-sample(nrow(MergedDT_beta),size=300000,replace=FALSE) # random number of pairs in order to fit the model (both significant and not significant)
rnd2<-sample(which(apply(apply(MergedDT_beta[rnd,],1, function (x) x==0),2,sum)==0),size=200000,replace=FALSE) # This is in order to have conditions that beta may be 0 
data.mash.bhat<-MergedDT_beta[rnd][rnd2]
data.mash.std<-MergedDT_std.err[rnd][rnd2]

# Random data equal number of SNPs per gene and distance to TSS 
#random_200K<-readRDS("~/myscratch/MacroMap/MacroMap/Analysis/eQTLs/Macromap_fds/analysis/mashR/files/final_dataset_200K.RDS")
#rnd<- which(MergedDT$gene_snp %in% unique(random_200K$gene_snp))

#data.mash.bhat<-MergedDT_beta[rnd]
#data.mash.std<-MergedDT_std.err[rnd]

data.mash = mash_set_data(as.matrix(data.mash.bhat[,-c(4)]),as.matrix(data.mash.std[,-c(4)])) 
Vhat = estimate_null_correlation_simple(data.mash)
data.mash.v =mash_update_data(data.mash,V=Vhat,ref=3)

#apply(apply(data.mash$Bhat,2, function (x) x==0),2, table) # I need to find a seed where there is no zero beta 
#apply(apply(data.mash$Shat,2, function (x) x==0),2, table) # I need to find a seed where there is no zero beta 

data = mash_set_data(as.matrix(sign_eQTLs_betas[,-c(4)]),as.matrix(sign_eQTLs_std.err[,-c(4)])) # load data that are significant (best snp per gene for every condition) 
data = mash_update_data(data,V=Vhat)

#U.pca = cov_pca(data,5)    # make data driven covariance matrices based on Bhat and  PCs 
#U.pca.z<-cov_pca_Z(data,24) # make data driven covariance matrices based on z score and  PCs (cov_pca_Z is my function)

data_to_test = mash_set_data(as.matrix(bhat_matrix_best_snp_gene[,-c(4)]),as.matrix(shat_matrix_best_snp_gene[,-c(4)])) # model cannot fit all these data points ? x
data_to_test = mash_update_data(data_to_test,V=Vhat,ref=3)

U.f = cov_flash(data_to_test)

U.c = cov_canonical(data_to_test)  # make canonical covariance matrices 

# extreme deconvolution 
U.ed = cov_ed(data_to_test, U.f) # refine covariance matrices 
#U.ed<-readRDS("~/myscratch/MacroMap/MacroMap/Analysis/eQTLs/Macromap_fds/analysis/mashR/files/flashR/U.ed.RDS")


m.ed = mash(data.mash.v, c(U.c,U.ed), outputlevel = 1) # fit model on order to estimate p (20000 couple of min ) # 200,000 around 3.5 hours 
m.ed.t<-mash(data_to_test,g=m.ed$fitted_g,fixg=TRUE) # fit the data to test 

saveRDS(m.ed,file=paste0(path,"m.ed.RDS"))

saveRDS(m.ed.t,file=paste0(path,"m.ed.t.RDS"))

saveRDS(U.ed,paste0(path,"U.ed.RDS"))

#########################################################################################################################################################################################
saveRDS(best_snp_gene.all.conditions,paste0(path,"best_snp_gene.all.conditions.RDS"))
saveRDS(bhat_matrix_best_snp_gene,paste0(path,"bhat_matrix_best_snp_gene.RDS"))
saveRDS(eQTLs_all_conditions,paste0(path,"eQTLs_all_conditions.RDS"))
saveRDS(shat_matrix_best_snp_gene,paste0(path,"shat_matrix_best_snp_gene.RDS"))
saveRDS(sign_eQTLs_betas,paste0(path,"sign_eQTLs_betas.RDS"))
saveRDS(sign_eQTLs_nominalpval,paste0(path,"sign_eQTLs_nominalpval.RDS"))
saveRDS(sign_eQTLs_std.err,paste0(path,"sign_eQTLs_std.err.RDS"))
saveRDS(conditions.order,paste0(path,"conditions.order.RDS"))

saveRDS(eQTLs_10FDR_250K,paste0(path,"eQTLs_10FDR_250K.RDS"))
saveRDS(eGenes_1,paste0(path,"eGenes_1.RDS"))
saveRDS(donors_1,paste0(path,"donors_1.RDS"))
saveRDS(labels_eQTLs,paste0(path,"labels_eQTLs.RDS"))
saveRDS(list_signif,paste0(path,"list_signif.RDS"))
saveRDS(list.df,paste0(path,"list.df.RDS"))
saveRDS(Vhat,paste0(path,"Vhat.RDS"))
