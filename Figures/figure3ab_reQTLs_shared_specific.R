library(mashr)
library(data.table)
library(reshape)
library(ggplot2)
library(RColorBrewer)

################################################################################################################################################################################################################################
# load data and mash model 
path_to_read_files<-"/nfs/users/nfs_n/np12/myscratch/MacroMap/Data/MacroMap_manuscript_figures/data/"
path_to_plot<-"/nfs/users/nfs_n/np12/myscratch/MacroMap/Data/MacroMap_manuscript_figures/plots/"

#path<-"/nfs/users/nfs_n/np12/myscratch/MacroMap/Analysis/eQTLs/Macromap_fds/analysis/mashR/analysis/analysisV2/8_flashR_best_eQTLs_random200K_commonbaseline/"
#path_to_plot<-"/nfs/users/nfs_n/np12/myscratch/MacroMap/Analysis/eQTLs/Macromap_fds/analysis/mashR/analysis/8_flashR_condition_by_condition/8_flashR_condition_by_condition_baseline/"
#path_to_save_files<-"/nfs/users/nfs_n/np12/myscratch/MacroMap/Analysis/eQTLs/Macromap_fds/analysis/mashR/analysis/8_flashR_condition_by_condition/8_flashR_condition_by_condition_baseline/"

eQTLs_all_conditions<-readRDS(paste0(path_to_read_files,"eQTLs_all_conditions.RDS"))
sign_eQTLs_betas<-readRDS(paste0(path_to_read_files,"sign_eQTLs_betas.RDS"))
sign_eQTLs_std.err<-readRDS(paste0(path_to_read_files,"sign_eQTLs_std.err.RDS"))
Vhat<-readRDS(paste0(path_to_read_files,"Vhat.RDS"))
rownames(sign_eQTLs_betas)<-read.table(paste0(path_to_read_files,"rownames_sign_eQTLs_betas.txt"),h=F)$V1
rownames(sign_eQTLs_std.err)<-read.table(paste0(path_to_read_files,"rownames_sign_eQTLs_betas.txt"),h=F)$V1

m.ed<-readRDS(paste0(path_to_read_files,"m.ed.RDS"))

data_to_test_all_eQTLs = mash_set_data(as.matrix(sign_eQTLs_betas[,-c(4)]),as.matrix(sign_eQTLs_std.err[,-c(4)])) # load data that are significant (best snp per gene for every condition) 
data_to_test_all_eQTLs = mash_update_data(data_to_test_all_eQTLs,V=Vhat,ref=3)

# Run this only once 
#m.ed.t_all_eQTLs<-mash(data_to_test_all_eQTLs,g=m.ed$fitted_g,fixg=TRUE) # fit the data to test 
#saveRDS(m.ed.t_all_eQTLs,file=paste0(path_to_read_files,"m.ed.t_all_eQTLs"))

# now load the model 
m.ed.t_all_eQTLs<-readRDS(paste0(path_to_read_files,"m.ed.t_all_eQTLs.RDS"))

cond_pcs<-read.table(paste0(path_to_read_files,"cond.eqtls.txt"))

conditions<-sapply(strsplit(as.character(cond_pcs$V1),"_1"),"[[",1)
pcs<-cond_pcs$V2
################################################################################################################################################################################################################################
# load signficant eQTLs per condition 
for (i in 1:length(conditions)) {
  cat("\nRead Input data",conditions[i],"\n");
  cond<-paste0("/nfs/users/nfs_n/np12/myscratch/MacroMap/Data/MacroMap_manuscript_figures/data/eQTLs_per_condition/",conditions[i],"_0.05_1MB_PC",pcs[i],".significant.stderr.txt")
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

rownames(m.ed.t_all_eQTLs$result$lfsr) <-rownames(sign_eQTLs_betas)

list_condition_by_condition_eQTLs_mashR<-lapply(list.df,function (x) m.ed.t_all_eQTLs$result$lfsr[rownames(m.ed.t_all_eQTLs$result$lfsr) %in% x$gene_snp,] )
list_condition_by_condition_eQTLs_mashR<-list_condition_by_condition_eQTLs_mashR[-c(3,4)]



response_QTLs<-NULL
for (i in 1:22) {
  response_QTLs[[i]]<-list_condition_by_condition_eQTLs_mashR[[i]][,i] <=0.05
} 

list.df2<-list.df[-c(3,4)]

shared<-NULL
response<-NULL

for (i in 1:22) {
  a<-list.df2[[i]]
  b<-response_QTLs[[i]]
  
  shared<-rbindlist(list(a[!b],shared),use.names=FALSE)
  response<-rbindlist(list(a[b],response),use.names=FALSE)
  
}

shared<-shared[!shared$Best_variant_in_cis_sLPS_6 %in% response$Best_variant_in_cis_sLPS_6,] # remove response from shared data 
shared<-shared[!duplicated(shared$gene_snp)]
response<-response[!duplicated(response$gene_snp)]
################################################################################################################################################################################################################################
# I tested 2 mash models first to define respone eQTLs and a second that ctrl are included as a condition 
# The second has been modelled with all conditions together 

# In the paper I use only the lfsr results which are based on the baseline model 
list_condition_by_condition_eQTLs_mashR_betas<-lapply(list.df,function (x) m.ed.t_all_eQTLs$result$PosteriorMean[rownames(m.ed.t_all_eQTLs$result$lfsr) %in% x$gene_snp,] )
list_condition_by_condition_eQTLs_mashR_betas<-list_condition_by_condition_eQTLs_mashR_betas[-c(3,4)]


shared_mashR_betas<-NULL
respone_mashR_betas<-NULL
cond.specific.mash.no.lfsr <-NULL
cond.specific.mash.no.lfsr.shared_per_con<-NULL

lfsr.mash    <- m.ed.t_all_eQTLs$result$lfsr # for the lsfr is better to test whether the response QTLs is significant compared to the baseline model 
colnames(lfsr.mash) <- unlist((strsplit(sapply(strsplit(colnames(lfsr.mash),"Beta_regression_","[[",2),"[[",2),"-")))
rownames(lfsr.mash)<-rownames(m.ed.t_all_eQTLs$result$lfsr) # This has been checked names are in same order 

# Identify how many reQTLs are shared or condition specific 
for (i in 1:22) {
  a<-list_condition_by_condition_eQTLs_mashR_betas[[i]][,i]
  b<-response_QTLs[[i]]
  names(a)<-names(b)
  shared_mashR_betas[[i]]<-a[!b]
  respone_mashR_betas[[i]]<-a[b]
  
  c<-which(rownames(lfsr.mash) %in% names(respone_mashR_betas[[i]]))
  
  lfsr.mash.f<-lfsr.mash[c,]
  sign.based_lfsr.f <- which(rowSums(lfsr.mash.f < 0.05) == 1)
  cond.specific.mash.no.lfsr[i]<-length(sign.based_lfsr.f)
  cond.specific.mash.no.lfsr.shared_per_con[[i]]<-table(factor(rowSums(lfsr.mash.f < 0.05),levels = 1:22))
}

m<-do.call(rbind,cond.specific.mash.no.lfsr.shared_per_con)
rownames(m)<-colnames(lfsr.mash)
mean(100*m[,1]/rowSums(m)) # 1.11% abstract paper quite rare 
#####################################################################################v###################################
m<-do.call(rbind,cond.specific.mash.no.lfsr.shared_per_con)
rownames(m)<-colnames(lfsr.mash)
m1<-m/rowSums(m)*100

a<-apply(m1,2,mean) # paper 13 states max sharing of reQTLs 

colors<-brewer.pal(8, "Set1")[c(1,2)]

pdf(paste0(path_to_plot,"Fig3b_Number_of_conditions_response_qtl_shared_boxplot_version.pdf"),useDingbats=FALSE,width = 8,height = 7)
par(mar = c(5, 6, 3, 5))
boxplot(m1,yaxt = "n",xaxt = "n",xlab="Number of stimulated conditions a reQTL is detected",col=colors[2],cex.lab=1.5)
axis(1, at=1:22, labels=1:22,cex.axis=0.75)
axis(side=4)
mtext("Detection rate", side = 4, line = 3,cex=1.5)
par(new = TRUE)
plot(cumsum(a), type="b",xaxt = "n",xlab="",col=colors[1],ylab="Cumulative percentage of reQTLs detected\nacross stimulated conditions",pch=19,cex=1.2,cex.lab=1.5)
abline(h=seq(0,100,by=10),col="gray",lty=3)
dev.off()

##########################################################################################################################################################################
m_ggplot<-melt(m1)
m_ggplot$X1<-factor(m_ggplot$X1,rev(levels(m_ggplot$X1)))
m_ggplot$X2<-factor(m_ggplot$X2)
m_ggplot$facet<-ifelse(m_ggplot$X2==1,FALSE,TRUE)

p <- ggplot(m_ggplot,aes(X2,X1)) + geom_tile(aes(fill=value),color = "white") +
  guides(fill=guide_colorbar("% of reQTLs")) +scale_fill_gradientn(colors=col_heat,guide="colorbar") 

p<-p +theme_bw(16) + theme(legend.justification = "center",legend.position="bottom") + xlab("Number of stimulated conditions a reQTL is detected") + ylab("Stimulated conditions")

p<-p + facet_grid(.~facet, space = "free", scales = "free_x") +theme(strip.text.x = element_blank())


pdf(paste0(path_to_plot,"/fig3a_Number_of_conditions_response_qtl_shared_ggplot_version.pdf"),useDingbats=FALSE)
print(p)
dev.off()
##########################################################################################################################################################################
#
