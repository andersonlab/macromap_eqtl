# plot examples of response QTLs that are also differentially expressed 

library(data.table)
library(RColorBrewer)
library(ggplot2)
library(DESeq2)
library(dplyr)

####################################################################################################################################################################################
# step 1 find the genes/snps to plot 

path_to_read_files<-"/nfs/users/nfs_n/np12/myscratch/MacroMap/Data/MacroMap_manuscript_figures/data/"
path_to_plot<-"/nfs/users/nfs_n/np12/myscratch/MacroMap/Data/MacroMap_manuscript_figures/plots/"
response_QTLs<-readRDS(file=paste0(path_to_read_files,"response_per_conditions.RDS"))

conditions<-c("CIL_24","CIL_6","Ctrl_24","Ctrl_6","IFNB_24","IFNB_6","IFNG_24","IFNG_6","IL4_24","IL4_6","LIL10_24","LIL10_6","P3C_24","P3C_6","R848_24",
              "R848_6","sLPS_24","sLPS_6","PIC_24","PIC_6","MBP_24","MBP_6","Prec_D0","Prec_D2")

conditions<-conditions[order(conditions)]
names(response_QTLs)<-conditions[-c(3,4)]

# load diff exprresion 
path='/nfs/users/nfs_n/np12/myscratch/MacroMap/Data/MacroMap_manuscript_figures/data/diff_exp/'
filelist = list.files(path= "/nfs/users/nfs_n/np12/myscratch/MacroMap/Data/MacroMap_manuscript_figures/data/diff_exp",pattern="^All.*.txt",full.names = F)
datalist = lapply(filelist, function(x)read.table(paste0(path,x), header=T)) 
names(datalist)[1:4]<-sapply(strsplit(filelist[c(1:4)],"_"), function (x) {paste(x[2],x[3],x[4],x[5],x[6],sep="_")} )
names(datalist)[5:24] <-sapply(strsplit(filelist[-c(1:4)],"_"), function (x) {paste(x[2],x[3],x[4],strsplit(x[5],".txt")[1],sep="_")} )
datalist<-datalist[-c(3,4)]

names(datalist)<-sapply(strsplit(names(datalist),"vs_"),"[[",2)
datalist<-datalist[order(names(datalist))]


response_diff_genes<-NULL
shared_diff_QTLs_genes<-NULL
up_down_regulation<-NULL
response_diff_genes_stats<-NULL
response_no_diff_genes_stats<-NULL

for (i in 1:length(response_QTLs)) {
  x<-response_QTLs[[i]]
  diff_exp_genes<-datalist[[i]]
  diff_exp_genes$padj[diff_exp_genes$padj==0]<-2.225074e-308
  
  response_QTLs_genes<-substr(names(x[which(x)]),1,15)
  shared_QTLs_genes<-substr(names(x[which(!x)]),1,15)
  
  #no_diff_exp_genes<-diff_exp_genes[diff_exp_genes$padj >0.5, ] # genes that have high p-values of not being diff exp for plotting 
  no_diff_exp_genes<-diff_exp_genes[ abs(diff_exp_genes$log2FoldChange) <= 1,] # include all to plot CTSA 
  diff_exp_genes <-diff_exp_genes[(abs(diff_exp_genes$log2FoldChange) > 1 & diff_exp_genes$padj <0.05),]
  
  response_diff_genes[[i]]<-table(response_QTLs_genes %in% diff_exp_genes$Gene)
  shared_diff_QTLs_genes[[i]]<-table(shared_QTLs_genes %in% diff_exp_genes$Gene)
  up_down_regulation[[i]]<-table(diff_exp_genes$log2FoldChange[diff_exp_genes$Gene %in% response_QTLs_genes] >0)
  response_diff_genes_stats[[i]]<-cbind(diff_exp_genes[diff_exp_genes$Gene %in% response_QTLs_genes,],reQTL=names(x[which(x)][which(response_QTLs_genes %in% diff_exp_genes$Gene)]))
  response_no_diff_genes_stats[[i]]<-cbind(no_diff_exp_genes[no_diff_exp_genes$Gene %in% response_QTLs_genes,],reQTL=names(x[which(x)][which(response_QTLs_genes %in% no_diff_exp_genes$Gene)]))
  
  
}

names(response_diff_genes)<-names(response_QTLs)
names(shared_diff_QTLs_genes)<-names(response_QTLs)
names(response_diff_genes_stats)<-names(response_QTLs)
names(response_no_diff_genes_stats)<-names(response_QTLs)

response_diff_genes_stats<-lapply(response_diff_genes_stats, function (x) x[order(x$log2FoldChange,x$padj),])
####################################################################################################################################################################################
## step 2 add lfsr 

path<-"/nfs/users/nfs_n/np12/myscratch/MacroMap/Data/MacroMap_manuscript_figures/data/"
sign_eQTLs_betas<-readRDS(paste0(path,"sign_eQTLs_betas.RDS"))
m.ed.t_all_eQTLs<-readRDS(paste0(path,"m.ed.t_all_eQTLs.RDS"))
rownames(m.ed.t_all_eQTLs$result$lfsr) <-rownames(sign_eQTLs_betas)

colnames(m.ed.t_all_eQTLs$result$lfsr)<-names(response_QTLs)
lfsr<-data.frame(m.ed.t_all_eQTLs$result$lfsr)
rownames(lfsr)<-read.table(paste0(path,"rownames_sign_eQTLs_betas.txt"),h=F)$V1

for (i  in 1:length(names(response_QTLs))) {
  reQTLs<-lfsr[rownames(lfsr) %in% response_diff_genes_stats[[i]]$reQTL,]
  reQTLs<-reQTLs[match(response_diff_genes_stats[[i]]$reQTL,rownames(reQTLs)),]
  response_diff_genes_stats[[i]]$lfsr<-reQTLs[,i]
  response_diff_genes_stats[[i]]<-response_diff_genes_stats[[i]][order(response_diff_genes_stats[[i]]$log2FoldChange,response_diff_genes_stats[[i]]$lfsr,decreasing = T),]
  
  reQTLs<-lfsr[rownames(lfsr) %in% response_no_diff_genes_stats[[i]]$reQTL,]
  reQTLs<-reQTLs[match(response_no_diff_genes_stats[[i]]$reQTL,rownames(reQTLs)),]
  response_no_diff_genes_stats[[i]]$lfsr<-reQTLs[,i]
  response_no_diff_genes_stats[[i]]<-response_no_diff_genes_stats[[i]][order(response_no_diff_genes_stats[[i]]$log2FoldChange,response_no_diff_genes_stats[[i]]$lfsr,decreasing = T),]
  
  
}

## step 3a. compare nominal betas between stimulation and ctrl to define gain/lost effects in reQTLs  
sign_eQTLs_betas<-readRDS(paste0(path,"sign_eQTLs_betas.RDS"))
rownames(sign_eQTLs_betas)<-read.table(paste0(path,"rownames_sign_eQTLs_betas.txt"),h=F)$V1
sign_eQTLs_std.err<-readRDS(paste0(path,"sign_eQTLs_std.err.RDS"))
rownames(sign_eQTLs_std.err)<-read.table(paste0(path,"rownames_sign_eQTLs_betas.txt"),h=F)$V1


sign_eQTLs_betas_comparison<-NULL
for (i in 1:length(response_QTLs)) {
  x<-response_diff_genes_stats[[i]]
  response_QTLs_per_cond<-x$reQTL
  
  if(i %in% c(1,2)) {
    sign_eQTLs_betas_comparison[[i]]<-data.frame(reQTL=response_QTLs_per_cond,data.frame(sign_eQTLs_betas[rownames(sign_eQTLs_betas) %in% response_QTLs_per_cond])[,c(i,3,4)],data.frame(sign_eQTLs_std.err[rownames(sign_eQTLs_std.err) %in% response_QTLs_per_cond])[,c(i,3)])
  } else {
    sign_eQTLs_betas_comparison[[i]]<-data.frame(reQTL=response_QTLs_per_cond,data.frame(sign_eQTLs_betas[rownames(sign_eQTLs_betas) %in% response_QTLs_per_cond])[,c(i+2,3,4)],data.frame(sign_eQTLs_std.err[rownames(sign_eQTLs_std.err) %in% response_QTLs_per_cond])[,c(i+2,3)])
  }
}

sign_eQTLs_betas_comparison_no_diff<-NULL
for (i in 1:length(response_QTLs)) {
  x<-response_no_diff_genes_stats[[i]]
  response_QTLs_per_cond<-x$reQTL
  
  if(i %in% c(1,2)) {
    sign_eQTLs_betas_comparison_no_diff[[i]]<-data.frame(reQTL=response_QTLs_per_cond,data.frame(sign_eQTLs_betas[rownames(sign_eQTLs_betas) %in% response_QTLs_per_cond])[,c(i,3,4)],data.frame(sign_eQTLs_std.err[rownames(sign_eQTLs_std.err) %in% response_QTLs_per_cond])[,c(i,3)])
  } else {
    sign_eQTLs_betas_comparison_no_diff[[i]]<-data.frame(reQTL=response_QTLs_per_cond,data.frame(sign_eQTLs_betas[rownames(sign_eQTLs_betas) %in% response_QTLs_per_cond])[,c(i+2,3,4)],data.frame(sign_eQTLs_std.err[rownames(sign_eQTLs_std.err) %in% response_QTLs_per_cond])[,c(i+2,3)])
  }
}
# combine them 
response_diff_genes_stats_betas<-NULL 
for (i in 1:length(response_QTLs)) {
  x<-response_diff_genes_stats[[i]]
  y<-sign_eQTLs_betas_comparison[[i]]
  response_diff_genes_stats_betas[[i]]<-x %>% left_join(y,by="reQTL")
}

z<-NULL
for (i in 1:length(response_QTLs)) {
  x<-response_diff_genes_stats_betas[[i]]
  
  a<-x %>% filter(log2FoldChange<0 & abs(x[11]) > abs(x[12])) %>% n_distinct("Gene")
  b<-x %>% filter(log2FoldChange<0 & abs(x[11]) < abs(x[12])) %>% n_distinct("Gene")
  c<-x %>% filter(log2FoldChange>0 & abs(x[11]) > abs(x[12])) %>% n_distinct("Gene")
  d<-x %>% filter(log2FoldChange>0 & abs(x[11]) < abs(x[12])) %>% n_distinct("Gene")
  z[[i]]<-data.frame(down_gain=a,down_lost=b,up_gain=c,up_lost=d)
}
stats_up_down<-t(do.call(rbind,(z)))
stats_up_down<-t(t(stats_up_down)/apply(stats_up_down,2, sum))

rowMeans(stats_up_down*100)
#down_gain down_lost   up_gain   up_lost 
#21.077277  3.138264 67.668374  8.116084 

# no diff exp reQTLs 
mean(sapply(sign_eQTLs_betas_comparison_no_diff,dim)[1,]/(sapply(sign_eQTLs_betas_comparison_no_diff,dim)[1,] + sapply(sign_eQTLs_betas_comparison,dim)[1,] ))

# results per condition just betas 
#lapply(sign_eQTLs_betas_comparison, function (x) table(factor(as.integer(abs(x[2]) > abs(x[3])),levels=c(0,1)))))

# SD into account corrects for power issues (e.g number of individuals in PIC_6)
#gain_lost_response<-t(do.call(rbind,lapply(sign_eQTLs_betas_comparison, function (x) table(factor(as.integer(abs(x[2]) > abs(x[3])- 1.96*x[5]),levels=c(0,1))))))

# check betas 
gain_lost_response_diff<-t(do.call(rbind,lapply(sign_eQTLs_betas_comparison, function (x) table(factor(as.integer(abs(x[2]) > abs(x[3])),levels=c(0,1)))))) 
gain_lost_response_nodiff<-t(do.call(rbind,lapply(sign_eQTLs_betas_comparison_no_diff, function (x) table(factor(as.integer(abs(x[2]) > abs(x[3])),levels=c(0,1)))))) 


range(100*(gain_lost_response_diff[1,]/ apply(gain_lost_response_diff,2,sum))) #  0.00000 14.28571 
mean(100*(gain_lost_response_diff[1,]/ apply(gain_lost_response_diff,2,sum))) # 11.25435 


####################################################################################################################################################################################
# 3b. compare betas stimulation vs ctrl_24 

eQTLs_betas_comparison_all<-NULL
for (i in 1:length(response_QTLs)) {
  eqtls_per_condition<-names(response_QTLs[[i]])
  x<-sign_eQTLs_betas[rownames(sign_eQTLs_betas) %in% eqtls_per_condition]
  
  if(i %in% c(1,2)) {
    eQTLs_betas_comparison_all[[i]]<- substr((eqtls_per_condition[abs(x[,1]) > abs(x[,3])]),1,15)
  } else {
    eQTLs_betas_comparison_all[[i]]<- substr((eqtls_per_condition[abs(x[,1]) > abs(x[,3])]),1,15)
  }
}

eQTLs_betas_comparison_reQTLsmashr<-NULL
for (i in 1:length(response_QTLs)) {
  x<-sign_eQTLs_betas_comparison[[i]]
  y<-sign_eQTLs_betas_comparison_no_diff[[i]]
  z<-rbind(x,y)
  eQTLs_betas_comparison_reQTLsmashr[[i]]<-substr(z$reQTL[abs(z[,2]) >abs(z[,3])],1,15)
}

length(unique(unlist(eQTLs_betas_comparison_reQTLsmashr))) # 2228 

## step 3 find examples and run them through the plotting eQTL script to extract SNPs/expression values 


# up regulated genes that has response QTL  ENSG00000136688.10_rs6724667_chr2_113032064 IL36G 
# down regulated genes that has a respone QTL ENSG00000105875.13_rs62481922_chr7_135206797 WDR91 
# reQTL gene no differential expression   ENSG00000257704.3_rs7248282_chr19_47265863  INAFM1 
# YPEL4 response in Ctrl                  ENSG00000166793.10_rs4939158_chr11_57790552  YPEL4


# manual run the above genes 

source("/nfs/users/nfs_n/np12/myscratch/MacroMap/Data/MacroMap_manuscript_figures/scripts/gencode_annotation.R")
source("/nfs/users/nfs_n/np12/myscratch/MacroMap/Data/MacroMap_manuscript_figures/scripts/ggplot_themes.R")

log10TPMs<-readRDS("/nfs/users/nfs_n/np12/myscratch/MacroMap/Analysis/QC/RDS_expression_covariates/Macromap_final_expression/protein_coding_lincRNA/RDS_log10_TPMs_filtd_macromap_fds.rds")
covariates<-readRDS("/nfs/users/nfs_n/np12/myscratch/MacroMap//Analysis/QC/RDS_expression_covariates/Macromap_final_expression/protein_coding_lincRNA/covariates_macromap_fds.rds")

log10TPMs_by_condition<-split(data.frame(t(log10TPMs[,-c(1:6)])),covariates$Stimulus_Hours)
log10TPMs_by_condition<-lapply(log10TPMs_by_condition,function (x) t(x))

col<-brewer.pal(8, "Dark2")


write.table('ENSG00000166793.10_rs4939158_chr11_57790552',file="/nfs/users/nfs_n/np12/myscratch/MacroMap/Data/MacroMap_manuscript_figures/scripts/extract_info_QTLtools/eQTL.txt",row.names=F,col.names=F,quote=F)

setwd("/nfs/users/nfs_n/np12/myscratch/MacroMap/Data/MacroMap_manuscript_figures/scripts/extract_info_QTLtools/")
system ("bash extract_phenotype_corrected.sh < eQTL.txt|bash")


gene<-ref("ENSG00000166793.10")
group<-"sLPS_6"
cond<-sapply(strsplit(group,"_"),"[[",1)
hours<-sapply(strsplit(group,"_"),"[[",2)

ens_id<-ens(gene)


gr1vsgr2<-readRDS(paste0("/nfs/users/nfs_n/np12/myscratch/MacroMap/Data/MacroMap_manuscript_figures/data/diff_exp/All_Ctrl_vs_",cond,".RDS"))


colors_all = colorRampPalette(brewer.pal(8, "Dark2"))(96) #  
colors_all_24<-as.character(colors_all[c(1,4,7,8,17,20,25,28,33,36,41,44,49,50,57,60,64,66,72,77,84,88,93,96)])
colorCond<-cbind.data.frame(conditions,colors_all_24,stringsAsFactors = FALSE)

gene_idx<-which(rownames(gr1vsgr2) == substr(ens_id,1,15))

fiss <- plotCounts(gr1vsgr2,gene_idx, intgroup = c("Stimulus","Hours"), returnData = TRUE)
fiss$log10_TPMs<-as.numeric(log10TPMs[ens_id,rownames(fiss)])
fiss$color<-colorCond$colors_all_24[match(paste0(fiss$Stimulus, "_", fiss$Hours),colorCond$conditions)]
fiss$Stimulus_Hours<-paste0(fiss$Stimulus, "_", fiss$Hours)
fiss$Expression<-residuals(lm(fiss$log10_TPMs ~ covariates$Library_prep[covariates$SampleID_RunID %in% rownames(fiss)]))
snp<-read.table(paste0("/nfs/users/nfs_n/np12/myscratch/MacroMap/Data/MacroMap_manuscript_figures/scripts/extract_info_QTLtools/snp_",group,".out.content.txt.gz"),h=T)

fiss$HipsciID<-covariates$HipsciID[which(covariates$SampleID_RunID %in% rownames(fiss) )]
fiss$Variant<-round(snp[,2][match(fiss$Hi,snp$sample)])

fiss<-fiss %>% filter(Stimulus_Hours %in% c(group,paste0("Ctrl_",hours))) %>% filter(!is.na(Variant)) 
p<- ggplot(fiss,aes(x=factor(round(Variant)),y=Expression,color=Stimulus_Hours)) + geom_boxplot(width=0.5,outlier.size = -1) +  
  facet_wrap(~Stimulus_Hours) + geom_point(size = 3,alpha=1, shape = 20,position = position_jitterdodge(jitter.width = 0.4)) + 
  geom_smooth(aes(group=Stimulus_Hours), method="lm",se=T) +
  scale_fill_manual(values = colorCond$colors_all_24[colorCond$conditions %in% unique(fiss$Stimulus_Hours)]) + 
  scale_color_manual(values = colorCond$colors_all_24[colorCond$conditions %in% unique(fiss$Stimulus_Hours)]) + xlab(names(snp)[2]) + ylab(paste0("Expression","\n",gene))

pdf(paste0(path_to_plot,"Fig2c.d.e_",gene,"_",group,"_boxplots.pdf"),useDingbats=FALSE)
print (p + theme_Publication(24) + theme(strip.background = element_blank(), strip.text = element_blank())) 
dev.off()
