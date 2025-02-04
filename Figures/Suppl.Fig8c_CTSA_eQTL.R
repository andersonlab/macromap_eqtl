library(data.table)
library(RColorBrewer)
library(ggplot2)
library(DESeq2)
library(dplyr)


source("/nfs/users/nfs_n/np12/myscratch/MacroMap/Data/MacroMap_manuscript_figures/scripts/gencode_annotation.R")
source("/nfs/users/nfs_n/np12/myscratch/MacroMap/Data/MacroMap_manuscript_figures/scripts/ggplot_themes.R")

path_to_plot<-"/nfs/users/nfs_n/np12/myscratch/MacroMap/Data/MacroMap_manuscript_figures/plots/"

log10TPMs<-readRDS("/nfs/users/nfs_n/np12/myscratch/MacroMap/Analysis/QC/RDS_expression_covariates/Macromap_final_expression/protein_coding_lincRNA/RDS_log10_TPMs_filtd_macromap_fds.rds")
covariates<-readRDS("/nfs/users/nfs_n/np12/myscratch/MacroMap//Analysis/QC/RDS_expression_covariates/Macromap_final_expression/protein_coding_lincRNA/covariates_macromap_fds.rds")

log10TPMs_by_condition<-split(data.frame(t(log10TPMs[,-c(1:6)])),covariates$Stimulus_Hours)
log10TPMs_by_condition<-lapply(log10TPMs_by_condition,function (x) t(x))

col<-brewer.pal(8, "Dark2")
conditions<-c("CIL_24","CIL_6","Ctrl_24","Ctrl_6","IFNB_24","IFNB_6","IFNG_24","IFNG_6","IL4_24","IL4_6","LIL10_24","LIL10_6","P3C_24","P3C_6","R848_24",
              "R848_6","sLPS_24","sLPS_6","PIC_24","PIC_6","MBP_24","MBP_6","Prec_D0","Prec_D2")

## CTSA CAD 

write.table('ENSG00000064601.16_rs3827066_chr20_45957384',file="/nfs/users/nfs_n/np12/myscratch/MacroMap/Data/MacroMap_manuscript_figures/scripts/extract_info_QTLtools/eQTL.txt",row.names=F,col.names=F,quote=F)

setwd("/nfs/users/nfs_n/np12/myscratch/MacroMap/Data/MacroMap_manuscript_figures/scripts/extract_info_QTLtools/")
system ("bash extract_phenotype_corrected.sh < eQTL.txt|bash")


gene<-ref("ENSG00000064601.16")
group<-"LIL10_24"
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

pdf(paste0(path_to_plot,"Suppl_Fig8c_",gene,"_",group,"_boxplots.pdf"),useDingbats=FALSE)
print (p + theme_Publication(24) + theme(strip.background = element_blank(), strip.text = element_blank())) 
dev.off()
