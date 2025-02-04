library(data.table)
library(ggplot2)
library(RColorBrewer)
library("ggrepel")
source(paste0("~/bin/Rscripts/","ggplot_themes.R"))

########################################################################################################################################################################
# load GTEx data 
temp = list.files(path = "/nfs/users/nfs_n/np12/myscratch/MacroMap/Data/GTEx_v8/permuted/",pattern="*.egenes.txt.gz",full.names = T)
myfiles = lapply(temp, fread)  

names(myfiles)<-sapply(strsplit(sapply(strsplit(temp,"[.]"),"[[",1),"permuted//"),"[[",2)

anno<-read.table("/nfs/users/nfs_n/np12/myscratch/MacroMap/Data/GTEx_v8/permuted/gencode/protein_cod_link_RNA.txt")

sign_per_tissue<- function (x) { 
  
  res<-table(x$qval[x$gene_id %in% anno$V5]<=0.05)
  res1<-cbind(res[2],res[2]/sum(res))
  return(res1)
}
########################################################################################################################################################################
# number of donors per condition 
donors_1<-readRDS("/nfs/users/nfs_n/np12/myscratch/MacroMap/Data/MacroMap_manuscript_figures/data/samples_per_condition.RDS")

gtex_new<-read.csv("/nfs/users/nfs_n/np12/myscratch/MacroMap/Data/MacroMap_manuscript_figures/data/gtex_stats.csv")
gtex_new<-gtex_new[gtex_new$Donors_geno_RNA >70,]
gtex_new<-gtex_new[!is.na(gtex_new$Tissue),]
gtex_new<-cbind(gtex_new,data.frame(tissue_id=names(myfiles),do.call(rbind,lapply(myfiles,sign_per_tissue))))

colors<-read.table("/nfs/users/nfs_n/np12/myscratch/MacroMap/Data/MacroMap_manuscript_figures/data/gtex_colors.txt",sep="\t",comment.char="&",h=T)
colors<-colors[match(gtex_new$tissue_id, colors$tissue_id),]
gtex_new<-gtex_new[match(gtex_new$tissue_id, colors$tissue_id),]
gtex_new$Color<-tolower(colors$color_hex)
gtex_new$eGenes<-as.numeric(as.character(gtex_new$eGenes))


pdf("/nfs/users/nfs_n/np12/myscratch/MacroMap/Data/MacroMap_manuscript_figures/plots/Fig2_support.gtex_plot_genes.pdf")
par(mar = c(5,6,2,2))
plot(gtex_new$Donors_geno_RNA,gtex_new$X1,col=gtex_new$Color,pch=20,cex=3,xlim=c(50,800),ylim=c(0,20000),ylab="Number of eGenes (5% FDR)",cex.lab=1.7,cex.axis=1.5,xlab="Sample size",font.lab=2 )
lm.gtex<-lm (X1 ~ Donors_geno_RNA,data=gtex_new)
abline(lm.gtex,lwd=4,col="grey50",lty=2)
points(mean(donors_1$V2),10170,col="red",pch=15,cex=2.5) # macromap 
dev.off()

pdf("/nfs/users/nfs_n/np12/myscratch/MacroMap/Data/MacroMap_manuscript_figures/plots/Fig2_support.gtex_plot_genes_precentages.pdf")
par(mar = c(5,7,2,2))
plot(gtex_new$Donors_geno_RNA,gtex_new$X2,col=gtex_new$Color,pch=20,cex=3,xlim=c(50,800),ylim=c(0,1),ylab="Number of eGenes (5% FDR)/\n Total expresssed",cex.lab=1.7,cex.axis=1.5,xlab="Sample size",font.lab=2 )
lm.gtex<-lm (X2 ~ Donors_geno_RNA,data=gtex_new)
abline(lm.gtex,lwd=4,col="grey50",lty=2)
points(mean(donors_1$V2),10170/14034,col="red",pch=18,cex=2.5) # macromap 
dev.off()


gtex_new$tissue_id<-as.character(gtex_new$tissue_id)
gtex_new$Tissue_id<-ifelse(gtex_new$tissue_id %in% gtex_new$tissue_id[c(12,21,22,31)], gtex_new$tissue_id,"")

d<-data.frame("","",mean(donors_1$V2),"","",10170/14034,"","","MacroMap",mean(donors_1$V2),10170/14034,"red","MacroMap")
colnames(d)<-colnames(gtex_new)
gtex_new<-rbind(gtex_new,d)

gtex_new$C<-"grey70"
gtex_new$C[which(!gtex_new$Tissue_id=="")]<-gtex_new$Color[which(!gtex_new$Tissue_id=="")]

p<-ggplot(gtex_new,aes(x = Donors_geno_RNA, y = X2))  + geom_point(size=4,color=gtex_new$C) + geom_smooth(method='lm',linetype="dashed",color="grey50",fill = "") + 
  geom_text_repel(aes(label= Tissue_id),color = factor(gtex_new$Color),bg.r=0,show.legend = FALSE,box.padding = 1,point.padding = 0.8,size=7,max.overlaps  = Inf) + xlim(c(0,750)) + ylim(0,1)  + 
  geom_point(aes(x=mean(donors_1$V2), y=10170/14034), colour="red",pch=18,cex=7) + xlab("Sample size") + ylab("Number of eGenes /\n Total expresssed") + 
  scale_x_continuous(breaks=seq(0, 700, by = 100))


pdf("/nfs/users/nfs_n/np12/myscratch/MacroMap/Data/MacroMap_manuscript_figures/plots/Fig2a.gtex_plot_genes_precentages_ggplot_version_paper.pdf",useDingbats=FALSE)
p + theme_Publication(20) + theme(strip.background = element_blank(), strip.text = element_blank())
dev.off()
########################################################################################################################################################################
# paper comparison 
gtex_similar<-gtex_new[gtex_new$Donors_geno_RNA >180 & gtex_new$Donors_geno_RNA  < 210,]
gtex_similar<-gtex_similar[gtex_similar$tissue_id!="MacroMap",]
macromap_eqtls<-read.table("/nfs/users/nfs_n/np12/myscratch/MacroMap/Analysis/eQTLs/Macromap_fds/analysis/eQTLs_per_TSS/1MB/cond.eqtls.txt")

macromap_eqtls$V4<-macromap_eqtls$V3/14034


# 14034 This is the correct number of tested genes. 14060 are with chrM included which I don't have genotypes 
########################################################################################################################################################################
