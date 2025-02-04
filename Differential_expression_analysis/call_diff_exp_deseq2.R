library("DESeq2")
library("biomaRt")
library("pheatmap")
library("RColorBrewer")
library("ggplot2")
library("limma")
library("BiocParallel")
library("edgeR")
library("RColorBrewer")
source("~/bin/Rscripts/lm.pval.R")
library("qvalue")
library("sva")
library("viridis")

volcano <- function (res,group1,group2,hours) 
{
  maintext<-ifelse(hours %in% c("6","24"),paste0(group1,"_vs_",group2,"\n",hours," hours"),paste0(group1,"_vs_",group2,"\n",hours))
  pdf(file=paste0("Volcano_",group1,"_vs_",group2,"_",hours,".pdf"))
  with(as.data.frame(res), plot(log2FoldChange, -log10(pvalue), pch=20, main=maintext, 
                                xlim=c(-max(log2FoldChange),max(log2FoldChange))),ylim=c(0,max(-log10(pvalue))),cex.axis=1.5,cex.main=1.5)
  # Add colored points: red if padj<0.05, orange of log2FC>1, green if both)
  with(subset(as.data.frame(res), padj<.05 ), points(log2FoldChange, -log10(pvalue), pch=20, col="red"))
  with(subset(as.data.frame(res), abs(log2FoldChange)>1), points(log2FoldChange, -log10(pvalue), pch=20, col="orange"))
  with(subset(as.data.frame(res), padj<.05 & abs(log2FoldChange)>1), points(log2FoldChange, -log10(pvalue), pch=20, col="blue"))
  legend("topleft",legend=c(paste0("5% FDR (",length(subset(as.data.frame(res),padj<.05 )$Gene),")"),
                            paste0("log2FC >1 not Sign (",length(subset(as.data.frame(res), abs(log2FoldChange)>1)$Gene),")"),
                            paste0("5%FDR & log2FC>1 (",length(subset(as.data.frame(res),padj<.05 & abs(log2FoldChange)>1)$Gene),")")),
         pch=c(20,20),col=c("red","orange","blue"),cex=1) #inset=c(-0.5,0) )
  dev.off()
}


register(MulticoreParam(5))
################################################################################################################################################################################################
# load data 
dds_expression<-readRDS("~/myscratch/MacroMap/Analysis/QC/RDS_expression_covariates/Macromap_final_expression/protein_coding_lincRNA/dds_expression_macromap_fds.rds")
covariates<-readRDS("~/myscratch/MacroMap/Analysis/QC/RDS_expression_covariates/Macromap_final_expression/protein_coding_lincRNA/covariates_macromap_fds.rds")

dds_expression$Stimulus_Hours<-as.factor(dds_expression$Stimulus_Hours)
dds_expression$Donor<-as.factor(dds_expression$Donor)


callDE<-function (dds_obj,group1,group2) 
{
  gr1vsgr2<-paste0("dds_comp_",group1,"_vs_",group2)
  #gr1vsgr2<-dds_obj[,dds_obj$Stimulus %in% c(group1,group2)]
  gr1vsgr2<-dds_obj[,grep(paste0(group1,"|",group2),colnames(dds_obj))]
  gr1vsgr2 <- estimateSizeFactors(gr1vsgr2)

  
  gr1vsgr2$Stimulus<-dropEmptyLevels(gr1vsgr2$Stimulus)
  gr1vsgr2$Stimulus_Hours<-dropEmptyLevels(gr1vsgr2$Stimulus_Hours)

  if (group2 %in% c("Prec_D0","Prec_D2")) {
    gr1vsgr2$Stimulus<-relevel(gr1vsgr2$Stimulus,"Ctrl") 
    } else {
    gr1vsgr2$Stimulus<-relevel(gr1vsgr2$Stimulus,ref=group1)
    gr1vsgr2$Hours<-relevel(as.factor(gr1vsgr2$Hours),"6","24")
  }
  

dat  <- counts(gr1vsgr2, normalized = TRUE)
idx  <- rowMeans(dat) > 1
dat  <- dat[idx, ]
mod  <- model.matrix(~ Stimulus_Hours , colData(gr1vsgr2))
mod0 <- model.matrix(~   1, colData(gr1vsgr2))
svseq <- svaseq(dat, mod, mod0, n.sv = 10)
colData(gr1vsgr2)<-cbind(colData(gr1vsgr2),data.frame(svseq$sv[,c(1:10)]))

if (group2 %in% c("Prec_D0","Prec_D2")) {
  gr1vsgr2<-DESeqDataSet(gr1vsgr2,design= ~ X1 + X2 + X3 + X4 + X5 + X6 + X7 + X8 + X9 + X10  +  Stimulus_Hours)
  gr1vsgr2<-DESeq(gr1vsgr2)
} else {
  gr1vsgr2<-DESeqDataSet(gr1vsgr2,design= ~ X1 + X2 + X3 + X4 + X5 + X6 + X7 + X8 + X9 + X10  +  Stimulus + Hours + Stimulus:Hours)
  gr1vsgr2<-DESeq(gr1vsgr2,test="LRT", reduced = ~ X1 + X2 + X3 + X4 + X5 + X6 + X7 + X8 + X9 + X10 +  Stimulus + Hours, parallel=TRUE)
}


genes<-readRDS("/lustre/scratch117/cellgen/team170/np12/MacroMap/Analysis/QC/RDS_expression_covariates/Macromap_final_expression/protein_coding_lincRNA/gene_expression.merged_macromap_fds.rds")
genes<-genes[match(rownames(gr1vsgr2),substr(rownames(genes),1,15)),]

saveRDS(gr1vsgr2,file=paste0("All_",group1,"_vs_",group2,".RDS"))

#resultsNames(gr1vsgr2)
if (group2 %in% c("Prec_D0","Prec_D2")) {
  res_6 <- results(gr1vsgr2)
  res_6$hgnc<-genes$hgnc_symbol
  res_6<-res_6[!is.na(res_6$padj),]
  res_6$Gene<-rownames(res_6)
  col_idx <- grep("Gene|hgnc", names(res_6))
  res_6 <- res_6[, c(col_idx, (1:ncol(res_6))[-col_idx])]
  write.table( res_6, file=paste0("All_",group1,"_vs_",group2,"_","6",".txt"),sep="\t",row.names=F,quote=F)
  
  res_6.sign<-res_6[which(res_6$padj <=0.05 & abs(res_6$log2FoldChange) >=1 ),]
  res_6.sign <-res_6.sign[ order(res_6.sign$pvalue), ] 
  write.table( res_6.sign, file=paste0("Diff_expressed_genes_at_5%_FDR_",group1,"_vs_",group2,"_","6",".txt"),sep="\t",row.names=F,quote=F)
  write.table( res_6.sign$Gene, file=paste0("Diff_expressed_genes_at_5%_FDR_",group1,"_vs_",group2,"_","6","_pathways.txt"),sep="\t",row.names=F,col.names=F,quote=F)
  volcano(res_6,group1,group2,"6") 
  
  } else {

# 6 hours IFNG vs Ctrl 
res_6 <- results(gr1vsgr2, name=paste0("Stimulus_",group2,"_vs_",group1), test="Wald", alpha=0.05)
res_6$hgnc<-genes$hgnc_symbol
res_6<-res_6[!is.na(res_6$padj),]
res_6$Gene<-rownames(res_6)
col_idx <- grep("Gene|hgnc", names(res_6))
res_6 <- res_6[, c(col_idx, (1:ncol(res_6))[-col_idx])]
write.table( res_6, file=paste0("All_",group1,"_vs_",group2,"_","6",".txt"),sep="\t",row.names=F,quote=F)

res_6.sign<-res_6[which(res_6$padj <=0.05 & abs(res_6$log2FoldChange) >=1 ),]
res_6.sign <-res_6.sign[ order(res_6.sign$pvalue), ] 
write.table( res_6.sign, file=paste0("Diff_expressed_genes_at_5%_FDR_",group1,"_vs_",group2,"_","6",".txt"),sep="\t",row.names=F,quote=F)
write.table( res_6.sign$Gene, file=paste0("Diff_expressed_genes_at_5%_FDR_",group1,"_vs_",group2,"_","6","_pathways.txt"),sep="\t",row.names=F,col.names=F,quote=F)
volcano(res_6,group1,group2,"6") 


# 24 hours IFNG vs Ctrl 
res_24 <- results(gr1vsgr2, contrast=list(c(paste0("Stimulus_",group2,"_vs_",group1),paste0("Stimulus",group2,".Hours24"))), test="Wald", alpha=0.05,parallel=TRUE)
res_24$hgnc<-genes$hgnc_symbol
res_24<-res_24[!is.na(res_24$padj),]
res_24$Gene<-rownames(res_24)
col_idx <- grep("Gene|hgnc", names(res_24))
res_24 <- res_24[, c(col_idx, (1:ncol(res_24))[-col_idx])]
write.table( res_24, file=paste0("All_",group1,"_vs_",group2,"_","24",".txt"),sep="\t",row.names=F,quote=F)

res_24.sign<-res_24[which(res_24$padj <=0.05 & abs(res_24$log2FoldChange) >=1 ),]
res_24.sign <-res_24.sign[ order(res_24.sign$pvalue), ] 
write.table( res_24.sign, file=paste0("Diff_expressed_genes_at_5%_FDR_",group1,"_vs_",group2,"_","24",".txt"),sep="\t",row.names=F,quote=F)
write.table( res_24.sign$Gene, file=paste0("Diff_expressed_genes_at_5%_FDR_",group1,"_vs_",group2,"_","24","_pathways.txt"),sep="\t",row.names=F,col.names=F,quote=F)
volcano(res_24,group1,group2,"24") 

# interaction term 
res_int <- results(gr1vsgr2, name=paste0("Stimulus",group2,".Hours24"), test="Wald", alpha=0.05)
res_int$hgnc<-genes$hgnc_symbol
res_int<-res_int[!is.na(res_int$padj),]
res_int$Gene<-rownames(res_int)
col_idx <- grep("Gene|hgnc", names(res_int))
res_int <- res_int[, c(col_idx, (1:ncol(res_int))[-col_idx])]
write.table( res_int, file=paste0("All_",group1,"_vs_",group2,"_","interactionTime",".txt"),sep="\t",row.names=F,quote=F)

res_int.sign<-res_int[which(res_int$padj <=0.05 & abs(res_int$log2FoldChange) >=1 ),]
res_int.sign <-res_int.sign[ order(res_int.sign$pvalue), ] 
write.table( res_int.sign, file=paste0("Diff_expressed_genes_at_5%_FDR_",group1,"_vs_",group2,"_","interactionTime",".txt"),sep="\t",row.names=F,quote=F)
write.table( res_int.sign$Gene, file=paste0("Diff_expressed_genes_at_5%_FDR_",group1,"_vs_",group2,"_","interactionTime","_pathways.txt"),sep="\t",row.names=F,col.names=F,quote=F)
volcano(res_int,group1,group2,"interaction_time") 


# heatmap of log2fold 
all.deg<-unique(c(rownames(res_6.sign),rownames(res_24.sign),rownames(res_int.sign)))

all.betas<-data.frame(res_6[all.deg,]$log2FoldChange,res_24[all.deg,]$log2FoldChange,res_int[all.deg,]$log2FoldChange)
colnames(all.betas)<-c(paste0(group1,"_vs_",group2,"_","6"),paste0(group1,"_vs_",group2,"_","24"),paste0(group1,"_vs_",group2,"_","interaction"))

mat<-all.betas
thr <- 3 
mat[mat < -thr] <- -thr
mat[mat > thr] <- thr

pdf(file=paste0("heatmap_",group1,"_vs_",group2,".pdf"))
pheatmap(mat, breaks=seq(from=-thr, to=thr, length=100),cluster_col=FALSE,color=magma(100),show_rownames=F)
dev.off()


##################################################################################################################################################################################
# log10 TPMs 
log10TPMs<-readRDS("/lustre/scratch117/cellgen/team170/np12/MacroMap/Analysis/QC/RDS_expression_covariates/Macromap_final_expression/protein_coding_lincRNA/RDS_log10_TPMs_filtd_macromap_fds.rds")
log10TPMs<-log10TPMs[match(rownames(gr1vsgr2),substr(rownames(log10TPMs),1,15)),]

conditions<-c("CIL_24","CIL_6","Ctrl_24","Ctrl_6","IFNB_24","IFNB_6","IFNG_24","IFNG_6","IL4_24","IL4_6","LIL10_24","LIL10_6","P3C_24","P3C_6","R848_24",
              "R848_6","sLPS_24","sLPS_6","PIC_24","PIC_6","MBP_24","MBP_6","Prec_D0","Prec_D2")
conditions<-conditions[order(conditions)]

colors_all = colorRampPalette(brewer.pal(8, "Dark2"))(96) #  
colors_all_24<-as.character(colors_all[c(1,4,7,8,17,20,25,28,33,36,41,44,49,50,57,60,64,66,72,77,84,88,93,96)])
colorCond<-cbind.data.frame(conditions,colors_all_24,stringsAsFactors = FALSE)

only_6<-res_6.sign[(!res_6.sign$Gene %in% c(res_24$Gene[res_24$padj <=0.05],res_int$Gene[res_int$padj <=0.05])),]
only_24<-res_24.sign[(!res_24.sign$Gene %in% c(res_6$Gene[res_6$padj <=0.05],res_int$Gene[res_int$padj <=0.05])),]
only_int<-res_int.sign[(!res_int.sign$Gene %in% c(res_24$Gene[res_24$padj <=0.05],res_6$Gene[res_int$padj <=0.05])),]

write.table( only_6, file=paste0("Diff_expressed_genes_at_5%_FDR_",group1,"_vs_",group2,"_","only_6",".txt"),sep="\t",row.names=F,quote=F)
write.table( only_24, file=paste0("Diff_expressed_genes_at_5%_FDR_",group1,"_vs_",group2,"_","only_24",".txt"),sep="\t",row.names=F,quote=F)
write.table( only_int, file=paste0("Diff_expressed_genes_at_5%_FDR_",group1,"_vs_",group2,"_","only_int",".txt"),sep="\t",row.names=F,quote=F)

write.table( only_6$Gene, file=paste0("Diff_expressed_genes_at_5%_FDR_",group1,"_vs_",group2,"_","only_6","_pathways.txt"),sep="\t",row.names=F,col.names=F,quote=F)
write.table( only_24$Gene, file=paste0("Diff_expressed_genes_at_5%_FDR_",group1,"_vs_",group2,"_","only_24","_pathways.txt"),sep="\t",row.names=F,col.names=F,quote=F)
write.table( only_int$Gene, file=paste0("Diff_expressed_genes_at_5%_FDR_",group1,"_vs_",group2,"_","only_int","_pathways.txt"),sep="\t",row.names=F,col.names=F,quote=F)



plot_boxplots<-function (gene,name,subname,path,order) 
{
gene_idx<-which(gene  ==  rownames(gr1vsgr2))

#fiss <- plotCounts(gr1vsgr2, order(res_6.sign$padj)[2], intgroup = c("Stimulus","Hours"), returnData = TRUE)
fiss <- plotCounts(gr1vsgr2,gene_idx, intgroup = c("Stimulus","Hours"), returnData = TRUE)
fiss$log10_TPMs<-as.numeric(log10TPMs[gene_idx,rownames(fiss)])
fiss$color<-colorCond$colors_all_24[match(paste0(fiss$Stimulus, "_", fiss$Hours),colorCond$conditions)]
fiss$Stimulus_Hours<-paste0(fiss$Stimulus, "_", fiss$Hours)
#ggplot(fiss,aes(x = Hours, y = log10_TPMs, color = Stimulus_Hours, group = Stimulus_Hours)) + geom_point() + stat_summary(fun.y=mean, geom="line")


p<-ggplot(fiss,aes(x = Hours, y = log10_TPMs, color = Stimulus_Hours, group =Stimulus_Hours))  + 
  geom_boxplot(width=0.5,outlier.size = -1) + 
  geom_point(size = 4, shape = 20,alpha=0.75, position = position_jitterdodge(jitter.width = 0.1)) +
  scale_fill_manual(values = colorCond$colors_all_24[colorCond$conditions %in% unique(fiss$Stimulus_Hours)]) + 
  scale_color_manual(values = colorCond$colors_all_24[colorCond$conditions %in% unique(fiss$Stimulus_Hours)]) +
  theme_classic() +  stat_summary(fun.y=median, geom="line",aes(group=Stimulus),size=1.5) + ggtitle(name) +
  theme(axis.title.x=element_text(size=22),axis.title.y=element_text(size=22),axis.text.x=element_text(size=22),
        axis.text.y=element_text(size=22),legend.title = element_text(size=15),legend.text = element_text(size=14),plot.title = element_text(hjust = 0.5,size=24)) 

p1<-ggplot(fiss,aes(x = Hours, y = log10_TPMs, color = Stimulus_Hours, group =Stimulus_Hours))  + 
  geom_boxplot(width=0.5,outlier.size = -1) + 
  geom_point(size = 4, shape = 20,alpha=0.75, position = position_jitterdodge(jitter.width = 0.1)) +
  scale_fill_manual(values = colorCond$colors_all_24[colorCond$conditions %in% unique(fiss$Stimulus_Hours)]) + 
  scale_color_manual(values = colorCond$colors_all_24[colorCond$conditions %in% unique(fiss$Stimulus_Hours)]) +
  theme_classic() + ggtitle(name) +
  theme(axis.title.x=element_text(size=22),axis.title.y=element_text(size=22),axis.text.x=element_text(size=22),
        axis.text.y=element_text(size=22),legend.title = element_text(size=15),legend.text = element_text(size=14),plot.title = element_text(hjust = 0.5,size=24)) 


pdf(file=paste0(path,"/",subname,"_boxplot_",group1,"_vs_",group2,"_",name,"_medians_order",order,".pdf"))
print(p)
dev.off()

pdf(file=paste0(path,"/",subname,"_boxplot_",group1,"_vs_",group2,"_",name,"_",order,".pdf"))
print(p1)
dev.off()

}


toplot<-c("res_6.sign","res_24.sign","res_int.sign","only_6","only_24","only_int")

for (k in 1:length(toplot)) {
  dir.create(toplot[k], showWarnings = FALSE)
  tmp<-get(toplot[k])
  if(dim (tmp)[1]==0) next
  for(i in seq(1:20))
  {
    plot_boxplots(tmp$Gene[i],ifelse(is.na(tmp$hgnc[i]),tmp$Gene[i],tmp$hgnc[i]),toplot[k],toplot[k],i)
  }
}
}
}


args = commandArgs(trailingOnly=TRUE) # input run args[1]==always controls a

dir_analysis<-paste0(getwd(),"/",args[1],"_vs_",args[2])
dir.create(file.path(dir_analysis), showWarnings = FALSE)
setwd(file.path(dir_analysis))

callDE(dds_expression,args[1],args[2])


########################################################################################################################################################################
