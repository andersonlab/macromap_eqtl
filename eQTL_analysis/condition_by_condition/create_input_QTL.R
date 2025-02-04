# Create fastQTL and QTLtools input expresion files

args <- commandArgs(trailingOnly = TRUE)

gene_expression<-read.table(args[1],header=T,check.names=FALSE)


header<-gene_expression[c(1:6)]
header$tss<-ifelse(header$Strand=="-",header$End-1,header$Start-1)

header.b<-NULL

for (i in 1: dim (header)[1])
{
  if(header[i,]$Strand=="+") {
    header.b[[i]]<-header[i,][c(1,7,3,2,6,5)]
    colnames(header.b[[i]])<-c("Chr","gene_start","gene_end","gene","length","strand")
  }
  else{
    header.b[[i]]<-header[i,][c(1,7,4,2,6,5)]
    colnames(header.b[[i]])<-c("Chr","gene_start","gene_end","gene","length","strand")
  }
}

header.c<-do.call(rbind.data.frame, header.b)

gene_expression_qtltools<-data.frame(header.c,gene_expression[-c(1:6)],check.names=FALSE)
gene_expression_fastqtl<-data.frame(header.c[1:4],gene_expression[-c(1:6)],check.names=FALSE)

# order by chr and position
chrOrder<-c(paste(1:22,sep=""),"X","Y","M")
gene_expression_qtltools$Chr<-factor(gene_expression_qtltools$Chr, levels=chrOrder)
gene_expression_fastqtl$Chr<-factor(gene_expression_fastqtl$Chr, levels=chrOrder)

# order by position
gene_expression_qtltools<-gene_expression_qtltools[order(gene_expression_qtltools$Chr,gene_expression_qtltools$gene_start),]
gene_expression_fastqtl<-gene_expression_fastqtl[order(gene_expression_fastqtl$Chr,gene_expression_fastqtl$gene_start),]

# Just the autosomes
if (args[2]=="autosomes") {
gene_expression_qtltools<-gene_expression_qtltools[as.numeric(gene_expression_qtltools$Chr) <= 22 ,]
gene_expression_fastqtl<-gene_expression_fastqtl[as.numeric(gene_expression_fastqtl$Chr) <= 22 ,]
}

if (args[2]=="autosomesX") {
gene_expression_qtltools<-gene_expression_qtltools[gene_expression_qtltools$Chr %in% c(1:22,"X") ,]
gene_expression_fastqtl<-gene_expression_fastqtl[gene_expression_fastqtl$Chr %in% c(1:22,"X") ,]
}


gene_expression_qtltools$Chr<-paste0("chr",gene_expression_qtltools$Chr)
gene_expression_fastqtl$Chr<-paste0("chr",gene_expression_fastqtl$Chr)

colnames(gene_expression_qtltools)[1]<-"#chr"
colnames(gene_expression_fastqtl)[1]<-"#chr"

if (args[3]=="fastQTL") {
  file_name<-paste0(strsplit(args[1],".txt"),"_fastQTL_input.txt")
  write.table(gene_expression_fastqtl,file=paste0(strsplit(args[1],".txt"),"_fastQTL_input.txt"),col.names=T, row.names=F,sep="\t",quote=F)
  } else {
    file_name<-paste0(strsplit(args[1],".txt"),"_QTLtools_input.txt")
    write.table(gene_expression_qtltools,file=paste0(strsplit(args[1],".txt"),"_QTLtools_input.txt"),col.names=T, row.names=F,sep="\t",quote=F)
}

system (paste("bgzip",file_name))
system(paste("tabix -p bed",paste0(file_name,".gz")))
#tabix -p bed $i\.fastQTL_format_final.txt.gz
#rm $i*tmp $i*header
