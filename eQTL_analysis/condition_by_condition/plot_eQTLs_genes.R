library(qvalue)
library(RColorBrewer)
library(data.table)
library(plyr)

COL = brewer.pal(9,"Set1");
args = commandArgs(trailingOnly=TRUE) # input run


tss<-c("50K","250K","500K","750K","1MB")
#tss<-c("500K")

PC=c('NO_PC',1,2,3,4,5,6,7,8,9,10,15,20,25,30,35,40,50,60,70,80)
axPC=c('-10',1,2,3,4,5,6,7,8,9,10,15,20,25,30,35,40,50,60,70,80)
nPC=length(PC)

#GENE
df.data.tmp <-NULL
df.data.tss <-NULL
df.data<-NULL

for (k in 1:length(tss)) {

    for (p in 1:nPC) {
    d = read.table(paste0("/nfs/users/nfs_n/np12/myscratch/MacroMap/Analysis/eQTLs/Macromap_fds/",args[1],"/analysis/results/",tss[k],"_", PC[p], ".GENE.txt.gz"),sep=" ", header=FALSE, stringsAsFactor=FALSE)
    cat(paste0("Reading file", paste0("/nfs/users/nfs_n/np12/myscratch/MacroMap/Analysis/eQTLs/Macromap_fds/",args[1],"/analysis/results/",tss[k],"_", PC[p], ".GENE.txt.gz"), "\n"))
    d = d[!is.na(d$V19),]
    d$q = qvalue(d$V19)$qvalue
    fdr001 = sum(d$q <= 0.01)
    fdr005 = sum(d$q <= 0.05)
    fdr01 = sum(d$q <= 0.10)

    df.data.tmp[[p]]<-data.frame(fdr001,fdr005,fdr01)
    }
  df.data[[k]]<-df.data.tmp
  }
df.data_50K<-do.call(rbind,df.data[[1]])
df.data_250K<-do.call(rbind,df.data[[2]])
df.data_500K<-do.call(rbind,df.data[[3]])
df.data_750K<-do.call(rbind,df.data[[4]])
df.data_1MB<-do.call(rbind,df.data[[5]])

df.data_50K$PC<-PC
df.data_250K$PC<-PC
df.data_500K$PC<-PC
df.data_750K$PC<-PC
df.data_1MB$PC<-PC

plot_all<-function (TSS) {

  max_no<-as.numeric(max(eval(parse(text=paste0("df.data_",'50K',"$fdr",'01'))),eval(parse(text=paste0("df.data_",'250K',"$fdr",'01'))),
              eval(parse(text=paste0("df.data_",'500K',"$fdr",'01'))),eval(parse(text=paste0("df.data_",'750K',"$fdr",'01'))),eval(parse(text=paste0("df.data_",'1MB',"$fdr",'01')))))
  max_rounded<-round_any(max_no,600,f=ceiling)


  plot(0, 0, xlim=c(-10, 85), ylim=c(0, max_rounded), type="n", xlab="#PCs", ylab="#eGenes", main=paste(args[1],"\n",paste0("TSS=",TSS)),xaxt='n',cex.axis=1.5,
       cex.lab=1.7,cex.main=1.7)
  abline(v=seq(0, 100, 10), col="lightgrey", lty=2)
  abline(h=seq(0, max_rounded, 100), col="lightgrey", lty=2)

points(axPC, eval(parse(text=paste0("df.data_",TSS,"$fdr",'001'))), type="o", col=COL[2], pch=20, lwd=2)
points(axPC, eval(parse(text=paste0("df.data_",TSS,"$fdr",'005'))), type="o", col=COL[3], pch=20, lwd=2)
points(axPC, eval(parse(text=paste0("df.data_",TSS,"$fdr",'01'))), type="o", col=COL[4], pch=20, lwd=2)

axis(side=1, at=axPC,labels=PC,cex.axis=1.5)
legend("topright", legend=c("FDR=0.01","FDR=0.05","FDR=0.10"), fill=COL[1:3], bg="white")

}

pdf(paste0("/nfs/users/nfs_n/np12/myscratch/MacroMap/Analysis/eQTLs/Macromap_fds/",args[1],"/analysis/results/",args[1],"_eqtl_numbers.pdf"),height = 8,width = 15)
par(mfrow=c(2,3),mar=c(5,5,4,2)+0.1)
plot_all("50K")
plot_all("250K")
plot_all("500K")
plot_all("750K")
plot_all("1MB")
dev.off()

pdf(paste0("/nfs/users/nfs_n/np12/myscratch/MacroMap/Analysis/eQTLs/Macromap_fds/",args[1],"/analysis/results/",args[1],"_50K_eqtl_numbers.pdf"),width = 12)
par(mar=c(5,5,4,2)+0.1)
plot_all("50K")
dev.off()

pdf(paste0("/nfs/users/nfs_n/np12/myscratch/MacroMap/Analysis/eQTLs/Macromap_fds/",args[1],"/analysis/results/",args[1],"_250K_eqtl_numbers.pdf"),width = 12)
par(mar=c(5,5,4,2)+0.1)
plot_all("250K")
dev.off()

pdf(paste0("/nfs/users/nfs_n/np12/myscratch/MacroMap/Analysis/eQTLs/Macromap_fds/",args[1],"/analysis/results/",args[1],"_500K_eqtl_numbers.pdf"),width = 12)
par(mar=c(5,5,4,2)+0.1)
plot_all("500K")
dev.off()

pdf(paste0("/nfs/users/nfs_n/np12/myscratch/MacroMap/Analysis/eQTLs/Macromap_fds/",args[1],"/analysis/results/",args[1],"_750K_eqtl_numbers.pdf"),width = 12)
par(mar=c(5,5,4,2)+0.1)
plot_all("750K")
dev.off()

pdf(paste0("/nfs/users/nfs_n/np12/myscratch/MacroMap/Analysis/eQTLs/Macromap_fds/",args[1],"/analysis/results/",args[1],"_1MB_eqtl_numbers.pdf"),width = 12)
par(mar=c(5,5,4,2)+0.1)
plot_all("1MB")
dev.off()

eGenes<-data.frame(df.data_1MB$PC[which.max(df.data_1MB$fdr005)],max(df.data_1MB$fdr005))
colnames(eGenes)<-c("PC","eGenes")
write.table(eGenes,file=paste0("/nfs/users/nfs_n/np12/myscratch/MacroMap/Analysis/eQTLs/Macromap_fds/",args[1],"/analysis/results/",args[1],"_1MB_eqtl_numbers.txt"),col.names = T,row.names = F,quote = F,sep="\t")


eGenes<-data.frame(df.data_750K$PC[which.max(df.data_750K$fdr005)],max(df.data_750K$fdr005))
colnames(eGenes)<-c("PC","eGenes")
write.table(eGenes,file=paste0("/nfs/users/nfs_n/np12/myscratch/MacroMap/Analysis/eQTLs/Macromap_fds/",args[1],"/analysis/results/",args[1],"_750K_eqtl_numbers.txt"),col.names = T,row.names = F,quote = F,sep="\t")

eGenes<-data.frame(df.data_500K$PC[which.max(df.data_500K$fdr005)],max(df.data_500K$fdr005))
colnames(eGenes)<-c("PC","eGenes")
write.table(eGenes,file=paste0("/nfs/users/nfs_n/np12/myscratch/MacroMap/Analysis/eQTLs/Macromap_fds/",args[1],"/analysis/results/",args[1],"_500K_eqtl_numbers.txt"),col.names = T,row.names = F,quote = F,sep="\t")

eGenes<-data.frame(df.data_50K$PC[which.max(df.data_50K$fdr005)],max(df.data_50K$fdr005))
colnames(eGenes)<-c("PC","eGenes")
write.table(eGenes,file=paste0("/nfs/users/nfs_n/np12/myscratch/MacroMap/Analysis/eQTLs/Macromap_fds/",args[1],"/analysis/results/",args[1],"_50K_eqtl_numbers.txt"),col.names = T,row.names = F,quote = F,sep="\t")

eGenes<-data.frame(df.data_250K$PC[which.max(df.data_250K$fdr005)],max(df.data_250K$fdr005))
colnames(eGenes)<-c("PC","eGenes")
write.table(eGenes,file=paste0("/nfs/users/nfs_n/np12/myscratch/MacroMap/Analysis/eQTLs/Macromap_fds/",args[1],"/analysis/results/",args[1],"_250K_eqtl_numbers.txt"),col.names = T,row.names = F,quote = F,sep="\t")
