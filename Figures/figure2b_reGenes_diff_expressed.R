library(data.table)
library(RColorBrewer)
library("colorspace")

# load eQTLs that are reQTL per conditions 
path_to_read_files<-"/nfs/users/nfs_n/np12/myscratch/MacroMap/Data/MacroMap_manuscript_figures/data/"
response_QTLs<-readRDS(file=paste0(path_to_read_files,"response_per_conditions.RDS"))

conditions<-c("CIL_24","CIL_6","Ctrl_24","Ctrl_6","IFNB_24","IFNB_6","IFNG_24","IFNG_6","IL4_24","IL4_6","LIL10_24","LIL10_6","P3C_24","P3C_6","R848_24",
              "R848_6","sLPS_24","sLPS_6","PIC_24","PIC_6","MBP_24","MBP_6","Prec_D0","Prec_D2")

conditions.order<-conditions[order(conditions)]
names(response_QTLs)<-conditions.order[-c(3,4)]

# load diff exprresion 
path='/nfs/users/nfs_n/np12/myscratch/MacroMap/Data/MacroMap_manuscript_figures/data/diff_exp/'
filelist = list.files(path= "/nfs/users/nfs_n/np12/myscratch/MacroMap/Data/MacroMap_manuscript_figures/data/diff_exp",pattern="^All.*.txt",full.names = F)
datalist = lapply(filelist, function(x)read.table(paste0(path,x), header=T)) 
names(datalist)[1:4]<-sapply(strsplit(filelist[c(1:4)],"_"), function (x) {paste(x[2],x[3],x[4],x[5],x[6],sep="_")} )
names(datalist)[5:24] <-sapply(strsplit(filelist[-c(1:4)],"_"), function (x) {paste(x[2],x[3],x[4],strsplit(x[5],".txt")[1],sep="_")} )
datalist<-datalist[-c(3,4)]

names(datalist)<-sapply(strsplit(names(datalist),"vs_"),"[[",2)
datalist<-datalist[order(names(datalist))]

# define up/down regulation/ diff expression and reGenes that are/are not diff exp 
response_diff_genes<-NULL
shared_diff_QTLs_genes<-NULL
up_down_regulation<-NULL
response_diff_genes_stats<-NULL

for (i in 1:length(response_QTLs)) {
  x<-response_QTLs[[i]]
  diff_exp_genes<-datalist[[i]]
  diff_exp_genes$padj[diff_exp_genes$padj==0]<-2.225074e-308
  
  response_QTLs_genes<-substr(names(x[which(x)]),1,15)
  shared_QTLs_genes<-substr(names(x[which(!x)]),1,15)
  
  diff_exp_genes <-diff_exp_genes[(abs(diff_exp_genes$log2FoldChange) > 1 & diff_exp_genes$padj <0.05),]
  
  response_diff_genes[[i]]<-table(response_QTLs_genes %in% diff_exp_genes$Gene)
  shared_diff_QTLs_genes[[i]]<-table(shared_QTLs_genes %in% diff_exp_genes$Gene)
  up_down_regulation[[i]]<-table(diff_exp_genes$log2FoldChange[diff_exp_genes$Gene %in% response_QTLs_genes] >0)
  response_diff_genes_stats[[i]]<-cbind(diff_exp_genes[diff_exp_genes$Gene %in% response_QTLs_genes,],reQTL=names(x[which(x)][which(response_QTLs_genes %in% diff_exp_genes$Gene)]))
}

names(response_diff_genes)<-names(response_QTLs)
names(shared_diff_QTLs_genes)<-names(response_QTLs)
names(response_diff_genes_stats)<-names(response_QTLs)

response_diff_genes_stats<-lapply(response_diff_genes_stats, function (x) x[order(x$padj),])


df_response_diff_genes<-data.frame(matrix(unlist(response_diff_genes), nrow=length(response_diff_genes), byrow=TRUE))
names(df_response_diff_genes)<-c("FALSE","TRUE")
rownames(df_response_diff_genes)<-names(response_QTLs)
df_response_diff_genes<-t(df_response_diff_genes)

df_up_down_regulation<-data.frame(matrix(unlist(up_down_regulation), nrow=length(up_down_regulation), byrow=TRUE))
names(df_up_down_regulation)<-c("FALSE","TRUE")
rownames(df_up_down_regulation)<-names(response_QTLs)
df_up_down_regulation<-t(df_up_down_regulation)


#################################################################################################################################################################################################################################
# plot them 
conditions<-conditions.order


sum_resp<-length(unique(sapply(strsplit(as.character(unlist(sapply( response_QTLs, function (x) as.character(names(x[x=="TRUE"]))))),"_"),"[[",1)))
response_QTLs<-sapply(response_QTLs, function (x) table(x))
colnames(response_QTLs)<-conditions[-c(3,4)]
response_QTLs<-response_QTLs[,order(response_QTLs[2,]/response_QTLs[1,])]
response_QTLs<-apply(response_QTLs,2,rev)
response_QTLs<-response_QTLs*-1

df_response_diff_genes<-df_response_diff_genes[,match(colnames(response_QTLs),colnames(df_response_diff_genes))]
df_response_diff_genes<-apply(df_response_diff_genes,2,rev)

df_up_down_regulation<-df_up_down_regulation[,match(colnames(response_QTLs),colnames(df_up_down_regulation))]
df_up_down_regulation<-apply(df_up_down_regulation,2,rev)

a<-rbind(df_up_down_regulation,df_response_diff_genes[2,])


col_heat<-sequential_hcl(8, palette = "Purp" )
col_brew<-c(brewer.pal(4, "Accent")[c(2,1)],brewer.pal(10, "Paired")[c(10,9)])

pdf("/nfs/users/nfs_n/np12/myscratch/MacroMap/Data/MacroMap_manuscript_figures/plots/Fig2b.response_eQTLs_being_diff_expressed.pdf")
par(mar = c(5,3.5,2,2))
par(mfrow=c(1,2))

mp<-barplot(response_QTLs,las=2,horiz=TRUE,col=c(col_heat[7],col_brew[2]),xlim=c(-4000,0),cex.names = 1.1,cex.main=1.3,main="response eQTLs",axisnames=FALSE,axes=F)
axis(1, at=seq(-4000,0,500), labels=seq(-4000,0, 500)*-1,las=3,cex.axis=1.3)
t<-round(apply(response_QTLs,2, function(x) x[1]/(x[1]+x[2])*100),1)
text(response_QTLs[1,]-1700,mp - 0.1, labels = paste0(t,"%"," (",response_QTLs[1,]*-1,") "), pos = 4,cex=1,col="black",font=2)

mp<-barplot(a,las=2,horiz=TRUE,col= c(col_heat[1],col_heat[4],col_heat[7]),xlim=c(0,1000),cex.names = 1.1,cex.main=1.3,main="Differentially expressed \nresponse eQTLs",axisnames=TRUE,axes=F)
axis(1, at=seq(0,1000,100), labels=seq(0, 1000, 100),las=3,cex.axis=1.3)
t<-round(apply(df_response_diff_genes,2, function(x) x[1]/(x[1]+x[2])*100),1)
text(df_response_diff_genes[1,]-50,mp - 0.1, labels =paste0(" (",df_response_diff_genes[1,],") ",t,"%"), pos = 4,cex=1,col="black",font=2)
dev.off()

round(mean(100*a[1,]/df_response_diff_genes[1,]),1) # 75.8
round(mean(100*a[2,]/df_response_diff_genes[1,]),1) # 24.2 

range(100*a[1,]/df_response_diff_genes[1,]) # 37.50000 91.01124
range(100*a[2,]/df_response_diff_genes[1,]) # 8.988764 62.500000




pdf("/nfs/users/nfs_n/np12/myscratch/MacroMap/Data/MacroMap_manuscript_figures/plots/Fig2_support.response_eQTLs_being_diff_expressed_no_up_down.pdf")
par(mar = c(5,3.5,2,2))
par(mfrow=c(1,2))

#barplot(response_QTLs,horiz=T,axes=F,las=1,axisnames=TRUE,las=2,cex.names = 1)
mp<-barplot(response_QTLs,las=2,horiz=TRUE,col=c(col_heat[5],col_brew[2]),xlim=c(-4000,0),cex.names = 1.1,cex.main=1.3,main="response eQTLs",axisnames=FALSE,axes=F)
axis(1, at=seq(-4000,0,500), labels=seq(-4000,0, 500)*-1,las=3,cex.axis=1.3)
t<-round(apply(response_QTLs,2, function(x) x[1]/(x[1]+x[2])*100),1)
text(response_QTLs[1,]-1700,mp - 0.1, labels = paste0(t,"%"," (",response_QTLs[1,]*-1,") "), pos = 4,cex=1,col="black",font=2)


#barplot(df_response_diff_genes,axes=F,horiz=T,axisnames=FALSE,xlim=c(-1000,0))
mp<-barplot(df_response_diff_genes,las=2,horiz=TRUE,col= c(col_heat[1],col_heat[5]),xlim=c(0,1000),cex.names = 1.1,cex.main=1.3,main="Differentially expressed \nresponse eQTLs",axisnames=TRUE,axes=F)
axis(1, at=seq(0,1000,100), labels=seq(0, 1000, 100),las=3,cex.axis=1.3)
t<-round(apply(df_response_diff_genes,2, function(x) x[1]/(x[1]+x[2])*100),1)
text(df_response_diff_genes[1,]-50,mp - 0.1, labels =paste0(" (",df_response_diff_genes[1,],") ",t,"%"), pos = 4,cex=1,col="black",font=2)
dev.off()




