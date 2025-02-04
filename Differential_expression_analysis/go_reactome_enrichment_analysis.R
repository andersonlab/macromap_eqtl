library("clusterProfiler")
library("org.Hs.eg.db")
library("UpSetR")
library("ReactomePA",lib.loc ="/nfs/users/nfs_n/np12/myscratch/tmp")

#BC Ctrl 24

args = commandArgs(trailingOnly=TRUE) # input run 

#args[1]<-"Ctrl"
#args[2]<-"CIL"
#args[3]<-"6"

dir_analysis<-paste0(getwd(),"/",args[1],"_vs_",args[2],"_",args[3])
dir.create(file.path(dir_analysis), showWarnings = FALSE)
setwd(file.path(dir_analysis))

#group1<-args[1]
#group2<-args[2]

callEnrichment<-function (group1,group2,hours) 
{

stimulus_com<-paste0(group1,"_vs_",group2)

  genes_path<-paste0("/nfs/users/nfs_n/np12/myscratch/MacroMap/Analysis/Diff_exp_analysis/Macromap_fds/",stimulus_com,"/Diff_expressed_genes_at_5%_FDR_",
                     stimulus_com,"_",args[3],"_pathways.txt")

  genes <-read.table(genes_path,h=F,stringsAsFactors = F)
back<-read.table("/lustre/scratch119/realdata/mdt3/teams/gaffney/np12/MacroMap/Analysis/Diff_exp_analysis/Macromap_fds/Ctrl_vs_CIL/All_back.txt",h=F,stringsAsFactors = F)

genes = bitr(genes$V1, fromType="ENSEMBL", toType="ENTREZID", OrgDb="org.Hs.eg.db")
back<-bitr(back$V1, fromType="ENSEMBL", toType="ENTREZID", OrgDb="org.Hs.eg.db")



ego2 <- enrichGO(gene         = genes$ENTREZID,
                 OrgDb         = org.Hs.eg.db,
                 ont           = "BP",
                 pAdjustMethod = "BH",
                 universe= back$ENTREZID,
                 pvalueCutoff  = 0.01,
                 qvalueCutoff  = 0.05)

ego3 <- setReadable(ego2, OrgDb = org.Hs.eg.db)
write.table(as.data.frame(ego3),file=paste0(stimulus_com,"_",args[3],".GO.txt"),quote = F,row.names = F,col.names = T,sep ="$")


reactome.enric <- enrichPathway(gene = genes$ENTREZID,
                         organism     = 'human',
                         pAdjustMethod = "BH",
                         universe= back$ENTREZID,
                         pvalueCutoff  = 0.01,
                         qvalueCutoff  = 0.05,
                         readable=TRUE)

write.table(as.data.frame(reactome.enric),file=paste0(stimulus_com,"_",args[3],".REAC.txt"),quote = F,row.names = F,col.names = T,sep="$")


}


callEnrichment(args[1],args[2],args[3])
