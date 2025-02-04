library(data.table)
gencode = fread("/nfs/users/nfs_n/np12/myscratch/MacroMap/Annotation/gencode.v27.annotation.ref.txt",header=F)
setnames(gencode, c("chr","Start","End","ensembl_id","gene_id"))
#gencode[,CHR:=gsub("chr","",CHR)]

ens = function(gene){
  if(length(gene) == 1){
    return(gencode[grep(paste0("^",gene,"$"),gencode$gene_id),]$ensembl_id)
  } else {
    return(sapply(gene, function(gg) gencode[grep(paste0("^",gg,"$"),gencode$gene_id),]$ensembl_id))
  }
}

ref = function(gene){
  if(length(gene) == 1){
    return(gencode[gencode$ensembl_id == gene,]$gene_id)
  } else {
    return(sapply(gene, function(gg) gencode[gencode$ensembl_id == gg,]$gene_id))
  }  
}

ref2 = function(gene){
  if(length(gene) == 1){
    return(gencode[substr(gencode$ensembl_id,1,15) == gene,]$gene_id)
  } else {
    return(sapply(gene, function(gg) gencode[substr(gencode$ensembl_id,1,15) == gg,]$gene_id))
  }
}

