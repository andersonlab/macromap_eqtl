args <- commandArgs(trailingOnly = TRUE)

a<-read.table(args[1],h=T,check.names=FALSE)
a1<-names(a)[-c(1)]
b<-read.table(args[2],h=T,check.names=FALSE)

b1<-names(b)[-c(1)]
a2<-a[!names(a) %in% a1[!a1 %in% b1]]
b2<-b[!names(b) %in% b1[!b1 %in% names(a2[-c(1)])]]
a3<-a2[-c(1)][match(names(b2)[-c(1)],names(a2)[-c(1)])]
#match(names(a3), names(b2))

write.table(data.frame(PCs=a$PCs,a3,check.names=F),file=args[3],sep="\t",quote=F,col.names=T,row.names=F)
write.table(b2,file=args[4],sep="\t",quote=F,col.names=T,row.names=F)
