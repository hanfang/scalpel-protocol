#!/usr/bin/env Rscript

indel=read.table("type_and_size.txt", header=FALSE)
colnames(indel)= c("type","size")
indel_30=indel[indel[,2]<=30,]
indel.table <- table(indel_30$type,factor(indel_30$size,lev=1:30)  )
pdf('indelsize_by_type.pdf', width=12, height=8)
barplot(indel.table, main="indel distribution within coding sequence (CDS)", xlab="", col=c("green","red"), legend = rownames(indel.table))
dev.off()