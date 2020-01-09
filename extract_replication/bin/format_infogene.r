#!/usr/bin/Rscript

library(data.table)
args = commandArgs(trailingOnly=TRUE)
DataGenes<-as.data.frame(fread(args[1], sep="\t", header=F))
DataGenes<-DataGenes[DataGenes$V3=="gene",]
DataGenes$V1<-gsub("chr","",DataGenes$V1)
DataGenes$GeneName<-sapply(strsplit(DataGenes$V9,split=";"),function(x)gsub("\"","",gsub(" ","",gsub("gene_name","",grep("gene_name",x, value=T)))))
DataGenes$GeneID<-sapply(strsplit(DataGenes$V9,split=";"),function(x)gsub("\"","",gsub(" ","",gsub("gene_id","",grep("gene_id",x, value=T)))))
DataGenes$PosGene<-1:nrow(DataGenes)
DataGenes<-DataGenes[,c(1,4,5,10, 11)]
names(DataGenes)<-c("CHR", "BEGIN", "END","GENE", "ID")
write.table(DataGenes, sep="\t", row.names=F, file=args[2], col.names=T)

