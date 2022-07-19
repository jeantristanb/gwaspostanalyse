#!/usr/bin/env Rscript

library("optparse")
library('qqman')
suppressMessages(require(data.table))

CumSumLogLik=function(pvalt) {
  pvalt=pvalt[order(pvalt)]
  logObs=-log10(pvalt)
  logExp=-log10(seq(1,length(pvalt))/(length(pvalt)))
  Coord=data.frame(logExp=logExp,logObs=logObs)
  return(Coord)
}

#SNP	CHR	BP

option_list = list(
  make_option(c("--list_files"), type="character", 
              help="file name contains all file with result column and header: Traits Type Site File FilePheno ", metavar="character"),
  make_option(c("--file_gwas"), type="character", 
              help="file name contains all file with result column and header: Traits Type Site File FilePheno ", metavar="character"),
  make_option(c("--maf"), type="numeric", default=0.01,
              help="output file name [default= %default]", metavar="character"),
  make_option(c( "--out_file"), type="character", default="out",
              help="output file name [default= %default]", metavar="character"),
  make_option(c("--head_freq"), type="character", default="A1FREQ",
              help="output file name [default= %default]", metavar="character"),
  make_option(c("--head_pval"), type="character", default="P_BOLT_LMM",
              help="output file name [default= %default]", metavar="character"),
  make_option(c("--head_rs"), type="character", default="SNP",
              help="output file name [default= %default]", metavar="character"),
  make_option(c("--head_chr"), type="character", default="CHR",
              help="output file name [default= %default]", metavar="character"),
  make_option(c("--head_bp"), type="character", default="BP",
              help="output file name [default= %default]", metavar="character"),
  make_option(c("--pheno"), type="character",
              help="output file name", metavar="character"),
  make_option(c("--type"), type="character",
              help="output file name", metavar="character"),
  make_option(c("--max_pval"), type="numeric", default = 0.001,
              help="output file name"),
  make_option(c("--type_out"), type="character", default="pdf",
              help="output file name [default= %default]", metavar="character")
);
args = commandArgs(trailingOnly=TRUE)
if(length(args)<2){
opt=list(maf=0.01,pheno="waist_hip_r_c_qc",out_file="waist_hip_r_c_qc.West.man.pdf",head_pval="P_BOLT_LMM",head_freq="A1FREQ",list_files="list_file.gwas",head_bp="BP",head_chr="CHR",head_rs="SNP",type="West",type_out="pdf")
}else{
opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);
}

if(is.null(opt[['list_files']])){
FileAll=data.frame(Pheno=strsplit(opt[['pheno']],split=',')[[1]], File=strsplit(opt[['file_gwas']],split=',')[[1]], Type=strsplit(opt[['type']], split=',')[[1]])
}else{
FileAll=read.table(opt[['list_files']], header=T, stringsAsFactors=F)
}

if(!is.null(opt[['pheno']])){
FileAll<-FileAll[FileAll$Pheno==opt[['pheno']] & FileAll$Type==opt[['type']],]
}
if(nrow(FileAll)!=1){
cat('not 1 line available in ', opt[['list_files']],"(pheno :",opt[['pheno']],")")
q()
}
print(FileAll$File)
ListeFile<-as.character(FileAll$File)
NameMod<-FileAll$Type
Maf<-opt[['maf']]
#SNP	CHR	BP	GENPOS
Pval<-opt[['head_pval']]
Freq<-opt[['head_freq']]
Rs<-opt[['head_rs']]
Chr<-opt[['head_chr']]
bp<-opt[['head_bp']]
if(opt[['type_out']]=="pdf")fctplot=pdf else if(opt[['type_out']]=="jpeg")fctplot=jpeg else if(opt[['type_out']]=="tiff")fctplot=tiff else q()
ListAll=list()
Data<-as.data.frame(fread(as.character(FileAll$File[1])))
if(!any(is.numeric(Data[!is.na(Data[,Pval]) & !is.infinite(Data[,Pval]),Pval]))){
Data[,Pval]<-as.numeric(Data[,Pval])
}
if(!any(is.numeric(Data[!is.na(Data[,Freq]) & !is.infinite(Data[,Freq]),Freq]))){
Data[,Freq]<-as.numeric(Data[,Freq])
}

Data<-Data[!is.infinite(rowSums(Data[, c(Pval,Freq)])) & Data[,Freq]>Maf & Data[,Freq]<1-Maf & Data[,Pval]<opt[['max_pval']],]
Data<-na.omit(Data[,c(Pval,Freq,Rs,Chr,bp)])

fctplot(opt[["out_file"]])
if(nrow(Data)>0){
manhattan(Data, chr = Chr, bp = bp, p = Pval, snp = Rs)
}else{
plot(c(1,10),c(1,10), xlab="", ylab="",main="", type="n")
text(5,5,"None result", cex=5)
}
dev.off()


