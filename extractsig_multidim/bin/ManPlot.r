#!/usr/bin/Rscript
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
  make_option(c("--file"), type="character", default="All_trait.in",
              help="file name contains all file with result column and header: Traits Type Site File FilePheno ", metavar="character"),
  make_option(c("--maf"), type="numeric", default=0.01,
              help="output file name [default= %default]", metavar="character"),
  make_option(c("--header"), type="character",
              help="output file name", metavar="character"),
  make_option(c("--listrs"), type="character",
              help="output file name", metavar="character"),
  make_option(c("--num"), type="character",
              help="output file name", metavar="character"),
  make_option(c("--max_pval"), type="numeric", default = 1,
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

#FileAll=read.table(opt[['list_files']], header=T, stringsAsFactors=F)
#if(!is.null(opt[['pheno']])){
header=opt [['header']]
Param<-matrix(unlist(strsplit(strsplit(header,split=',')[[1]],split=':')),ncol=2, byrow=T)
Param[,1]<-toupper(Param[,1])
File<-opt[['file']]
Maf<-as.numeric(opt[['maf']])
MaxPval<-as.numeric(opt[['max_pval']])
listpos=read.table(opt[['listrs']], header=F)

#SNP	CHR	BP	GENPOS
#"RSID:rs,CHRO:chr,POS:ps,A1:allele1,A2:allele0,BETA:beta,SE:se,PVAL:p_wald,FREQA1:af,TYPE:rin,INFO:LogNewAcr,REFSIG:1"
Pval<-Param[Param[,1]=='PVAL',2]
if(any(Param[,1]=='FREQA1'))Freq<-Param[Param[,1]=='FREQA1',2] else Freq<-NA
Rs<-Param[Param[,1]=='RSID',2]
Chr<-Param[Param[,1]=='CHRO',2]
bp<-Param[Param[,1]=='POS',2]
out<-paste(opt[['num']],'_',Param[Param[,1]=='TYPE',2],"_", Param[Param[,1]=='INFO',2],sep='')

#if(opt[['type_out']]=="pdf")fctplot=pdf else if(opt[['type_out']]=="jpeg")fctplot=jpeg else if(opt[['type_out']]=="tiff")fctplot=tiff else q()
fctplot=eval(parse(text = opt[['type_out']]))
if( opt[['type_out']]=='pdf'){
width=7*2
}else{

width=480*2
}
ListAll=list()
Data<-as.data.frame(fread(File,header=T))
listrs<-Data[paste(Data[,Chr],Data[,bp]) %in% paste(listpos[,1], listpos[,2]),Rs]
if(!any(is.numeric(Data[!is.na(Data[,Pval]) & !is.infinite(Data[,Pval]),Pval]))){
Data[,Pval]<-as.numeric(Data[,Pval])
}
balisefreq=T
if(!is.na(Freq)){
Data[,Freq]<-as.numeric(Data[,Freq])
balisefreq<- Data[,Freq]>Maf & Data[,Freq]<1-Maf
}

Data<-Data[balisefreq & !is.infinite(Data[, Pval]) & !is.na(Data[, Pval]),]
ValPval=CumSumLogLik(Data[,Pval])
write.table(ValPval, file=paste(out,'qq',sep='.'), row.names=F, col.names=F, sep='\t')

Data<-Data[Data[,Pval]<MaxPval,]
Data<-na.omit(Data[,c(Pval,Rs,Chr,bp)])

fctplot(paste(out,opt[['type_out']],sep='.'), width=width)
if(nrow(Data)>0){
manhattan(Data, chr = Chr, bp = bp, p = Pval, snp = Rs, highlight=listrs)
}else{
plot(c(1,10),c(1,10), xlab="", ylab="",main="", type="n")
text(5,5,"None result", cex=5)
}
dev.off()
