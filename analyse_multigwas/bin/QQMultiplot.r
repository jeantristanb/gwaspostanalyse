#!/usr/bin/Rscript
library("optparse")
suppressMessages(require(data.table))


CumSumLogLik=function(pvalt) {
  pvalt=pvalt[order(pvalt)]
  logObs=-log10(pvalt)
  logExp=-log10(seq(1,length(pvalt))/(length(pvalt)))
  Coord=data.frame(logExp=logExp,logObs=logObs)
  return(Coord)
}


option_list = list(
  make_option(c("--list_files"), type="character", default="All_trait.in",
              help="file name contains all file with result column and header: Traits Type Site File FilePheno ", metavar="character"),
  make_option(c("--maf"), type="numeric", default=0.01,
              help="output file name [default= %default]", metavar="character"),
  make_option(c("--out_file"), type="character", default="out",
              help="output file name [default= %default]", metavar="character"),
  make_option(c("--head_freq"), type="character", default="A1FREQ",
              help="output file name [default= %default]", metavar="character"),
  make_option(c("--head_pval"), type="character", default="P_BOLT_LMM",
              help="output file name [default= %default]", metavar="character"),
  make_option(c("--pheno"), type="character",
              help="output file name", metavar="character"),
  make_option(c("--type_out"), type="character", default="jpeg",
              help="output file name [default= %default]", metavar="character"),
  make_option(c("--nThread"), type="integer", default=3,
              help="output file name [default= %default]")
);
args = commandArgs(trailingOnly=TRUE)
if(length(args)<2){
opt=list(maf=0.01,pheno="subcutaneous_fat_qc",out_file="subcutaneous_fat_qc.qq.jpeg",head_pval="P_BOLT_LMM",head_freq="A1FREQ",list_files="list_file.gwas",type_out="jpeg")

}else{
opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);
}

FileAll=read.table(opt[['list_files']], header=T, stringsAsFactors=F)
if(!is.null(opt[['pheno']])){
FileAll<-FileAll[FileAll$Pheno==opt[['pheno']],]
}
if(nrow(FileAll)==0){
cat('not lines available in ', opt[['list_files']],"(pheno :",opt[['pheno']],")")
q()
}
ListeFile<-FileAll$File
NameMod<-FileAll$Type
Maf<-opt[['maf']]

Pval<-opt[['head_pval']]
Freq<-opt[['head_freq']]

CmtFile<-1
ListResLog<-list()
MafMax<-1-Maf
for(File in ListeFile){
print(File)
Data<-as.data.frame(fread(File,header=T,nThread=opt[['nThread']]), )
if(!any(is.numeric(Data[!is.na(Data[,Pval]) & !is.infinite(Data[,Pval]),Pval]))){
Data[,Pval]<-as.numeric(Data[,Pval])
}
if(!any(is.numeric(Data[!is.na(Data[,Freq]) & !is.infinite(Data[,Freq]),Freq]))){
Data[,Freq]<-as.numeric(Data[,Freq])
}
Data<-na.omit(Data[ ,c(Pval,Freq)])
Data<-Data[!is.infinite(rowSums(Data[,c(Pval,Freq)])) & Data[,Freq]>Maf & Data[,Freq]<MafMax,]
if(nrow(Data)>0){
ListResLog[[NameMod[CmtFile]]]<-CumSumLogLik(Data[,Pval])
}
CmtFile<-CmtFile+1
}

if(opt[['type_out']]=="pdf")fctplot=pdf else if(opt[['type_out']]=="jpeg")fctplot=jpeg else if(opt[['type_out']]=="tiff")fctplot=tiff else q()

### density
fctplot(opt[["out_file"]])
RangeXLim<-range(sapply(ListResLog,function(x)return(range(x[,1]))))
RangeYLim<-range(sapply(ListResLog,function(x)return(range(x[,2]))))
plot(RangeXLim, RangeYLim, type="n", xlab="Exp", ylab="Obs")
Cmt<-1
for(Col in names(ListResLog)){
lines(ListResLog[[Col]][,1], ListResLog[[Col]][,2], col=Cmt, lwd=2)
Cmt<-Cmt+1
}
abline(a=0,b=1,lty=2, lwd=2)
legend("topleft", legend=names(ListResLog), col=1:Cmt, lty=1, lwd=2, bty="n")
dev.off()

