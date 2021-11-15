#!/usr/bin/Rscript

listpackage=c("data.table", "plyr", "optparse", "qqman", "knitr", "kableExtra", "ggplot2")
namepack=rownames(installed.packages())
for(package in listpackage[!c(listpackage %in% namepack)]){
install.packages(package, lib=.libPaths()[2], repos='http://cran.us.r-project.org')
}
library(data.table)
library(plyr)
library("optparse")
library('qqman')
library('knitr')
library(kableExtra)
library("ggplot2")


as.characterspe<-function(x,round){
y<-rep(NA, length(x))
bal=!is.na(x) & x>10**-round
y[bal]<-as.character(format(round(x[bal],round), scientific=F))
y[!bal]<-as.character(format(x[!bal], scientific=T,digits=round))
y
}


GetMergeFile<-function(ListFile){

Cmt<-1
for(File in ListFile){
if(length(readLines(File))>1){
Data<-read.csv(File,stringsAsFactors=F)
if(Cmt==1){
DataF<-Data
}else{
DataF<-rbind.fill(DataF,Data)
}
Cmt<-Cmt+1
}
}
if(Cmt==1)return(as.data.frame(matrix(nrow=0, ncol=3)))
else return(DataF)
}

GetStatValue=function(value){
value=value[!is.na(value)]
if(length(unique(value))>5){
#return(data.frame(N=length(value), Mean=mean(value), Sd=sd(value), Med=media(value),Min=min(value), Max=max(value)))
return(data.frame(Param=c("N","Mean", "Sd","Med","Min","Max"),Val=c(length(value), mean(value), sd(value), median(value),min(value),max(value))))
}else{
TmpData<-as.data.frame(table(value))
names(TmpData)<-c("Param", "Val")
TmpData<-rbind(data.frame(Param="N", Val=length(value)) ,TmpData)
return(TmpData)
}
}

option_list = list(
  make_option(c( "--list_files"), type="character",
              help="file name contains all file with result column and header: Traits Type Site File FilePheno ", metavar="character"),
  make_option(c("--list_files_othertraits"), type="character",
              help="file name contains all file with result column and header: Traits Type Site File FilePheno ", default=""),
  make_option(c("--maf"), type="numeric", default=0.01,
              help="output file name [default= %default]", metavar="character"),
  make_option(c("--out"), type="character",
              help="output file name [default= %default]", metavar="character"),
  make_option(c("-c", "--covariates"), type="character",
              help="covariates to statistical analysis [default= %default]", default="", metavar="character"),
  make_option(c("--max_pval"), type="numeric", default = 0.001,
              help="output file name"),
  make_option(c("--wind_merge"), type="numeric", default = 200000,
              help="output file name"),
  make_option(c("--max_pval_rep"), type="numeric", default = 0.001,
              help="output file name"),
  make_option(c("--gwas_cat"), type="character",
              help="output file name"),
  make_option(c("--genes_info"), type="character",
              help="output file name"),
  make_option(c("--all_file"), type="character",
              help="output file name"),
  make_option(c("--size_win"), type="numeric", default = 25000 ,
              help="output file name"),
  make_option(c("--head_freq"), type="character", default="A1FREQ",
              help="output file name [default= %default]", metavar="character"),
  make_option(c("--head_pval"), type="character", default="P_BOLT_LMM",
              help="output file name [default= %default]", metavar="character"),
  make_option(c("--head_beta"), type="character", default="BETA",
              help="output file name [default= %default]", metavar="character"),
  make_option(c("--head_rs"), type="character", default="SNP",
              help="output file name [default= %default]", metavar="character"),
  make_option(c("--head_se"), type="character", default="SE",
              help="output file name [default= %default]", metavar="character"),
  make_option(c("--type_out"), type="character", default="pdf",
              help="output file name [default= %default]", metavar="character")
);

args = commandArgs(trailingOnly=TRUE)
if(length(args)<2){
opt=list(maf=0.01,gwas_cat="/dataE/AWIGenGWAS/shared/ResultGWAS/Ressource/GWAS_Catalog_V37.tsv",max_pval=1.0E-6,max_pval_rep=1.0E-5,size_win=25000,list_files_othertraits="list_file_other.gwas",all_file="/home/jeantristan/Travail/GWAS/Analyse_pipeline/.tmplist_file_analyse_cholesterol_1_qc",list_files="list_file.gwas",wind_merge=250000,head_rs="SNP",head_beta="BETA",head_se="SE",head_pval="P_BOLT_LMM",covariates="age,sex")
}else{
opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);
}

FreqHead=opt[['head_freq']]
PvalHead=opt[['head_pval']]
BetHead=opt[['head_beta']]
SetHead=opt[['head_se']]


ListForExcel=list()
Param=data.frame(param=c(), value=c())
AllFileToMerge<-readLines(opt[["all_file"]])
Pheno=AllFileToMerge[1]
if(is.null(opt[['out']])){
out=Pheno
}else{
out=opt[[out]]
}
covariate_list=as.vector(unlist(strsplit(opt[['covariates']], split=",")))
Param<-data.frame(Param=names(opt), value=as.vector(unlist(opt)))
ListForExcel[["param"]]=Param
if(!is.null(opt[['list_files']])){
FileInfo=read.table(opt[['list_files']] ,header=T, stringsAsFactors=F)
FileInfo<-FileInfo[FileInfo$Pheno==Pheno,]
ListForExcel[["file_info"]]=FileInfo
if("PhenoFile" %in% names(FileInfo)){
Cmt<-1
ListDataPheno<-list()
for(PhenoFile in FileInfo$PhenoFile){
DataPhenoTmp<-read.table(PhenoFile, header=T, stringsAsFactors=F)
DataPhenoTmp<-na.omit(DataPhenoTmp[ ,c("FID","IID", Pheno,covariate_list)])
if("Plink" %in% names(FileInfo)){
print("deleted ind")
PlinkFam<-read.table(paste(FileInfo$Plink[Cmt],".fam",sep=""), header=F)
DataPhenoTmp<-merge(DataPhenoTmp,PlinkFam[,c(1,2)], by=c(1,2))
}
ListDataPheno[[FileInfo$Type[Cmt]]]<-DataPhenoTmp
ResPheno<-GetStatValue(ListDataPheno[[FileInfo$Type[Cmt]]][,Pheno])
names(ResPheno)[2]<-FileInfo$Type[Cmt]
if(Cmt==1)ResPhenoF<-ResPheno
else ResPhenoF<-merge(ResPhenoF,ResPheno,all=T)
Cmt<-Cmt+1
ListForExcel[["trait_distribution"]]<-ResPhenoF
}
}else{
ResPhenoF=as.data.frame(matrix(nrow=0, col=3))
}
}else{
FileInfo=as.data.frame(matrix(nrow=0, col=3))
ResPhenoF=as.data.frame(matrix(nrow=0, col=3))
ListDataPheno=list()
}

NSig<-GetMergeFile(grep("_nsig.csv", AllFileToMerge, value=T))
if(nrow(NSig)>0){
ToCompute<-c("All",grep("Alpha",names(NSig) ,value=T))
NSig<-NSig[order(NSig$Chro),c("Chro","Type", ToCompute)]
DataFNSig<-aggregate(.~Type,data=NSig[,c("Type", ToCompute)], FUN=sum)
DataFNSig$Chro="All"
DataFNSigF<-rbind(NSig,DataFNSig[,c("Chro","Type", ToCompute)])
for(x in ToCompute[-1])DataFNSig[,x]<-paste(as.characterspe(DataFNSig[,x]/DataFNSig[,"All"]*100,1),DataFNSig[,x])
ListForExcel[["distsig_all"]]=DataFNSigF
}else{
DataFNSig<-as.data.frame(matrix(nrow=0, col=3))
}

dir.create("figure/")
QQPloFigI<-c(grep(".qq.pdf", AllFileToMerge,value=T), grep(".qq.jpeg", AllFileToMerge,value=T))
cat(QQPloFigI,"\n")
Ext<-paste(".",rev(strsplit(QQPloFigI, split=".",fixed=T)[[1]])[1],sep="")
cat(Ext,"\n")
QQPloFig<-paste(gsub(".","-",gsub(Ext,"",paste("figure/",basename(gsub("_","-",QQPloFigI)), sep="")), fixed=T),Ext,sep="")
file.copy(QQPloFigI,QQPloFig)
ManPloFigI<-grep(".man.pdf", AllFileToMerge,value=T)

ManPloFig<-c()
for(ManPloI in ManPloFigI){
ManPlo<-paste(gsub(".","-",gsub(".pdf","",paste("figure/",basename(gsub("_","-",ManPloI)), sep="")), fixed=T),".pdf",sep="")
file.copy(ManPloI,ManPlo)
ManPloFig<-c(ManPloFig,ManPlo)
}
FileRnw="MergeAll.Rnw"


RsResume=GetMergeFile(grep("_rsresume.csv", AllFileToMerge, value=T))
RsDetail=GetMergeFile(grep("_rsgwasdet.csv", AllFileToMerge, value=T))
ListForExcel[["rsgwasdet"]]=RsDetail
ListForExcel[["rsresume"]]=RsResume

WindResume=GetMergeFile(grep("_windresume.csv", AllFileToMerge, value=T))
WindDetail=GetMergeFile(grep("_windgwasdet.csv", AllFileToMerge, value=T))
ListForExcel[["windgwasdet"]]=WindDetail
ListForExcel[["windresume"]]=WindResume

OtherTrait=GetMergeFile(grep("_othertrait.csv", AllFileToMerge, value=T))
ListForExcel[["othertrait"]]=OtherTrait

knit(FileRnw)

file.copy(gsub(".Rnw",".tex",basename(FileRnw)), paste(Pheno,".tex",sep=""),overwrite=T)
library(openxlsx)
options(java.parameters = "-Xmx8000m")
write.xlsx(ListForExcel,file= paste(Pheno,".xlsx",sep=""))






