#!/usr/bin/Rscript

library(data.table)
library('knitr')
library("optparse")
library('qqman')
library(kableExtra)

as.characterspe<-function(x,round){
y<-rep(NA, length(x))
bal=!is.na(x) & x>10**-round
y[bal]<-as.character(format(round(x[bal],round), scientific=F))
y[!bal]<-as.character(format(x[!bal], scientific=T,digits=round))
y
}



option_list = list(
  make_option("--list_pdf", type="character",
              help="file gwas contains resultat ", metavar="character"),
  make_option("--man_fig", type="character",
              help="file gwas contains resultat ", metavar="character"),
  make_option("--qq_fig", type="character",
              help="file gwas contains resultat ", metavar="character"),
  make_option("--lamb_fig", type="character",
              help="file gwas contains resultat ", metavar="character"),
  make_option("--lamb_val", type="character",
              help="file gwas contains resultat ", metavar="character"),
  make_option("--rs_info", type="character",
              help="file gwas contains resultat ", metavar="character"),
  make_option("--csv_res", type="character",
              help="file gwas contains resultat ", metavar="character"),
  make_option("--csv_out", type="character",
              help="file gwas contains resultat ", metavar="character")
)

opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser)


listpath=strsplit(Sys.getenv('PATH'), split=':')[[1]]

filekniti='do_report.Rnw'
if(!file.exists(filekniti)){
fileknit=""
listpath=strsplit(Sys.getenv('PATH'), split=':')[[1]]
for(path in listpath){
if(file.exists(paste(path,filekniti,sep="/"))){
fileknit=paste(path,filekniti,sep="/")
}
}
}else{
fileknit=filekniti
}
DirPWD=getwd()
#listpdf=strsplit(opt[['list_pdf']],split=',')[[1]]
listpdf=readLines(opt[['list_pdf']])
if(length(listpdf)>0){
PdfMat=cbind(matrix(unlist(strsplit(gsub('.pdf$','',listpdf), split='_')),nrow=length(listpdf),byrow=T),data.frame(FileName=paste(DirPWD,'/',listpdf,sep='')))
PdfMat<-PdfMat[,-1]
names(PdfMat)<-c('Chro', 'Pos', 'Info','Type', 'path')
}else{
PdfMat<-as.data.frame(matrix(nrow=0,ncol=5))
names(PdfMat)<-c('Chro', 'Pos', 'Info','Type', 'path')
}

## qqplot
FileQQ=paste(DirPWD,opt[['qq_fig']],sep='/')
DataLamb=read.table(opt[['lamb_val']],header=T)
FileLamb=paste(DirPWD,opt)
DataCSVI<-read.csv(opt[['csv_res']])
infors=read.table(opt[['rs_info']] ,header=F)
DataCSV<-merge(DataCSVI, infors, by.x=c('chr','bp'), by.y=c("V1","V2"))
DataCSV<-DataCSV[,c("V3","V4",names(DataCSVI))]

ListFigMan=paste(PWD,unlist(strsplit(opt[['man_fig']],split=',')),sep='/')

headCSV<-c("begin windows", "end windows", strsplit(readLines(opt[['csv_res']], 1),split=',')[[1]])
DataCSVPrint<-DataCSV
names(DataCSVPrint)<-headCSV
write.csv(DataCSVPrint,row.names=F, file=opt[['csv_out']])
knit(fileknit)
system('pdflatex do_report')
system('pdflatex do_report')
#library(openxlsx)
#options(java.parameters = "-Xmx8000m")
#write.xlsx(listout_excel,file= paste("resume_tab.xlsx",sep=""))



