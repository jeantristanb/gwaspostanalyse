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
  make_option("--csv_res", type="character",
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
listpdf=strsplit(opt[['list_pdf']],split=',')[[1]]
PdfMat=cbind(matrix(unlist(strsplit(gsub('.pdf$','',listpdf), split='_')),nrow=length(listpdf),byrow=T),data.frame(FileName=paste(DirPWD,'/',listpdf,sep='')))
PdfMat<-PdfMat[,-1]
names(PdfMat)<-c('Chro', 'Pos', 'Info','Type', 'path')
print(PdfMat)
DataCSV<-read.csv(opt[['csv_res']])
headCSV<-strsplit(readLines(opt[['csv_res']], 1),split=',')[[1]]
knit(fileknit)
system('pdflatex do_report')
system('pdflatex do_report')
#library(openxlsx)
#options(java.parameters = "-Xmx8000m")
#write.xlsx(listout_excel,file= paste("resume_tab.xlsx",sep=""))



