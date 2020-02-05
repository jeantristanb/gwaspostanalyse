#!/usr/bin/Rscript
library(plyr)
args = commandArgs(trailingOnly=TRUE)
if(length(args)!=2){
print('argument number false for mergefile')
print(args)
q(2)
}

listfile=strsplit(args[1], split=',')[[1]]
CmtRs<-1
Cmt<-1
listres=list()
for(File in listfile){
DataI<-read.table(File, header=T, sep='\t')
Head=strsplit(readLines(File, 1),split='\t')[[1]]
names(DataI)<-Head
if('rsid' %in% names(DataI)){
Data<-DataI[,!c(names(DataI)=='rsid')]
if(CmtRs==1)DataRs<-DataI[,c('chr','bp','rsid')]
else DataRs<-merge(DataRs,unique(DataI[,c('chr','bp','rsid')]), all=T, by=c('chr','bp'))
CmtRs<-CmtRs+1
}else{
Data<-DataI
}
Type=as.character(unique(Data$Type))
if(!(Type %in% names(listres)))listres[[Type]]=Data
else listres[[Type]]=merge(listres[[Type]] ,Data, by=c('Type','chr','bp'),all=T)
Cmt<-Cmt+1
}
print(names(listres))
Cmt<-1
for(Type in names(listres)){
if(Cmt==1)DataF<-listres[[Type]]
else DataF<-rbind.fill(DataF, listres[[Type]])
Cmt<-Cmt+1
}
if(CmtRs>1){
DataRsF<-cbind(DataRs[,c(1,2)],rsid=apply(DataRs[,-c(1,2)],1, function(x){x<-x[!is.na(x)];paste(gsub(' ','',unique(x)),collapse=';')}))
DataF2<-merge(DataRsF,DataF,  by=c('chr','bp'),all=T)
head(DataF2)
DataF2<-DataF2[,c('Type', 'rsid', 'chr','bp',names(DataF2)[!(names(DataF2) %in% c('Type', 'rsid', 'chr','bp'))])]
}else{
DataF2<-DataF
}
write.csv(DataF2, file=paste(args[2],'.csv',sep=''),row.names=F, quote=F)

