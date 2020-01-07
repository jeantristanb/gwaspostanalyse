#!/usr/bin/Rscript
args = commandArgs(trailingOnly=TRUE)
listfile=strsplit(args[1],split=',')[[1]]
Out=args[2]
#DataAllBloc<-read.table('All_block.blocks.det', head=T)
Cmt=1
for(file in listfile){
Data<-read.table(file, header=T)
if(Cmt==1)DataAllBloc<-Data
else DataAllBloc<-rbind(DataAllBloc,Data)
Cmt<-Cmt+1
}

Cmt<-1
for(Chro in unique(DataAllBloc$CHR)){
DataAllBlocChro<-DataAllBloc[DataAllBloc$CHR==Chro,]
DataAllBlocChro[1,'BP1']<-1
DataAllBlocChro[-nrow(DataAllBlocChro),'BP2']<-DataAllBlocChro[2:nrow(DataAllBlocChro),'BP1']-1
if(Cmt==1)DataAllBlocN<-DataAllBlocChro
else DataAllBlocN<-rbind(DataAllBlocN,DataAllBlocChro)
Cmt<-Cmt+1
}
DataAllBlocN$BP1<-format(as.integer(DataAllBlocN$BP1),  drop0trailing = TRUE, scientific = FALSE)
DataAllBlocN$BP2<-format(as.integer(DataAllBlocN$BP2),  drop0trailing = TRUE, scientific = FALSE)
write.table(DataAllBlocN, file=Out, row.names=F, quote=F, col.names=T, sep='\t')

