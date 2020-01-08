MakeGwasCatInfo<-function(DataI, Chro,Pos,Info){
Cmt<-1
listhead=strsplit(Info,split=",")[[1]]
for(head in listhead){
DataHead<-aggregate(as.formula(paste(head,"~",paste(Chro,Pos,sep="+"),sep="")), data=unique(DataI[,c(Chro,Pos,head)]),FUN=function(x)paste(unique(x), collapse=";"))
if(Cmt==1)DataF<-data
else DataF<-merge(DataF,data,all=T, by=c(Chro,Pos))
Cmt<-Cmt+1
}
DataF
}
