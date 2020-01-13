MakeGwasCatInfo<-function(DataI, Chro,Pos,Info){
Cmt<-1
listhead=strsplit(Info,split=",")[[1]]
for(head in listhead){
DataHead<-aggregate(as.formula(paste(head,"~",paste(Chro,Pos,sep="+"),sep="")), data=unique(DataI[,c(Chro,Pos,head)]),FUN=function(x)paste(unique(x), collapse=";"))
if(Cmt==1)DataF<-DataHead
else DataF<-merge(DataF,DataHead,all=T, by=c(Chro,Pos))
Cmt<-Cmt+1
}
tmpnbpubli<-aggregate(as.formula(paste(head,"~",paste(Chro,Pos,sep="+"),sep="")), data=DataI, FUN=length)
names(tmpnbpubli)<-c(Chro,Pos,'nbpubli')
DataF<-merge(DataF,tmpnbpubli)
DataF
}
as.characterspe<-function(x,round){
y<-rep(NA, length(x))
bal=!is.na(x) & x>10**-round
y[bal]<-as.character(format(round(x[bal],round), scientific=F))
y[!bal]<-as.character(format(x[!bal], scientific=T,digits=round))
y
}

PlotWindLocusZoom<-function(Chro,BeginBlock, EndBlock, WindSize,DataGWAS, ChroGW,PosGW, RsGW,PvalGW,DataGC,ChroGC, PosGC, DataClump, ListRs, DataGene,ChroGE, BeginGE, EndGE, NameGene, DataBlock){

t_col <- function(color, percent = 50, name = NULL) {
#      color = color name, percent = % transparency, name = an optional name for the color
## Get RGB values for named color
rgb.val <- col2rgb(color)
## Make new color using input color as base and alpha set by transparency
t.col2 <- rgb(rgb.val[1], rgb.val[2], rgb.val[3],max = 255,alpha = (100 - percent) * 255 / 100,names = name)
## Save the color
invisible(t.col2)
}
ScaleBp=1000000
ColI="gray88";ColInitial=t_col(ColI,50)
ColC="red3";ColClump=t_col(ColC,75)
ColB='grey';ColBox=t_col(ColB,75)
ColB2='red';ColBox2=t_col(ColB2,75)
ColG='black';ColGWAS=t_col(ColG,50)
##defined lim
xlimI=c(max(1,BeginBlock-WindSize),EndBlock+WindSize)

## select ataGWAS
DataGWAS2<-DataGWAS[DataGWAS[,ChroGW]==Chro & DataGWAS[, PosGW]>=xlimI[1] & DataGWAS[, PosGW]<=xlimI[2],]
DataGWAS2$LogPval<- -log10(DataGWAS2[,PvalGW])
DataGWAS2$bp2<-DataGWAS2[, PosGW]/ScaleBp

DataClump2<-DataClump[DataClump$CHR==Chro & (DataClump$BP>=BeginBlock & DataClump$BP<=EndBlock) ,]
DataClump2$BP2<-DataClump2$BP/ScaleBp
BaliseRsClump<-DataGWAS2[,RsGW] %in% DataClump2$SNP
## BaliseRsClump
ListClump<-as.character(DataClump2[,"SP2"])
ListClump<-gsub('\\([0-9]\\)','',unlist(strsplit(ListClump,split=',')))


DataGWAS2$col<-ColInitial
DataGWAS2$bg<-ColI
DataGWAS2$pch=21
DataGWAS2$text<-NA
DataGWAS2$col[BaliseRsClump]<-ColC
DataGWAS2$bg[BaliseRsClump]<-ColC
DataGWAS2$pch[BaliseRsClump]<-22
DataGWAS2$text[BaliseRsClump]<-DataGWAS2[BaliseRsClump, RsGW]

ylim<-range(DataGWAS2$LogPval);ylim=c(0,ylim[2]*1.1)
xlim<-xlimI/ScaleBp

BaliseRsClump<-DataGWAS2[,RsGW] %in% ListClump
DataGWAS2$bg[BaliseRsClump]<-ColClump
DataGWAS2$col[BaliseRsClump]<-ColClump
DataGWAS2$bg[BaliseRsClump]<-ColClump
DataGWAS2$col[BaliseRsClump]<-ColClump

DataGWASCat2<-DataGC[DataGC[,ChroGC]==Chro,]
DataGWASCat2$PosBegin372<-DataGWASCat2[,PosGC]/ScaleBp
DataGWASCat2$col<-'grey'

ListRsCat<-gsub(' ','', strsplit(as.character(ListRs),split=',')[[1]])
DataGWASCat2$col[DataGWASCat2[,PosGC] %in% DataGWAS2[DataGWAS2[, RsGW] %in% ListRsCat,PosGW]]<-'red'

PosGWASCat<-DataGWAS2[,PosGW] %in% DataGWASCat2[,PosGC]
DataGWAS2$bg[PosGWASCat]<-ColGWAS
DataGWAS2$col[PosGWASCat]<-ColG
DataGWAS2$text[PosGWASCat]<-DataGWAS2$RSID[PosGWASCat]

DataGene$BEGIN2<-DataGene[,BeginGE]/ScaleBp
DataGene$END2<-DataGene[,EndGE]/ScaleBp
DataGene$CHR<-DataGene[,ChroGE]

Balise1<-DataGene$CHR==Chro & ((DataGene$END2>=xlim[1] & DataGene$END2<=xlim[2]) | (DataGene$BEGIN2>=xlim[1] & DataGene$BEGIN2<=xlim[2]))
Balise2<-DataGene$CHR==Chro & ((xlim[1]>=DataGene$BEGIN2 & xlim[1]<=DataGene$END2) |(xlim[2]>=DataGene$BEGIN2 & xlim[2]<=DataGene$END2))
DataGene2<-DataGene[Balise1 | Balise2,]


DataBlock2<- DataBlock[DataBlock$CHR==Chro,]
DataBlock2$BP12<-DataBlock2$BP1/ScaleBp
DataBlock2$BP22<-DataBlock2$BP2/ScaleBp
nb=2+nrow(DataGene2)
nf <- layout( matrix(c(1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,2,2,2,2,2),4,5,byrow=T), respect=T)
par(mar=c(1,4,4,1))
plot(xlim, ylim,type='n', xaxt='n', ylab='-log10(pval)')
points(DataGWAS2$bp2, DataGWAS2$LogPval, pch=DataGWAS2$pch, col=DataGWAS2$col, bg=DataGWAS2$bg)
text(DataGWAS2$bp2[!is.na(DataGWAS2$RSID)], DataGWAS2$LogPval[!is.na(DataGWAS2$RSID)]*1.05,  DataGWAS2$text[!is.na(DataGWAS2$RSID)], cex=0.8, col='red')
par(mar=c(5,4,0,1))

plot(xlim, c(0,nb), type='n', yaxt='n',xlab=paste('chr', Chro,' (Mb)', sep=''), ylab='')
rect(BeginBlock/ScaleBp,nb-.5, EndBlock/ScaleBp,nb-0.1, col=ColB2, bg=ColBox2)
rect(DataBlock2$BP12,nb-0.5, DataBlock2$BP22,nb-0.1, col=ColBox)
arrows(DataGWASCat2$PosBegin372,nb-1,DataGWASCat2$PosBegin372,nb-0.6,code=2,lwd=0.5, length=0.1, col=DataGWASCat2$col)
if(nrow(DataGene2)>0){
for(CmtGene in 1:nrow(DataGene2)){
lines(c(DataGene2[CmtGene,'BEGIN2'],DataGene2[CmtGene,'END2'] ),c(CmtGene-0.5, CmtGene-0.5))
xpostxt<-DataGene2[CmtGene,'BEGIN2']+abs(DataGene2[CmtGene,'BEGIN2']-DataGene2[CmtGene,'END2'])/2
if(xpostxt<xlim[1])xpostxt<-xlim[1]
if(xpostxt>xlim[2])xpostxt<-DataGene2[CmtGene,'BEGIN2']
if(xpostxt<xlim[1])xpostxt<-xlim[1]
text(c(xpostxt,xpostxt), c(CmtGene,CmtGene)-0.25, DataGene2[CmtGene,'GENE'])
}
}
}

strspl=function(x, sep=";"){
x<-as.character(x)
gsub(" $", "", gsub("^ ", "",unlist(strsplit(x, split=sep))))
}

