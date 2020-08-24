MakeGwasCatInfo<-function(DataI, Chro,Pos,Info){
if(!is.null(Info)){
Cmt<-1
if(length(Info)==1){
listhead=strsplit(Info,split=",")[[1]]
}else{
listhead=Info
}
for(head in listhead){
DataHead<-aggregate(as.formula(paste(head,"~",paste(Chro,Pos,sep="+"),sep="")), data=unique(DataI[,c(Chro,Pos,head)]),FUN=function(x)paste(unique(x), collapse=";"))
if(Cmt==1)DataF<-DataHead
else DataF<-merge(DataF,DataHead,all=T, by=c(Chro,Pos))
Cmt<-Cmt+1
}

tmpnbpubli<-aggregate(as.formula(paste(head,"~",paste(Chro,Pos,sep="+"),sep="")), data=DataI, FUN=length)
names(tmpnbpubli)<-c(Chro,Pos,'nbpubli')
DataF<-merge(DataF,tmpnbpubli)
}else{
DataI$Nb<-1
head='Nb'
tmpnbpubli<-aggregate(as.formula(paste(head,"~",paste(Chro,Pos,sep="+"),sep="")), data=DataI, FUN=length)
names(tmpnbpubli)<-c(Chro,Pos,'nbpubli')
DataF<-tmpnbpubli
}
DataF
}
as.characterspe<-function(x,round){
y<-rep(NA, length(x))
bal=!is.na(x) & x>10**-round
y[bal]<-as.character(format(round(x[bal],round), scientific=F))
y[!bal]<-as.character(format(x[!bal], scientific=T,digits=round))
y
}

#
#
#
BuildWind<-function(DataGWASCat,ChroGC,PosGC, InfoGC,DataResBlockI, DataGWAS,PvalGW,BetaGW, RsGW, A1GW, A0GW,SeGW, AfGW,listInfoGC){
DataResBlock<-merge(DataResBlockI,DataGWAS, by.x='RsClump', by.y=RsGW)
DataGWASAll<-merge(DataResBlock,DataGC, by.x=c("Chro", "BPInfo"), by.y=c(ChroGC,PosGC), all.x=T)
DataGWASAll<-DataGWASAll[!is.na(DataGWASAll$nbpubli),]

## merge  between DataGWAS and Data used
Cmt<-1
for(Col in c(RsGC,listInfoGC, 'RsClump')){
tmp<-aggregate(as.formula(paste(Col,"~Chro+BeginBlock+EndBlock")), data=DataGWASAll, FUN=function(x)paste(unique(strspl(x)), collapse=", "))
if(Cmt==1)GWASWind<-tmp
else GWASWind<-merge(GWASWind, tmp,all=T,by=c("Chro","BeginBlock","EndBlock"))
Cmt<-Cmt+1
}
DataGWASAll$Wind<-paste(DataGWASAll$Chro,DataGWASAll$BeginBlock,DataGWASAll$EndBlock)
Cmt<-1
for(Wind in unique(DataGWASAll$Wind)){
SubData<-DataGWASAll[DataGWASAll$Wind==Wind,]
PClump<-min(SubData$PClump)
Rslump<-SubData$RsClump[which.min(SubData$PClump)]
Chrolump=unique(SubData$Chro)
if(!is.null(AfGW))AfClump<-SubData[which.min(SubData$PClump), AfGW] else AfClump=NA
if(!is.null(BetaGW))BetaClump<-SubData[which.min(SubData$PClump),BetaGW] else BetaClump=NA
if(!is.null(A1GW ))A1Clump<-SubData[which.min(SubData$PClump), A1GW] else A1Clump=NA
if(!is.null(A0GW ))A0Clump<-SubData[which.min(SubData$PClump), A0GW] else A0Clump=NA
if(!is.null(SeGW))SeClump<-SubData[which.min(SubData$PClump), SeGW] else SeClump=NA
if(!is.null(PvalGW))PWClump<-SubData[which.min(SubData$PClump), PvalGW] else PWClump<-NA
if(Cmt==1)DataPclumpMin<-data.frame(Wind=Wind,MinPClump=PClump,MinPRsClump=Rslump, MinAfClump=AfClump, MinBetaClump=BetaClump, MinA1Clump=A1Clump, MinA0Clump=A0Clump, MinSeClump=SeClump, PClump=PWClump)
else DataPclumpMin<-rbind(DataPclumpMin ,data.frame(Wind=Wind,MinPClump=PClump,MinPRsClump=Rslump, MinAfClump=AfClump, MinBetaClump=BetaClump, MinA1Clump=A1Clump, MinA0Clump=A0Clump, MinSeClump=SeClump, PClump=PWClump))
Cmt<-Cmt+1
}
GWASWind$Wind<-paste(GWASWind$Chro,GWASWind$BeginBlock,GWASWind$EndBlock)
DataAllWithInfo<-merge(GWASWind,DataPclumpMin, by="Wind")
DataAllWithInfo$p.adjust.bonf=p.adjust(DataAllWithInfo$MinPClump, 'bonferroni')
DataAllWithInfo
}

PlotWindLocusZoom<-function(Chro,BeginBlock, EndBlock, WindSize,DataGWAS, ChroGW,PosGW, RsGW,PvalGW,DataGC,ChroGC, PosGC, DataClump, ListRs, DataGene,ChroGE, BeginGE, EndGE, NameGene, DataBlock){


t_col <- function(color, percent = 50, name = NULL) {
  #      color = color name
  #    percent = % transparency
  #       name = an optional name for the color

## Get RGB values for named color
rgb.val <- col2rgb(color)

## Make new color using input color as base and alpha set by transparency
t.col <- rgb(rgb.val[1], rgb.val[2], rgb.val[3],
             max = 255,
             alpha = (100 - percent) * 255 / 100,
             names = name)

## Save the color
invisible(t.col)
}
## END
ScaleBp=1000000
ColI="blue";ColInitial<-t_col(ColI,70, "bluenorm")
ColC="red3";ColClump=t_col(ColC,70)
ColB='grey';ColBox=t_col(ColB,50)
ColB2='red';ColBox2=t_col(ColB2,75)
ColG='black';ColGWAS=t_col(ColG,50)
##defined lim
xlimI=c(max(1,BeginBlock-WindSize),EndBlock+WindSize)

## select ataGWAS
DataGWAS2<-DataGWAS[DataGWAS[,ChroGW]==Chro & DataGWAS[, PosGW]>=xlimI[1] & DataGWAS[, PosGW]<=xlimI[2],]
DataGWAS2$LogPval<- -log10(DataGWAS2[,PvalGW])
DataGWAS2$bp2<-DataGWAS2[, PosGW]/ScaleBp

DataGWASCat2<-DataGC[DataGC[,ChroGC]==Chro,]
DataGWASCat2$PosBegin372<-DataGWASCat2[,PosGC]/ScaleBp
DataGWASCat2$col<-'grey'
ListRsCat<-gsub(' ','', strsplit(as.character(ListRs),split=',')[[1]])
DataGWASCat2$col[DataGWASCat2[,PosGC] %in% DataGWAS2[DataGWAS2[, RsGW] %in% ListRsCat,PosGW]]<-'red'

DataClump2<-DataClump[DataClump$CHR==Chro & (DataClump$BP>=BeginBlock & DataClump$BP<=EndBlock) ,]
DataClump2$BP2<-DataClump2$BP/ScaleBp
## BaliseRsClump
ListClump<-as.character(DataClump2[,"SP2"])
ListClump<-gsub('\\([0-9]\\)','',unlist(strsplit(ListClump,split=',')))



DataGWAS2$col<-ColInitial
DataGWAS2$bg<-ColInitial
DataGWAS2$pch=21
DataGWAS2$text<-NA
DataGWAS2$coltext<-NA




BaliseRsClump<-DataGWAS2[,RsGW] %in% DataClump2$SNP
#DataGWAS2$col[BaliseRsClump]<-ColC
DataGWAS2$bg[BaliseRsClump & DataGWAS2$col!=ColGWAS]<-ColC
DataGWAS2$cex<-1#((DataGWAS2$LogPval-min(DataGWAS2$LogPval))/(max(DataGWAS2$LogPval)-min(DataGWAS2$LogPval))+0.1)
DataGWAS2$pch[BaliseRsClump]<-22
DataGWAS2$text[BaliseRsClump]<-as.character(DataGWAS2[BaliseRsClump, RsGW])
DataGWAS2$coltext[BaliseRsClump]<-ColC

ylim<-range(DataGWAS2$LogPval);ylim=c(0,ylim[2]*1.1)
xlim<-xlimI/ScaleBp

BaliseRsOtherCl<-DataGWAS2[,RsGW] %in% ListClump
DataGWAS2$bg[BaliseRsOtherCl]<-ColClump
DataGWAS2$col[BaliseRsOtherCl & DataGWAS2$col!=ColGWAS]<-ColClump

PosGWASCat<-DataGWAS2[,PosGW] %in% DataGWASCat2[,PosGC]
DataGWAS2$col[PosGWASCat]<-ColGWAS
DataGWAS2$text[PosGWASCat]<-as.character(DataGWAS2[PosGWASCat, RsGW])
DataGWAS2$coltext[PosGWASCat]<-'black'
DataGWAS2$pch[PosGWASCat]<-23


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
#c(bottom, left, top, right)
nf <- layout( matrix(c(1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,2,2,2,2,2),4,5,byrow=T), respect=T)
par(mar=c(0,4,4,1))
plot(xlim, ylim,type='n', xaxt='n', ylab='-log10(pval)', bty='l')
points(DataGWAS2$bp2, DataGWAS2$LogPval, pch=DataGWAS2$pch, col=DataGWAS2$col, bg=DataGWAS2$bg, cex=DataGWAS2$cex)
text(DataGWAS2$bp2[!is.na(DataGWAS2[,RsGW])], DataGWAS2$LogPval[!is.na(DataGWAS2[,RsGW])]*1.05,  DataGWAS2$text[!is.na(DataGWAS2[,RsGW])], cex=0.8, col=DataGWAS2$coltext)
par(mar=c(5,4,0,1))

plot(xlim, c(0,nb+1), type='n', yaxt='n',xlab=paste('chr', Chro,' (Mb)', sep=''), ylab='')
rect(BeginBlock/ScaleBp,nb+1-.9, EndBlock/ScaleBp,nb+1-0.1, col=ColB2, bg=ColBox2)
rect(DataBlock2$BP12,nb+1-0.9, DataBlock2$BP22,nb+1-0.1, col=ColBox)
arrows(DataGWASCat2$PosBegin372,nb+1-1.3,DataGWASCat2$PosBegin372,nb+1-0.9,code=2,lwd=0.5, length=0.1, col='black')
if(nrow(DataGene2)>0){
for(CmtGene in 1:nrow(DataGene2)){
lines(c(DataGene2[CmtGene,'BEGIN2'],DataGene2[CmtGene,'END2'] ),c(CmtGene-0.5, CmtGene-0.5))
xpostxt<-DataGene2[CmtGene,'BEGIN2']+abs(DataGene2[CmtGene,'BEGIN2']-DataGene2[CmtGene,'END2'])/2
if(xpostxt<xlim[1])xpostxt<-xlim[1]
if(xpostxt>xlim[2])xpostxt<-DataGene2[CmtGene,'BEGIN2']
if(xpostxt<xlim[1])xpostxt<-xlim[1]
text(c(xpostxt,xpostxt), c(CmtGene,CmtGene)+0.1, DataGene2[CmtGene,'GENE'], cex=0.75)
}
}
}

strspl=function(x, sep=";"){
x<-as.character(x)
gsub(" $", "", gsub("^ ", "",unlist(strsplit(x, split=sep))))
}

