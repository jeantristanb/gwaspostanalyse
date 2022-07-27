#!/usr/bin/Rscript

as.characterspe<-function(x,round){
y<-rep(NA, length(x))
bal=!is.na(x) & x>10**-round
y[bal]<-as.character(format(round(x[bal],round), scientific=F))
y[!bal]<-as.character(format(x[!bal], scientific=T,digits=round))
y
}


library(data.table)
library("optparse")
library('qqman')
WriteCSV<-function(Data, Ent, Type, row.name=F){
write.csv(Data, paste(Ent,Type,".csv",sep=""))
}

CumSumLogLik=function(pvalt) {
  pvalt=pvalt[order(pvalt)]
  logObs=-log10(pvalt)
  logExp=-log10(seq(1,length(pvalt))/(length(pvalt)))
  Coord=data.frame(logExp=logExp,logObs=logObs)
  return(Coord)
}
#SNP    CHR     BP

option_list = list(
  make_option(c( "--list_files"), type="character",
              help="file name contains all file with result column and header: Traits Type Site File FilePheno ", metavar="character"),
  make_option(c( "--list_files_othertraits"), type="character",
              help="file name contains all file with result column and header: Traits Type Site File FilePheno ", metavar="character"),
  make_option(c("--maf"), type="numeric", default=0.01,
              help="output file name [default= %default]", metavar="character"),
  make_option(c("--out"), type="character", default="out",
              help="output file name [default= %default]", metavar="character"),
  make_option(c( "--head_freq"), type="character", default="A1FREQ",
              help="output file name [default= %default]", metavar="character"),
  make_option(c("--head_pval"), type="character", default="P_BOLT_LMM",
              help="output file name [default= %default]", metavar="character"),
  make_option(c("--head_beta"), type="character", default="BETA",
              help="output file name [default= %default]", metavar="character"),
  make_option(c("--head_rs"), type="character", default="SNP",
              help="output file name [default= %default]", metavar="character"),
  make_option(c( "--head_se"), type="character", default="SE",
              help="output file name [default= %default]", metavar="character"),
  make_option(c( "--head_A1"), type="character",
              help="output file name [default= %default]", metavar="character"),
  make_option(c( "--head_A2"), type="character",
              help="output file name [default= %default]", metavar="character"),
  make_option(c( "--head_chr"), type="character", default="CHR",
              help="output file name [default= %default]", metavar="character"),
  make_option(c("--head_bp"), type="character", default="BP",
              help="output file name [default= %default]", metavar="character"),
  make_option(c( "--head_chr_gc"), type="character", default="CHR",
              help="output file name [default= %default]", metavar="character"),
  make_option(c("--head_bp_gc"), type="character", default="BP",
              help="output file name [default= %default]", metavar="character"),
  make_option(c("--pheno"), type="character",
              help="output file name", metavar="character"),
  make_option(c("--max_pval"), type="numeric", default = 0.001,
              help="output file name"),
  make_option(c("--wind_merge"), type="numeric", default = 200000,
              help="output file name"),
  make_option(c( "--max_pval_rep"), type="numeric", default = 0.001,
              help="output file name"),
  make_option(c("--num_chr"), type="character", default = "1",
              help="output file name"),
  make_option(c("--gwas_cat"), type="character", 
              help="output file name"),
  make_option(c("--genes_info"), type="character", 
              help="output file name"),
  make_option(c( "--size_win"), type="numeric", default = 25000 ,
              help="output file name"),
  make_option(c( "--type_out"), type="character", default="pdf",
              help="output file name [default= %default]", metavar="character")
);



MergeChr<-function(x){
x<-x[!is.na(x)]
if(length(x)==0)return("-")
paste(unique(x), collapse="; ")
}
args = commandArgs(trailingOnly=TRUE)
if(length(args)<2){
#opt=list(list_files="/home/jeantristan/Travail/GWAS/VImputed2/TransfertAnaly/list_file.gwas", maf=0.01, out="Test/test", head_freq="A1FREQ", head_rs="SNP", head_chr="CHR", head_bp="BP", pheno="bmi_qc", max_pval=1*10**-5, max_pval_rep=1*10**-4, type_out="pdf", num_chr="1", head_pval="P_BOLT_LMM", head_A1="ALLELE1", head_A2="ALLELE0", gwas_cat="/dataE/AWIGenGWAS/shared/ResultGWAS/Ressource/GWAS_Catalog_V37.tsv", size_win=25000, genes_info="/dataE/AWIGenGWAS/shared/ResultGWAS/Ressource/gencode.v19.genes", wind_merge=200000)
#opt=list(maf=0.01,pheno="ldl_qc",out="ldl_qc_21",head_pval="P_BOLT_LMM",head_freq="A1FREQ",list_files="list_file.gwas",head_bp="BP",head_chr="CHR",head_rs="SNP",head_beta="BETA",head_se="SE",head_A1="ALLELE1",head_A2="ALLELE0",max_pval=1.0E-6,max_pval_rep=1.0E-5,num_chr="21",gwas_cat="/dataE/AWIGenGWAS/shared/ResultGWAS/Ressource/GWAS_Catalog_V37.tsv",genes_info="/dataE/AWIGenGWAS/shared/ResultGWAS/Ressource/gencode.v19.genes",size_win=25000,list_files_othertraits="list_file_other.gwas",wind_merge=250000)
opt=list(maf=0.01,pheno="cholesterol_1_qc",out="cholesterol_1_qc_11",head_pval="P_BOLT_LMM",head_freq="A1FREQ",list_files="list_file.gwas",head_bp="BP",head_chr="CHR",head_rs="SNP",head_beta="BETA",head_se="SE",head_A1="ALLELE1",head_A2="ALLELE0",max_pval=1.0E-6,max_pval_rep=1.0E-5,num_chr=11,gwas_cat="/dataE/AWIGenGWAS/shared/ResultGWAS/Ressource/GWAS_Catalog_V37.tsv",genes_info="/dataE/AWIGenGWAS/shared/ResultGWAS/Ressource/gencode.v19.genes",size_win=25000,list_files_othertraits="list_file_other.gwas",wind_merge=250000)
}else{
opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);
}
CheckHeader<-function(x, data, headname){
if(!any(x %in% names(data))){
cat('header ',x,'not found in ', headname, '\n')
q(status=2)
}
}

## keept argument
ChroHead=opt[['head_chr']]
RsHead=opt[['head_rs']]
BpHead=opt[['head_bp']]
A1Head=opt[['head_A1']]
A2Head=opt[['head_A2']]

FreqHead=opt[['head_freq']]
PvalHead=opt[['head_pval']]
BetHead=opt[['head_beta']]
SetHead=opt[['head_se']]

ChroHeadGc=opt[['head_chr_gc']]
BpHeadGc=opt[['head_bp_gc']]

ListHeadCommon=c(ChroHead, RsHead, BpHead,A1Head,A2Head)
ListHeadCommon<-ListHeadCommon[!is.null(ListHeadCommon)]
ListHeadType=c(FreqHead,BetHead,SetHead,PvalHead)

Alpha=opt[['max_pval']]
AlphaRep=opt[['max_pval_rep']]

##
Maf<-opt[['maf']]
ChroNum=opt[['num_chr']]
EntOut=opt[['out']]
SizeWind<-opt[['size_win']]

FileAll=read.table(opt[['list_files']], header=T, stringsAsFactors=F)
FileAll<-FileAll[FileAll$Pheno==opt[['pheno']],]
ListRes=list()
Cmt<-1

if(!is.null(opt[['gwas_cat']])){
DataGWASCat<-read.table(opt[['gwas_cat']] ,sep="\t",header=T, stringsAsFactors=F)
CheckHeader(ChroHeadGc, DataGWASCat,'gwas catalog');CheckHeader(BpHeadGc, DataGWASCat,'gwas catalog')
DataGWASCat<-DataGWASCat[as.character(DataGWASCat[,ChroHeadGc])==ChroNum ,]
DataGWASCat$PosCat<-1:nrow(DataGWASCat)
}
if(!is.null(opt[['genes_info']])){
DataGeneInfo<-unique(read.table(opt[['genes_info']] ,sep="\t",header=T, stringsAsFactors=F))
DataGeneInfo<-DataGeneInfo[as.character(DataGeneInfo$CHR)==ChroNum,]
DataGeneInfo$PosGene<-1:nrow(DataGeneInfo)
}
for(File in FileAll$File){
Data<-fread(basename(File))
Data<-as.data.frame(Data[as.character(Data[[ChroHead]])==ChroNum])
if(!any(is.numeric(Data[!is.na(Data[,PvalHead]) & !is.infinite(Data[,PvalHead]),PvalHead]))){
Data[,PvalHead]<-as.numeric(Data[,PvalHead])
}
if(!any(is.numeric(Data[!is.na(Data[,FreqHead]) & !is.infinite(Data[,FreqHead]),FreqHead]))){
Data[,FreqHead]<-as.numeric(Data[,FreqHead])
}

Maf2<-1-Maf
Data<-Data[!is.infinite(rowSums(Data[, c(PvalHead,FreqHead)])) & Data[,FreqHead]>=Maf & Data[,FreqHead]<=Maf2 & !is.na(Data[,FreqHead]),]
if(nrow(Data)>0){
ListRes[[FileAll$Type[Cmt]]]<-Data
Cmt<-Cmt+1
}
}
GetPercSig<-function(Data, PvalHead, listalpha){
ListValue=nrow(Data)
for(alpha in listalpha)ListValue<-c(ListValue ,length(Data[Data[,PvalHead]<=alpha,PvalHead]))
tmp=as.data.frame(matrix(ListValue, nrow=1))
names(tmp)<-c("All", paste("Alpha", listalpha))
tmp
}
## compute percentage
Cmt<-1
for(Nam in names(ListRes)){
ResSig<-GetPercSig(ListRes[[Nam]], PvalHead, c(0.05, Alpha, AlphaRep))
ResSig$Type=Nam
ResSig$Chro=ChroNum
if(Cmt==1)ResSigF<-ResSig
else ResSigF<-rbind(ResSigF ,ResSig)
Cmt<-Cmt+1
}
WriteCSV(ResSigF, EntOut, "_nsig")

## select significatif
ListRsAlpha1<-c()
for(Nam in names(ListRes)){
ListRsAlpha1<-c(ListRsAlpha1,ListRes[[Nam]][ListRes[[Nam]][,PvalHead]<Alpha,RsHead])
}
ListRsAlpha1<-unique(ListRsAlpha1)

## Used PosBegin, PosEnd, NumCat
#Data1<-DataFByRsTmp;Data2<-DataGWASCat[,c("BEGIN", "END", "PosCat")];SizeWind<-100000
CompareWindowsChro<-function(PosDeb, PosFin,Groupe,ListPosDeb, ListPosFin, ListGroup, Lim){
PosDeb<-PosDeb-Lim;PosFin<-PosFin+Lim
tmp<-which(((PosDeb>=ListPosDeb & PosDeb<=ListPosFin) | (PosFin>=ListPosDeb & PosFin<=ListPosFin) | (ListPosDeb>=PosDeb & ListPosDeb<=PosFin) | (ListPosFin>=PosDeb & ListPosFin<=PosFin)))
if(length(tmp)==0)return(NA)
ListGroup[tmp]
}


MergeAllChro<-function(Data1, Data2, SizeWind){
for(Cmt in 1:nrow(Data1)){
print(head(Data1))
print(head(Data2))
Wind<-CompareWindowsChro(Data1[Cmt,1], Data1[Cmt,2], Data1[Cmt,3], Data2[,1], Data2[,2],Data2[,3],SizeWind)
Tmp<-data.frame(Pos1=rep(Data1[Cmt,3], length(Wind)), Pos2=Wind)
if(Cmt==1)DataF<-Tmp
else DataF<-rbind(DataF,Tmp)
}
DataF
}


if(length(ListRsAlpha1)>0){
Cmt<-1
for(Nam in names(ListRes)){
Data<-ListRes[[Nam]][ListRes[[Nam]][,RsHead] %in% ListRsAlpha1,c(ListHeadCommon, ListHeadType)]
names(Data)[names(Data) %in% ListHeadType]<-paste(names(Data)[names(Data) %in% ListHeadType],"_",Nam,sep="")
if(Cmt==1)DataFByRs<-Data
else DataFByRs<-merge(DataFByRs,Data, all=T)
Cmt<-Cmt+1
}
if(!is.null(opt[['gwas_cat']])){
DataFByRs$PosGWAS<-1:nrow(DataFByRs)
DataFByRsTmp<-data.frame(BEGIN=DataFByRs[,BpHead]-100, END=DataFByRs[,BpHead]+100, DataFByRs$PosGWAS)
Res<-MergeAllChro(DataFByRsTmp, DataGWASCat[,c(BpHeadGc, BpHeadGc, "PosCat")], SizeWind)
names(Res)<-c("PosGWAS", "PosCat");Res<-Res[!is.na(Res[,2]),]
DataFByRsFGWAS<-merge(DataFByRs,Res,by="PosGWAS",all.x=T)
DataFByRsFGWAS<-merge(DataFByRsFGWAS,DataGWASCat,by="PosCat",all.x=T, suffixes = c("",".cat"),)
DataFByRsFGWAS<-DataFByRsFGWAS[,!(names(DataFByRsFGWAS) %in% c("PosGWAS", "PosCat"))]
WriteCSV(DataFByRsFGWAS, EntOut, "_rsgwasdet")
AggMapped<-aggregate(as.formula(paste("MAPPED_TRAIT~",BpHead,sep="")), data=DataFByRsFGWAS, FUN=MergeChr,na.action="na.pass")
AggRef<-aggregate(as.formula(paste("Ref~",BpHead,sep="")), data=DataFByRsFGWAS, FUN=MergeChr,na.action="na.pass")
DataFByRsF<-merge(merge(DataFByRs,AggMapped,by=BpHead,all.x=T), AggRef, by=BpHead, all.x=T)
}else{
DataFByRsF<-DataFByRs
}
if(!is.null(opt[['genes_info']])){
DataFByRsFT<-DataFByRsF
DataFByRsFT$PosGWAS<-1:nrow(DataFByRsFT)
DataFByRsTmp<-data.frame(BEGIN=DataFByRsFT[,BpHead]-100, END=DataFByRsFT[,BpHead]+100,DataFByRsFT$PosGWAS)
Res<-MergeAllChro(DataFByRsTmp, DataGeneInfo[,c("BEGIN", "END", "PosGene")], SizeWind)
names(Res)<-c("PosGWAS", "PosGene")
DataGeneInfo$InfoTmp<-DataGeneInfo$GENE #paste(DataGeneInfo$GENE," (",DataGeneInfo$BEGIN,"-",DataGeneInfo$END,")",sep="")
TmpRes<-merge(Res,DataGeneInfo[,c("PosGene","InfoTmp")],all.x=T)
GeneMapped<-aggregate(InfoTmp~PosGWAS, data=TmpRes, FUN=MergeChr,na.action="na.pass")
names(GeneMapped)<-c("PosGWAS", "Genes")
DataFByRsFT<-merge(DataFByRsFT,GeneMapped,by="PosGWAS",all.x=T)
DataFByRsF<-DataFByRsFT[,!c(names(DataFByRsFT) %in% c("PosGWAS", "PosGene"))]
}
}else{
Head<-c(ListHeadCommon, paste(rep(ListHeadType,length(names(ListRes))),"_",names(ListRes),sep=""))
DataFByRs<-data.frame(matrix(nrow=0, ncol=length(Head)));names(DataFByRs)<-Head
if(!is.null(opt[['gwas_cat']])){
TmpGWAS<-data.frame(matrix(nrow=0, ncol=ncol(DataGWASCat)));names(TmpGWAS)<-names(DataGWASCat)
names(TmpGWAS)[names(TmpGWAS) %in% names(DataFByRs)]<-paste(names(TmpGWAS)[names(TmpGWAS) %in% names(DataFByRs)],".cat",sep="")
DataFByRsF<-cbind(DataFByRs,TmpGWAS)
WriteCSV(DataFByRsF[,!(names(DataFByRsF) %in% c("PosGWAS", "PosCat"))], EntOut, "_gwasdet")
Head<-c(ListHeadCommon, paste(rep(ListHeadType,length(names(ListRes))),"_",names(ListRes),sep=""), "MAPPED_TRAIT", "Ref")
DataFByRsF<-data.frame(matrix(nrow=0, ncol=length(Head)));names(DataFByRsF)<-Head
}
if(!is.null(opt[['genes_info']])){
DataFByRsF<-data.frame(matrix(nrow=0, ncol=ncol(DataFByRsF)+1));names(DataFByRsF)<-c(Head,"Genes")
}
}
WriteCSV(DataFByRsF,EntOut, "_rsresume" )

DoWindowsPval<-function(TraitRes, Pvalue=10**-4, sizeWindows=100000, GetRsGroup=F, Rs="SNP", Chr="CHR", Pos="BP", Freq="A1FREQ", Beta="BETA", Pval="P_BOLT_LMM"){
#TraitRes<-ListDataRsBest; Pvalue<-1; sizeWindows<-opt[['wind_merge']]; Rs=RsHead;Chr=ChroHead; Pos=BpHead; Freq=FreqHead; Beta=BetHead; Pval=PvalHead
Cmt<-1
ListBestRs<-list()
for(Pop in names(TraitRes)){
ResPop<-TraitRes[[Pop]]
ResPop<-ResPop[,c(Rs,Chr,Pos, Freq, Beta, Pval)]
names(ResPop)[4:ncol(ResPop)]<-paste(names(ResPop)[4:ncol(ResPop)],"_",Pop,sep="")
if(Cmt==1)ResF<-ResPop
else ResF<-merge(ResF , ResPop, by=c(Rs,Chr,Pos),all=T)
Cmt<-Cmt+1
}
names(ResF)<-gsub("-","_", names(ResF))
print(head(ResF))
ListPva<-grep(Pval, names(ResF), value=T)
if(length(ListPva)>1){
ResF2<-ResF[apply(ResF[,ListPva],1, function(x)length(x[!is.na(x) & x<Pvalue])>0),]
}else {
ResF2<-ResF[!is.na(ResF[,ListPva]) & ResF[,ListPva]<Pvalue,]
}
ResF2<-ResF2[order(ResF2[,Chr], ResF2[,Pos]),]
CmtG<-1
for(Chro in unique(ResF2[,Chr])){
ResF2Chr<-ResF2[ResF2[,Chr]==Chro,]
if(nrow(ResF2Chr)>0){
if(nrow(ResF2Chr)==1){
ResF2Chr$Groupe<-paste(Chro,"1",sep="_")
}else{
IsDiff<-c(F,ResF2Chr[2:length(ResF2Chr[,Pos]), Pos]-ResF2Chr[1:(length(ResF2Chr[,Pos])-1),Pos]>sizeWindows)
Tmp<-c(1,which(IsDiff),length(IsDiff)+1)
ResF2Chr$Groupe<-paste(Chro,rep(1:(length(Tmp)-1),Tmp[2:length(Tmp)]-Tmp[1:(length(Tmp)-1)]),sep="_")
}}else{
ResF2Chr<-cbind(ResF2Chr,Groupe=matrix(nrow=0,ncol=1)) #$Groupe<-paste(Chro,1,sep="_")
}
if(CmtG==1)ResF2WithGr<-ResF2Chr
else ResF2WithGr<-rbind(ResF2WithGr ,ResF2Chr)
CmtG<-CmtG+1
}
if(CmtG==1)return(list())
MinGr<-aggregate(as.formula(paste(Pos,"~Groupe")),data=ResF2WithGr,min);names(MinGr)<-c("Groupe","Min_bp")
ChroGr<-aggregate(as.formula(paste(Chr,"~Groupe",sep="")),data=ResF2WithGr,unique);names(ChroGr)<-c("Groupe","Chro")
MaxGr<-aggregate(as.formula(paste(Pos,"~Groupe")),data=ResF2WithGr,max);names(MaxGr)<-c("Groupe","Max_bp")
Npos<-aggregate(as.formula(paste(Pos,"~Groupe")),data=ResF2WithGr,length);names(Npos)<-c("Groupe","nb_pos")
AllGr<-merge(merge(merge(MinGr,MaxGr), Npos), ChroGr)
for(PvalName in ListPva){
Pval<-aggregate(as.formula(paste(PvalName,"~Groupe",sep="")), data=ResF2WithGr, FUN=function(x)min(x, na.rm=T),na.action=NULL)
AllGr<-merge(AllGr,Pval,all=T)
}
if(GetRsGroup){
return(list(AllGr=AllGr,RsGroup=ResF2WithGr))
}else return(AllGr)
}


ListRs<-c()
for(Nam in names(ListRes)){
ListRs<-c(ListRs,ListRes[[Nam]][ListRes[[Nam]][,PvalHead]<AlphaRep,RsHead])
}
ListDataRsBest<-c()
for(Name in names(ListRes))ListDataRsBest[[Name]]<-ListRes[[Name]][ListRes[[Name]][,PvalHead]<AlphaRep,]

suppressWarnings(BestSolWind1<-DoWindowsPval(ListDataRsBest, 1, opt[['wind_merge']], Rs=RsHead,Chr=ChroHead, Pos=BpHead, Freq=FreqHead, Beta=BetHead, Pval=PvalHead))
if(length(BestSolWind1)>0){

BestListRs2<-c()
for(File in names(ListRes)){
TmpF<-ListRes[[File]]
for(CmtBestSolWind1 in 1:nrow(BestSolWind1)){
BestListRs2<-c(BestListRs2,TmpF[(BestSolWind1$Min_bp[CmtBestSolWind1]-SizeWind)<TmpF[, BpHead] & TmpF[,BpHead]<=(BestSolWind1$Max_bp[CmtBestSolWind1]+SizeWind),RsHead])
}
}
BestListRs2<-unique(BestListRs2)
ListDataRsBestF2<-list()
for(File in names(ListRes))ListDataRsBestF2[[File]]<-ListRes[[File]][ListRes[[File]][,RsHead] %in% BestListRs2 ,]
suppressWarnings(BestSolWindTmp<-DoWindowsPval(ListDataRsBestF2, 1, opt[['wind_merge']],GetRsGroup=T ,Rs=RsHead,Chr=ChroHead, Pos=BpHead, Freq=FreqHead, Beta=BetHead, Pval=PvalHead))

BestSolWind<-BestSolWindTmp[['AllGr']]
RsAndInfo<-BestSolWindTmp$RsGroup
#ListPval<-grep(Pval ,names(BestSolWind), value=T)
if(!is.null(opt[['gwas_cat']])){
BestSolWind$Num<-1:nrow(BestSolWind)

for(Cmt in 1:nrow(BestSolWind)){
ResRs<-CompareWindowsChro(BestSolWind$Min_bp[Cmt]-100, BestSolWind$Max_bp[Cmt]+100,BestSolWind$Num[Cmt], DataGWASCat$BEGIN, DataGWASCat$END, DataGWASCat$PosCat, Lim=opt[['size_win']])
ListPOSGWAS<-rep(BestSolWind$Num[Cmt], length(ResRs))
Tmp<-data.frame(PosGWAS=ListPOSGWAS, PosCat=ResRs)
if(Cmt==1){
CorrRsGwasCatGen<-Tmp
}else{
CorrRsGwasCatGen<-rbind(CorrRsGwasCatGen,Tmp)
}
}
BestSolWindGWAS<-merge(BestSolWind,CorrRsGwasCatGen,by.x="Num", by.y="PosGWAS",all.x=T)
BestSolWindGWAS<-merge(BestSolWindGWAS, DataGWASCat,by.x="PosCat", by.y="PosCat", all.x=T)
WriteCSV(BestSolWindGWAS, EntOut, "_windgwasdet")
AggMapped<-aggregate(MAPPED_TRAIT~Min_bp+Max_bp, data=BestSolWindGWAS, FUN=MergeChr,na.action="na.pass")
AggPubli<-aggregate(Ref~Min_bp+Max_bp, data=BestSolWindGWAS, FUN=MergeChr,na.action="na.pass")
BestSolWind<-merge(merge(BestSolWind,AggMapped, by=c("Min_bp","Max_bp"),all.x=T), AggPubli, by=c("Min_bp","Max_bp"),all.x=T)
}else{
writeLines("NA", file=paste(EntOut, "_windgwasdet.csv",sep=""))
}

if(!is.null(opt[['genes_info']])){
for(Cmt in 1:nrow(BestSolWind)){
ResRsGenes<-CompareWindowsChro(BestSolWind$Min_bp[Cmt]-100, BestSolWind$Max_bp[Cmt]+100,BestSolWind$Num[Cmt], DataGeneInfo$BEGIN, DataGeneInfo$END, DataGeneInfo$PosGene, Lim=opt[['size_win']])
ListPOSGWASGenes<-rep(BestSolWind$Num[Cmt], length(ResRsGenes))
GeneTmp<-data.frame(PosGWAS=ListPOSGWASGenes, PosCat=ResRsGenes)
if(Cmt==1){
CorrRsGwasGenes<-GeneTmp
}else{
CorrRsGwasGenes<-rbind(CorrRsGwasGenes,GeneTmp)
}
}
BestSolWindgenes<-merge(BestSolWind,CorrRsGwasGenes,by.x="Num", by.y="PosGWAS",all=T)
BestSolWindgenes<-merge(BestSolWindgenes,DataGeneInfo,by.x="PosCat", by.y="PosGene", all.x=T)
AggGene<-aggregate(GENE~Min_bp+Max_bp, data=BestSolWindgenes, FUN=MergeChr,na.action="na.pass")
BestSolWind<-merge(BestSolWind,AggGene ,by=c("Min_bp","Max_bp"),all.x=T)
}
WriteCSV(BestSolWind,EntOut, "_windresume" )
}else{
writeLines("NA", con=paste(EntOut, "_windresume.csv",sep=""))
writeLines("NA", con=paste(EntOut, "_windgwasdet.csv",sep=""))
}

if(!is.null(opt[['list_files_othertraits']])){
library(plyr)
if(length(ListRsAlpha1)>0){
FileAllOther=read.table(opt[['list_files_othertraits']], header=T, stringsAsFactors=F)
FileAllOther=FileAllOther[FileAllOther$Pheno!=opt[['pheno']], ]
CmtTrait=1
NameModOther2<-c()
for(OtherTraits  in unique(FileAllOther$Pheno)){
ListeFileOther<-FileAllOther$File[FileAllOther$Pheno==OtherTraits]
NameModOther<-FileAllOther$Type[FileAllOther$Pheno==OtherTraits]
ListResOther=list()
BestRsOther<-data.frame(BestRs=ListRsAlpha1)
CmtFile<-1
for(File in ListeFileOther){
Data<-as.data.frame(fread(File,header=T))
Data<-Data[Data[,RsHead] %in% ListRsAlpha1,]
#if(nrow(Data)>0){
Res<-data.frame(Rs=Data[,RsHead] ,paste(as.characterspe(Data[,BetHead],2),as.characterspe(Data[,FreqHead],2),as.characterspe(Data[,PvalHead],3),sep=", "), Data[,PvalHead])
names(Res)[2:3]<-c(NameModOther[CmtFile], paste('pval_', NameModOther[CmtFile],sep=''))
NameModOther2<-c(NameModOther2,NameModOther[CmtFile])
BestRsOther<-merge(BestRsOther,Res,by=1, all.x=T)
#}
CmtFile<-CmtFile+1
}
RsAndInfo2<-RsAndInfo[,c(RsHead,ChroHead,BpHead)]
BestRsOther<-merge(RsAndInfo2, BestRsOther,by.y="BestRs", by.x=RsHead)
BestRsOther$Pheno<-OtherTraits
if(CmtTrait==1)BestRsOtherF<-BestRsOther
else BestRsOtherF<-rbind.fill(BestRsOtherF,BestRsOther)
CmtTrait<-CmtTrait+1
}
if(CmtTrait>1){
BestRsOtherF<-BestRsOtherF[apply(BestRsOtherF[,grep('pval_', names(BestRsOtherF))],1,function(x)any(x[!is.na(x) & is.finite(x)]<AlphaRep)),c(RsHead, unique(NameModOther2), "Pheno")]
BestRsOtherF<-merge(DataFByRsF[, c(RsHead,ChroHead,BpHead)],BestRsOtherF, all.y=T)
WriteCSV(BestRsOtherF,EntOut, "_othertrait")
}else writeLines("None",paste(EntOut, "_othertrait.csv",sep=""))
}else{
writeLines("None",paste(EntOut, "_othertrait.csv",sep=""))
}
}




