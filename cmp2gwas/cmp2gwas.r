library(data.table)
library(gridExtra)
library(ggplot2)
library(optparse)
option_list = list(
  make_option(c("--data_f1"), type="character",
              help="data files with phenotype and FID", metavar="character"),
  make_option(c("--data_f2"), type="character",
              help="data files with phenotype and FID", metavar="character"),
  make_option(c("--chr_f1"), type="character",
              help="data files with phenotype and FID", metavar="character"),
  make_option(c("--chr_f2"), type="character",
              help="data files with phenotype and FID", metavar="character"),
  make_option(c("--bp_f1"), type="character",
              help="data files with phenotype and FID", metavar="character"),
  make_option(c("--bp_f2"), type="character",
              help="data files with phenotype and FID", metavar="character"),
  make_option(c("--a0_f1"), type="character",
              help="data files with phenotype and FID", metavar="character"),
  make_option(c("--a0_f2"), type="character",
              help="data files with phenotype and FID", metavar="character"),
  make_option(c("--a1_f1"), type="character",
              help="data files with phenotype and FID", metavar="character"),
  make_option(c("--a1_f2"), type="character",
              help="data files with phenotype and FID", metavar="character"),
  make_option(c("--beta_f1"), type="character",
              help="data files with phenotype and FID", metavar="character"),
  make_option(c("--beta_f2"), type="character",
              help="data files with phenotype and FID", metavar="character"),
  make_option(c("--se_f1"), type="character",
              help="data files with phenotype and FID", metavar="character"),
  make_option(c("--se_f2"), type="character",
              help="data files with phenotype and FID", metavar="character"),
  make_option(c("--af_f1"), type="character",
              help="data files with phenotype and FID", metavar="character"),
  make_option(c("--af_f2"), type="character",
              help="data files with phenotype and FID", metavar="character"),
  make_option(c("--head_f1"), type="character",
              help="data files with phenotype and FID", metavar="character"),
  make_option(c("--head_f2"), type="character",
              help="data files with phenotype and FID", metavar="character"),
  make_option(c("--out"), type="character", default="out",
              help="output file name [default= %default]", metavar="character")
);
#File<-opt[['data_f1']];listinfo=opt;listhead<-grep('data_f1',grep('_f1',names(opt),value=T), invert=T, value=T)
readfile<-function(File, listinfo, listhead){
Data<-fread(File)
#newname<-names(listinfo[listhead])
#value<-unlist(listinfo[listhead])
initialname<-c()
newname<-c()
for(val in listhead){
if(!is.null(listinfo[[val]])){
initialname<-c(initialname,listinfo[[val]])
newname<-c(newname,val)
}
}
Data<-Data[ ,as.vector(initialname),with = FALSE]
names(Data)<-newname
return(Data)
}



args = commandArgs(trailingOnly=TRUE)
if(length(args)==0){
#rsid	chr	bp	a1	a0	beta	se	p	freqA1	n
opt=list(
data_f1='/home/jeantristan/Data/DataGWASCKD/UACR/Format/Teumeretal2019/formatted_20180202-UACR_overall-AA-nstud_1-SumMac_400.rsid.format.awigen',
chr_f1='chr',
bp_f1='bp',
beta_f1='beta',
a1_f1='a1',
a0_f1='a0',
se_f1='se',
af_f1='freqA1',
p_f1='p',
#chr	rs	ps	n_miss	allele1	allele0	af	beta	se	logl_H1	l_remle	p_wald
data_f2='/home/jeantristan/Travail/GWAS/GWAS_CKD/ImputedDataV3/Result/agesexpca/Res/LogNewAcr/LogNewAcr_All_agesexpca_20190120.imp.stat',
chr_f2='chr',
bp_f2='ps',
beta_f2='beta',
a1_f2='allele1',
a0_f2='allele0',
se_f2='se',
af_f2='af',
p_f2='p_wald'
head_f1='Teumer'
head_f2='Awigen'
)
}else{
opt = OptionParser(option_list=option_list);
opt = parse_args(opt);
}

Head1=opt[['head_f1']]
Head2=opt[['head_f1']]
GetCommon<-function(ChrPos1, ChrPos2,Head1, Head2){
NbPos1<-length(ChrPos1)
NbPos2<-length(ChrPos2)
NbComm<-length(ChrPos1[ChrPos1 %in% ChrPos2])
ResumeNbBase<-data.frame(Type=c(paste('Just ',Head1,sep=''), 'Common', paste('Just ',Head2,sep='')), Count=c(NbPos1-NbComm, NbComm,NbPos2-NbComm))
ResumeNbBase$Perc<-ResumeNbBase$Count/sum(ResumeNbBase$Count)*100
ResumeNbBase
}

Data1<-readfile(opt[['data_f1']], opt,grep('data_f1',grep('_f1',names(opt),value=T), invert=T, value=T))
Data2<-readfile(opt[['data_f2']], opt,grep('data_f2',grep('_f2',names(opt),value=T), invert=T, value=T))


DataAll=merge(Data1,Data2, by.x=c('chr_f1', 'bp_f1'),  by.y=c('chr_f2', 'bp_f2'), all=T)

## exchange b and af 
DataAll$af_f2_old<-DataAll$af_f2
DataAll$beta_f2_old<-DataAll$beta_f2
DataAll$a1_f2_old<-DataAll$a1_f2
DataAll$a0_f2_old<-DataAll$a0_f2
balise<-!is.na(DataAll$a1_f1) & !is.na(DataAll$a1_f2) & DataAll$a1_f1!=DataAll$a1_f2

DataAll$af_f2[balise]<- 1 - DataAll$af_f2_old[balise]
DataAll$beta_f2[balise]<- - DataAll$beta_f2_old[balise]
DataAll$a1_f2[balise]<-DataAll$a0_f2_old[balise]
DataAll$a0_f2[balise]<-DataAll$a1_f2_old[balise]

## 
## co
Data1$ChrPos<-paste(Data1$chr_f1, Data1$bp_f1)
Data2$ChrPos<-paste(Data2$chr_f2, Data2$bp_f2)
DataAll$ChrPos<-paste(DataAll$chr_f1, DataAll$bp_f1)


ResumeNbBase<-GetCommon(Data1$ChrPos,Data2$ChrPos,Head1, Head2)
## 
plotcommon<-ggplot(data=ResumeNbBase, aes(x=Type, y=Perc)) +
  geom_bar(stat="identity", fill="steelblue")+geom_text(aes(label=Count), vjust=1.6, color="white", size=3.5)+
  theme_minimal()

##
Cmt<-1
methodr2="spearman"
for(alpha in c(1,0.05,0.01,10**-3,10**-4,10**-5,10**-6,10**-7)){
Balise<-(DataAll$p_f1<alpha |DataAll$p_f2<alpha) #& DataAll$a1_f2== DataAll$a1_f1 & DataAll$a0_f2== DataAll$a0_f1
DataAllSig<-DataAll[Balise,]
Nbbaselim1<-GetCommon(DataAllSig$ChrPos[!is.na(DataAllSig$p_f1) & DataAllSig$p_f1<alpha],DataAllSig$ChrPos[!is.na(DataAllSig$p_f2)],Head1, Head2)
Nbbaselim2<-GetCommon(DataAllSig$ChrPos[!is.na(DataAllSig$p_f2) & DataAllSig$p_f2<alpha],DataAllSig$ChrPos[!is.na(DataAllSig$p_f1)],Head1, Head2)

## plot
for(alpha2 in c(1,0.05,0.01,10**-3,10**-4)){
Bal1SigCmp1<-!is.na(DataAllSig$p_f1) & DataAllSig$p_f1<alpha & !is.na(DataAllSig$p_f2)  & DataAllSig$a1_f2== DataAllSig$a1_f1 & DataAllSig$a0_f2== DataAllSig$a0_f1 & DataAllSig$p_f2<alpha2
corbeta1<-cor(DataAllSig$beta_f1[Bal1SigCmp1], DataAllSig$beta_f2[Bal1SigCmp1], method=methodr2)
Bal1SigCmp2<-!is.na(DataAllSig$p_f1) & DataAllSig$p_f2<alpha & !is.na(DataAllSig$p_f2)  & DataAllSig$a1_f2== DataAllSig$a1_f1 & DataAllSig$a0_f2== DataAllSig$a0_f1  & DataAllSig$p_f1<alpha2
corbeta2<-cor(DataAllSig$beta_f1[Bal1SigCmp2], DataAllSig$beta_f2[Bal1SigCmp2], method=methodr2)
#plot(DataAllSig$beta_f1[Bal1SigCmp2], DataAllSig$beta_f2[Bal1SigCmp2])
PercSamSig2<-sum((diag(table(DataAllSig$beta_f1[Bal1SigCmp2]<0, DataAllSig$beta_f2[Bal1SigCmp2]<0)/length(DataAllSig$beta_f2[Bal1SigCmp2])))*100)
PercSamSig1<-sum((diag(table(DataAllSig$beta_f1[Bal1SigCmp1]<0, DataAllSig$beta_f2[Bal1SigCmp1]<0)/length(DataAllSig$beta_f2[Bal1SigCmp1])))*100)
NbVal1<-length(which(Bal1SigCmp1))
NbVal2<-length(which(Bal1SigCmp2))

corbeta1abs<-cor(abs(DataAllSig$beta_f1[Bal1SigCmp1]), abs(DataAllSig$beta_f2[Bal1SigCmp1]), method=methodr2)
corbeta2abs<-cor(abs(DataAllSig$beta_f1[Bal1SigCmp2]), abs(DataAllSig$beta_f2[Bal1SigCmp2]), method=methodr2)

Stat1<-data.frame(type=Head1, alpha=alpha, nbnofound=Nbbaselim1$Count[1], nbFound=Nbbaselim1$Count[2], r2Other=corbeta1, r2otherabs=corbeta1abs,PercSameSig=PercSamSig1, alpha2=alpha2,countr=NbVal1)
Stat2<-data.frame(type=Head2, alpha=alpha, nbnofound=Nbbaselim2$Count[1], nbFound=Nbbaselim2$Count[2], r2Other=corbeta2, r2otherabs=corbeta2abs, PercSameSig=PercSamSig2, alpha2=alpha2, countr=NbVal2)
Stat<-rbind(Stat1,Stat2)
## comparison of freq and beta

if(Cmt==1)StatAlpha<-Stat
else StatAlpha<-rbind(StatAlpha,Stat)
Cmt<-Cmt+1
}
}
StatAlpha$PercNotFound<-StatAlpha$nbnofound/(StatAlpha$nbnofound+StatAlpha$nbFound)*100
StatAlpha$Count<-StatAlpha$nbnofound+StatAlpha$nbFound
StatAlpha$logalpha<- - log10(StatAlpha$alpha)

commonalpha<-ggplot(data=unique(StatAlpha[,c('logalpha','PercNotFound','type', 'Count')]), aes(x=logalpha, y=PercNotFound,fill=type)) +  scale_fill_brewer(palette="Paired")+
  geom_bar(stat="identity",position=position_dodge())+geom_text(aes(label=Count), vjust=1.6, color="black", size=3.5)+
  labs(title="% not found in other data set", 
         x="-log10(alpha)", y = "Percentage not found")+
  theme_minimal()


##
StatAlpha2<-StatAlpha
StatAlpha$countr[StatAlpha$countr>100]<-""
StatAlpha2$alpha2<-as.character(StatAlpha2$alpha2)
percsamesig1<-ggplot(data=StatAlpha2[StatAlpha2$type==Head1 ,], aes(x=logalpha, y=PercSameSig,fill=alpha2)) +  scale_fill_brewer(palette="Paired")+
  geom_bar(stat="identity",position=position_dodge())+geom_text(aes(label=countr), vjust=1.6, color="black", size=3.5)+
  labs(title=paste("% beta same sign : ", Head1,sep=""),
         x="-log10(alpha)", y = "Percentage not found")+
  theme_minimal()

percsamesig2<-ggplot(data=StatAlpha2[StatAlpha2$type==Head2 ,], aes(x=logalpha, y=PercSameSig,fill=alpha2)) +  scale_fill_brewer(palette="Paired")+
  geom_bar(stat="identity",position=position_dodge())+geom_text(aes(label=countr), vjust=1.6, color="black", size=3.5)+
  labs(title=paste("% beta same sign : ",Head2,sep=''),
         x="-log10(alpha)", y = "% beta same sign")+
  theme_minimal()



tiff(paste(opt[['out']],'.tiff',sep=''))
grid.arrange(plotcommon, commonalpha, percsamesig1,percsamesig2,nrow=2, ncol=2)
dev.off()

write.csv(StatAlpha2,file=paste(opt[['out']],'.csv',sep=''), row.names=F)
write.table(DataAll, file=paste(opt[['out']],'.merge',sep=''), row.names=F)
