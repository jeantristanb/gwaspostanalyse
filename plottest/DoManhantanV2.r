library(qqman)
library(data.table)
#https://genome.sph.umich.edu/wiki/Code_Sample:_Generating_QQ_Plots_in_R
source('/home/jeantristan/Travail/GWAS/AllVar/pca_prune_50_10_0.1/Analyse/PlotForAbstract/QQTest.r')
source('/home/jeantristan/Travail/GWAS/AllVar/pca_prune_50_10_0.1/Analyse/PlotForAbstract/manhattanModif.R')
args = commandArgs(trailingOnly=TRUE)
Trait=args[1]

DataInfoRs=read.table("InfoRs",header=T)
listres=list()
#Trait="og_friedewald_qc"
#for(Trait in c("og_friedewald_qc","og_hdl_qc","og_triglycerides_qc","og_cholesterol_1_qc")){
listsnp=as.character(DataInfoRs[!is.na(DataInfoRs[,Trait]) & DataInfoRs[,Trait]=="Y","Signals"])

File=grep(".imp.bolt",dir(paste("/dataE/AWIGenGWAS/shared/ResultGWAS/ImputedRes_20190108/Res/",Trait,sep=""),"_All_",full.names=T),value=T)
Data<-fread(File)
Data<-as.data.frame(Data[Data[['A1FREQ']]>0.01 & Data[['A1FREQ']]<0.99,])
Data<-Data[Data$P_BOLT_LMM<10**-3,]
#SNP	CHR	BP	GENPOS	ALLELE1	ALLELE0	A1FREQ	INFO	BETA	SE	P_BOLT_LMM_INF	P_BOLT_LMM
highlight=list()
highlight[['red']]=as.character(DataInfoRs[ DataInfoRs$InfoFreq=="African" & !is.na(DataInfoRs[,Trait]) & DataInfoRs[,Trait]=="Y" ,'Signals'])
highlight[['green3']]=as.character(DataInfoRs[ DataInfoRs$InfoFreq!="African" & !is.na(DataInfoRs[,Trait]) & DataInfoRs[,Trait]=="Y" ,'Signals'])
tiff(paste(Trait,"_All.man2.tiff",sep=""),width = 480*2, height = 480)
manhattan2(Data, chr = "CHR", bp = "BP", p = "P_BOLT_LMM", snp = "SNP", highlight=highlight, annotateSnp=listsnp)
dev.off()
#}
