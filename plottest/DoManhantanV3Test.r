library(qqman)
library(data.table)
#https://genome.sph.umich.edu/wiki/Code_Sample:_Generating_QQ_Plots_in_R
source('QQTest.r')
source('manhattanModif.R')
args = commandArgs(trailingOnly=TRUE)
Trait=args[1]

DataInfoRs=read.csv("InfoRs2",header=T,stringsAsFactors=F)
listres=list()
Trait="og_friedewald_qc"
#for(Trait in c("og_friedewald_qc","og_hdl_qc","og_triglycerides_qc","og_cholesterol_1_qc")){
listsnp=as.character(DataInfoRs[!is.na(DataInfoRs[,Trait]) & DataInfoRs[,Trait]=="Y","Signals"])

File=grep(".imp.bolt",dir(paste("/dataE/AWIGenGWAS/shared/ResultGWAS/ImputedRes_20190108/Res/",Trait,sep=""),"_All_",full.names=T),value=T)
Data<-fread(File)
Data<-as.data.frame(Data[Data[['A1FREQ']]>0.01 & Data[['A1FREQ']]<0.99,])
#Data<-Data[Data$P_BOLT_LMM<10**-3,]
#SNP	CHR	BP	GENPOS	ALLELE1	ALLELE0	A1FREQ	INFO	BETA	SE	P_BOLT_LMM_INF	P_BOLT_LMM
highlight=list()
highlight[['red']]=as.character(DataInfoRs[ DataInfoRs$Comments=="Novel signal" & DataInfoRs[,Trait]=="Y" ,'Signals'])
highlight[['blue']]=as.character(DataInfoRs[ DataInfoRs$Comments=="Novel lead" & DataInfoRs[,Trait]=="Y" ,'Signals'])
highlight[['green3']]=as.character(DataInfoRs[ DataInfoRs$Comments=="Known SNP"  & DataInfoRs[,Trait]=="Y" ,'Signals'])
Max<-max(-log10(Data[,"P_BOLT_LMM"])*1.1)
DataInfoRsTrait<-DataInfoRs[ DataInfoRs[,Trait]=="Y" ,]
DataInfoRsTraitUni<-DataInfoRsTrait[c(DataInfoRsTrait[1:(nrow(DataInfoRsTrait)-1),"Loci"]!=DataInfoRsTrait[2:(nrow(DataInfoRsTrait)),"Loci"],T), ]
DataSNP<-Data[Data[,"SNP"] %in% DataInfoRsTraitUni[,"Signals"],]
merge(DataSNP,DataInfoRsTraitUni, by.x="SNP",by.y="Signals")
DataInfoRsTraitUni<-DataInfoRsTraitUni[,c("Signals","Loci")]
#DataInfoRsTraitUni$Loci[DataInfoRsTraitUni$Loci=="PCSK9"]<-"PCSK9 CELSR2"
#DataInfoRsTraitUni$Loci[DataInfoRsTraitUni$Loci=="CELSR2"]<-""
names(DataInfoRsTraitUni)<-c("SNP", "Genes")
DataInfoRsTraitUni$Genes[grep("intergenic",DataInfoRsTraitUni$Genes)]<-""
source("manhattanModif.R")
tiff(paste(Trait,"_All.man3.tiff",sep=""),width = 480*2, height = 480)
manhattan2(Data, chr = "CHR", bp = "BP", p = "P_BOLT_LMM", snp = "SNP", highlight=highlight, annotateSnp=listsnp, suggestiveline=-log10(5e-5), genes=DataInfoRsTraitUni)
dev.off()
#}
