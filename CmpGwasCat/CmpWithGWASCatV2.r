library(data.table)

computedher<-function(beta, se, af,N){
#https://journals.plos.org/plosone/article/file?type=supplementary&id=info:doi/10.1371/journal.pone.0120758.s001
#https://www.netflix.com/watch/80205392?trackId=155573560
maf<-af
maf[!is.na(af) & af>0.5]<- 1 - maf[!is.na(af) & af>0.5]
ba<-!is.na(beta) & !is.na(se) & !is.na(maf) & !is.na(N)
a<-rep(NA, length(beta))
b<-rep(NA, length(beta))
a<-2*(beta[ba]**2)*(maf[ba]*(1-maf[ba]))
b<-2*(se[ba]**2)*N[ba]*maf[ba]*(1-maf[ba])
res<-rep(NA, length(beta))
res[ba]<-a[ba]/(a[ba]+b[ba])
res
}
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

transform_gwascat<-function(GwasCat){
NamesGWAS<-names(GwasCat)
GwasCat$risk.allele.cat<-sapply(strsplit(as.character(GwasCat$STRONGEST.SNP.RISK.ALLELE),split='-'), function(x)return(x[2]))
GwasCat$risk.allele.af<-as.numeric( GwasCat$RISK.ALLELE.FREQUENCY)
IC<-t(sapply(strsplit(gsub("[", "", sapply(strsplit(GwasCat$X95..CI..TEXT., split=']',fixed=T),function(x)x[1]),fixed=T), split="-"),function(x){
if(length(x)==2){return(c(as.numeric(x[1]),as.numeric(x[2])))
}else{
return(c(-as.numeric(x[2]),as.numeric(x[3])))
}
}))
IC<-data.frame(lower.cat=IC[,1], upper.cat=IC[,2])
GwasCat<-cbind(GwasCat,IC)
GwasCat$beta.cat<-GwasCat$OR.or.BETA
GwasCatF<-GwasCat
GwasCatF$nsample.cat<-sapply(strsplit(GwasCatF$INITIAL.SAMPLE.SIZE,split="[ ]"),function(x)sum(as.integer(gsub(",", "",grep("[0-9]", x,value=T)))))
GwasCatF$sd.cat<-(GwasCatF$upper.cat - GwasCatF$beta.cat)/1.96
GwasCatF$sd.cat2<-(GwasCatF$upper.cat - GwasCatF$beta.cat)/1.96*sqrt(GwasCatF$nsample.cat)
GwasCatF$Z.gwascat<-GwasCatF$beta.cat/(GwasCatF$sd.cat)
GwasCatF$h2.cat=computedher(GwasCatF$beta.cat, GwasCatF$sd.cat, GwasCatF$risk.allele.af,GwasCatF$nsample.cat)

GwasCatF$Qc<-T
GwasCatF$Qc<-GwasCatF$Qc & apply(GwasCatF[,c('beta.cat','sd.cat','risk.allele.af', 'nsample.cat','lower.cat','upper.cat')], 1,function(x)all(!is.na(x)))
GwasCatF$Qc<-GwasCatF$Qc & GwasCat$risk.allele.af<1 & GwasCat$risk.allele.af>0 & GwasCat$beta.cat>GwasCat$lower.cat & GwasCat$beta.cat<GwasCat$upper.cat
GwasCatF[,c(NamesGWAS, 'risk.allele.cat','beta.cat','lower.cat','upper.cat','nsample.cat','sd.cat','Z.gwascat','risk.allele.af', 'Qc','h2.cat')]
}

GetDataWithGC<-function(GwasCat, File, chrhead, bphead,betahead,sehead,a1head, a2head,afhead,N, chrheadcat="Chro37", bpheadcat="PosBegin37", head=""){
require(data.table)
DataGWAS<-fread(File)
DataGWAS[["ChrPos"]]<-paste(DataGWAS[[chrhead]],DataGWAS[[bphead]])
GwasCat[,'ChrPos']<-paste(GwasCat[,"Chro37"], GwasCat[,"PosBegin37"])
DataGWAS<-as.data.frame(DataGWAS[DataGWAS[['ChrPos']] %in% as.character(GwasCat[,'ChrPos'])])
GwasCatData<-merge(DataGWAS,GwasCat, by="ChrPos", all=T)
GwasCatData$Qc<-GwasCatData$Qc & !is.na(GwasCatData[,a1head]) & (GwasCatData[,a1head]==GwasCatData$risk.allele.cat | GwasCatData[,a2head]==GwasCatData$risk.allele.cat)
BaliseChange<-GwasCatData$Qc & GwasCatData[,a1head]!=GwasCatData$risk.allele.cat
GwasCatData$risk.allele.af.old<-GwasCatData$risk.allele.af
GwasCatData$risk.allele.af[BaliseChange]<- 1 - GwasCatData$risk.allele.af.old[BaliseChange]
GwasCatData$risk.allele[BaliseChange]<-GwasCatData$allele0[BaliseChange]
GwasCatData$beta.cat.old<-GwasCatData$beta.cat
GwasCatData$beta.cat[BaliseChange]<- -GwasCatData$beta.cat.old[BaliseChange]
GwasCatData$Z.gwascat.old<- GwasCatData$Z.gwascat
GwasCatData$Z.gwascat[BaliseChange]<- -GwasCatData$Z.gwascat.old[BaliseChange]
GwasCatData$upper.cat.old<- GwasCatData$upper.cat
GwasCatData$lower.cat.old<- GwasCatData$lower.cat
GwasCatData$lower.cat[BaliseChange]<- -GwasCatData$upper.cat.old[BaliseChange]
GwasCatData$upper.cat[BaliseChange]<- -GwasCatData$lower.cat.old[BaliseChange]
GwasCatData$h2<-computedher(GwasCatData[,betahead], GwasCatData[,sehead], GwasCatData[,afhead],rep(N, nrow(GwasCatData)))
GwasCatData$Z<-GwasCatData[,betahead]/GwasCatData[,sehead]
#names(DataGWAS)[!(names(DataGWAS) %in% ChrPos)]<-paste(names(DataGWAS)[!(names(DataGWAS) %in% ChrPos)])
names(GwasCatData)[names(GwasCatData) %in% c(names(DataGWAS),'h2','Z')]<-paste(names(GwasCatData)[names(GwasCatData) %in% c(names(DataGWAS),'h2','Z')], head,sep='')
return(GwasCatData[,grep(".old",names(GwasCatData),invert=T,value=T)])
}

plotfreq<-function(dataall,freq1, freq2,xlab='Awigen frequencies', ylab='GWAS cat frequencies'){
nbcat=20
height=1
perctrans<-90
plot(c(0,1), c(0,1), type='n', cex=0.5, xlab=xlab, ylab='Gwas catalog frequencies', xlim=c(0,1), ylim=c(0,1))
ht<-hist(dataall[,freq1], nbcat, plot=F)
ht$counts<-ht$counts/sum(ht$counts)*height
ht2<-hist(dataall[,freq2], nbcat, plot=F)
ht2$counts<-ht2$counts/sum(ht2$counts)*height
rect(rep(0,length(ht2$counts)),(1:length(ht2$counts))/nbcat,
           ht2$counts,(1:length(ht2$counts))/nbcat+1/nbcat, col=t_col('red', perctrans), bg=t_col('red', perctrans))
plot(ht, add=T, col=t_col("orange"))
points(dataall[,freq1], dataall[,freq2], pch=22, bg=t_col("blue") ,col=t_col("blue"), cex=0.5)
abline(a=0,b=1, lty=2 ,col=t_col('red'), lwd=2)
}

plotZ<-function(dataall,Z1, Z2, xlab='Z (GWAS cat)', ylab='Z (AWIGEN)'){
r2<-cor(dataall[,Z1],dataall[,Z2], method='spearman')
r2abs<-cor(abs(dataall[,Z1]), abs(dataall[,Z2]), method='spearman')
plot(dataall[,Z1], dataall[,Z2], pch=22, cex=0.5,bg=t_col("blue") ,col=t_col("blue"), xlab=xlab, ylab=ylab)
text(min(dataall[,Z1])-min(dataall[,Z1])*0.1,max(dataall[,Z2]), paste(paste("r2 :",round(r2,2)), paste("\n         r2 (abs) :",round(r2abs,2)), sep=""))
abline(h=0, col='red', lty=2)
abline(v=0, col='red', lty=2)
}



GwasCatI<-read.table('~/Data/GWASCat/GWASCat_27_March_2020_ckd.tsv', header=T, sep='\t', stringsAsFactors=F)
GwasCat<-transform_gwascat(GwasCatI)
GwasCat$PosCat<-1:nrow(GwasCat)

Pheno="MDRDNoEth";ListMappedTrait<-grep('glomerular filtration rate', GwasCat$MAPPED_TRAIT,value=T)
Pheno="CKDEPINoEth";ListMappedTrait<-grep('glomerular filtration rate', GwasCat$MAPPED_TRAIT,value=T)
chrheadcat="Chro37";bpheadcat="PosBegin37"
chrhead<-'chr';bphead='ps';betahead<-'beta';sehead='se';N=11000;a1head<-"allele1";a2head="allele0";afhead<-'af'
File<-paste('/home/jeantristan/Travail/GWAS/GWAS_CKD/ImputedDataV4/Result/agesexhivdmhtnpcabmi/Res/',Pheno,'/',Pheno,'_All_agesexhivdmhtnpcabmi_20190120.imp.stat',sep="")
DataCKDEPI<-GetDataWithGC(GwasCat, File, chrhead, bphead,betahead,sehead,a1head, a2head,afhead,N, chrheadcat="Chro37", bpheadcat="PosBegin37", head='.ckdepi')
svg(paste('difffrequencies_gwasawigen_',Pheno,'.svg',sep=''))
plotfreq(DataCKDEPI[!is.na(DataCKDEPI$af.ckdepi) & DataCKDEPI$Qc & (DataCKDEPI$MAPPED_TRAIT %in% ListMappedTrait),],'af.ckdepi', 'risk.allele.af',xlab='Awigen frequencies', ylab='GWAS cat frequencies')
dev.off()
svg(paste('diffZ_gwasawigen_',Pheno,'.svg',sep=''))
plotZ(DataCKDEPI[!is.na(DataCKDEPI$af.ckdepi) & DataCKDEPI$Qc & (DataCKDEPI$MAPPED_TRAIT %in% ListMappedTrait),],"Z.gwascat", "Z.ckdepi", xlab='Z (GWAS cat)', ylab='Z (AWIGEN)')
dev.off()

svg(paste('diffh2_gwasawigen_',Pheno,'.svg',sep=''))
plotZ(DataCKDEPI[!is.na(DataCKDEPI$af.ckdepi) & DataCKDEPI$Qc & (DataCKDEPI$MAPPED_TRAIT %in% ListMappedTrait),],"h2.cat", "h2.ckdepi", xlab='h2 (GWAS cat)', ylab='h2 (AWIGEN)')
dev.off()

Pheno="LogNewAcr";ListMappedTraitUACR<-unique(grep('albumin', GwasCat$MAPPED_TRAIT,value=T));ListMappedTraitUACR<-ListMappedTrait[ListMappedTrait!='serum non-albumin protein measurement']
File<-paste('/home/jeantristan/Travail/GWAS/GWAS_CKD/ImputedDataV4/Result/agesexhivdmhtnpcabmi/Res/',Pheno,'/',Pheno,'_All_agesexhivdmhtnpcabmi_20190120.imp.stat',sep="")
DataUacr<-GetDataWithGC(GwasCat, File, chrhead, bphead,betahead,sehead,a1head, a2head,afhead,N, chrheadcat="Chro37", bpheadcat="PosBegin37",head='.uacr')


svg(paste('difffrequencies_gwasawigen_',Pheno,'.svg',sep=''))
plotfreq(DataUacr[DataUacr$Qc & (DataUacr$MAPPED_TRAIT %in% ListMappedTraitUACR),],'af.uacr', 'risk.allele.af',xlab='Awigen frequencies', ylab='GWAS cat frequencies')
dev.off()

svg(paste('diffh2_gwasawigen_',Pheno,'.svg',sep=''))
plotZ(DataUacr[!is.na(DataUacr$af.uacr) & DataUacr$Qc & (DataUacr$MAPPED_TRAIT %in% ListMappedTraitUACR),],"h2.cat", "h2.uacr", xlab='h2 (GWAS cat)', ylab='h2 (AWIGEN)')
dev.off()

svg(paste('diffZ_gwasawigen_',Pheno,'.svg',sep=''))
plotZ(DataUacr[!is.na(DataUacr$af.uacr) & DataUacr$Qc & (DataUacr$MAPPED_TRAIT %in% ListMappedTraitUACR),],"Z.gwascat", "Z.uacr", xlab='Z (GWAS cat)', ylab='Z (AWIGEN)')
dev.off()


All<-merge(DataCKDEPI, DataUacr[,c('PosCat', grep('uacr',names(DataUacr), value=T))], by='PosCat')
All$Use<-'NoPASS'
All$Use[!is.na(All$af.uacr) & All$Qc]<-'Other'
All$Use[!is.na(All$af.uacr) & All$MAPPED_TRAIT %in% ListMappedTraitUACR]<-'UACR'
All$Use[!is.na(All$af.uacr) & All$MAPPED_TRAIT %in% ListMappedTrait]<-'eGFR'
write.csv(All, row.names=F, file='resumeall.csv')
