library(data.table)
library(graphics)
library("optparse")

computedher<-function(beta, se, af,N){
#https://journals.plos.org/plosone/article/file?type=supplementary&id=info:doi/10.1371/journal.pone.0120758.s001
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

transform_gwascat<-function(GwasCat, chrheadcat, bpheadcat){
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
GwasCatF$Z.cat<-GwasCatF$beta.cat/(GwasCatF$sd.cat)
GwasCatF$h2.cat=computedher(GwasCatF$beta.cat, GwasCatF$sd.cat, GwasCatF$risk.allele.af,GwasCatF$nsample.cat)

GwasCatF$Qc<-T
GwasCatF$Qc<-GwasCatF$Qc & apply(GwasCatF[,c('beta.cat','sd.cat','risk.allele.af', 'nsample.cat','lower.cat','upper.cat')], 1,function(x)all(!is.na(x)))
GwasCatF$Qc<-GwasCatF$Qc & GwasCat$risk.allele.af<1 & GwasCat$risk.allele.af>0 & GwasCat$beta.cat>GwasCat$lower.cat & GwasCat$beta.cat<GwasCat$upper.cat
GwasCatF[,c(NamesGWAS, 'risk.allele.cat','beta.cat','lower.cat','upper.cat','nsample.cat','sd.cat','Z.cat','risk.allele.af', 'Qc','h2.cat')]
}

GetDataWithGC<-function(GwasCat, File, chrhead, bphead,betahead,sehead,a1head, a2head,afhead,N, chrheadcat="Chro37", bpheadcat="PosBegin37", head=""){
require(data.table)
DataGWAS<-fread(File)
DataGWAS[["ChrPos"]]<-paste(DataGWAS[[chrhead]],DataGWAS[[bphead]])
GwasCat[,'ChrPos']<-paste(GwasCat[,chrheadcat], GwasCat[,bpheadcat])
DataGWAS<-as.data.frame(DataGWAS[DataGWAS[['ChrPos']] %in% as.character(GwasCat[,'ChrPos'])])
GwasCatData<-merge(DataGWAS,GwasCat, by="ChrPos", all=T)
GwasCatData$Qc<-GwasCatData$Qc & !is.na(GwasCatData[,a1head]) & (GwasCatData[,a1head]==GwasCatData$risk.allele.cat | GwasCatData[,a2head]==GwasCatData$risk.allele.cat)
BaliseChange<-GwasCatData$Qc & GwasCatData[,a1head]!=GwasCatData$risk.allele.cat
GwasCatData$risk.allele.af.old<-GwasCatData$risk.allele.af
GwasCatData$risk.allele.af[BaliseChange]<- 1 - GwasCatData$risk.allele.af.old[BaliseChange]
GwasCatData$risk.allele[BaliseChange]<-GwasCatData$allele0[BaliseChange]
GwasCatData$beta.cat.old<-GwasCatData$beta.cat
GwasCatData$beta.cat[BaliseChange]<- -GwasCatData$beta.cat.old[BaliseChange]
GwasCatData$Z.cat.old<- GwasCatData$Z.cat
GwasCatData$Z.cat[BaliseChange]<- -GwasCatData$Z.cat.old[BaliseChange]
GwasCatData$upper.cat.old<- GwasCatData$upper.cat
GwasCatData$lower.cat.old<- GwasCatData$lower.cat
GwasCatData$lower.cat[BaliseChange]<- -GwasCatData$upper.cat.old[BaliseChange]
GwasCatData$upper.cat[BaliseChange]<- -GwasCatData$lower.cat.old[BaliseChange]
if(!is.null(betahead) & !is.null(sehead) & !is.null(afhead)){
GwasCatData$h2<-computedher(GwasCatData[,betahead], GwasCatData[,sehead], GwasCatData[,afhead],rep(N, nrow(GwasCatData)))
GwasCatData$Z<-GwasCatData[,betahead]/GwasCatData[,sehead]
}else{
GwasCatData$h2<-NA
GwasCatData$Z<-NA
}
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
#points(dataall[,freq1], dataall[,freq2], pch=22, bg=t_col("blue") ,col=t_col("blue"), cex=0.5)
points(dataall[,freq1], dataall[,freq2], pch=22, ,col=t_col("blue",90), cex=0.2)
rect(rep(0,length(ht2$counts)),(1:length(ht2$counts))/nbcat,
           ht2$counts,(1:length(ht2$counts))/nbcat+1/nbcat, col=t_col('red', perctrans), bg=t_col('red', perctrans))
plot(ht, add=T, col=t_col("orange"))
abline(a=0,b=1, lty=2 ,col=t_col('red'), lwd=2)
}

#dataall<-DataInfo[!is.na(DataInfo[,headfreq]) & DataInfo$Qc,];freq1<-headfreq;freq2<-'risk.allele.af';xlab="";ylab=""

plotfreq<-function(dataall,freq1, freq2,cex_pt,alpha_pt,xlab='Awigen frequencies', ylab='GWAS cat frequencies'){
nbcat=20
height=1
perctrans<-90
ht<-hist(dataall[,freq1], nbcat, plot=F)
ht$counts<-ht$counts/sum(ht$counts)*height
ht2<-hist(dataall[,freq2], nbcat, plot=F)
ht2$counts<-ht2$counts/sum(ht2$counts)*height
layout.matrix <- matrix(c(2, 1, 0, 3), nrow = 2, ncol = 2)

layout(mat = layout.matrix,
       heights = c(max(ht$counts)*1.2, 1), # Heights of the two rows
       widths = c(2, max(ht2$counts)*2.4)) # Widths of the two columns
#layout.show(3)

#c(bottom, left, top, right)
par(mar=c(5, 4, 1, 1))
plot(c(0,1), c(0,1), type='n', cex=0.5, xlab=xlab, ylab=ylab, xlim=c(0,1), ylim=c(0,1), bty='n')
points(dataall[,freq1], dataall[,freq2], pch=22, ,col=t_col("blue",95), bg=t_col("blue",95),cex=cex_pt)
abline(a=0,b=1, lty=2 ,col=t_col('red'), lwd=2)
par(mar=c(0,4,1,0))
plot(ht, col=t_col("orange"), xlab="", ylab="", xlim=c(0,1),  xaxt='n',main="")
par(mar=c(5,0,1,1))
plot(c(0, max(ht2$counts)), c(0,1), type='n', cex=0.5, xlab='', ylab='', xlim=c(0,max(ht2$counts)), ylim=c(0,1), yaxt='n',bty='n')
rect(rep(0,length(ht2$counts)),(1:length(ht2$counts))/nbcat,
           ht2$counts,(1:length(ht2$counts))/nbcat-1/nbcat, col=t_col('red', perctrans), bg=t_col('red', perctrans))
}



plotZ<-function(dataall,Z1, Z2, xlab='Z (GWAS cat)', ylab='Z (AWIGEN)'){
r2<-cor(dataall[,Z1],dataall[,Z2], method='spearman')
r2abs<-cor(abs(dataall[,Z1]), abs(dataall[,Z2]), method='spearman')
plot(dataall[,Z1], dataall[,Z2], pch=22, cex=0.5,bg=t_col("blue") ,col=t_col("blue"), xlab=xlab, ylab=ylab)
text(min(dataall[,Z1])-min(dataall[,Z1])*0.1,max(dataall[,Z2]), paste(paste("r2 :",round(r2,2)), paste("\n         r2 (abs) :",round(r2abs,2)), sep=""))
abline(h=0, col='red', lty=2)
abline(v=0, col='red', lty=2)
}

option_list = list(
  make_option(c("--gwascat"), type="character", default=NULL, 
              help="gwas catalog input", metavar="character"),
  make_option(c("--gc_chr"), type="character", default=NULL, 
              help="gwas catalog input", metavar="character"),
  make_option(c("--in_bp"), type="character", default=NULL, 
              help="gwas catalog input", metavar="character"),
  make_option(c("--in_chr"), type="character", default=NULL, 
              help="gwas catalog input", metavar="character"),
  make_option(c("--in_af"), type="character", default=NULL, 
              help="gwas catalog input", metavar="character"),
  make_option(c("--in_beta"), type="character", default=NULL, 
              help="gwas catalog input", metavar="character"),
  make_option(c("--in_se"), type="character", default=NULL, 
              help="gwas catalog input", metavar="character"),
  make_option(c("--in_a1"), type="character", default=NULL, 
              help="gwas catalog input", metavar="character"),
  make_option(c("--gwas"), type="character", default=NULL, 
              help="gwas catalog input", metavar="character"),
  make_option(c("--in_a2"), type="character", default=NULL, 
              help="gwas catalog input", metavar="character"),
  make_option(c("--gc_bp"), type="character", default=NULL, 
              help="gwas catalog input", metavar="character"),
  make_option(c("--cex_pt"), type="numeric", default=0.15,
              help="gwas catalog input", metavar="character"),
  make_option(c("--alpha_pt"), type="numeric", default=90,
              help="gwas catalog input", metavar="character"),
  make_option(c("--n"), type="numeric", default=11000,
              help="gwas catalog input", metavar="numeric"),
  make_option(c("--head"), type="character", default="", 
              help="head", metavar="character"),
  make_option(c("--out"), type="character", default="out", 
              help="output file name [default= %default]", metavar="character")
); 

checkhead<-function(x, head){
if(is.null(x)){
cat(head, ' : null\n exit')
q()
}
return(x)
}


args = commandArgs(trailingOnly=TRUE)
if(length(args)==0){
#GwasCatI<-read.table('gwascat/GWASCat_27_March_2020_diab.tsv', header=T, sep='\t', stringsAsFactors=F);File='Afr_GC_1000G_form.freq';chrheadcat="Chro37";bpheadcat="PosBegin37";chrhead<-'CHR';bphead='BP';betahead<-NULL;sehead=NULL;N=NULL;a1head<-"A1";a2head="A2";afhead<-'MAF'
headout=""
cex_pt=0.15
alpha_pt=0.9
Out='test'
GwasCatI<-read.table("../gwascatalog/egfr_pheno.tsv", header=T, sep='\t', stringsAsFactors=F);chrheadcat="ChroNewRs";bpheadcat="PosNewRs"
chrhead='chr';bphead='ps';betahead='beta';sehead='se';N=11000;a1head="allele1";a2head="allele0";afhead='af'
File="/home/jeantristan/Travail/GWAS/GWAS_CKD/ImputedDataV4/Result/agesexhivdmhtnpcabmi/Res/LogNewAcr/LogNewAcr_All_agesexhivdmhtnpcabmi_20190120.imp.stat"



}else{
opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);
GwasCatI<-read.table(checkhead(opt[['gwascat']], 'gwascat'), header=T, sep='\t', stringsAsFactors=F)
File<-opt[['gwas']];chrheadcat=checkhead(opt[['gc_chr']], 'gc_chr');bpheadcat=checkhead(opt[['gc_bp']], 'gc_bp')
chrhead<-opt[['in_chr']];bphead=opt[['in_bp']];betahead<-opt[['in_beta']];sehead=opt[['in_se']];N=opt[['n']];a1head<-opt[['in_a1']];a2head=opt[['in_a2']];afhead<-opt[['in_af']];headout=checkhead(opt[['head']], 'head')
alpha_pt=checkhead(opt[['alpha_pt']], 'alpha_pt')
Out=checkhead(opt[['out']], 'out')
cex_pt=checkhead(opt[['cex_pt']], 'cex_pt')
}
if(!is.null(betahead) & !is.null(sehead) & !is.null(afhead)){
checkhead(betahead, 'in_beta');checkhead(sehead, 'in_se');checkhead(afhead, 'in_af');checkhead(a1head, 'in_a1');checkhead(a2head, 'in_a2');checkhead(bphead, 'in_bp');checkhead(chrhead, 'in_chr')
baliseZ<-T
}else{
baliseZ<-F
}
if(any(!(c(bpheadcat, chrheadcat) %in% names(GwasCatI)))){
paste(chrheadcat, 'or ', chrheadcat , ' not found in ', paste(names(GwasCatI),collapse=','))
q()
}
GwasCat<-transform_gwascat(GwasCatI, chrheadcat,bpheadcat)
GwasCat$PosCat<-1:nrow(GwasCat)

## data 
DataInfo<-GetDataWithGC(GwasCat,File, chrhead, bphead,betahead,sehead,a1head, a2head,afhead,N,chrheadcat= chrheadcat, bpheadcat=bpheadcat, head="")
svg(paste(Out,'_cmpfreq.svg',sep=''))
plotfreq(DataInfo[!is.na(DataInfo[,afhead]) & DataInfo$Qc,],afhead, 'risk.allele.af',cex_pt=cex_pt,alpha_pt=alpha_pt,xlab=headout, ylab='GWAS Catalog')
dev.off()

## plot of 
#c(bottom, left, top, right)’
#the four sides of the plot.  The default is ‘c(5, 4, 4, 2) +
plothtfreq<-function(dataall, freqhead, col, ylim=NA, plot=T){
nbcat=20
height<-1
perctrans<-90
ht<-hist(dataall[,freqhead], nbcat, plot=F)
ht$counts<-ht$counts/sum(ht$counts)*height *100
par(mar=c(4,4,1,0))
if(plot)plot(ht, col=t_col(col), xlab=xlab, ylab=ylab, xlim=c(0,1),  main="", ylim=ylim)
return(ht)
}
freqhead<-afhead;dataall<-DataInfo[!is.na(DataInfo[,afhead]) & DataInfo$Qc,];ylab="% of SNPs";xlab="Frequency";col="orange"
htaf<-plothtfreq(dataall, freqhead, col, NA, F)
freqhead<-'risk.allele.af';dataall<-DataInfo[!is.na(DataInfo[,afhead]) & DataInfo$Qc,];ylab="% of SNPs";xlab="Frequency";col="red"
htcat<-plothtfreq(dataall, freqhead, col, NA, F)
ylim<-range(htaf$counts, htcat$counts)
ylim[1]<-0

freqhead<-afhead;dataall<-DataInfo[!is.na(DataInfo[,afhead]) & DataInfo$Qc,];ylab="% of SNPs";xlab="Frequency";col="orange"
svg(paste(Out,'_histfreq_',afhead,'.svg',sep=''),width = 7*1.5, height = 7)
plothtfreq(dataall, freqhead, col, ylim, T)
dev.off()

freqhead<-'risk.allele.af';dataall<-DataInfo[!is.na(DataInfo[,afhead]) & DataInfo$Qc,];ylab="% of SNPs";xlab="Frequency";col="red"
svg(paste(Out,'_histfreq_gwascat.svg',sep=''), width = 7*1.5, height = 7)
plothtfreq(dataall, freqhead, col, ylim, T)
dev.off()



names(DataInfo)
table(is.na(DataInfo$h2), is.na(DataInfo$h2.cat))
table(is.na(DataInfo$Z), is.na(DataInfo$Z.cat))
if(baliseZ){
svg(paste(Out,'_cmpZ.svg',sep=''))
plotZ(DataInfo[!is.na(DataInfo$af) & DataInfo$Qc ,],"Z", "Z.cat", xlab=headout, ylab='GWAS Catalog')
dev.off()

svg(paste(Out,'_cmph2.svg',sep=''))
plotZ(DataInfo[!is.na(DataInfo$af) & DataInfo$Qc,],"h2", "h2.cat", xlab=headout, ylab='GWAS Catalog')
dev.off()
}
write.csv(DataInfo, row.names=F, quote=T, file=paste(Out,"_resume.csv",sep=''))

