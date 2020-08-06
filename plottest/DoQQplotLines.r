library(data.table)
#library(qqman)
#https://genome.sph.umich.edu/wiki/Code_Sample:_Generating_QQ_Plots_in_R
#source('QQTest.r')


listres=list()
for(Trait in c("og_friedewald_qc","og_hdl_qc","og_triglycerides_qc","og_cholesterol_1_qc")){
File=grep(".imp.bolt",dir(paste("/dataE/AWIGenGWAS/shared/ResultGWAS/ImputedRes_20190108/Res/",Trait,sep=""),"_All_",full.names=T),value=T)
Data<-fread(File)
Data<-Data[Data[['A1FREQ']]>0.01 & Data[['A1FREQ']]<0.99,]
pvector<-Data[['P_BOLT_LMM']]
listres[[Trait]]=data.frame(o=-log10(sort(pvector, decreasing = FALSE)), e=-log10(ppoints(length(pvector))))
#listres[[Trait]]=Data[['P_BOLT_LMM']]
}
ylab = expression(Observed ~ ~-log[10](italic(p)))
xlab = expression(Expected ~ ~-log[10](italic(p)))
#xlim=range(sapply(listres, function(x)range(x[,1], na.rm=T)))
#ylim=range(sapply(listres, function(x)range(x[,2], na.rm=T)))
names(listres)<-c("ldl","hdl","triglycerides", "cholesterol")
tiff("CmpQQPlotLines.tiff")
RangeXLim<-range(sapply(listres,function(x)return(range(x[,2]))))
RangeYLim<-range(sapply(listres,function(x)return(range(x[,1]))))
plot(RangeXLim, RangeYLim, type="n", xlab=xlab, ylab=ylab, bty="n")
Cmt<-1
for(Col in names(listres)){
lines(listres[[Col]][,2], listres[[Col]][,1], col=Cmt, lwd=2)
Cmt<-Cmt+1
}
abline(0, 1, col = "red", lwd=3, lty=2)
legend("topleft", legend=names(listres), col=1:Cmt, lty=1, lwd=2, bty="n")
dev.off()


