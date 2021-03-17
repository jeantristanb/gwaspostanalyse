#!/usr/bin/env Rscript

estlambda<- function(data, plot=FALSE, proportion=1.0,
                        method="regression", filter=TRUE, df=1,... ) {
        data <- data[which(!is.na(data))]
        if (proportion>1.0 || proportion<=0)
                stop("proportion argument should be greater then zero and less than or equal to one")

        ntp <- round( proportion * length(data) )
        if ( ntp<1 ) stop("no valid measurements")
        if ( ntp==1 ) {
                warning(paste("One measurement, lambda = 1 returned"))
                return(list(estimate=1.0, se=999.99))
        }
        if ( ntp<10 ) warning(paste("number of points is too small:", ntp))
        if ( min(data)<0 ) stop("data argument has values <0")
        if ( max(data)<=1 ) {
#               lt16 <- (data < 1.e-16)
#               if (any(lt16)) {
#                       warning(paste("Some probabilities < 1.e-16; set to 1.e-16"))
#                       data[lt16] <- 1.e-16
#               }
                data <- qchisq(data, 1, lower.tail=FALSE)
        }
        if (filter)
        {
                data[which(abs(data)<1e-8)] <- NA
        }
        data <- sort(data)
        ppoi <- ppoints(data)
        ppoi <- sort(qchisq(ppoi, df=df, lower.tail=FALSE))
        out <- list()
        if (method=="regression") {
                s <- summary( lm(data~0+ppoi) )$coeff
                out$estimate <- s[1,1]
                out$se <- s[1,2]
        } else if (method=="median") {
                out$estimate <- median(data, na.rm=TRUE)/qchisq(0.5, df)
                out$se <- NA
        } else if (method=="KS") {
                limits <- c(0.5, 100)
                out$estimate <- estLambdaKS(data, limits=limits, df=df)
                if ( abs(out$estimate-limits[1])<1e-4 || abs(out$estimate-limits[2])<1e-4 )
                        warning("using method='KS' lambda too close to limits, use other method")
                out$se <- NA
        } else {
                stop("'method' should be either 'regression' or 'median'!")
        }

        if (plot) {
                lim <- c(0, max(data, ppoi,na.rm=TRUE))
#               plot(ppoi,data,xlim=lim,ylim=lim,xlab="Expected",ylab="Observed", ...)
                oldmargins <- par()$mar
                par(mar=oldmargins + 0.2)
                plot(ppoi, data,
                     xlab=expression("Expected " ~ chi^2),
                     ylab=expression("Observed " ~ chi^2),
                     ...)
                abline(a=0, b=1)
                abline(a=0, b=out$estimate, col="red")
                par(mar=oldmargins)
        }

        out
}

library(data.table)
library("optparse")
#https://genome.sph.umich.edu/wiki/Code_Sample:_Generating_QQ_Plots_in_R
getscriptdir<-function(){
initial.options <- commandArgs(trailingOnly = FALSE)
file.arg.name <- "--file="
script.name <- sub(file.arg.name, "", initial.options[grep(file.arg.name, initial.options)])
script.basename <- dirname(script.name)
return(script.basename)
}
FileSc<-'QQTest.r' 
if(!file.exists(FileSc)){
FileSc<-paste(getscriptdir() ,FileSc,sep='/')
}
source(FileSc)
option_list = list(
  make_option(c("--info"), type="character", default=NULL, 
           help="file contains pheno file by line", metavar="character"),
  make_option(c("--pval_head"), type="character", default=NULL, 
           help="file contains pheno file by line", metavar="character"),
  make_option(c("--freq_head"), type="character", default=NULL, 
           help="file contains pheno file by line", metavar="character"),
  make_option(c("--maf"), type="numeric", default=NULL, 
           help="file contains pheno file by line", metavar="character"),
  make_option(c("--test"), type="character", default="F", 
           help="file contains pheno file by line", metavar="character"),
  make_option(c("--metasoft"), type="character", default="F", 
           help="file contains pheno file by line", metavar="character"),
  make_option(c("--type_plot"), type="character", default="pdf", 
           help="file contains pheno file by line", metavar="character"),
  make_option(c("--out"), type="character", default="out.txt", 
           help="output file name [default= %default]", metavar="character")
); 


args = commandArgs(trailingOnly=TRUE)
balfreq=T
balisemetasoft=F
if(length(args)>0){
opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);
FileInfo=opt[["info"]]
headfreq=opt[['freq_head']]
headpval=opt[['pval_head']]
maf=opt[['maf']]
typeplot=opt[['type_plot']]
outplot=opt[['out']]
if(opt[['metasoft']]=="T")balisemetasoft=T
test=F
print(is.null(opt[["test"]]))
if(!is.null(opt[["test"]]) & opt[["test"]]=="T"){
test=T
}
}else{
FileInfo="/home/jeantristan/Travail/GWAS/lipid_gwas/smartpca/Analyse/QQPlot/VMars21/infoall.tsv"
headpval="P_BOLT_LMM"
headfreq="A1FREQ"
maf=0.01
typeplot="pdf"
test=T
outplot="test"
}
if(is.null(typeplot)){
cat('type plot does not exist')
q()
}
fctplot=get(typeplot)
listres=list()

InfoRed<-read.table(FileInfo, stringsAsFactors=F)
if(is.null(headpval)){
cat("pval header not found")
q(2)
}
if(is.null(maf) | is.null(headfreq)){
balfreq=F
}else{
if(maf>0.5 | maf< 0){
cat("maf > 0.5 or maf<0")
q(2)
}
}
Cmt<-1
listtrait2<-c()
for(Trait in InfoRed[,1]){
if(test)Data<-fread(InfoRed[Cmt,2], nrows=10000)
else Data<-fread(InfoRed[Cmt,2])
if(balfreq){
Data<-Data[!is.na(Data[[headpval]]) & !is.na(Data[[headfreq]]) & Data[[headfreq]]>maf & Data[[headfreq]]<(1-maf),]
}
if(balisemetasoft)names(Data)<-c(names(Data)[-1], 'None')
listres[[Trait]]=Data[[headpval]]
#lbdaa<-estlambda(listres[[Trait]][,headpval])
lbdaareg<-estlambda(listres[[Trait]])
lbdaamedi<-median(qchisq(listres[[Trait]], df=1, lower.tail=FALSE)) / qchisq(0.5, 1)
tmpdf<-data.frame(var=Trait,file=InfoRed[Cmt,2],lb_estreg=lbdaareg$estimate, ld_sereg=lbdaareg$se, lb_estmed=lbdaamedi)
if(Cmt==1)lamda<-tmpdf
else lamda<-rbind(lamda,tmpdf)
Cmt<-Cmt+1
Trait2<-paste(Trait, ' (',round(lbdaareg$estimate,2),')',sep='')
listtrait2<-c(listtrait2, Trait2)
}
#xlim=range(sapply(listres, function(x)range(x[,1], na.rm=T)))
#ylim=range(sapply(listres, function(x)range(x[,2], na.rm=T)))
fctplot(paste(outplot, typeplot,sep='.'))
names(listres)<-listtrait2
qqunif.plot(listres,  aspect = "fill")
dev.off()
write.csv(lamda, row.names=F, file=paste(outplot, 'csv',sep='.'))


