#!/usr/bin/Rscript
PlotLambda<-function(MatOr, Head=1, HeadVal="lambda",HeadLower="lower", HeadUpper="upper", xlab=HeadVal){
library(RColorBrewer)
YLim<-c(-1,nrow(MatOr)+length(unique(MatOr[,Head])))
XLim<-range(MatOr[,c(HeadLower, HeadUpper)] ,na.rm=T);XLim[2]<-XLim[2]*1.15
ListePch<-21:25
qual_col_pals = brewer.pal.info[brewer.pal.info$category == 'qual',]
col_vector = unlist(mapply(brewer.pal, qual_col_pals$maxcolors, rownames(qual_col_pals)))
par(mar=c(4, 6, 1, 1) +0.1)
plot(XLim,YLim,xlim=XLim, ylim=YLim,type='n', bty="n", yaxt='n', main='', xlab=xlab,ylab="")
abline(v=1, lty=2,lwd=2)
Cmt2<-1;Cmt<-1
for(Lab in unique(MatOr[,Head])){
MatOrLab<-MatOr[MatOr[,Head]==Lab,];NbRow<-nrow(MatOrLab)
sapply(1:NbRow,function(x)lines(c(MatOrLab[x,HeadLower],MatOrLab[x,HeadUpper]), rep(Cmt2+x-1,2), col=col_vector[Cmt], lty=3, lwd=1.5))
points(MatOrLab[,HeadVal] ,Cmt2:(Cmt2+NbRow-1), pch=ListePch[1:NbRow], col=col_vector[Cmt],bg=col_vector[Cmt], pty=2)
mtext(Lab, side = 2,  at =Cmt2+0.5, cex=0.75, las=2)
Cmt2<-Cmt2+NbRow+1;Cmt<-Cmt+1
}
}

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
#		lt16 <- (data < 1.e-16)
#		if (any(lt16)) {
#			warning(paste("Some probabilities < 1.e-16; set to 1.e-16"))
#			data[lt16] <- 1.e-16
#		}
                data <- qchisq(data, 1, lower.tail=FALSE)
        }
        if (filter)
        {
                data[which(abs(data)<1e-8)] <- NA
        }
        data <- sort(data)
        ppoi <- ppoints(data)
        ppoi <- sort(qchisq(ppoi, df=df, lower.tail=FALSE))
        data <- data[1:ntp]
        ppoi <- ppoi[1:ntp]
#	s <- summary(lm(data~offset(ppoi)))$coeff
#       bug fix thanks to Franz Quehenberger

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
#		plot(ppoi,data,xlim=lim,ylim=lim,xlab="Expected",ylab="Observed", ...)
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
suppressMessages(require(data.table))
library(RColorBrewer)
library("optparse")

option_list = list(
  make_option(c("--listfile"), type="character",
              help="file name contains all file with result column and header: Traits Type Site File FilePheno ", metavar="character"),
  make_option(c("--type"), type="character",
              help="file name contains all file with result column and header: Traits Type Site File FilePheno ", metavar="character"),
  make_option(c("--out"), type="character", default=0.01,
              help="output file name [default= %default]", metavar="character")

)
opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);

listfilename=strsplit(opt[['listfile']],split=',')[[1]]
listfile=list()
listlamb=list()
CmtF=1
for(File in listfilename){
head=gsub("\\.qq$","",basename(File))
tmp<-as.data.frame(fread(File, header=F))
tmp<-na.omit(tmp)
tmp<-tmp[!is.infinite(tmp[,1]) & !is.infinite(tmp[,2]),]
listfile[[head]]=tmp
lam=estlambda(10**(-tmp[,2]))
lam2<-data.frame(head=head, lambda=lam$estimate ,se=lam$se)
if(CmtF==1)DataFra<-lam2
else DataFra<-rbind(DataFra,lam2)
CmtF<-CmtF+1
}
listhead=names(listfile)
DataFra$lower=DataFra$lambda - 1.96 * DataFra$se
DataFra$upper=DataFra$lambda + 1.96 * DataFra$se

write.table(DataFra,row.names=F, col.names=T, file=paste(opt[['out']],'.lamda' ,sep=''),quote=F, sep='\t')



xlim=range(unlist(sapply(listhead, function(x)range(listfile[[x]][,1], na.rm=T))))
ylim=range(unlist(sapply(listhead, function(x)range(listfile[[x]][,2], na.rm=T))))

color = grDevices::colors()[grep('gr(a|e)y', grDevices::colors(), invert = T)]
col=sample(color, length(listhead))
fctplot=eval(parse(text = opt[['type']]))
out=paste(opt[['out']],  opt[['type']],sep='.')
fctplot(paste(opt[['out']],"_lambda.",opt[['type']],sep=''))
PlotLambda(DataFra)
dev.off()

fctplot(out)
plot(xlim, ylim, type='n', xlab='Expected', ylab='Observed', main='',xlim=xlim, ylim=ylim)
for(Cmt in 1:length(listhead)){
lines(listfile[[listhead[Cmt]]][,1], listfile[[listhead[Cmt]]][,2],col=col[Cmt])
}
legend('topleft', legend=listhead, col=col, lty=1,bty='n')
dev.off()

