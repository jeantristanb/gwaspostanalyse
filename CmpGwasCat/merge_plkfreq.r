library("optparse")
library(data.table)

option_list = list(
  make_option(c("--freq_plk"), type="character", default=NULL,
              help="gwas catalog input", metavar="character"),
  make_option(c("--out"), type="character", default=NULL,
              help="gwas catalog input", metavar="character"),
  make_option(c("--bim_plk"), type="character", default=NULL,
              help="gwas catalog input", metavar="character"))

opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);

datafrq<-fread(opt[['freq_plk']])
databim<-fread(opt[['bim_plk']], header=F)
head(databim)
databim<-databim[,c("V1","V2","V4"), with = FALSE]#with = FALSE]
head(databim)
names(databim)<-c("CHR", "SNP","BP")
allFreq<-merge(datafrq, databim, by=c("CHR", "SNP"),all=F)
write.table(allFreq, row.names=F, col.names=T, sep="\t", quote=F,file=opt[["out"]])

