library("optparse")

option_list = list(
  make_option(c("--gwascat"), type="character", default=NULL,
              help="gwas catalog input", metavar="character"),
  make_option(c("--gc_chr"), type="character", default=NULL,
              help="gwas catalog input", metavar="character"),
  make_option(c("--gc_bp"), type="character", default=NULL,
              help="gwas catalog input", metavar="character"),
  make_option(c("--gc_rs"), type="character", default=NULL,
              help="gwas catalog input", metavar="character"),
  make_option(c("--sep"), type="character", default="\t",
              help="gwas catalog input", metavar="character"),
  make_option(c("--out"), type="character", default=NULL,
              help="gwas catalog input", metavar="character")
)

opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);
data<-read.table(opt[['gwascat']], header=T, sep=opt[['sep']], stringsAsFactors=F)
cat(c(opt[["gc_chr"]], opt[["gc_bp"]], opt[["gc_bp"]]))
chropos<-unique(data[,c(opt[["gc_chr"]], opt[["gc_bp"]], opt[["gc_bp"]], opt[["gc_rs"]])])
write.table(chropos, row.names=F, col.names=F, file=opt[["out"]], quote=F, sep='\t')

