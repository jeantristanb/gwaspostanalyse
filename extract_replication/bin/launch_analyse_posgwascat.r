#!/usr/bin/Rscript

library(data.table)
library('knitr')
library("optparse")
library('qqman')

#    analyse_posgwascat.r --chro_header_info ${params.head_chr_gwascat} --bp_header_info ${params.head_bp_gwascat} --gwas_cat $infogwas --file_gwas $gwas --chro_header_gwas ${params.head_chr}  --bp_header_gwas ${params.head_bp} --rs_header_gwas ${params.head_rs} --pval_header_gwas ${params.head_pval}
#    analyse_posgwascat.r --chro_header_gwascat ${params.head_chr_gwascat} --bp_header_gwascat ${params.head_bp_gwascat} --gwas_cat $infogwas --file_gwas $gwas --chro_header_gwas ${params.head_chr}  --bp_header_gwas ${params.head_bp} --rs_header_gwas ${params.head_rs} --pval_header_gwas ${params.head_pval}

option_list = list(
  make_option("--gwas_file", type="character",
              help="file gwas contains resultat ", metavar="character"),
  make_option("--bp_gwas", type="character",
              help="bp header for gwas file", metavar="character"),
  make_option("--chro_gwas", type="character",
              help="chro header for gwas file", metavar="character"),
  make_option("--pval_gwas", type="character",
              help="pvalue header for gwas file", metavar="character"),
  make_option("--beta_gwas", type="character",
              help="beta header for gwas file", metavar="character"),
  make_option("--se_gwas", type="character",
              help="beta header for gwas file", metavar="character"),
  make_option("--af_gwas", type="character",
              help="beta header for gwas file", metavar="character"),
  make_option("--rs_gwas", type="character",
              help="beta header for gwas file", metavar="character"),
  make_option("--gwas_cat", type="character",
              help="file gwas contains resultat ", metavar="character"),
  make_option("--chro_gwascat", type="character",
              help="file gwas contains resultat ", metavar="character"),
  make_option("--bp_gwascat", type="character",
              help="file gwas contains resultat ", metavar="character"),
  make_option("--af_gwascat", type="character",
              help="file gwas contains resultat ", metavar="character"),
  make_option("--beta_gwascat", type="character",
              help="file gwas contains resultat ", metavar="character"),
  make_option("--print_gwascat", type="character",
              help="file gwas contains resultat ", metavar="character"),
  make_option("--pval_gwascat", type="character",
              help="file gwas contains resultat ", metavar="character"),
  make_option("--info_gwascat", type="character",
              help="file gwas contains resultat ", metavar="character"),
  make_option("--threshpval_gwascat", type="double",
              help="file gwas contains resultat ", default=5*10**-8),
  make_option("--threshpval", type="double",
              help="file gwas contains resultat ", default=0.05)
)

opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser)
listcheck=c("gwas_cat", "gwas_file", 'bp_gwas', 'chro_gwas', 'chro_gwascat', 'bp_gwascat', 'pval_gwas')
#opt=list(gwas_file="/home/jeantristan/Travail/git/test_gwaspostanalyse/work/3d/cc5dd2aee5080bd5a275d1bf36f6e9/out.sub_gwas", chro_gwas="chr",bp_gwas="bp", pval_gwas="PVALUE_RE2", gwas_cat='~/Data/GWASCat/GWASCat_141019_ckd.tsv', chro_gwascat='Chro37', bp_gwascat='PosBegin37', pval_gwascat='P.VALUE',threshpval_gwascat=5*10**-8, pval_gwas='PVALUE_RE2')
for(arg in listcheck){
if(is.null(opt[[arg]])){
cat('args ',arg, ' not found \nexit\n')
q(2)
}
}
ChroGC=opt[['chro_gwascat']]
PosGC=opt[['bp_gwascat']]
GWASCat=opt[['gwas_cat']] 
InfoGC=opt[['MakeGwasCatInfo']]

ChroGW=opt[['chro_gwas']]
BPGW=opt[['bp_gwas']]
Pval=opt[['pval_gwas']]
GwasFile=opt[['gwas_file']]
PvalGC=opt[['pval_gwascat']]
LimPvalGC=opt[['threshpval_gwascat']]
LimPval=opt[['threshpval']]


filekniti='analyse_posgwascat.Rnw'
if(!file.exists(filekniti)){
fileknit=""
listpath=strsplit(Sys.getenv('PATH'), split=':')[[1]]
for(path in listpath){
if(file.exists(paste(path,filekniti,sep="/"))){
fileknit=paste(path,filekniti,sep="/")
}
}
}else{
fileknit=filekniti
}

filelibi='libanalyse_gwas.r'
if(!file.exists(filelibi)){
filelib=""
listpath=strsplit(Sys.getenv('PATH'), split=':')[[1]]
for(path in listpath){
if(file.exists(paste(path,filelibi,sep="/"))){
filelib=paste(path,filelibi,sep="/")
}
}
}else{
filelib=filelibi
}

if(filelib==""){
cat('not found ', filelib, ,'in path\n')
q(2)
}


source(filelib)


if(fileknit==""){
cat('not found ', fileknit, ,'in path\n')
q(2)
}

DataGWAS=as.data.frame(fread(GwasFile, header=T))
DataGWASCat<-read.table(GWASCat,sep='\t', header=T)
LimPval=opt[['threshpval']]


knit(fileknit)
