#!/usr/bin/Rscript

library(data.table)
library('knitr')
library("optparse")
library('qqman')
library(kableExtra)
#    analyse_posgwascat.r --chro_header_info ${params.head_chr_gwascat} --bp_header_info ${params.head_bp_gwascat} --gwas_cat $infogwas --file_gwas $gwas --chro_header_gwas ${params.head_chr}  --bp_header_gwas ${params.head_bp} --rs_header_gwas ${params.head_rs} --pval_header_gwas ${params.head_pval}
#    analyse_posgwascat.r --chro_header_gwascat ${params.head_chr_gwascat} --bp_header_gwascat ${params.head_bp_gwascat} --gwas_cat $infogwas --file_gwas $gwas --chro_header_gwas ${params.head_chr}  --bp_header_gwas ${params.head_bp} --rs_header_gwas ${params.head_rs} --pval_header_gwas ${params.head_pval}

option_list = list(
  make_option("--gwas_file", type="character",
              help="file gwas contains resultat ", metavar="character"),
  make_option("--res_block", type="character",
              help="bp header for gwas file", metavar="character"),
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
  make_option("--rs_gwascat", type="character",
              help="beta header for gwas file", metavar="character"),
  make_option("--gwas_cat", type="character",
              help="file gwas contains resultat ", metavar="character"),
  make_option("--chro_gwascat", type="character",
              help="file gwas contains resultat ", metavar="character"),
  make_option("--bp_gwascat", type="character",
              help="file gwas contains resultat ", metavar="character"),
  make_option("--af_gwascat", type="character",
              help="file gwas contains resultat ", metavar="character"),
  make_option("--size_win_kb", type="double",
              help="file gwas contains resultat "),
  make_option("--beta_gwascat", type="character",
              help="file gwas contains resultat ", metavar="character"),
  make_option("--print_gwascat", type="character",
              help="file gwas contains resultat ", metavar="character", default=NA),
  make_option("--pval_gwascat", type="character",
              help="file gwas contains resultat ", metavar="character", default=NA),
  make_option("--info_gene", type="character",
              help="file gwas contains resultat ", metavar="character"),
  make_option("--threshpval_gwascat", type="double",
              help="file gwas contains resultat ", default=5*10**-8),
  make_option("--haploblocks", type="character",
              help="file gwas contains resultat "),
  make_option("--clump", type="character",
              help="file gwas contains resultat "),
  make_option("--threshpval", type="double",
              help="file gwas contains resultat ", default=0.05),
  make_option("--a1_gwas", type="character",
              help="beta header for gwas file", metavar="character"),
  make_option("--a0_gwas", type="character",
              help="beta header for gwas file", metavar="character")
)

opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser)
listcheck=c("gwas_cat", "gwas_file", 'bp_gwas', 'chro_gwas', 'chro_gwascat', 'bp_gwascat', 'pval_gwas', 'print_gwascat','size_win_kb')
#opt=list(gwas_file="/home/jeantristan/Travail/git/test_gwaspostanalyse/work/3d/cc5dd2aee5080bd5a275d1bf36f6e9/out.sub_gwas", chro_gwas="chr",bp_gwas="bp", pval_gwas="PVALUE_RE2", gwas_cat='~/Data/GWASCat/GWASCat_141019_ckd.tsv', chro_gwascat='Chro37', bp_gwascat='PosBegin37', pval_gwascat='P.VALUE',threshpval_gwascat=5*10**-8, pval_gwas='PVALUE_RE2')
for(arg in listcheck){
if(is.null(opt[[arg]])){
cat('args ',arg, ' not found \nexit\n')
q(2)
}
}
GWASCat=opt[['gwas_cat']] 
ChroGC=opt[['chro_gwascat']]
PosGC=opt[['bp_gwascat']]
InfoGC=opt[['print_gwascat']]
PvalGC=opt[['pval_gwascat']]
RsGC=opt[['rs_gwascat']]

if(is.null(PvalGC) | length(PvalGC)==0 |PvalGC=='NA' | PvalGC=='na' | PvalGC=='GC' | is.na(PvalGC)){
PvalGC=NULL
}
if(is.null(InfoGC) | length(InfoGC)==0 |InfoGC=='NA' | InfoGC=='na' | InfoGC=='GC' | is.na(InfoGC)){
InfoGC=NULL
}else{
InfoGC=unlist(strsplit(InfoGC, split=','))
}

#ChroGC='Chro37';PosGC='PosBegin37';InfoGC='DISEASE.TRAIT,REPORTED.GENE.S.,MAPPED_GENE,INITIAL.SAMPLE.SIZE';#;RsGC='RSID'
#RsGC='SNPS'
#PvalGW="PVALUE_RE2";ChroGW='chr';PosGW='bp';RsGW="RSID";AfGW=NULL;BetaGW=NULL;A1GW=NULL;A0GW=NULL;SeGW=NULL;SeGW=NULL;

GwasFile=opt[['gwas_file']]
ChroGW=opt[['chro_gwas']]
PosGW=opt[['bp_gwas']]
SeGW=opt[['se_gwas']]
BetaGW=opt[['beta_gwas']]
RsGW=opt[['rs_gwas']]
A1GW=opt[['a1_gwas']];
A0GW=opt[['a0_gwas']]
PvalGW=opt[['pval_gwas']]
AfGW=opt[['af_gwas']]
LimPvalGC=opt[['threshpval_gwascat']]
LimPval=opt[['threshpval']]

InfoGene=opt[['info_gene']]
ChroGE="CHR";BeginGE="BEGIN";EndGE="END";NameGene="GENE"
Haploblock=opt[['haploblocks']]

filekniti='analyse_posgwascat_byblock.Rnw'
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
listout_excel=list()

if(fileknit==""){
cat('not found ', fileknit, ,'in path\n')
q(2)
}



DataGWAS=as.data.frame(fread(GwasFile, header=T))
DataGWASCat<-read.table(GWASCat,sep='\t', header=T)
DataGene=read.table(InfoGene ,sep='\t', header=T)
DataHaplo=read.table(Haploblock,sep='\t', header=T)
SizeWind=opt[['size_win_kb']]*1000
LimPval=opt[['threshpval']]
DataClump=read.table(opt[['clump']], header=T)
DataResBlockI=read.table(opt[['res_block']], header=T)
DirPWD=getwd()

GCHeadTmp<-c(ChroGC,PosGC, PvalGC,InfoGC, RsGC)
GCHeadTmpNF<-GCHeadTmp[!(GCHeadTmp %in% names(DataGWASCat))]
if(length(GCHeadTmpNF)>0){
cat('not found', GCHeadTmpNF, 'in info\n')
q(2)
}

GWHeadTmp<-c(ChroGW,PosGW, SeGW, BetaGW, RsGW, PvalGW,AfGW)
GWHeadTmp<-GWHeadTmp[!(GWHeadTmp %in% names(DataGWAS))]
if(length(GWHeadTmp)>0){
cat('not found', GWHeadTmp, 'in gwas')
q(2)
}

knit(fileknit)
library(openxlsx)
options(java.parameters = "-Xmx8000m")
if(length(listout_excel)==0){
listout_excel[['res_windld']]=matrix('nothing found',nrow=1,ncol=1)
dir.create('figure/')
writeLines(c('no result'),con='figure/none')
}
write.xlsx(listout_excel,file= paste("resume_tab.xlsx",sep=""))


system('pdflatex analyse_posgwascat_byblock')
system('pdflatex analyse_posgwascat_byblock')
