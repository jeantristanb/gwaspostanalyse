source('/home/jeantristan/Travail/git/gwaspostanalyse/extract_replication/bin/libanalyse_gwas.r')
library(data.table)
library("optparse")
#library('qqman')

Test=T
BestSol=F
ChroBP=F
if(Test==F){
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
              help="beta header for gwas file", metavar="character"),
  make_option("--out", type="character",default="out",
              help="beta header for gwas file", metavar="character"),
)





opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser)

GWASCat=opt[['gwas_cat']]
ChroGC=opt[['chro_gwascat']]
PosGC=opt[['bp_gwascat']]
InfoGC=opt[['print_gwascat']]
PvalGC=opt[['pval_gwascat']]
RsGC=opt[['rs_gwascat']]


FileClump=opt[['clump']]
FileGWAS=opt[['gwas_file']]
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
FileResBlocI=opt[['res_block']]
WindSize=opt[['size_win_kb']]*1000
InfoGene=opt[['info_gene']]
HaploBlock=opt[['haploblocks']]
FileBlock<-opt[['clump']]
InfoGC=opt[['print_gwascat']]
ChroGE="CHR";BeginGE="BEGIN";EndGE="END";NameGene="GENE"
OutPatt=opt[['out']]
if(!is.null(opt[['chr_plot']]) & !is.null(opt[['bp_plot']])){
ChrPlot=opt[['chr_plot']]
BpPlot=as.integer(opt[['bp_plot']])
ChroBPPlot=T
}
if(!is.null(opt[[bestsol]]))BestSol=T
}else{
FileResBlocI<-'Old/4ef048faf214607018dd5a234cf49f/out.clump.ldbloc.detail'
FileClump='Old/4ef048faf214607018dd5a234cf49f/out.plk.clumped'
FileGWAS<-'Old/4ef048faf214607018dd5a234cf49f/out.sub_gwas'
GWASCat<-'Old/4ef048faf214607018dd5a234cf49f/GWASCat_27_March_2020_ckd.tsv'
InfoGene='Old/4ef048faf214607018dd5a234cf49f/gene_info.gene'
HaploBlock='Old/4ef048faf214607018dd5a234cf49f/block_all.blocks'
WindSize<-10000

ChroGW<-"chr"
PosGW<-"ps"
RsGW="rs"
AfGW<-'af'
BetaGW<-'beta'
A1GW<-'allele1'
A0GW<-'allele0'
SeGW<-'se'
PvalGW="p_wald"
ChroGE="CHR"
BeginGE="BEGIN"
EndGE="END"
NameGene="GENE"
InfoGC="DISEASE.TRAIT,REPORTED.GENE.S.,MAPPED_GENE,INITIAL.SAMPLE.SIZE"
ChroGC="Chro37";PosGC="PosBegin37";RsGC="SNPS"
OutPatt="test"
}
DataResBlockI<-read.table(FileResBlocI, header=T)
DataClump<-read.table(FileClump, header=T)
DataGWAS<-read.table(FileGWAS,header=T)
DataGWASCat=read.table(GWASCat, header=T)
DataGene<-read.table(InfoGene, header=T)
DataBlock=read.table(HaploBlock,sep='\t', header=T)
listInfoGC<-strsplit(InfoGC, split=',')[[1]]

#ListRs="rs73041161"
DataGC<-MakeGwasCatInfo(DataGWASCat,ChroGC,PosGC, paste(c(InfoGC, RsGC),collapse=','))
DataWind<-BuildWind(DataGWASCat,ChroGC,PosGC, InfoGC,DataResBlockI,  DataGWAS,PvalGW,BetaGW, RsGW, A1GW, A0GW,SeGW, AfGW,listInfoGC)
#if()BestSol<-DataWind[order(DataWind$PClump)[1],]
if(BestSol)AllSolAn<-DataWind[order(DataWind$PClump)[1],]
else if(ChroBPPlot){
#ChrPlot=opt[['chr_plot']]
#BpPlot=as.integer(opt[['bp_plot']])
#ChroBPPlot=T
AllSolAn<-DataWind[DataWind$Chro==ChrPlot & (DataWind$BeginBlock-WindSize)<BpPlot & (DataWind$EndBlock+WindSize)>BpPlot,]
}else {
AllSolAn<-DataWind[DataWind$p.adjust.bonf<0.05,]
}

for(CmtWind in 1:nrow(AllSolAn)){
SolAn<-AllSolAn[CmtWind,]
Chro<-SolAn$Chro;BeginBlock<-SolAn$BeginBlock;EndBlock<-SolAn$EndBlock;ListRs=as.character(SolAn[,'MinPRsClump'])
pdf(paste(OutPatt,'_', Chro,'_',BeginBlock,'_',EndBlock,'_',WindSize,'pdf',sep=''))
PlotWindLocusZoom(SolAn$Chro,SolAn$BeginBlock, SolAn$EndBlock, WindSize,DataGWAS, ChroGW,PosGW, RsGW,PvalGW,DataGC,ChroGC, PosGC, DataClump, ListRs, DataGene,ChroGE, BeginGE, EndGE, NameGene, DataBlock)
dev.off()
}
