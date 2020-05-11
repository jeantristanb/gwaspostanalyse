# Transferability for gwas result

## Installation
* python3 : matplotlib,pandas,numpy
* R : data.table
* nextflow
* plink
* pdflatex

## Parameters
###

### blocs for transferabilties
* if `haploblocks` file not gave by user, pipeline build with plink file(s) and plink blocs with 
* `haploblocks` : contains blocs to analyse
* build block :
 * `plkref_haploblocks` : list of bed file separate by comma to build blocs (one or more )
 * `plk_othopt_haploblocks` : other option for plink 
 * create folder _blocs_ where blocs list in one file
 * used `--blocks` option to build plink option

### Transferability
research by position, windows, clump and block tranferability 
* `haploblocks` : 
* `file_gwas` : file gwas
 * `head_pval` : ok
 * `head_freq` : ?
 * `head_bp` : bp 
 * `head_chr` : chr
 * `head_rs` : rs
 * `head_beta` : ?
 * `head_se`  : ?
 * `head_A1` :  effector allele for gwas file
 * head_A2` : Other allele for gwas file
* `gwas_cat` : position to check transferability, must be separate by tabulation
 * `head_chro_gwascat` : header for chr of gwas catalog
 * `head_bp_gwascat` : header for bp of gwas 
 * `head_pval_gwascat` : header for bp of gwas 
 * `head_rs_gwascat` : header for rs of gwas catalog
 * `info_gwascat` : headers of gwas catalog to print separate by a comma 
* `size_wind` : position around each pos of gwas cat
* `threshpval` : threshold pvalue after adjustment [ default : 0.05]
* `genes_file` : gene file for annotation, 
  * must be contains : "CHR"	"BEGIN"	"END"	"GENE"	"ID"
  * if not present will be download to :
  * `gene_file_ftp` ftp to download a file of genes and format [ default : "ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_19/gencode.v19.annotation.gtf.gz" ]


### other parameters :
  * `plink_bin` : binary for plink [default : plink ]
  * `mem_plink` : memory for plink [default : '10G']
  * `cpu_plink` : cpu number for plink [default :2 ]


