# Transferability for gwas result

## Installation
* python3 : matplotlib,pandas,numpy
* R : 
* nextflow
* plink

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
 * `head_freq`
 * `head_bp` : bp 
 * `head_chr` : chr
 * `head_rs` : rs
 * `head_beta`
 * `head_se` 
 * `head_A1`
 * head_A2`
* `gwas_cat` : position to check transferability, must be separate by tabulation
 * head_chr_gwascat
 * head_bp_gwascat
* `size_wind` : position around each pos of gwas cat

### other parameters :
  * `plink_bin` : binary for plink [default : plink ]
  * `mem_plink` : memory for plink [default : '10G']
  * `cpu_plink` : cpu number for plink [default :2 ]


