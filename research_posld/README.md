# research chip in ld with 
Objectif of pipeline is to identifie chip (or other) in ld with significant position 
## strategie clump 
 * extract pos in ld between pos  in list infoand other position
## arguments
* `list_info` :  file contains chromosome, position to search
 * `chro_header_info` : chro header
 * `-bp_header_info` : bp header
* `file_gwas` : file gwas
 * `chro_header_gwas` : chro file of gwas
 * `rs_header_gwas` : 
 * `bp_header_info` : header for bp
 * `pval_header_gwas` : header of pvalue
* `out` : out header
* parametre ld :
 * `size_win_kb` : windows for clump and ld
 * `bfile`
 * `minpvalue` : min pvalue
* parametre plink :
 * `keep` : option keep of plink
 * `cpus` : cpu number for plink 
 * `binplinks` : directory of plink

