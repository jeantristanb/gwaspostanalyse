#!/usr/bin/env nextflow
/*
 * Authors       :
 *
 *
 *      Jean-Tristan Brandenburg
 *
 *  On behalf of the H3ABionet Consortium
 *  2015-2018
 *
 *
 * Description  : Nextflow pipeline for Wits GWAS.
 *
 */

//---- General definitions --------------------------------------------------//

import java.nio.file.Paths



def helps = [ 'help' : 'help' ]
allowed_params = ["cut_maf", "output_dir", "pb_around_rs", "mem_req", "work_dir","mem_req","big_time", "output","nb_cpu" ]
allowed_params_blocks = ["haploblocks", "plkref_haploblocks", "plk_othopt_haploblocks"]
allowed_params_other=["max_forks", "strandreport", "manifest", "idpat", "accessKey", "access-key", "secretKey", "secret-key","region", "AMI","maxInstances","instance-type", "instanceType", "bootStorageSize", "boot-storage-size", "max-instances", "sharedStorageMount", "shared-storage-mount", "scripts"]
allowed_params_head = ["head_pval", "head_freq", "head_bp", "head_chr", "head_rs", "head_beta", "head_se", "head_A1", "head_A2","max_pval", "max_pval_rep", "gwas_cat", "size_win", "genes_info","files_othertraits", "wind_merge"]
allowed_params+=allowed_params_head
allowed_params+=allowed_params_other
params.each { parm ->
  if (! allowed_params.contains(parm.key)) {
    println "\nUnknown parameter : Check parameter <$parm>\n";
  }
}


def params_help = new LinkedHashMap(helps)
params.queue      = 'batch'
params.work_dir   = "$PWD"
params.input_listfiles = "${params.work_dir}/list_files.input"
params.output_dir = "${params.work_dir}/output"
params.output = "replication"
params.cut_maf = 0.01

params.pb_around_rs = 25000

params.mem_req="8G"
params.big_time="100H"

params.head_pval = "P_BOLT_LMM"
params.head_freq = "A1FREQ"
params.head_bp = "BP"
params.head_chr = "CHR"
params.head_rs = "SNP"
params.out = ""
params.data = ""
params.head_beta="BETA"
params.head_se="SE"
params.head_A1="ALLELE1"
params.head_A2="ALLELE0"

params.wind_merge=100000

params.max_pval_rep=10**-6
params.size_win=25000
params.nb_cpu = 3

params.gwas_cat="/dataE/AWIGenGWAS/shared/ResultGWAS/Ressource/GWAS_Catalog_V37.tsv"
// haploblocks information
params.haploblocks=""
params.plkref_haploblocks=""
params.plk_othopt_haploblocks=""

//params plink
params.plink_bin='plink'
params.mem_plink='10G'
params.cpu_plink=2

 
//BedFileI=/dataE/AWIGenGWAS/shared/ResultGWAS/Ressource/locuszoom/data/1000G/genotypes/2014-10-14/EUR/chr$Chro
//plink -bfile $BedFileI --blocks no-pheno-req --out $DirOut/chr$Chro"_block"  &
if(params.haploblocks==''){
if(params.plkref_haploblocks==''){

}else{
listfile=params.plkref_haploblocks
ch_plkred_haplo=Channel.fromPath(listfile.split(',').collect{it+'.bed'}).merge(Channel.fromPath(listfile.split(',').collect{it+'.bim'})).merge(Channel.fromPath(listfile.split(',').collect{it+'.fam'}))
process ExtractBlocksBedFile{
    memory params.mem_plink
    time params.big_time
    cpus params.cpu_plink
    input :
     set file(bed), file(bim), file(fam) from ch_plkred_haplo
    output :
     // chr1_blocks.blocks.det
     file("${out}.blocks.det") into ch_block_bed
    script :
     plkhead=bed.baseName 
     out=plkhead
     """
     ${params.plink_bin} -bfile $plkhead --blocks no-pheno-req --out $out  --threads ${params.cpu_plink} ${params.plk_othopt_haploblocks}
     """
}
ch_block_merg=ch_block_bed.collect()
process MergeBlocsBedFile{
     input :
       file(allfile) from ch_block_merg
     publishDir "${params.output_dir}/blocs/", overwrite:true, mode:'copy'
     output :
        file(out) into ch_haploblocks
     script :
        out='block_all.blocks' 
        allfilejoin=allfile.join(",")
        """
        mergformat_haploblock.r $allfilejoin $out 
        """
}
}
}else{
ch_haploblocks=Channel.fromPath(params.haploblocks)
}




