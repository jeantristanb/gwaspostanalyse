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
allowed_params = ["cut_maf", "output_dir", "pb_around_rs", "mem_req", "work_dir","mem_req","big_time", "output","nb_cpu" , "input_dir","input_pat", "file_gwas", "gwas_cat", "site_wind", "r2_clump","min_pval_clump"]
allowed_params_blocks = ["haploblocks", "plkref_haploblocks", "plk_othopt_haploblocks"]
allowed_params_other=["max_forks", "strandreport", "manifest", "idpat", "accessKey", "access-key", "secretKey", "secret-key","region", "AMI","maxInstances","instance-type", "instanceType", "bootStorageSize", "boot-storage-size", "max-instances", "sharedStorageMount", "shared-storage-mount", "scripts"]
allowed_params_headinfo=["head_chr_gwascat", "head_bp_gwascat"]
allowed_params_head = ["head_pval", "head_freq", "head_bp", "head_chr", "head_rs", "head_beta", "head_se", "head_A1", "head_A2"]
allowed_params+=allowed_params_head
allowed_params+=allowed_params_other
allowed_params+=allowed_params_blocks
allowed_params+=allowed_params_headinfo
params.each { parm ->
  if (! allowed_params.contains(parm.key)) {
    println "\nUnknown parameter : Check parameter <$parm>\n";
  }
}

checker = { fn ->
   if (fn.exists())
       return fn;
    else
       error("\n\n------\nError in your config\nFile $fn does not exist\n\n---\n")
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

params.r2_clump=0.1
params.min_pval_clump=0.1
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


bed = Paths.get(params.input_dir,"${params.input_pat}.bed").toString()
bim = Paths.get(params.input_dir,"${params.input_pat}.bim").toString()
fam = Paths.get(params.input_dir,"${params.input_pat}.fam").toString()
raw_src_ch= Channel.create()

Channel
    .from(file(bed),file(bim),file(fam))
    .buffer(size:3)
    .map { a -> [checker(a[0]), checker(a[1]), checker(a[2])] }
    .set { raw_src_ch }


fileinfogwas_bed=Channel.fromPath(params.gwas_cat)

process SubBedFile{
    memory params.mem_plink
    time params.big_time
    cpus params.cpu_plink
   input :
     file(file_info) from fileinfogwas_bed
     set file(bed), file(bim), file(fam) from raw_src_ch
   output :
       set file("${out}.bed"), file("${out}.bim"), file("${out}.fam") into sub_plkfile
   script :
     plk=bed.baseName
     out=plk+"_sub"
     """
     subsample_plk.py --out $out --chro_header_info ${params.head_chr_gwascat} --bp_header_info ${params.head_bp_gwascat} --list_info ${file_info}  --bfile $plk --cpus ${params.cpu_plink} --wind_size ${params.size_win}
     """
}

filegwas=Channel.fromPath(params.file_gwas)
fileinfogwas_rep=Channel.fromPath(params.gwas_cat)


process ComputedReplication{
  input :
     set file(bed), file(bim), file(fam) from sub_plkfile
     file(gwas) from filegwas
     file(file_info) from fileinfogwas_rep
     file(haploblock) from ch_haploblocks
   output :
   script :
     plk=bed.baseName
     //allowed_params_head = ["head_pval", "head_freq", "head_bp", "head_chr", "head_rs", "head_beta", "head_se", "head_A1", "head_A2"]
     """ 
     extract_posclum.py --bfile $plk --chro_header_info ${params.head_chr_gwascat} --bp_header_info ${params.head_bp_gwascat} --list_info ${file_info}  --bfile $plk --cpus ${params.cpu_plink} --wind_size ${params.size_win} --out ${params.output} --chro_header_gwas ${params.head_chr}  --bp_header_gwas ${params.head_bp} --rs_header_gwas ${params.head_rs} --pval_header_gwas ${params.head_pval} --minpval ${params.min_pval_clump} --r2 ${params.r2_clump} --binplink ${params.plink_bin} --file_gwas $gwas --file_block_ld $haploblock
     """ 

}
//echo "python3 /home/jeantristan/Travail/GWAS/PythonScript/extract_posclum.py --list_info $FileCat --file_gwas $FileGWAS --wind_size $WindSize --out $Out --chro_header_info $ChrInfo --bp_header_info $BpInfo --chro_header_gwas $ChrGWAS --bp_header_gwas $BpGWAS --rs_header_gwas $RsGWAS --pval_header_gwas $PvalGWAS --bfile $bfile --maxpval $pval --r2 $R2 --file_block_ld All_block.blocks.extented.det" >> $BashFile




