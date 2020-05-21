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
allowed_params = ["cut_maf", "output_dir", "pb_around_rs", "mem_req", "work_dir","mem_req","big_time", "output","nb_cpu" , "input_dir","input_pat", "file_gwas", "gwas_cat", "site_wind", "r2_clump","min_pval_clump", "size_win_kb"]
allowed_params_blocks = ["haploblocks", "plkref_haploblocks", "plk_othopt_haploblocks"]
allowed_params_other=["max_forks", "strandreport", "manifest", "idpat", "accessKey", "access-key", "secretKey", "secret-key","region", "AMI","maxInstances","instance-type", "instanceType", "bootStorageSize", "boot-storage-size", "max-instances", "sharedStorageMount", "shared-storage-mount", "scripts"]
allowed_params_headinfo=["head_chro_gwascat", "head_bp_gwascat", "head_pval_gwascat"]
allowed_params_head = ["head_pval", "head_freq", "head_bp", "head_chr", "head_rs", "head_beta", "head_se", "head_A1", "head_A0"]
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


params.mem_req="8G"
params.big_time="100H"

params.head_pval = "P_BOLT_LMM"
params.head_freq = ""
params.head_bp = "BP"
params.head_chr = "CHR"
params.head_rs = "SNP"
params.out = ""
params.data = ""
params.head_beta=""
params.head_se=""
params.head_A1="ALLELE1"
params.head_A0="ALLELE0"


params.max_pval_rep=10**-6
params.size_win_kb=25
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
params.genes_file=""
params.gene_file_ftp="ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_19/gencode.v19.annotation.gtf.gz"

//params gwas cat 
params.head_bp_gwascat="PosBegin37"
params.head_chro_gwascat="Chro37"
params.head_pval_gwascat=""
params.head_rs_gwascat="SNPS"


params.info_gwascat="DISEASE.TRAIT,REPORTED.GENE.S.,MAPPED_GENE,INITIAL.SAMPLE.SIZE"
params.threshold_pval_gwascat=1
params.threshpval=0.05

 
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
        file(out) into (ch_haploblocks, ch_haploblocks_byblock, ch_haploblocks_bypos, ch_haploblocks_clump)
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
ch_haploblocks_byblock=Channel.fromPath(params.haploblocks)
ch_haploblocks_bypos=Channel.fromPath(params.haploblocks)
ch_haploblocks_clump=Channel.fromPath(params.haploblocks)
}
if(params.genes_file==""){
process GetGeneFile{
     publishDir "${params.output_dir}/genes/", overwrite:true, mode:'copy'
     output:
       file(out) into (geneinfo_ch_pos, geneinfo_ch_byblock, geneinfo_ch_byclump) 
     script :
       out='gene_info.gene'
       """
       wget -O tmpfile.gz ${params.gene_file_ftp}
       gunzip tmpfile.gz
       format_infogene.r tmpfile $out  
       """

}

}else{
geneinfo_ch_pos=Channel.fromPath(params.genes_file)
geneinfo_ch_clump=Channel.fromPath(params.genes_file)
geneinfo_ch_byblock=Channel.fromPath(params.genes_file)
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

filegwassub=Channel.fromPath(params.file_gwas)
process SubBedFile{
    memory params.mem_plink
    time params.big_time
    cpus params.cpu_plink
   input :
     file(file_info) from fileinfogwas_bed
     set file(bed), file(bim), file(fam) from raw_src_ch
     file(gwas) from filegwassub
   output :
       set file("${out}.bed"), file("${out}.bim"), file("${out}.fam") into sub_plkfile
   script :
     plk=bed.baseName
     out=plk+"_sub"
     """
     subsample_plk.py --out $out --chro_header_info ${params.head_chro_gwascat} --bp_header_info ${params.head_bp_gwascat} --list_info ${file_info}  --bfile $plk --cpus ${params.cpu_plink} --size_win_kb ${params.size_win_kb} --chro_header_gwas ${params.head_chr}  --bp_header_gwas ${params.head_bp} --rs_header_gwas ${params.head_rs}  --file_gwas $gwas --a1_header_gwas ${params.head_A1} --a0_header_gwas ${params.head_A0}

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
   publishDir "${params.output_dir}/tmpfile/", overwrite:true, mode:'copy'
   output :
      file("${out}.clump.bypos") into clump_betweenpos
      file("${out}.clump.byrs") into res_byclump
      file("${out}.in.list_info") into res_bypos
      file("${out}.sub_gwas") into (file_sub_gwas_bypos,file_sub_gwas_byblock, file_sub_gwas_byclump)
      file("${out}.clump.ldbloc.detail") into res_byblock
      file("${out}.plk.clumped") into (clump_res_byblock,clump_res_pos, clump_res_clump)
   script :
     plk=bed.baseName
     out=params.output
     """ 
     extract_posclum.py --bfile $plk --chro_header_info ${params.head_chro_gwascat} --bp_header_info ${params.head_bp_gwascat} --list_info ${file_info}  --bfile $plk --cpus ${params.cpu_plink} --size_win_kb ${params.size_win_kb} --out ${params.output} --chro_header_gwas ${params.head_chr}  --bp_header_gwas ${params.head_bp} --rs_header_gwas ${params.head_rs} --pval_header_gwas ${params.head_pval} --minpval ${params.min_pval_clump} --r2 ${params.r2_clump} --binplink ${params.plink_bin} --file_gwas $gwas --file_block_ld $haploblock
     """ 
}

fileinfogwas_ana_byclump=Channel.fromPath(params.gwas_cat)
filegwas_anclump=Channel.fromPath(params.file_gwas)
process AnalyzeByClump{
  input :
     file(infogwas) from fileinfogwas_ana_byclump
     file(file_res) from res_byclump
     file(gwas) from file_sub_gwas_byclump
     file(geneinfo) from geneinfo_ch_clump
     file(haploblocks) from ch_haploblocks_clump
     file(clump) from clump_res_clump
  publishDir "${params.output_dir}/clump/", overwrite:true, mode:'copy'
  output :
    file("$out")
    file("$outxlxs")
    file('figure/*')
  script :
    pvalgwascat =  (params.head_pval_gwascat!='') ? " --threshpval_gwascat ${params.threshold_pval_gwascat}   --pval_gwascat ${params.head_pval_gwascat}" : ""
    out=params.output+"_byclump.pdf"
    outxlxs=params.output+"_byclump.xlsx"
    betagwascat= (params.head_beta!='') ? " --beta_gwas ${params.head_beta}" : ""
    freqgwascat= (params.head_freq!='') ? " --af_gwas ${params.head_freq}" : ""
    segwascat= (params.head_freq!='') ? " --se_gwas ${params.head_freq}" : ""
    """
    launch_analyse_posgwascat_byclump.r --res_clump $file_res --chro_gwascat ${params.head_chro_gwascat} --bp_gwascat ${params.head_bp_gwascat} --gwas_cat $infogwas --gwas_file $gwas --chro_gwas ${params.head_chr}  --bp_gwas ${params.head_bp} --rs_gwas ${params.head_rs} $pvalgwascat --pval_gwas ${params.head_pval} --threshpval ${params.threshpval} --print_gwascat ${params.info_gwascat} --info_gene $geneinfo --haploblocks $haploblocks --clump $clump --size_win_kb ${params.size_win_kb} --rs_gwascat ${params.head_rs_gwascat} --a1_gwas ${params.head_A1} --a0_gwas ${params.head_A0}  $betagwascat $freqgwascat $segwascat
    mv analyse_posgwascat_byclump.pdf $out
    mv resume_tab.xlsx $outxlxs
    """
}





fileinfogwas_ana_bypos=Channel.fromPath(params.gwas_cat)
filegwas_anpos=Channel.fromPath(params.file_gwas)
process AnalyzeByPos{
  input :
     file(infogwas) from fileinfogwas_ana_bypos 
     file(file_res) from res_bypos
     file(gwas) from file_sub_gwas_bypos
     file(geneinfo) from geneinfo_ch_pos
     file(haploblocks) from ch_haploblocks_bypos
     file(clump) from clump_res_pos
  publishDir "${params.output_dir}/bypos/", overwrite:true, mode:'copy'
  output :
    file("$out")
    file("$outxlxs")
    file('figure/*')
  script :
    pvalgwascat =  (params.head_pval_gwascat!='') ? " --threshpval_gwascat ${params.threshold_pval_gwascat}   --pval_gwascat ${params.head_pval_gwascat}" : ""
    out=params.output+"_bypos.pdf"
    outxlxs=params.output+"_bypos.xlsx"
    betagwascat= (params.head_beta!='') ? " --beta_gwas ${params.head_beta}" : ""
    freqgwascat= (params.head_freq!='') ? " --af_gwas ${params.head_freq}" : ""
    segwascat= (params.head_freq!='') ? " --se_gwas ${params.head_freq}" : ""
    """
    launch_analyse_posgwascat_bypos.r --chro_gwascat ${params.head_chro_gwascat} --bp_gwascat ${params.head_bp_gwascat} --gwas_cat $infogwas --gwas_file $gwas --chro_gwas ${params.head_chr}  --bp_gwas ${params.head_bp} --rs_gwas ${params.head_rs} $pvalgwascat --pval_gwas ${params.head_pval} --threshpval ${params.threshpval} --print_gwascat ${params.info_gwascat} --info_gene $geneinfo --haploblocks $haploblocks --clump $clump --size_win_kb ${params.size_win_kb}  $betagwascat $freqgwascat $segwascat
    mv analyse_posgwascat_bypos.pdf $out
    mv resume_tab.xlsx $outxlxs
    """
}

fileinfogwas_ana_byblock=Channel.fromPath(params.gwas_cat)
filegwas_byblock=Channel.fromPath(params.file_gwas)
process AnalyzeByBlock{
  input :
     file(infogwas) from fileinfogwas_ana_byblock
     file(file_res) from res_byblock
     file(gwas) from file_sub_gwas_byblock
     file(geneinfo) from geneinfo_ch_byblock
     file(haploblocks) from ch_haploblocks_byblock
     file(clump) from clump_res_byblock
  publishDir "${params.output_dir}/byblock/", overwrite:true, mode:'copy'
  output :
    file("$out")
    file("$outxlxs")
    file('figure/*')
  script :
    pvalgwascat =  (params.head_pval_gwascat!='') ? " --pval_gwascat ${params.head_pval_gwascat} --threshpval_gwascat ${params.threshold_pval_gwascat}" : ""
    betagwascat= (params.head_beta!='') ? " --beta_gwas ${params.head_beta}" : ""
    freqgwascat= (params.head_freq!='') ? " --af_gwas ${params.head_freq}" : ""
    segwascat= (params.head_freq!='') ? " --se_gwas ${params.head_freq}" : ""
    out=params.output+"_byblock.pdf"
    outxlxs=params.output+"_byblock.xlsx"
    """
    launch_analyse_posgwascat_byblock.r --res_block $file_res --chro_gwascat ${params.head_chro_gwascat} --bp_gwascat ${params.head_bp_gwascat} --gwas_cat $infogwas --gwas_file $gwas --chro_gwas ${params.head_chr}  --bp_gwas ${params.head_bp} --rs_gwas ${params.head_rs} $pvalgwascat --pval_gwas ${params.head_pval} --threshpval ${params.threshpval} --print_gwascat ${params.info_gwascat} --info_gene $geneinfo --haploblocks $haploblocks --clump $clump --size_win_kb ${params.size_win_kb} --rs_gwascat ${params.head_rs_gwascat} --a1_gwas ${params.head_A1} --a0_gwas ${params.head_A0} $betagwascat $freqgwascat $segwascat
    mv analyse_posgwascat_byblock.pdf $out
    mv resume_tab.xlsx $outxlxs
    """
}



