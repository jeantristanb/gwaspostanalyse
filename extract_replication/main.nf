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

//params gwas cat 
params.head_bp_gwascat="Chro37"
params.head_chro_gwascat="Pos37"
params.head_pval_gwascat="P.VALUE"
//"SNPS"	"DATE.ADDED.TO.CATALOG"	"PUBMEDID"	"FIRST.AUTHOR"	"DATE"	"JOURNAL"	"LINK"	"STUDY"	"DISEASE.TRAIT"	"INITIAL.SAMPLE.SIZE"	"REPLICATION.SAMPLE.SIZE"	"REGION"	"CHR_ID"	"CHR_POS"	"REPORTED.GENE.S."	"MAPPED_GENE"	"UPSTREAM_GENE_ID"	"DOWNSTREAM_GENE_ID"	"SNP_GENE_IDS"	"UPSTREAM_GENE_DISTANCE"	"DOWNSTREAM_GENE_DISTANCE"	"STRONGEST.SNP.RISK.ALLELE"	"MERGED"	"SNP_ID_CURRENT"	"CONTEXT"	"INTERGENIC"	"RISK.ALLELE.FREQUENCY"	"P.VALUE"	"PVALUE_MLOG"	"P.VALUE..TEXT."	"OR.or.BETA"	"X95..CI..TEXT."	"PLATFORM..SNPS.PASSING.QC."	"CNV"	"MAPPED_TRAIT"	"MAPPED_TRAIT_URI"	"STUDY.ACCESSION"	"GENOTYPING.TECHNOLOGY"	"PosBegin38"	"PosEnd38"	"Chr38"	"Ref"	"Alt"	"Chro37"	"PosBegin37"	"PosEnd37"

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
     subsample_plk.py --out $out --chro_header_info ${params.head_chr_gwascat} --bp_header_info ${params.head_bp_gwascat} --list_info ${file_info}  --bfile $plk --cpus ${params.cpu_plink} --size_win_kb ${params.size_win_kb}
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
      file("${out}.clump.bypos") into clump_betweenpos
      file("${out}.clump.byrs") into clump_usingrs
      file("${out}.in.list_info") into res_bypos
      file("${out}.sub_gwas") into file_sub_gwas 
      file("${out}.clump.ldbloc.detail") into clump_ldbloc
   script :
     plk=bed.baseName
     out=params.output
     """ 
     extract_posclum.py --bfile $plk --chro_header_info ${params.head_chr_gwascat} --bp_header_info ${params.head_bp_gwascat} --list_info ${file_info}  --bfile $plk --cpus ${params.cpu_plink} --size_win_kb ${params.size_win_kb} --out ${params.output} --chro_header_gwas ${params.head_chr}  --bp_header_gwas ${params.head_bp} --rs_header_gwas ${params.head_rs} --pval_header_gwas ${params.head_pval} --minpval ${params.min_pval_clump} --r2 ${params.r2_clump} --binplink ${params.plink_bin} --file_gwas $gwas --file_block_ld $haploblock
     """ 
}

fileinfogwas_ana_bypos=Channel.fromPath(params.gwas_cat)
filegwas_anpos=Channel.fromPath(params.file_gwas)
process AnalyzeByPos{
  input :
     file(infogwas) from fileinfogwas_ana_bypos 
     file(file_res) from res_bypos
     file(gwas) from file_sub_gwas
  script :
    pvalgwascat =  (params.head_pval_gwascat!='') ? " --pval_gwascat ${params.head_pval_gwascat} --threshpval_gwascat ${params.threshold_pval_gwascat}" : ""

    """
    launch_analyse_posgwascat.r --chro_gwascat ${params.head_chr_gwascat} --bp_gwascat ${params.head_bp_gwascat} --gwas_cat $infogwas --gwas_file $gwas --chro_gwas ${params.head_chr}  --bp_gwas ${params.head_bp} --rs_gwas ${params.head_rs} $pvalgwascat --pval_gwas ${params.head_pval} --threshpval ${params.threshpval} --print_gwascat ${params.info_gwascat}

    """
}
//echo "python3 /home/jeantristan/Travail/GWAS/PythonScript/extract_posclum.py --list_info $FileCat --file_gwas $FileGWAS --wind_size $WindSize --out $Out --chro_header_info $ChrInfo --bp_header_info $BpInfo --chro_header_gwas $ChrGWAS --bp_header_gwas $BpGWAS --rs_header_gwas $RsGWAS --pval_header_gwas $PvalGWAS --bfile $bfile --maxpval $pval --r2 $R2 --file_block_ld All_block.blocks.extented.det" >> $BashFile




