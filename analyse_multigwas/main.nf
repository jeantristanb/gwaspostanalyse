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
allowed_params = ["pheno", "covariates","input_listfiles", "cut_maf", "output_dir", "pb_around_rs", "mem_req", "work_dir","mem_req","big_time", "output","nb_cpu" ]
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
params.cut_maf = 0.01
params.pb_around_rs = 25000
params.mem_req="8G"
params.big_time="100H"
params.pheno=""
params.files_othertraits=""
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
params.max_pval=10**-6
params.max_pval_rep=10**-5
params.gwas_cat="/dataE/AWIGenGWAS/shared/ResultGWAS/Ressource/GWAS_Catalog_V37.tsv"
params.size_win=25000
params.genes_info="/dataE/AWIGenGWAS/shared/ResultGWAS/Ressource/gencode.v19.genes"
params.nb_cpu = 3



if(params.pheno=="" || params.input_listfiles==""){
println "pheno or input_listfiles not found"
System.exit(-1)
}

list_file_qq=Channel.create()
list_file_manh=Channel.create()
list_file_lchro=Channel.create()
list_file_analys=Channel.create()

Channel.fromPath(params.input_listfiles).separate(list_file_qq,list_file_manh, list_file_lchro, list_file_analys) { a -> [a,a,a,a] }
pheno_label_ch_qq = Channel.from(params.pheno.split(","))
process drawQQPlot {
    memory params.mem_req
    time params.big_time
    cpus params.nb_cpu
    input:
      file(list_file) from list_file_qq
    each this_pheno from pheno_label_ch_qq
    output:
      set this_pheno, file (output) into report_pca_pheno
    script:
      output="${this_pheno}.qq.jpeg"
    """
    QQMultiplot.r --maf ${params.cut_maf}  --pheno ${this_pheno} --out_file $output  --head_pval ${params.head_pval} --head_freq ${params.head_freq} --list_files $list_file --type_out jpeg --nThread ${params.nb_cpu}
    """
}

def GetTypeSex(File,ListTrait){
   File theInfoFile = new File( File )
   AllLines=theInfoFile.readLines()
   Head=AllLines[0].split()
   PosType = Head.findIndexOf { it ==~ /Type/ }  
   PosTrait = Head.findIndexOf { it ==~ /Pheno/ }
   ListRes=[]
   for(line in AllLines){
     SplLines=line.split()
     if(SplLines.size() >0 && ListTrait.contains(SplLines[PosTrait])){
        ListRes+=[[SplLines[PosTrait], SplLines[PosType]]]
     }
   }
   return ListRes
}
ListOfTraitType=GetTypeSex(params.input_listfiles,params.pheno.split(","))
list_traits_type=Channel.from(ListOfTraitType)

process drawManPlot {
    memory params.mem_req
    time params.big_time
    cpus params.nb_cpu
    input:
      file(list_file) from list_file_manh
    each  pheno_trait from list_traits_type
    output:
      set this_pheno,file (output) into report_pca_man
    script:
      this_pheno = pheno_trait[0]
      type       = pheno_trait[1]
      output="${this_pheno}.${type}.man.pdf"
    """
    ManPlot.r --maf ${params.cut_maf}  --pheno ${this_pheno} --out_file $output  --head_pval ${params.head_pval} --head_freq ${params.head_freq} --list_files $list_file --head_bp  ${params.head_bp} --head_chr ${params.head_chr} --head_rs ${params.head_rs} --type $type
    """
}

process getListeChro {
   memory params.mem_req
   time params.big_time
   input :
     file(list_file) from list_file_lchro
   output :
       stdout into listchro
   script:
    """
    head -2 $list_file > .tmp
    GetListChro.py --list_file .tmp --head_chr $params.head_chr
    """
}

check2=Channel.create()
listchro2=listchro.flatMap { list_str -> list_str.split() }.tap ( check2)
pheno_label_ch_analy = Channel.from(params.pheno.split(","))
if(params.files_othertraits){
list_file_other=Channel.fromPath(params.files_othertraits)
}else{
list_file_other=file('NO_FILE')
}

process ComputeStat{
   memory params.mem_req
   time params.big_time
   cpus params.nb_cpu
   input :
       file(list_file)  from list_file_analys
       file(other_file) from list_file_other
   each pheno from pheno_label_ch_analy 
   each chro from  listchro2
   output :
       //set pheno,file("${out}_nsig.csv"), file("${out}_rsgwasdet.csv"), file("${out}_rsresume.csv"),file("${out}_windgwasdet.csv"), file("${out}_windresume.csv"),file("${out}_othertrait.csv") into stats_multi
       set pheno, file("${out}*.csv") into stats_multi
   script :
    out="${pheno}_${chro}"
    other_trait       =  (params.files_othertraits) ? "--list_files_othertraits $other_file" : ""
    """
    AnalyseAll.r --maf ${params.cut_maf}  --pheno ${pheno} --out $out  --head_pval ${params.head_pval} --head_freq ${params.head_freq} --list_files $list_file --head_bp  ${params.head_bp} --head_chr ${params.head_chr} --head_rs ${params.head_rs}  --head_beta ${params.head_beta} --head_se ${params.head_se} --head_A1 ${params.head_A1} --head_A2 ${params.head_A2} --max_pval ${params.max_pval} --max_pval_rep ${params.max_pval_rep} --num_chr $chro --gwas_cat ${params.gwas_cat} --genes_info ${params.genes_info}  --size_win ${params.size_win} $other_trait --wind_merge ${params.wind_merge}
    """
}


report_pca_mangt=report_pca_man.groupTuple()
stats_multi_mer=stats_multi.groupTuple().join(report_pca_mangt).join(report_pca_pheno)
multiphenofile=Channel.from(1..params.pheno.split(",").size()).combine(Channel.fromPath(params.input_listfiles))
//multiphenofile=Channel.fromPath(params.input_listfiles)


if(params.covariates)covar="--covariates ${params.covariates}"
else covar=""

if(params.files_othertraits)othertraits="--list_files_othertraits ${params.files_othertraits}"
else othertraits=""

process MergeStat{
   memory params.mem_req
   time params.big_time
   input :
      val(listfile) from stats_multi_mer
      set val, file(fileinput) from multiphenofile
   publishDir "${params.output_dir}/${params.out}${pheno}", overwrite:true ,mode:'copy'
   output :
     set file("${pheno}.xlsx"), file("figure/*") ,file("${pheno}.tex"), file("${pheno}.pdf") into stats_multi_merg
   script :
     strs=listfile.flatten()
     pheno=strs[0]
     strs=strs.join("\n")
     list="${params.work_dir}/.tmplist_file_analyse_${pheno}"
     new File(list).withWriter{ it << strs } 
     """
     launc_MergeAll.r  --maf ${params.cut_maf} --gwas_cat ${params.gwas_cat} --max_pval ${params.max_pval} --max_pval_rep ${params.max_pval_rep} --size_win ${params.size_win} $othertraits --all_file $list --list_files $fileinput --wind_merge ${params.wind_merge} --head_rs ${params.head_rs}  --head_beta ${params.head_beta} --head_se ${params.head_se} --head_pval ${params.head_pval} $covar --head_freq ${params.head_freq}
     pdflatex $pheno >& /dev/null
     pdflatex $pheno 
     """
}








