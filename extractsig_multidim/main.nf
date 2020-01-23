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
allowed_params = ["cut_maf", "output_dir", "pb_around_rs", "mem_req", "work_dir","mem_req","big_time", "output","nb_cpu" , "input_dir","input_pat", "file_gwas", "gwas_cat", "pval_thresh"]
allowed_params_blocks = ["haploblocks", "plkref_haploblocks", "plk_othopt_haploblocks"]
allowed_params_other=["max_forks", "strandreport", "manifest", "idpat", "accessKey", "access-key", "secretKey", "secret-key","region", "AMI","maxInstances","instance-type", "instanceType", "bootStorageSize", "boot-storage-size", "max-instances", "sharedStorageMount", "shared-storage-mount", "scripts"]
allowed_params_headinfo=["head_chr_gwascat", "head_bp_gwascat"]
#allowed_params_head = ["head_pval", "head_freq", "head_bp", "head_chr", "head_rs", "head_beta", "head_se", "head_A1", "head_A2"]
#allowed_params+=allowed_params_head
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

#params.head_pval = "P_BOLT_LMM"
#params.head_freq = "A1FREQ"
#params.head_bp = "BP"
#params.head_chr = "CHR"
#params.head_rs = "SNP"
params.data = ""
#params.head_beta="BETA"
#params.head_se="SE"
#params.head_A1="ALLELE1"
#params.head_A2="ALLELE0"


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
#params.head_bp_gwascat="Chro37"
#params.head_chro_gwascat="Pos37"
#params.head_pval_gwascat="P.VALUE"
#params.head_rs_gwascat="SNPS"

params.info_gwascat="DISEASE.TRAIT,REPORTED.GENE.S.,MAPPED_GENE,INITIAL.SAMPLE.SIZE"
params.threshold_pval_gwascat=1
params.threshpval=0.05

 
//BedFileI=/dataE/AWIGenGWAS/shared/ResultGWAS/Ressource/locuszoom/data/1000G/genotypes/2014-10-14/EUR/chr$Chro
//plink -bfile $BedFileI --blocks no-pheno-req --out $DirOut/chr$Chro"_block"  &

params.r2_clump=0.1
params.min_pval_clump=0.1

def configfile_analysis(file){
   sep=','
   File theInfoFile = new File( file )
   if( !theInfoFile.exists() ) {
      println "File "+ file+" does not exist"
      exit(0)
   }
   def lines = theInfoFile.readLines()
   def SplH=lines[0].split(sep)
   resInfo=[]
   resFile=[]
   infoFile=[]
   lines.remove(0)
   PosCmp=-1
   CmtL=0
   for(line in lines){
       def splLine=line.split(sep)
       def cmtelem=0
       def SubRes=[]
       while (cmtelem < SplH.size()){
            if(SplH[cmtelem].toUpperCase()=='FILE'){
               resFile.add(splLine[cmtelem])
            }
            if(SplH[cmtelem].toUpperCase()=='INFO'){
               resFile.add(splLine[cmtelem])
            }
            if(splLine[cmtelem]!='NA' && splLine[cmtelem]!='' && SplH[cmtelem]!='FILE' && SplH[cmtelem]!='INFO'){
                 SubRes.add(SplH[cmtelem].toUpperCase()+':'+splLine[cmtelem])
            }
            cmtelem+=1
       }
       resInfo.add(SubRes.join(','))
       CmtL+=1
   }
 return([resFile,resInfo, NumRef])
}

info_file=configfile_analysis(params.file_config)


liste_filesi_ch=Channel.fromPath(info_file[0]).merge(Channel.from(info_file[1]))

process ExtractSig{
    memory ma_mem_req
    input :
      set file(file_assoc), val(info_file) from liste_filesi_ch
    output :
      file(filers) into (listrs_sig)
    script :
       newfile_assoc=file_assoc+".modif"
       """
       extract_rssig.py  --input_file $file_assoc --out_file $newfile_assoc --info_file $info_file --threshold ${params.pval_thresh}
       """
}


