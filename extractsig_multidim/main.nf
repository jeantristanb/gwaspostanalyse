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
allowed_params = ["cut_maf", "output_dir", "pb_around_rs", "mem_req", "work_dir","mem_req","big_time", "output","nb_cpu" , "input_dir","input_pat", "file_gwas", "pval_thresh", "plot_locuszoom"]
//allowed_params_blocks = ["haploblocks", "plkref_haploblocks", "plk_othopt_haploblocks", 'pval_thresh', "around_rs"]
allowed_params_other=["max_forks", "strandreport", "manifest", "idpat", "accessKey", "access-key", "secretKey", "secret-key","region", "AMI","maxInstances","instance-type", "instanceType", "bootStorageSize", "boot-storage-size", "max-instances", "sharedStorageMount", "shared-storage-mount", "scripts"]
allowed_params_headinfo=["head_chr_gwascat", "head_bp_gwascat", 'gwas_cat']
//allowed_params_head = ["head_pval", "head_freq", "head_bp", "head_chr", "head_rs", "head_beta", "head_se", "head_A1", "head_A2"]
//allowed_params+=allowed_params_head
allowed_params+=allowed_params_other
//allowed_params+=allowed_params_blocks
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
params.cut_maf = 0.00


params.mem_req="8G"
params.big_time="100H"

params.data = ""
params.around_rs=250000


params.max_pval_rep=10**-6
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
//params.head_bp_gwascat="Chro37"
//params.head_chro_gwascat="Pos37"
//params.head_pval_gwascat="P.VALUE"
//params.head_rs_gwascat="SNPS"

params.info_gwascat="DISEASE.TRAIT,REPORTED.GENE.S.,MAPPED_GENE,INITIAL.SAMPLE.SIZE"
params.threshold_pval_gwascat=1
params.pval_thresh=5*10**-8

params.loczm_bin  = ""
params.loczm_pop = "AFR"
params.loczm_build = "hg19"
params.loczm_source ="1000G_March2012"
params.loczm_gwascat = ""
params.plot_locuszoom = 1


 
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
   InfoFile=[]
   TypeFile=[]
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
               InfoFile.add(splLine[cmtelem])
            }
            if(SplH[cmtelem].toUpperCase()=='TYPE'){
               TypeFile.add(splLine[cmtelem])
            }
            if(splLine[cmtelem]!='NA' && splLine[cmtelem].toUpperCase()!='' && SplH[cmtelem].toUpperCase()!='FILE'){
                 SubRes.add(SplH[cmtelem].toUpperCase()+':'+splLine[cmtelem])
            }
            cmtelem+=1
       }
       resInfo.add(SubRes.join(','))
       CmtL+=1
   }
 return([resFile,resInfo, InfoFile, TypeFile])
}

info_file=configfile_analysis(params.file_config)


liste_filesi_ch=Channel.fromPath(info_file[0]).merge(Channel.from(info_file[1])).merge(Channel.from(1..info_file[0].size()))

process ExtractRsSig{
    input :
      set file(file_assoc), val(head_file), num from liste_filesi_ch
    output :
      file(filers) into (listrs_sig)

    script :
       filers=file_assoc+"_${num}.rs"
       """
       extract_rssigv2.py  --input_file $file_assoc --out_file $filers --info_file \"$head_file\" --threshold ${params.pval_thresh} --cut_maf  ${params.cut_maf}
       """
}


listrs_sig_col=listrs_sig.collect()
process MergeSigRs{
   input : 
      file(allfilers) from listrs_sig_col
   output :
      file(allrs) into (allrs_sig,allrs_sig_plot)
   script :
      allrs='allrs.rs'
      allfilersjoin=allfilers.join(" ")
      """
      cat $allfilersjoin|sort|uniq > $allrs
      """
}

process DefineWind{
   input :
       file(allrs) from allrs_sig_plot
   output:
       stdout into listvalposchro
       file('rs_infowind') into rs_infowind
   script :
      """
      define_windsig.py --file_rs $allrs --around ${params.around_rs} --out rs_infowind
      """
}
liste_filesi_ch_2=Channel.fromPath(info_file[0]).merge(Channel.from(info_file[1])).merge(Channel.from(info_file[2])).merge(Channel.from(1..info_file[0].size()))
liste_filesi_ch_3=liste_filesi_ch_2.combine(allrs_sig)

process ExtractInfoSig{
       input :
         set file(file_assoc), val(head_file), val(info_file), val(num),file(allrs) from liste_filesi_ch_3
       output :
         file(fileout) into listres_sig
       script:
        fileout="${file_assoc}_${num}.sub_file"
        """
        extract_possig.py  --input_file $file_assoc --out_file $fileout --head_info \"$head_file\" --info_file $allrs 
        """
}
asso_sig_col=listres_sig.collect()

process MergeInfoSig{
       input :
         file(filesub) from asso_sig_col
       publishDir params.output_dir, overwrite:true, mode:'copy'
       output :
          file(filemerge) into resumcsv
       script : 
         filemerge=params.output+'.csv'
         listfile=filesub.join(',')
         """
         mergefile_sig.r $listfile  ${params.output}
         """ 
}

poschro_ch = Channel.create()
check = Channel.create()
listvalposchro.flatMap { list_str -> list_str.split() }.tap ( check) .set { poschro_ch}
poschro_chF=Channel.fromPath(info_file[0]).merge(Channel.from(info_file[1])).merge(Channel.from(info_file[2])).merge(Channel.from(info_file[3])).combine(poschro_ch)




if(params.plot_locuszoom==1){
process FormatForLocusZoom{
    input :
      set file(file_assoc), val(head_file), val(info_file), val(typefile),val(poschro) from poschro_chF
    output :
      set file(filenamers), val(poschro),file(fileout), val(info_file), val(typefile) into chforgwascat
     script :
     fileout="formatforgwascat.out"
     filenamers="tmp.rs"
     """
     extract_position_forlocuszoom.py --input_file $file_assoc --out_file $fileout --head_info \"$head_file\" --chropos $poschro --around_rs ${params.around_rs} &> $filenamers
     """
}

chforgwascat_F=chforgwascat.filter { it[0] !='NA\n' }


if(params.loczm_gwascat!=""){
loczm_gwascat=" --gwas-cat ${params.loczm_gwascat}"
}else{
loczm_gwascat=""
}

process PlotLocusZoom{
  input :
    set file(filers), val(poschro),file(filegwas), val(info_file), val(typefile) from chforgwascat_F
  publishDir "${params.output_dir}/figure/$outdir", overwrite:true, mode:'copy'
     errorStrategy 'ignore'
  output:
   file("${out}.svg")
   file("${out}.pdf" )  into locuszoom_type
  script :
     //filetmp=filegwas.baseName.toString()
     //outdir=filetmp.take(filetmp.toString().lastIndexOf('.'))
     //rs=filers.readLines()[0]
     outdir=filegwas.baseName
     //rs=rs.replace('\n','')
     poschro=poschro.replace(':','_')
     out=outdir+'_'+poschro+'_'+info_file+"_"+typefile
 """
   rs=`head -1 $filers|awk '{print \$1}'`
   ${params.loczm_bin} --epacts  $filegwas --delim tab --refsnp  \$rs --flank ${params.around_rs} --pop ${params.loczm_pop} --build ${params.loczm_build} --source ${params.loczm_source} $loczm_gwascat --svg  -p out --no-date
  mv */*.svg $out".svg"
  mv */*.pdf $out".pdf"
 """
}
locuszoom_col=locuszoom_type.collect()
}else{
locuszoom_col=file('FILE_NO')
}
process DoReport{
  input :
    file(locuszoom) from locuszoom_col     
    file(csv) from resumcsv
    file(rswind) from rs_infowind
  publishDir "${params.output_dir}/", overwrite:true, mode:'copy'
  output :
    set file("$outpdf"), file("$outtex"), file("$outcsv")
  script :
   lzm=params.plot_locuszoom==1 ? locuszoom.join('\n') : ""
   tmpfile=params.work_dir+'/.tempfile'
   File writ = new File("${params.work_dir}/.tempfile")
   writ.write(lzm)
   outpdf="${params.output}.pdf"
   outtex="${params.output}.tex"
   outcsv="${params.output}_short.csv"
   """
   launch_doreport.r --list_pdf $tmpfile --csv_res $csv --rs_info $rswind --csv_out  $outcsv
   mv do_report.pdf $outpdf
   mv do_report.tex $outtex
   """
}
