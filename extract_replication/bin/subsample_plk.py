#!/usr/bin/env python3

import sys
import os
import argparse

def write_list_info(file_inf, ChroHeadInf,BpHeadInf, Wind, Out):
   print("---- begin : read pos to search in "+args.list_info+"----")
   read=open(file_inf)
   write=open(Out+".tmp",'w')
   head=read.readline()
   spl=head.replace('\n','').replace('"','').split('\t')
   poschro=spl.index(ChroHeadInf)
   posbp=spl.index(BpHeadInf)
   dic_list_inf={}
   list_inf={}
   for line in read :
      spl=line.split('\t')
      newbegin=int(spl[posbp])-Wind-100
      if newbegin<1 :
         newbegin=1
      res="\t".join([str(x).replace('"','') for x in [spl[poschro],int(newbegin),int(int(spl[posbp])+Wind+100), spl[poschro]+":"+spl[posbp]]])
      write.write(res+'\n')
   print("---- end : read pos to search in "+args.list_info+"----")
   write.close()
   return Out+".tmp"

def parseArguments():
    parser = argparse.ArgumentParser(description='extract annotation for specific position')
    parser.add_argument('--chro_header_info',type=str,help="header of chro in info file", required=True)
    parser.add_argument('--bp_header_info',type=str,help="header of pval in info files", required=True)
    parser.add_argument('--list_info',type=str,required=True, help="file contains rs chro pos used to do size_win_kbs")
    parser.add_argument('--out', type=str,help="header for out")
    parser.add_argument('--size_win_kb',type=float,required=True,help="size_win_kbs size in kb")
    parser.add_argument('--bin_plink',type=str,required=False,help="size_win_kbs size in kb", default="plink")
    parser.add_argument('--keep',type=str,required=False,help="size_win_kbs size in kb")
    parser.add_argument('--bfile',type=str,help="bfile to defined clump in gwas file", required=True)
    parser.add_argument('--cpus',type=str,help="cpus for plink", required=False, default=1)
    parser.add_argument('--file_gwas', type=str, help='gwas file for change rs with chr and bp')
    parser.add_argument('--chro_header_gwas', type=str, help='gwas file for change rs with chr and bp')
    parser.add_argument('--bp_header_gwas', type=str, help='gwas file for change rs with chr and bp')
    parser.add_argument('--rs_header_gwas', type=str, help='gwas file for change rs with chr and bp')
    parser.add_argument('--a1_header_gwas', type=str, help='gwas file for change rs with chr and bp')
    parser.add_argument('--a0_header_gwas', type=str, help='gwas file for change rs with chr and bp')
    args = parser.parse_args()
    args = parser.parse_args()
    return args

args=parseArguments()


ChroHeadInf=args.chro_header_info
BpHeadInf=args.bp_header_info
size_win_kb=args.size_win_kb

filebed = write_list_info(args.list_info, ChroHeadInf,BpHeadInf, size_win_kb*1000, args.out)
#plink -bfile $bfileI --keep $FileInd --extract range $FileCat".range" --make-bed --out $bfile  --keep-allele-order &

Cmd=args.bin_plink+" -bfile "+args.bfile+" --extract range " +filebed+ " --make-bed --out "+args.out+ ".tmp --keep-allele-order --threads "+args.cpus

if args.keep :
 Cmd+=" --keep "+args.keep

os.system(Cmd)

#--chro_header_gwas ${params.head_chr}  --bp_header_gwas ${params.head_bp} --rs_header_gwas ${params.head_rs}  --file_gwas $gwas
## read bim file to define chr, rs
readbim=open(args.out+'.tmp.bim')
dicbim={}
#1	rs555500075	0	10352	TA	T
for line in readbim :
    splline=line.split()
    if splline[0] not in dicbim :
      dicbim[splline[0]]={} 
    dicbim[splline[0]][splline[3]]=[splline[1],splline[4],splline[5]] 
readbim.close()
filerseq=open('tmp_eq.rs','w')
#filerseq.write('NewID\tOldID\n')
readgwas=open(args.file_gwas)
head=readgwas.readline().split()
posrs=head.index(args.rs_header_gwas)
posbp=head.index(args.bp_header_gwas)
poschro=head.index(args.chro_header_gwas)
if args.a1_header_gwas and args.a0_header_gwas:
  nbchange=0
  posa0=head.index(args.a0_header_gwas)
  posa1=head.index(args.a1_header_gwas)
  for line in  readgwas :
     spl=line.split()
     if (spl[poschro] in dicbim) and (spl[posbp] in dicbim[spl[poschro]]) :
      info=dicbim[spl[poschro]][spl[posbp]]
      if (spl[posa0]==info[1] and spl[posa1]==info[2]) or (spl[posa0]==info[2] and spl[posa1]==info[1]) :
        if info[0]!=spl[posrs] :
          filerseq.write(spl[posrs]+'\t'+info[0]+'\n')
          nbchange+=1
else :
  nbchange=0
  for line in  readgwas :
     spl=line.split()
     if (spl[poschro] in dicbim) and (spl[posbp] in dicbim[spl[poschro]]) :
      info=dicbim[spl[poschro]][spl[posbp]]
      if info[0]!=spl[posrs] :
         filerseq.write(spl[posrs]+'\t'+info[0]+'\n') 
         nbchange+=1
filerseq.close()

if nbchange>0 :
     Cmd=args.bin_plink+" -bfile "+args.out+".tmp --make-bed --out "+args.out+ " --keep-allele-order --threads "+args.cpus+" --update-name  tmp_eq.rs 1 2   "
     os.system(Cmd)
else :
    os.system('mv '+args.out+'.tmp.bim '+ args.out+'.bim') 
    os.system('mv '+args.out+'.tmp.bed '+ args.out+'.bed') 
    os.system('mv '+args.out+'.tmp.fam '+ args.out+'.fam') 

