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
    parser.add_argument('--list_info',type=str,required=True, help="file contains rs chro pos used to do windows")
    parser.add_argument('--out', type=str,help="header for out")
    parser.add_argument('--wind_size',type=float,required=True,help="windows size in kb")
    parser.add_argument('--bin_plink',type=str,required=False,help="windows size in kb", default="plink")
    parser.add_argument('--keep',type=str,required=False,help="windows size in kb")
    parser.add_argument('--bfile',type=str,help="bfile to defined clump in gwas file", required=True)
    parser.add_argument('--cpus',type=str,help="cpus for plink", required=False, default=1)
    args = parser.parse_args()
    return args

args=parseArguments()


ChroHeadInf=args.chro_header_info
BpHeadInf=args.bp_header_info
window=args.wind_size

filebed = write_list_info(args.list_info, ChroHeadInf,BpHeadInf, window*1000, args.out)
#plink -bfile $bfileI --keep $FileInd --extract range $FileCat".range" --make-bed --out $bfile  --keep-allele-order &

Cmd=args.bin_plink+" -bfile "+args.bfile+" --extract range " +filebed+ " --make-bed --out "+args.out+ " --keep-allele-order --threads "+args.cpus

if args.keep :
 Cmd+=" --keep "+args.keep

os.system(Cmd)


