#!/usr/bin/env python3

import sys
import pandas as pd
import argparse
import numpy as np

def parseArguments():
    parser = argparse.ArgumentParser(description='extract list chro of list file')
    parser.add_argument('--list_file',type=str,required=True)
    parser.add_argument('--head_chr',type=str,required=True,help="File with phenotype and covariate data")
    args = parser.parse_args()
    return args

args = parseArguments()

ReadFile=open(args.list_file)
Head=ReadFile.readline().replace('\n','').split()
PosFile=Head.index("File")
#echo -e "Pheno\tType\tFile"> list_file.gwas
PosPheno=Head.index("Pheno")
ListChro=set([])
for line in ReadFile :
   SplL=line.replace("\n","").split()
   if len(SplL)>0 :
       ReadGWAS=open(SplL[PosFile])
       Head2=ReadGWAS.readline().replace('\n','').split()
       PosHeadChr=Head2.index(args.head_chr)
       for lireGWAS in ReadGWAS :
           ListChro.add(lireGWAS.split()[PosHeadChr])
       ReadGWAS.close()

for x in ListChro :
   print(x)

