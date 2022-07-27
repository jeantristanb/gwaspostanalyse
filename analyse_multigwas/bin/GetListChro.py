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

listfile=args.list_file.split(',')
ListChro=set([])
for filea in  listfile :
 ReadGWAS=open(filea)
 Head2=ReadGWAS.readline().replace('\n','').split()
 PosHeadChr=Head2.index(args.head_chr)
 for lireGWAS in ReadGWAS :
    ListChro.add(lireGWAS.split()[PosHeadChr])
 ReadGWAS.close()

for x in ListChro :
   print(x)

