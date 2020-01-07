#!/usr/bin/env python3

import sys
import pandas as pd
import argparse
import numpy as np

def parseArguments():
    parser = argparse.ArgumentParser(description='extract list chro of list file')
    parser.add_argument('--list_file',type=str,required=True)
    parser.add_argument('--head_chr',type=str,required=True,help="File with phenotype and covariate data")
    parser.add_argument('--head_bp',type=str,required=True,help="File with phenotype and covariate data")
    parser.add_argument('--head_freq',type=str,required=True,help="File with phenotype and covariate data")
    parser.add_argument('--head_rs',type=str,required=True,help="File with phenotype and covariate data")
    parser.add_argument('--head_pv',type=str,required=True,help="File with phenotype and covariate data")
    parser.add_argument('--pval_max',type=float,required=True,help="File with phenotype and covariate data")
    parser.add_argument('--freq_min',type=float,required=True,help="File with phenotype and covariate data")
    parser.add_argument('--pheno', type=str,required=True,help="File with phenotype and covariate data")
    args = parser.parse_args()
    return args

def ReadFile2(File,  HeadFreq, MinFreq):
    Data=pd.read_csv(File,delim_whitespace=True)
    Data=Data[(Data[HeadFreq]>MinFreq) & (Data[HeadFreq]<(1-MinFreq))] 
    return(Data)

def GetWindows() :



#File="/dataE/AWIGenGWAS/shared/ResultGWAS/Imputed/Boltlmm/Res/pulse_average_qc/pulse-average-qc_West_age-sex_231018.imp.bolt"
#SNP	CHR	BP	GENPOS	ALLELE1	ALLELE0	A1FREQ	INFO	BETA	SE	P_BOLT_LMM_INF	P_BOLT_LMM
#Data=ReadFile2(File, "A1FREQ", 0.01)
#Tmp=Data.loc[DataSig,]
Test=True
if Test :
  list_file="/home/jeantristan/Travail/GWAS/VImputed2/TransfertAnaly/list_file.gwas"
  head_freq="A1FREQ"
  head_rs="SNP"
  head_chr="CHR"
  head_bp="BP"
  head_pv="P_BOLT_LMM"
  freq_min=0.01
  pval_max = 10**-5
  head_beta="BETA"
  head_a1="ALLELE1"
  head_a2="ALLELE0"
  pheno="bp_dia_average_qc"
else :
  args = parseArguments()
  list_file = args.list_file
  head_freq=args.head_freq
  head_rs=args.head_rs
  head_chr=args.head_chr
  head_pv=args.head_pv
  head_beta=args.head_beta
  freq_min=args.freq_min
  pval_max=args.pval_max
  pheno=args.pheno

ReadFile=open(list_file)
Head=ReadFile.readline().replace('\n','').split()
PosFile=Head.index("File")
PosPheno=Head.index("Pheno")
PosType=Head.index("Type")
ListData={}
ListSig=[]

for line in ReadFile :
   SplL=line.replace("\n","").split()
   if SplL[PosPheno]==pheno  :
      Data=ReadFile2(SplL[PosFile],  head_freq, freq_min)
      #Data=Data.set_index(head_rs)
      ListData[SplL[PosType]]=Data
      ListSig+=Data[Data[head_pv]<pval_max][head_rs].tolist()

ListSig=list(set(ListSig))
## merge data by
Cmt=1
for MD in ListData.keys() :
   Data=ListData[MD][ListData[MD][head_rs].isin(ListSig)]
   Data=Data[[head_rs,head_chr,head_bp,head_a1,head_a2,head_freq,head_pv,head_beta]]
   Data=Data.rename(columns={head_freq: "freq_"+MD, head_pv: 'Pv_'+MD, head_beta:'B_'+MD})  
   if Cmt==1 :
      DataF=Data#pd.merge(left, right, how='inner'
   else :
      DataF=pd.merge(DataF,Data, how='outer')
   Cmt+=1







