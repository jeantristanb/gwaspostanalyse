#!/usr/bin/env python3

import sys
import os
import argparse

def parseArguments():
    parser = argparse.ArgumentParser(description='transform file and header')
    parser.add_argument('--file_rs',type=str,help="input file",required=True)
    parser.add_argument('--around',type=int,help="input file",required=True)
    args = parser.parse_args()
    return args

#rsID,Chro,Pos,A1,A2,Beta,Se,Pval,N,freqA1,direction,Imputed,Sep,File,IsRefFile
args=parseArguments()
dicl={}
readin=open(args.file_rs)
around=args.around
for r in readin:
   rs=r.replace('\n','').split()
   if rs[0] not in dicl :
      dicl[rs[0]]={}
   dicl[rs[0]][rs[1]]=[]

for chro in dicl.keys() :
    dicchro=dicl[chro]
    listposchro=[int(x) for x in dicchro.keys()]
    for pos in dicchro: 
         posint=int(pos)
         dicl[chro][pos]=[str(x) for x in listposchro if abs(x-posint)<=around ]

for chro in dicl.keys() :
    balise=True
    while balise :
      dicchro=dicl[chro]
      listnb=[len(dicl[chro][x]) for x in dicchro]
      listkey=[x for x in dicchro]
      key=listkey[listnb.index(max(listnb))]
      print(chro+":"+key)
      listkey=list(dicl[chro][key])
      for x in listkey :
         if x in  dicl[chro] :
           del dicl[chro][x] 
      if len(dicl[chro].keys())==0:
         balise=False
