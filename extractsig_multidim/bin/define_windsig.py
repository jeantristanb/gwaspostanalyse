#!/usr/bin/env python3

import sys
import os
import argparse

def parseArguments():
    parser = argparse.ArgumentParser(description='transform file and header')
    parser.add_argument('--file_rs',type=str,help="input file",required=True)
    parser.add_argument('--around',type=int,help="input file",required=True)
    parser.add_argument('--out',type=int,help="output file",required=True)
    args = parser.parse_args()
    return args

#rsID,Chro,Pos,A1,A2,Beta,Se,Pval,N,freqA1,direction,Imputed,Sep,File,IsRefFile
args=parseArguments()
dicl={}
dicchr={}
readin=open(args.file_rs)
around=args.around
for r in readin:
   rs=r.replace('\n','').split()
   ## contains all pos
   if rs[0] not in dicchr :
      dicchr[rs[0]]=set([])
   dicchr[rs[0]].add(int(rs[1]))
   if rs[2]=='1' :
      if rs[0] not in dicl :
         dicl[rs[0]]={}
      dicl[rs[0]][rs[1]]=[]

for chro in dicl.keys() :
    dicchro=dicl[chro]
    listposchro=list(dicchr[chro])
    for pos in dicchro: 
         posint=int(pos)
         dicl[chro][pos]=[str(x) for x in listposchro if abs(x-posint)<=around ]

Write=open(args.out, 'w')
for chro in dicl.keys() :
    balise=True
    dicchro=dicl[chro]
    ## we contrasted pos => del all pos where that at least another position not in save (alteady significant in other dataset)       
    ## listpos good in dic
    listkey=[x for x in dicchro] 
    for pos in listkey :
        nbnotgood=len([x for x in dicchro[pos] if x not in listkey])
        if nbnotgood>0:
            del dicchro[pos]
    while balise :
      listkey=[x for x in dicchro] 
      listnb=[len(dicl[chro][x]) for x in dicchro]
      ## listpos good in dic
      listkey=[x for x in dicchro] 
      key=listkey[listnb.index(max(listnb))]
      print(chro+":"+key)
      Write.write("\n".join([chro+'\t'+str(x) for x in dicchro[key]])+'\n')
      listkey=list(dicl[chro][key])
      minpos=min(dicl[chro][key])
      maxpos=max(dicl[chro][key])
      for x in listkey :
         if x in  dicl[chro] :
           del dicl[chro][x] 
      if len(dicl[chro].keys())==0:
         balise=False
