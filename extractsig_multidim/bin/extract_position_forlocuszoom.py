#!/usr/bin/env python3

import sys
import os
import argparse

def GetSep(Sep):
   ListOfSep=["\t"," ",","]
   if len(Sep)>2 :
      Sep=Sep.upper()[:3]
      if Sep=='COM' :
         Sep=','
      elif Sep=='TAB' :
         Sep='\t'
      elif Sep=='WHI' :
         Sep=' '
   if Sep not in ListOfSep :
      return None
   return Sep


#extract_rssig.py  --input_file $file_assoc --out_file $newfile_assoc --info_file $info_file --threshold ${params.pval_thresh}

def parseArguments():
    parser = argparse.ArgumentParser(description='transform file and header')
    parser.add_argument('--input_file',type=str,help="input file",required=True)
    parser.add_argument('--head_info',type=str,help="input file",required=True)
#--chropos $poschro --around_rs ${params.around_rs}
    parser.add_argument('--out_file',type=str,help="output file",required=True)
    parser.add_argument("--chropos", help="", type=str,required=True)
    parser.add_argument("--around_rs", help="", type=int,required=True)
    args = parser.parse_args()
    return args

#rsID,Chro,Pos,A1,A2,Beta,Se,Pval,N,freqA1,direction,Imputed,Sep,File,IsRefFile
args=parseArguments()
infohead=args.head_info.split(",")
l_infohead=[x.split(":")[0] for x in infohead]
l_filehead=[x.split(":")[1] for x in infohead]



#namers=l_filehead[l_infohead.index('rsID')]
if "SEP" not in l_infohead :
   sep=None
else :
   sep=GetSep(l_filehead[l_infohead.index('SEP')])

if 'CHRO' in l_infohead :
   namechro=l_filehead[l_infohead.index('CHRO')]
else :
   namechro='NA'

if 'POS' in l_infohead :
   namepos=l_filehead[l_infohead.index('POS')]
else :
   namepos='NA'

if 'RSID' in l_infohead:
   namersid=l_filehead[l_infohead.index('RSID')]
else :
   namersid='NA'
#rsID,Chro,Pos,A1,A2,Beta,Se,Pval,N,freqA1,direction,Imputed,Sep,File,Info,Type

if namechro=='NA' or namepos=='NA' :
   sys.exit('not found chro head '+namechro+' or pos head '+namepos)

ChroPos=args.chropos
ChroPos2=ChroPos.split(':')
ChroCheck=ChroPos2[0]
PosCheck=int(ChroPos2[1])
PosMinCheck=PosCheck-args.around_rs
PosMaxCheck=PosCheck+args.around_rs

HeadToSave=['Pval','freqA1']
NewHead=["PVALUE","MAF"]

HeadToSave1=[x.upper() for x in HeadToSave]

ReadFile=open(args.input_file)
Header=ReadFile.readline().replace('\n','').split(sep)
if namechro not in Header or namepos not in Header:
   sys.exit('not found chro head '+namechro+' or pos head '+namepos+' or pval head '+'in file '+args.input_file) 

PosChro=Header.index(namechro)
PosPos=Header.index(namepos)
if namersid!='NA':
  PosRsID=Header.index(namersid)


PosToSave=[]
NewHeader=[]
#chropo
#       small.columns=["#CHROM","BEGIN","END","MARKER_ID","PVALUE"]
NewHeader+=['#CHROM', 'BEGIN', "END", "MARKER_ID"]

CmtHead=0
for x in HeadToSave1 :
   if x in l_infohead:
      nametmp=l_filehead[l_infohead.index(x)]
      if nametmp!='NA' :
           PosToSave.append(Header.index(nametmp))
           NewHeader.append(NewHead[CmtHead])
   CmtHead+=1


WriteRes=open(args.out_file, 'w')
WriteRes.write("\t".join(NewHeader)+'\n')
MinDistPos=1000000000
RsToUse='NA'
for line in ReadFile :
     Spl=line.split(sep)
     PosLine=int(Spl[PosPos])
     if ChroCheck==Spl[PosChro] and PosLine>=PosMinCheck and PosLine<=PosMaxCheck:
          if namersid=='NA' :
            RSID=Chro+':'+str(PosLine)
          else :
            RSID=Spl[PosRsID]
          if abs(PosCheck-PosLine)<MinDistPos :
            MinDistPos=abs(PosCheck-PosLine)
            RsToUse = RSID
          WriteRes.write(Spl[PosChro]+'\t'+Spl[PosPos]+'\t'+Spl[PosPos]+'\t'+ RSID+'\t'+"\t ".join([Spl[x] for x in PosToSave])+'\n')
print(RsToUse)
