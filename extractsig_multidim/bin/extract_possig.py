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
def GetInfo(File):
    Read=open(File)
    ListPos=set([])
    for Line in Read:
       Line=Line.replace('\n','')
       Spl=Line.split('\t')
       ListPos.add(Spl[0]+':'+Spl[1])
    return ListPos

def parseArguments():
    parser = argparse.ArgumentParser(description='transform file and header')
    parser.add_argument('--input_file',type=str,help="input file",required=True)
    parser.add_argument('--head_info',type=str,help="input file",required=True)
    parser.add_argument('--out_file',type=str,help="output file",required=True)
    parser.add_argument("--info_file", help="", type=str,required=True)
    args = parser.parse_args()
    return args

#rsID,Chro,Pos,A1,A2,Beta,Se,Pval,N,freqA1,direction,Imputed,Sep,File,IsRefFile
args=parseArguments()
infohead=args.head_info.split(",")
l_infohead=[x.split(":")[0].upper() for x in infohead]
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

if 'INFO' in l_infohead :
   nameInfo=l_filehead[l_infohead.index('INFO')]
else :
   nameInfo='NA'

if 'TYPE' in l_infohead :
   nameType=l_filehead[l_infohead.index('TYPE')]+'\t'
else :
   nameType=''


if namechro=='NA' or namepos=='NA' :
   sys.exit('not found chro head '+namechro+' or pos head '+namepos)

ListPos=GetInfo(args.info_file)

HeadToSave=['A1','A2','Beta','Se','Pval','N','freqA1']
HeadToSave1=[x.upper() for x in HeadToSave]

ReadFile=open(args.input_file)
Header=ReadFile.readline().replace('\n','').split(sep)
if namechro not in Header or namepos not in Header:
   sys.exit('not found chro head '+namechro+' or pos head '+namepos+' or pval head '+'in file '+args.input_file) 

PosChro=Header.index(namechro)
PosPos=Header.index(namepos)

PosToSave=[PosChro, PosPos]
NewHeader=[]
if nameType!='NA' :
   NewHeader.append('Type')
NewHeader+=['chr', 'bp']
if namersid!='NA':
   PosToSave.append(Header.index(namersid))
   NewHeader.append('rsid')

CmtHead=0
for x in HeadToSave1 :
   if x in l_infohead:
      nametmp=l_filehead[l_infohead.index(x)]
      if nametmp!='NA' :
           PosToSave.append(Header.index(nametmp))
           NewHeader.append(HeadToSave[CmtHead]+' ('+nameInfo+')')
   CmtHead+=1

if 'OTHER' in l_infohead :
   allnameType=l_filehead[l_infohead.index('OTHER')]
   if allnameType.upper()!='NA':
      allnameType=allnameType.split(';')
      for x in allnameType :
        PosToSave.append(Header.index(x))
        NewHeader.append(x)


WriteRes=open(args.out_file, 'w')
WriteRes.write("\t".join(NewHeader)+'\n')
for line in ReadFile :
     Spl=line.split(sep)
     if Spl[PosChro]+':'+Spl[PosPos] in ListPos :
           WriteRes.write(nameType+"\t ".join([Spl[x] for x in PosToSave])+'\n')

