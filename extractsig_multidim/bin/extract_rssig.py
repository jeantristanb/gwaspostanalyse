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
    parser.add_argument('--out_file',type=str,help="output file",required=True)
    parser.add_argument("--info_file", help="", type=str,required=True)
    parser.add_argument("--cut_maf", help="", type=float, default=0.0)
    parser.add_argument("--threshold", type=float,help="threshold to extract rs")
    args = parser.parse_args()
    return args

def CheckFreqNull(x) :
   return True 

def CheckFreq(x) :
   fltfreq=float(x[PosFreq])
   if fltfreq>=cut_maf_min and fltfreq<=cut_maf_max :
     return True
   else :
     return False
#rsID,Chro,Pos,A1,A2,Beta,Se,Pval,N,freqA1,direction,Imputed,Sep,File,IsRefFile
args=parseArguments()
infohead=args.info_file.split(",")
l_infohead=[x.split(":")[0] for x in infohead]
l_filehead=[x.split(":")[1] for x in infohead]
threshold=args.threshold

#nomrs=l_filehead[l_infohead.index('rsID')]
if "SEP" not in l_infohead :
   sep=None
else :
   sep=GetSep(l_filehead[l_infohead.index('SEP')])

if 'CHRO' in l_infohead :
   nomchro=l_filehead[l_infohead.index('CHRO')]
else :
   nomchro='NA'

if 'POS' in l_infohead :
   nompos=l_filehead[l_infohead.index('POS')]
else :
   nompos='NA'

if 'FREQA1' in l_infohead :
   nomfreq=l_filehead[l_infohead.index('FREQA1')]
else :
   nomfreq='NA'
     

if 'PVAL' in l_infohead :
    nompval=l_filehead[l_infohead.index('PVAL')]
else :
    nompval='NA'

if nomchro=='NA' or nompos=='NA' or nompval=='NA' :
   sys.exit('not found chro head '+nomchro+' or pos head '+nompos+' or pval head '+nompval) 
ReadFile=open(args.input_file)
Header=ReadFile.readline().replace('\n','').split(sep)
if nomchro not in Header or nompos not in Header or nompval not in Header :
   print(Header) 
   sys.exit('not found chro head '+nomchro+' or pos head '+nompos+' or pval head '+nompval+' in file '+args.input_file) 
PosChro=Header.index(nomchro)
PosPos=Header.index(nompos)
PosPval=Header.index(nompval)
if args.cut_maf>0 and nomfreq!='NA' :
   PosFreq=Header.index(nomfreq)
   FreqFctCheck=CheckFreq
   cut_maf_min=args.cut_maf
   cut_maf_max=1-args.cut_maf
else :
   FreqFctCheck=CheckFreqNull
   
 
WriteRes=open(args.out_file, 'w')

for line in ReadFile :
     Spl=line.split(sep)
     StrPval=Spl[PosPval]
     if StrPval!='NA' and StrPval!='Na' and StrPval!='na' and StrPval!='nan':
       FloatPval=float(Spl[PosPval])
       if FloatPval <= threshold and FreqFctCheck(Spl): 
          WriteRes.write(Spl[PosChro]+'\t'+Spl[PosPos]+'\n') 
