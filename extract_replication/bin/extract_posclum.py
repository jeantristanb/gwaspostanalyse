#!/usr/bin/env python3
import sys
import os
import pandas as pd
import argparse
import numpy as np
### Extract of file Chromosome, Position, Added Windows
### dic_list_inf : retur dictonnary by chromosome  with widos
### list_inf :return a dictionnary to save result
def CheckPosInLPos(listpos, pos):
   for lpos in listpos :
     if pos>=lpos[0] and pos<=lpos[1] :
        return True
   return False

def read_list_info(file_inf, ChroHeadInf,BpHeadInf, Wind):
   print("---- begin : read pos to search in "+args.list_info+"----")
   read=open(file_inf)
   head=read.readline()
   spl=head.replace('"','').split('\t')
   poschro=spl.index(ChroHeadInf)
   posbp=spl.index(BpHeadInf)
   dic_list_inf={}
   list_inf={}
   for line in read :
      spl=line.replace('"','').split('\t')
      if spl[poschro] not in dic_list_inf :
         dic_list_inf[spl[poschro]]=[]
         list_inf[spl[poschro]]={}
      dic_list_inf[spl[poschro]].append([int(spl[posbp])-Wind, int(spl[posbp])+Wind, int(spl[posbp])])
      list_inf[spl[poschro]][int(spl[posbp])]=None
   for chro in dic_list_inf.keys():
      dic_list_inf[chro].sort()
   print("---- end : read pos to search in "+args.list_info+"----")
   return (dic_list_inf,list_inf)
      
def WriteGWASPlk(file_gwas,infopos, fileplkgwas, filesubgwas,args, infors, ResLD=None) :
   print("---- begin : format gwas for plk ----")
   readgwas=open(file_gwas)
   headgwas=readgwas.readline()
   splhead=headgwas.replace('\n','').split()
   ChroPosGwas=splhead.index(args.chro_header_gwas)
   BpPosGwas=splhead.index(args.bp_header_gwas)
   RsPosGwas=splhead.index(args.rs_header_gwas)
   PvPosGwas=splhead.index(args.pval_header_gwas)
   cmtline=0
   writegwasplk=open(fileplkgwas,'w')
   writegwassub=open(filesubgwas,'w')
   writegwasplk.write("SNP\tP\n")
   listchro=infopos.keys()
   infogwas={}
   writegwas_linfo=open(args.out+'.in.list_info','w')
   writegwas_linfo.write(headgwas)
   writegwassub.write(headgwas)
   for line in readgwas:
     spl=line.split()
     chro=spl[ChroPosGwas]
     bp=int(spl[BpPosGwas])
     rs=spl[RsPosGwas]
     Pv=spl[PvPosGwas]
     if (spl[ChroPosGwas] in listchro) :
       if bp in infors[chro] :
          infors[chro][bp]=[rs,spl[PvPosGwas]]
          writegwas_linfo.write(line)
       if ResLD :
         if chro not in infogwas :
            infogwas[chro]={}
         if CheckPosInLPos(infopos[chro], float(bp)) or CheckPosInLPos(ResLD[chro], float(bp)): 
           infogwas[chro][bp]=[rs,Pv]
           writegwasplk.write(rs+"\t"+Pv+"\n")
           writegwassub.write(line)
           cmtline+=1
       elif CheckPosInLPos(infopos[spl[ChroPosGwas]], float(bp)) :
         if chro not in infogwas :
            infogwas[chro]={}
         infogwas[chro][bp]=[rs,Pv]
         writegwasplk.write(spl[RsPosGwas]+"\t"+Pv+"\n")
         writegwassub.write(line)
         cmtline+=1
   readgwas.close()
   writegwasplk.close()
   writegwassub.close()
   print("---- end : format gwas for plk ----")
   return infogwas

## What done : used position to analyse, to search in data_clump what is position most closest
## write final 
## input :
### info_pos : list of pos to anlyse
### data_clum : data have been clumped
### window : windows around each position 
###

def GetResByPos(infopos,data_clump, windows_size_kb, out) :
   print("---- begin :merge clump and dataI by pos ----")
   Cmt=0
   listpos=[]
   listposi=[]
   listchro=[]
   for chro in infopos.keys() :
     data_clump_chr=data_clump[data_clump['CHR']==str(chro)]
     for info in infopos[chro] :
      pos=info[2]
      lmin=[abs(x-pos) for x in data_clump_chr['BP']]
      if len(lmin)>0 and min(lmin)<windows_size_kb*1000 :
        minlmin=min(lmin)
        posmin=lmin.index(minlmin)
        Tmp=data_clump_chr.iloc[posmin]
        listpos.append(pos)
        listchro.append(chro)
        listposi.append(Tmp['BP'])
        Cmt+=1
   dat2=pd.DataFrame({'CHR':listchro,'BP':listposi,'BPI':listpos})
   dat2=dat2.merge(data_clump)
   dat2.to_csv(out+".clump.bypos",sep="\t",header=True,index=False,na_rep="NA")
   print("---- end :merge clump and dataI by pos ----")


def GetResByRs(infopos,data_clump, windows_size_kb, infors,out):
  cmtpos=0
  print("---- begin :merge clump and dataI by rs ----")
  listpos=[]
  listposi=[]
  listchro=[]
  listrsi=[]
  for chro in infopos.keys() :
     data_clump_chr=data_clump[data_clump['CHR']==str(chro)]
     for info in infopos[chro] :
        pos=info[2]
        if infors[chro][pos] :
           rs=infors[chro][pos][0]+"("
           cmtpos=0
           found=False
           for x in data_clump_chr['SP2'] :
             if rs in x :
               Tmp=data_clump_chr.iloc[cmtpos]
               listpos.append(pos)
               listchro.append(chro)
               listposi.append(Tmp['BP'])
               listrsi.append(infors[chro][pos][0])
               found =True
             cmtpos+=1
           if found ==False :
             print('not found '+rs+' for pos '+str(chro)+' '+str(pos))
        else :
           print('not found rs '+str(chro)+' '+str(pos))
  dat2=pd.DataFrame({'CHR':listchro,'BP':listposi,'BPI':listpos, 'RSI':listrsi})
  dat2=dat2.merge(data_clump)
  dat2.to_csv(args.out+".clump.byrs",sep="\t",header=True,index=False,na_rep="NA")
  print("---- end :merge clump and dataI by rs ----")

###
def GetLDInfo(infopos, ldfile):
   readfile=open(ldfile)
   header=readfile.readline()
   dicldwind={}
   for line in readfile :
      spl=line.replace('\n','').split()
      chro=spl[0]
      posbeg=int(spl[1])
      posend=int(spl[2])
      if chro in infopos :
        listpos=infopos[chro]
        #windld=[posbeg,posend,'']
        posgood=''
        listpos_a=[]
        for pos in listpos:
            if pos[2]>=posbeg and pos[2]<=posend : 
               posgood+=str(pos[2])+';' 
               listpos_a.append(pos)
        if posgood!='':
           if chro not in dicldwind : 
              dicldwind[chro]=[]
           infold=[posbeg,posend,posgood, listpos_a]
           #pos[3]=infold
           dicldwind[chro].append(infold)
   return dicldwind 
def GetResByLDWind(ld_info, data_clump, infors,out):
    ## we extracted all ld_block with info_pos
  cmtpos=0
  print("---- begin :merge clump and dataI by  ----")
  Header="Chro\tBeginBlock\tEndBlock\tBPClump\tPClump\tRsClump\tBPInfo"
  Write=open(out+".clump.ldbloc.detail", 'w')
  Write.write(Header+'\n')
  for chro in ld_info.keys() :
     data_clump_chr=data_clump[data_clump['CHR']==str(chro)]
     for info in ld_info[chro] :
        posbegin=int(info[0])
        posend=int(info[1])
        for index, row in data_clump_chr.iterrows():
          if int(row['BP'])>= posbegin and int(row['BP'])<=posend :
            chaine=chro+"\t"+str(posbegin)+"\t"+str(posend)+"\t"+str(row['BP'])+"\t"+str(row['P'])+"\t"+infors[chro][row['BP']][0]+"\t"
            AllChaine=""
            for Snp in info[3] :
               AllChaine+=chaine+"\t"+str(Snp[2])+"\n"
            Write.write(AllChaine) 
  


 
 
def parseArguments():
    parser = argparse.ArgumentParser(description='extract annotation for specific position')
    parser.add_argument('--list_info',type=str,required=True, help="file contains rs chro pos used to do windows")
    parser.add_argument('--file_gwas',type=str,required=True, help="file contains GWAS result need result")
    parser.add_argument('--size_win_kb',type=float,required=True,help="windows size in kb")
    parser.add_argument('--out', type=str,help="header for out")
    parser.add_argument('--chro_header_info',type=str,help="header of chro in info file", required=True)
    parser.add_argument('--bp_header_info',type=str,help="header of pval in info files", required=True)
    
    parser.add_argument('--chro_header_gwas',type=str,help="header of rs in gwas file", required=True)
    parser.add_argument('--rs_header_gwas',type=str,help="header of rs in gwas file", required=True)
    parser.add_argument('--pval_header_gwas',type=str,help="header of pval in gwas file", required=True)
    parser.add_argument('--bp_header_gwas',type=str,help="header of bp in GWAS files", required=True)
    parser.add_argument('--binplinks',type=str,help="binary of plink", default="plink")
    parser.add_argument('--minpval',type=str,help="clump option", default="0.1")
    parser.add_argument('--r2',type=str,help="clump option ", default="0.1")
    parser.add_argument('--cpus',type=str,help="clump option ", default="2")
    parser.add_argument('--bfile',type=str,help="bfile to defined clump in gwas file", required=True)
    parser.add_argument('--keep',type=str,help="option keep of plink ", default=None)
    parser.add_argument('--file_block_ld', type=str, help="block ld plink file obtained with output --block : .blocks.det", default=None, required=True)
    parser.add_argument('--resume',type=int,help="n header in inp files", default=0)
    args = parser.parse_args()
    return args

print("---- begin programs ----")
args = parseArguments()
read_infopos=open(args.list_info)
ChroHeadInf=args.chro_header_info
BpHeadInf=args.bp_header_info
windows_size_kb=args.size_win_kb
if args.resume == 0 :
   renew=True
else :
   renew = False
clean=False
### 
(infopos,infors)  = read_list_info(args.list_info, ChroHeadInf,BpHeadInf, windows_size_kb*1000)

### format for plink
fileplkgwas=args.out+".formarplk"
filesubgwas=args.out+".sub_gwas"
#if os.path.isfile(fileplkgwas)==False or renew:
if args.file_block_ld :
   ResLd=GetLDInfo(infopos, args.file_block_ld)    
else :
   ResLd=None
infogwas=WriteGWASPlk(args.file_gwas, infopos,fileplkgwas,filesubgwas, args, infors, ResLd) 

fileoutpltmp=args.out+".plk"
fileoutpl=fileoutpltmp+".clumped"
keep=""
if args.keep :
   " --keep "+args.keep
if os.path.isfile(fileoutpl)==False or renew:
   print("---- plink launch ----")
   os.system(args.binplinks+" -bfile "+args.bfile+" --clump "+fileplkgwas+" --clump-p1 "+ args.minpval+" --clump-p2 1 "+" --clump-kb "+str(windows_size_kb)+"  --clump-r2 "+args.r2+" --out "+fileoutpltmp+" "+keep+" --threads "+args.cpus)
   print("---- end plink ----")
else :
    print("not plink using file "+fileoutpl)

if os.path.isfile(fileoutpl)==False: 
   sys.exit("doesn't found "+fileoutpl+"\t"+"exit")

## read clum result
data_clump=pd.read_csv(fileoutpl,delim_whitespace=True)
data_clump['CHR']=data_clump['CHR'].astype(str)


GetResByPos(infopos,data_clump, windows_size_kb, args.out)
GetResByRs(infopos,data_clump, windows_size_kb, infors,args.out)

if args.file_block_ld :
   GetResByLDWind(ResLd, data_clump, infogwas,args.out)






#if clean==True:
#  os.remove(fileplkgwas)
#  os.remove(fileoutpl)
