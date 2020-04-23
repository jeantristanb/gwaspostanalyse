#!/usr/bin/env python3
import re
import os
import numpy as np
import argparse

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
   head=read.readline().replace('\n','')
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
      if spl[posbp]!='NA':
        postmp=int(float(spl[posbp]))
        dic_list_inf[spl[poschro]].append([postmp-Wind,postmp+Wind, postmp])
        list_inf[spl[poschro]][postmp]=None
   for chro in dic_list_inf.keys():
      dic_list_inf[chro].sort()
   print("---- end : read pos to search in "+args.list_info+"----")
   return (dic_list_inf,list_inf)

def read_bim_spe(filebim, infopos, infors):
    rbim=open(filebim)
    dbim={}
    listrs=[]
    for line in rbim:
      spline=line.replace('\n','').split()
      if (spline[0] in infopos) and CheckPosInLPos(infopos[spline[0]], int(spline[3])):
         if spline[0] not in dbim:
           dbim[spline[0]]=[]
         dbim[spline[0]].append(int(float(spline[3])))
    return dbim



def parseArguments():
    parser = argparse.ArgumentParser(description='extract annotation for specific position')
    parser.add_argument('--list_info',type=str,required=True, help="file contains rs chro pos used to do windows")
    parser.add_argument('--file_gwas',type=str,required=True, help="file contains GWAS result need result")
    parser.add_argument('--size_win_kb',type=float,required=False,help="windows size in kb", default=100)
    parser.add_argument('--out', type=str,help="header for out", default='test')
    parser.add_argument('--chro_header_info',type=str,help="header of chro in info file", required=True)
    parser.add_argument('--bim_search',type=str,help="header of chro in info file", required=True)
    parser.add_argument('--bp_header_info',type=str,help="header of pval in info files", required=True)
    parser.add_argument('--chro_header_gwas',type=str,help="header of rs in gwas file", required=True)
    parser.add_argument('--rs_header_gwas',type=str,help="header of rs in gwas file", required=True)
    parser.add_argument('--pval_header_gwas',type=str,help="header of pval in gwas file", required=True)
    parser.add_argument('--bp_header_gwas',type=str,help="header of bp in GWAS files", required=True)
    parser.add_argument('--binplinks',type=str,help="binary of plink", default="plink")
    parser.add_argument('--minpval',type=str,help="clump option", default="0.1")
    parser.add_argument('--r2',type=str,help="clump option ", default="0.1")
    parser.add_argument('--cpus',type=str,help="cpus number for plink", default="2")
    parser.add_argument('--bfile',type=str,help="bfile to defined clump in gwas file", required=True)
    parser.add_argument('--keep',type=str,help="option keep of plink ", default=None)
    parser.add_argument('--fam',type=str,help="option keep of plink ", default=None)
    parser.add_argument('--resume',type=int,help="n header in inp files", default=0)
    args = parser.parse_args()
    return args

print("---- begin programs ----")
args = parseArguments()
read_infopos=open(args.list_info)
ChroHeadInf=args.chro_header_info
BpHeadInf=args.bp_header_info
windows_size_kb=args.size_win_kb

filesubgwas=args.out+".sub_gwas"
ResLd=None
## extraction of information position of interest
(infopos,infors)  = read_list_info(args.list_info, ChroHeadInf,BpHeadInf, windows_size_kb*1000)
biminfo=read_bim_spe(args.bim_search, infopos, infors)

infogwas=args.out+".infogwas"
chaineld=[]
for chro in infopos:
  chaineld+=[chro+"\t"+str(x[2])+"\t"+str(x[2])+"\t"+chro+":"+str(x[2]) for x in infopos[chro]]
for chro in biminfo: 
  chaineld+=[chro+"\t"+str(x)+"\t"+str(x)+"\t"+chro+":"+str(x) for x in biminfo[chro]]
winf=open(infogwas, 'w')
winf.write("\n".join(set(chaineld)))
winf.close()
cmd=args.binplinks+" -bfile "+args.bfile+" --threads "+args.cpus+"  --r2 --extract range "+ infogwas+" -out "+args.out+" --maf 0.00001 "+"--ld-window-kb "+str(windows_size_kb) + "   --ld-window-r2 0 --ld-window 99999 "#+" --blocks-min-maf " #+" --ld-snps range "+infogwas)
#print(cmd)
#os.system(cmd)
#os.system(args.binplinks+" -bfile "+args.bfile+" --threads "+args.cpus+"  --r2  -out "+args.out+" --maf 0.00001 "+"--ld-window-kb "+str(windows_size_kb)+" --ld-snps range "+infogwas)

##
infoposchr={}
for chro in infopos:
 infoposchr[chro]=[x[2] for x in infopos[chro]]
print('keys')

readld=open(args.out+'.ld')
tmp=readld.readline()
resld={}
for line in readld :
  line2=re.sub('[ ]+', ' ', line)
  line2=re.sub('^[ ]+', '', line2).replace('\n','')
  spl=line2.split()
  chro=spl[0]
  Pos1=int(spl[1])
  Pos2=int(spl[4])
  ## case where solution
  #if  ((Pos1 in infoposchr[chro]) or (Pos2 in infoposchr[chro]))==False :
  #  print('a '+chro+' '+str(Pos1)+' '+str(Pos2))
  #if  ((Pos1 in biminfo[chro]) or (Pos2 in biminfo[chro])):
  #  print('b '+chro+' '+str(Pos1)+' '+str(Pos2))
  if ((chro in infoposchr) and (chro in biminfo)) and ((Pos1 in infoposchr[chro]) or (Pos2 in infoposchr[chro])) and ((Pos1 in biminfo[chro]) or (Pos2 in biminfo[chro])):
   ##
   R2=float(spl[6])
   if Pos1 in infoposchr[chro]: 
      PosI=Pos1 
      PosB=Pos2
      Rs1=spl[2]
      Rs2=spl[5]
   else :
      PosI=Pos2 
      PosB=Pos1
      Rs2=spl[2]
      Rs1=spl[5]
   if chro not in resld:
     resld[chro]={} 
   if PosI not in resld[chro]:
     resld[chro][PosI]={}
   if PosB in resld[chro][PosI]:
      print(str(PosI)+' '+str(PosB)+' found already')
   resld[chro][PosI][PosB]=[Rs1,Rs2,R2]

#for Chro in resld :
# for Pos in resld[Chro]:
#print(infopos.keys())
#print(resld.keys())
writeall=open(args.out+".res", 'w')
def takesecond(x):
   return x[1]
def takethird(x):
   return x[2]
headfile="chro\tpos\tInBim\tcloses\tclose2\tnearr2\tr2\tnearr22\tr22"
for chro in infopos :
  for posrang in infopos[chro]:
    pos = posrang[2]
    resall=[chro,pos]
    if chro in biminfo and pos in biminfo[chro]:
       IsInBim="T"
    else :
       IsInBim="F"
    resall.append(IsInBim)
    if pos not in resld[chro] :
      print("not found "+' '+chro+' '+str(pos))
      resall+=['NA','NA','NA','NA','NA',"NA"]
    else :
      ##search closest position
      respos=[[x,abs(x-pos), resld[chro][pos][x][2]] for x in resld[chro][pos]]
      respos.sort(key=takesecond)
      cp1="NA"
      cp2="NA"
      if len(respos)>0 :
        cp1=respos[0][0]
      if len(respos)>1 :
        cp2=respos[1][0]
      resall+=[cp1,cp2]
      respos.sort(key=takethird)
      cp1="NA"
      cp2="NA"
      r1="NA"
      r2="NA"
      if len(respos)>0 :
        cp1=respos[0][0]
        r1=respos[0][2]
      if len(respos)>1 :
        cp2=respos[1][0]
        r2=respos[1][2]
      resall+=[cp1,r1,cp2,r2]
    writeall.write("\t".join(resall)+'\n')
## read file bim contains position of chip
#print(infspos)
#print(infors)
## write gwas information
#(infogwas,otherpval)=WriteGWASPlk(args.file_gwas, infopos,fileplkgwas,filesubgwas, args, infors, ResLd)
##Clump

##



