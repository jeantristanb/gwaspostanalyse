#!/usr/bin/env python3

import gzip
import sys
import argparse

def spl_gz(z,sep=None) :
   return z.decode('utf-8').replace('\n','').split(sep)

def spl_nogz(z,sep=None) :
   return z.replace('\n','').split(sep)


def openf(File) :
  def is_gz_file(filepath):
      with open(filepath, 'rb') as test_f:
        return test_f.read(2) == b'\x1f\x8b'
  def checkexists(path_to_file) :
    if exists(path_to_file)==False :
     sys.exist('file '+ path_to_file+' doesn t exist')
  balisegz=False
  if is_gz_file(File) :
    readf=gzip.open(File)
    balisegz=True
  else :
    readf=open(File)
  spl=spl_nogz
  if balisegz:
     spl=spl_gz
  return (readf, balisegz,spl)


def parseArguments():
    parser = argparse.ArgumentParser(description='extract list chro of list file')
    parser.add_argument('--list_file',type=str,required=True)
    parser.add_argument('--head_chr',type=str,required=True,help="File with phenotype and covariate data")
    args = parser.parse_args()
    return args

args = parseArguments()

listfile=args.list_file.split(',')
ListChro=set([])

for File in  listfile :
       (ReadGWAS, balisegz, spl)=openf(File)
       for line in ReadGWAS :
            Head2=spl(ReadGWAS.readline())#.replace('\n','').split()
            if args.head_chr in Head2 :
               PosHeadChr=Head2.index(args.head_chr)
               nelem=len(Head2)
               break
       for lireGWAS in ReadGWAS :
           tmp=spl(lireGWAS)
           if len(tmp) == nelem and (args.head_chr not in tmp):
              ListChro.add(tmp[PosHeadChr])
       ReadGWAS.close()

for x in ListChro :
   print(x)

