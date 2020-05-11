#!/usr/bin/env python3

import sys
import os
import argparse

## algorithm: selected on file position 
def ExtractPosFile(File, headchro, headbp, headpval,threshold, sep=None):
    listchro={}
    readfile=open(File)
    headli=readfile.readline()
    spline=headlin.split(sep)
    poschro=spline.index(headchro)
    posbp=spline.index(headbp)
    pospval=spline.index(headpval)
    for line in readfile :
      spline=line.split(sep)
      if float(spline[pospval]) <= float(threshold) :
          pos=int(spline[posbp])
          if spline[poschro] not in listchro:
             listchro[spline[poschro]]={}
          listchro[spline[poschro]][pos]=spline
    readfile.close()
    return listchro

def definewind(lispos, limwind):
  lispos=list(set(listpos))
  lispos.sort()
  beginwind=max(min(listpos)-limwind,1)
  endwind=min(listpos)+limwind
  listwind=[]
  for x in listpos:
    if x>=endwind:  
      listwind.append([beginwind,endwind, listposwind]) 
      beginwind=x-1
      endwind=x+limwind
      listposwind=[]
    else :
      endwind+=limwind
    listposwind.append(x)
  if len(listposwind)>0 :
      listwind.append([beginwind,endwind, listposwind]) 
  return listposwind
   

