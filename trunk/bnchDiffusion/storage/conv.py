import numpy as np
import matplotlib.pyplot as plt

#read in data
flux1={}
flux2={}
pos={}
ks={}
dim={}

#load data
ins = ['1','2','4','8','16','32']
for i in ins:
  flux1[i]={}
  flux2[i]={}
  pos[i]=[]
  ks[i]=0.0
  dim[i]=1
  inFile = file(i+'.out','r')
  for line in inFile:
    if line.startswith('DIM'): dim[i]=int(line.split(',')[-1].strip())
    elif line.startswith('k'): ks[i]=float(line.split(',')[-1].strip())
    else:
      newpos = line.split(' | ')[0].split(',')
      pos[i].append((float(newpos[0]),float(newpos[1])))
      flux1[i][pos[i][-1]]=float(line.split(' | ')[1].split(',')[0])
      flux2[i][pos[i][-1]]=float(line.split(' | ')[1].split(',')[1])
