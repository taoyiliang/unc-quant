import numpy as np
import matplotlib.pyplot as plt


def chooseRegion(pair):
  x=pair[0]
  y=pair[1]
  #choose x
  rx=ry=0
  if x<50.0: rx=0
  elif x<100.0: rx=1
  elif x<150.0: rx=2
  elif x<=200.: rx=3
  else: print 'X not in bounds:',x

  if y<50.0: ry=0
  elif y<100.0: ry=1
  elif y<150.0: ry=2
  elif y<=200.: ry=3
  else: print 'Y not in bounds:',y

  return rx,ry

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
    elif line.startswith('k'): ks[i]=1.0+float(line.split(',')[-1].strip().split('+')[1])
    else:
      newpos = line.split(' | ')[0].split(',')
      pos[i].append((float(newpos[0]),float(newpos[1])))
      flux1[i][pos[i][-1]]=float(line.split(' | ')[1].split(',')[0])
      flux2[i][pos[i][-1]]=float(line.split(' | ')[1].split(',')[1])

#sort flux into regions
regFlux1={}
regFlux2={}
for i in ins:
  regFlux1[i]={}
  regFlux2[i]={}
  for x in range(4):
    for y in range(4):
      regFlux1[i][(x,y)]=[]
      regFlux2[i][(x,y)]=[]
  for key in flux1[i].keys():
    rx,ry = chooseRegion(key)
    regFlux1[i][(rx,ry)].append(flux1[i][key])
    regFlux2[i][(rx,ry)].append(flux2[i][key])

#turn into averages
for i in ins:
  for x in range(4):
    for y in range(4):
      regFlux1[i][(x,y)]=np.average(regFlux1[i][(x,y)])
      regFlux2[i][(x,y)]=np.average(regFlux2[i][(x,y)])

#do convergence for fluxes, k
gamK = []
gamF1 = []
gamF2 = []
imax = ins[-1]
for g,i in enumerate(ins[1:-1]):
  j=g+1
  im = ins[j-1]
  print j,i,im,
  gamK.append(abs(ks[imax]-ks[i])/abs(ks[imax]-ks[im]))
  num1=0.0
  den1=0.0
  num2=0.0
  den2=0.0
  for x in range(4):
    for y in range(4):
      num1+=(regFlux1[imax][(x,y)]-regFlux1[i][(x,y)])**2
      num2+=(regFlux2[imax][(x,y)]-regFlux2[i][(x,y)])**2
      den1+=(regFlux1[imax][(x,y)]-regFlux1[im][(x,y)])**2
      den2+=(regFlux2[imax][(x,y)]-regFlux2[im][(x,y)])**2
  num1=np.sqrt(num1)/16.
  num2=np.sqrt(num2)/16.
  den1=np.sqrt(den1)/16.
  den2=np.sqrt(den2)/16.
  print '|',num1,den1,'|',num2,den2
  gamF1.append(num1/den1)
  gamF2.append(num2/den2)
print 'f1',gamF1
print 'f2',gamF2
print 'k',gamK
