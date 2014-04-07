import matplotlib.pyplot as plt
import cPickle as pk
import numpy as np

#case = 'simple'
#case = 'source'
#case = '1dsup'
case = '2d_5v'

douni = True
#douni = False

#donorm = True
donorm = False

cases=['2-2-2-2-2',
       '4-2-2-2-2',
       '2-4-2-2-2',
       '2-2-4-2-2',
       '2-2-2-4-2',
       '2-2-2-2-4',
       '4-4-2-2-2',
       '6-6-2-2-2']

if donorm:
##### NORMALS
  plt.figure()

#get Monte Carlo run
  ctrs=[]
  bins=[]
  for line in file(case+'_MC_normal.csv','r'):
    if line[:3] in ['Cen','Ave','2nd','Var']:
      continue
    ctrs.append(float(line.split(',')[0]))
    bins.append(float(line.split(',')[1]))
  plt.plot(ctrs,bins,'k-',label='MC')

  for order in cases:
#for order in [2]:#,4,8,16]:
    ct,bn=pk.load(file(case+'_SC'+order+'_normal.pk','r'))
    plt.plot(ct,bn,label='SC'+order)

  plt.legend()
  plt.title(case+' Problem, Normal Unc')
  plt.xlabel('Solution Value')
  plt.ylabel('Solution Frequency')
#plt.axis([0.5,2.5,0,2.1]) #source
  plt.axis([0.9,1.1,0,30]) #1d
#plt.axis([0.6,2,0,2]) #poly

if douni:
##### UNIFORM
  plt.figure()

#get Monte Carlo run
  ctrs=[]
  bins=[]
  for line in file(case+'_MC_uniform.csv','r'):
    if line[:3] in ['Cen','Ave','2nd','Var']:
      continue
    ctrs.append(float(line.split(',')[0]))
    bins.append(float(line.split(',')[1]))
  plt.plot(ctrs,bins,'ko',label='MC')

  for order in cases:#['2-2','2-4','4-2','4-4']:#,4,8,16,32]:
  #for order in [2]:#,4,8,16]:
    ct,bn=pk.load(file(case+'_SC'+order+'_uniform.pk','r'))
    plt.plot(ct,bn,label='SC'+order)

  plt.legend()
  plt.title(case+' Problem, Bivariate Uniform Unc')
  plt.xlabel('Solution Value')
  plt.ylabel('Solution Frequency')
#plt.axis([0.5,2.5,0,2.1]) #source
  #plt.axis([0.96,1.05,0,100]) #poly
plt.show()

