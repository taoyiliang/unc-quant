import matplotlib.pyplot as plt
import cPickle as pk
import numpy as np

case = '2d'

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

#for order in [2,4,8,16]:
for order in [2,8,16,32]:
  ct,bn=pk.load(file(case+'_SC'+str(order)+'_normal.pk','r'))
  plt.plot(ct,bn,label='SC'+str(order))

plt.legend()
plt.title(case+' Problem, Normal Unc')
plt.xlabel('Solution Value')
plt.ylabel('Solution Frequency')
#plt.axis([0.5,2.5,0,2.1]) #source
#plt.axis([0.6,2,0,2]) #1d
plt.axis([0.9,1.3,0,25]) #poly

plt.show()

