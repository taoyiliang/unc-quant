from plot1 import addPlot
from plotMC import addPlot as pltMC
from plotHDMR import addPlot as addHDMR
import matplotlib.pyplot as plt
import numpy as np

mom=2

for line in file('ref.inp','r'):
  if line.startswith('N'):N=int(line.strip().split('=')[1])
  elif line.startswith('ref'):
    ref=line.strip().split('=')[1].split(',')
    ref[0]=int(ref[0])
    ref[1]=float(ref[1])
    ref[2]=float(ref[2])

slnplot=plt.figure()
errplot=plt.figure()
plt.figure(errplot.number)

for line in file('list.inp','r'):
  fname = line.strip()
  ary=fname.split('.')[0].split('_')
  if ary[0]=='MC':
    pltMC(fname,'MC',ref=ref,r=mom,slnfig=slnplot)
  elif ary[0] in ['HC','TD']:
    ary=ary[0]#.remove(ary[1])
    addPlot(fname,ary,ref=ref,r=mom,slnfig=slnplot)
  elif ary[0]=='hdmr':
    namelist = fname.split('.')[0].split('_')
    name = '_'.join([namelist[0],namelist[1],namelist[-1]])
    addHDMR(fname,name,ref=ref,r=mom,slnfig=slnplot)


plt.title(r'Error in $<R^%i>$; Attenuation, $N$=%i'%(mom,N))
plt.xlabel(r'PDE Solves $\eta$')
plt.ylabel('Rel. Error')
plt.legend(loc=3)
#plt.axis([1,1e4,1e-10,1e-2])

plt.figure(slnplot.number)
#plt.plot([1,1e5],[ref[1],ref[1]],'k:')
plt.title(r'Solution for $<R^%i>$; Attenuation, $N$=%i'%(mom,N))
plt.xlabel(r'PDE Solves $\eta$')
plt.ylabel(r'$<R>$')
plt.legend(loc=4)
#plt.axis([1,int(1e5),2.5e-2,4.5e-2])
plt.gca().set_xscale('log')


plt.show()
