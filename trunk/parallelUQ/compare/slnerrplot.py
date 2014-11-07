from plot1 import addPlot
from plotMC import addPlot as pltMC
from plotHDMR import addPlot as addHDMR
import matplotlib.pyplot as plt
import numpy as np

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
    pltMC(fname,'MC',ref=ref,slnfig=slnplot)
  elif ary[0] in ['HC','TD']:
    ary.remove(ary[1])
    print ary
    addPlot(fname,'_'.join(ary),ref=ref,slnfig=slnplot)
  elif ary[0]=='hdmr':
    addHDMR(fname,fname.split('.')[0],ref=ref,slnfig=slnplot)

#addPlot(MC,'MC',ref=ref)
#for h in HCs:
#  addPlot(h,'HC_iso',ref=ref,slnfig=slnplot)
#for h in TDs:
#  addPlot(h,'TD_iso',ref=ref,slnfig=slnplot)

#for h in anis:
#  addPlot(h,h.split('_')[0]+'_aniso',ref=ref)#+h.split('.')[0].split('_')[1],ref=ref)
#pltMC(MC,'MC',ref=ref,slnfig=slnplot)

#xs=np.linspace(2,32000,3)
#ysMC = 1e-2/np.sqrt(xs)
#ysHC1 = 1e-0*xs**(-3)
#ysHC2 = 5e-5*xs**(-0.5)
#plt.loglog(xs,ysMC,'k:',label=r'$c\ \eta^{-1/2}$')
#plt.loglog(xs,ysHC1,'k-.',label=r'$c\ \eta^{-4}$')
#plt.loglog(xs,ysHC2,'k:')#,label=r'$C_2\eta^{-1/2}$')
plt.title(r'Error in $<k>$; $N$=%i'%N)
plt.xlabel(r'PDE Solves $\eta$')
plt.ylabel('Rel. Error')
plt.legend(loc=3)
#plt.axis([1,1e4,1e-10,1e-2])

plt.figure(slnplot.number)
plt.title(r'Solution for $\rho$; $N=$%i'%N)
plt.xlabel(r'PDE Solves $\eta$')
plt.ylabel(r'$\rho$')
plt.legend(loc=4)
plt.gca().set_xscale('log')
plt.axis([1,int(1e5),0.998,1.012])


plt.show()
