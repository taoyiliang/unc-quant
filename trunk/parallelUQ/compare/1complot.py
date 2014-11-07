from plot1 import addPlot
from plotMC import addPlot as pltMC
from plotHDMR import addPlot as addHDMR
import matplotlib.pyplot as plt
import numpy as np

MC = 'MC_h5_N14.samples'
HCs = []#'HC_h5_N5_iso.moments']
TDs = ['TD_h5_N14.moments']
hdmr = ['hdmr_TD_N14_H1.out',
        'hdmr_TD_N14_H2.out']
#HCs = ['HC_h1_iso.moments',
#       'HC_h3_iso.moments',
#       'HC_h5_iso.moments']

#N=float(300003)
#ref= [0,3.0159136527772754e+05/N,3.0332781549076462e+05/N]
for line in file('ref.inp','r'):
  ref=line.strip().split(',')
  ref[0]=int(ref[0])
  ref[1]=float(ref[1])
  ref[2]=float(ref[2])
  break

pltMC(MC,'MC',ref=ref)
for h in HCs:
  addPlot(h,'HC_iso',ref=ref)#+h.split('.')[0].split('_')[1],ref=ref)
for h in TDs:
  addPlot(h,'TD_iso',ref=ref)#+h.split('.')[0].split('_')[1],ref=ref)
for h in hdmr:
  addHDMR(h,h.split('.')[0],ref=ref)
#for h in anis:
#  addPlot(h,h.split('_')[0]+'_aniso',ref=ref)#+h.split('.')[0].split('_')[1],ref=ref)
#pltMC(MC,'MC',ref=ref)

xs=np.linspace(1,100000,3)
ysMC = 5e-1/np.sqrt(xs)
plt.loglog(xs,ysMC,'k:',label=r'$5\ \eta^{-1/2}$')

plt.title(r'Error in $\left<k\right>$; $N$=14')
plt.xlabel(r'PDE Solves $\eta$ (~$L$)')
plt.ylabel('Rel. Error')
plt.legend(loc=3)
#plt.axis([1,1e5,1e-6,1e-2])
plt.show()
