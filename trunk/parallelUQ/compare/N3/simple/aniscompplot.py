from plot1 import addPlot
from plotMC import addPlot as pltMC
import matplotlib.pyplot as plt
import numpy as np

MC = 'MC_h5_iso.samples'
plots= ['HC_h5_iso.moments',
        'HC_h5_1-1-2-1-1.moments',
        'HC_h5_1-1-4-4-8.moments',
        'HC_h5_1-1-2-2-4.moments',
        'HC_h5_8-8-4-4-1.moments']

ref=[0,1.0013454858344246,1.8647257169268627e-04]

for p in plots:
  addPlot(p,p.split('_')[2].split('.')[0],ref=ref)
pltMC(MC,'MC',ref=ref)

xs=np.linspace(2,32000,3)
ysMC = 1e-2/np.sqrt(xs)
ysHC1 = 1e-0*xs**(-3)
#ysHC2 = 5e-5*xs**(-0.5)
#plt.loglog(xs,ysMC,'k:',label=r'$c\ \eta^{-1/2}$')
#plt.loglog(xs,ysHC1,'k-.',label=r'$c\ \eta^{-4}$')
#plt.loglog(xs,ysHC2,'k:')#,label=r'$C_2\eta^{-1/2}$')
plt.title(r'Error in $k$; $N$=5, $h\sim$0.018')
plt.xlabel(r'PDE Solves $\eta$')
plt.ylabel('Rel. Error')
plt.legend(loc=2)
plt.axis([1,1e4,1e-9,1e-2])
plt.show()
