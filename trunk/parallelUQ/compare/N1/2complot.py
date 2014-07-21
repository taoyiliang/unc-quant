from plot1 import addPlot
from plotMC import addPlot as pltMC
import matplotlib.pyplot as plt
import numpy as np

r=2

MC = 'MC_h5_N1_iso.samples'
HCs = ['HC_h5_N1_iso.moments']
TDs = ['TD_h5_N1_iso.moments']
#HCs = ['HC_h1_iso.moments',
#       'HC_h3_iso.moments',
#       'HC_h5_iso.moments']

ref=[0,1.0010574538553680e+00,7.1917735918303194e-05]

for h in HCs:
  addPlot(h,'HC_iso',ref=ref,r=r)#+h.split('.')[0].split('_')[1],ref=ref)
for t in TDs:
  addPlot(t,'TD_iso',ref=ref,r=r)#+h.split('.')[0].split('_')[1],ref=ref)
pltMC(MC,'MC',ref=ref)

#xs=np.linspace(2,32000,3)
#ysMC = 1e-2/np.sqrt(xs)
#ysHC1 = 1e-0*xs**(-3)
#ysHC2 = 5e-5*xs**(-0.5)
#plt.loglog(xs,ysMC,'k:',label=r'$c\ \eta^{-1/2}$')
#plt.loglog(xs,ysHC1,'k-.',label=r'$c\ \eta^{-4}$')
#plt.loglog(xs,ysHC2,'k:')#,label=r'$C_2\eta^{-1/2}$')
plt.title(r'Error in var($k$); $N$=1, $h\sim$0.018')
plt.xlabel(r'PDE Solves $\eta$')
plt.ylabel('Rel. Error')
plt.legend(loc=3)
plt.axis([1,1e4,1e-12,1e-0])
plt.show()
