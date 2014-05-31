from plot1 import addPlot
from plotMC import addPlot as pltMC
import matplotlib.pyplot as plt
import numpy as np

r=2
MC = 'MC_h5_iso.samples'
isos = ['HC_h5_iso.moments',
        'TD_h5_iso.moments']
anis = ['HC_h5_aniso.moments',
        'TD_h5_aniso.moments']
#HCs = ['HC_h1_iso.moments',
#       'HC_h3_iso.moments',
#       'HC_h5_iso.moments']

ref=[0,1.0025998595112884,1.7000448118653644e-04]

#addPlot(MC,'MC',ref=ref)
for h in isos:
  addPlot(h,h.split('_')[0]+'_iso',ref=ref,r=r)
for h in anis:
  addPlot(h,h.split('_')[0]+'_aniso',ref=ref,r=r)
pltMC(MC,'MC',ref=ref,r=r)

xs=np.linspace(2,32000,3)
ysMC = 1e-0/np.sqrt(xs)
ysHC1 = 5e-0*xs**(-2)
#ysHC2 = 5e-5*xs**(-0.5)
#plt.loglog(xs,ysMC,'k:',label=r'$c\ \eta^{-1/2}$')
#plt.loglog(xs,ysHC1,'k-.',label=r'$c\ \eta^{-2}$')
#plt.loglog(xs,ysHC2,'k:')#,label=r'$C_2\eta^{-1/2}$')
plt.title(r'Error in var($k$); $N$=3, $h\sim$0.018')
plt.xlabel(r'PDE Solves $\eta$')
plt.ylabel('Rel. Error')
plt.legend(loc=3)
plt.axis([10,2e4,1e-7,1e-0])
plt.show()
