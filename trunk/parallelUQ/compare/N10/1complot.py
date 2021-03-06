from plot1 import addPlot
from plotMC import addPlot as pltMC
import matplotlib.pyplot as plt
import numpy as np

MC = 'MC_h5_iso.samples'
HCs = ['HC_h5_iso.moments']
anis = ['HC_h5_1124244884.moments',
        'HC_h5_1234244554.moments']
#HCs = ['HC_h1_iso.moments',
#       'HC_h3_iso.moments',
#       'HC_h5_iso.moments']

ref=[0,9.3045730355252931e-01,1.2261150319231584e-03]

#addPlot(MC,'MC',ref=ref)
for h in HCs:
  addPlot(h,'HC_iso',ref=ref)#+h.split('.')[0].split('_')[1],ref=ref)
#for h in TDs:
#  addPlot(h,'TD_iso',ref=ref)#+h.split('.')[0].split('_')[1],ref=ref)
for h in anis:
  addPlot(h,h.split('_')[0]+'_'+h.split('_')[-1].split('.')[0],ref=ref)#+h.split('.')[0].split('_')[1],ref=ref)
pltMC(MC,'MC',ref=ref)

#xs=np.linspace(2,32000,3)
#ysMC = 1e-2/np.sqrt(xs)
#ysHC1 = 1e-0*xs**(-3)
#ysHC2 = 5e-5*xs**(-0.5)
#plt.loglog(xs,ysMC,'k:',label=r'$c\ \eta^{-1/2}$')
#plt.loglog(xs,ysHC1,'k-.',label=r'$c\ \eta^{-4}$')
#plt.loglog(xs,ysHC2,'k:')#,label=r'$C_2\eta^{-1/2}$')
plt.title(r'Error in $k$; $N$=5, $h\sim$0.018')
plt.xlabel(r'PDE Solves $\eta$')
plt.ylabel('Rel. Error')
plt.legend(loc=2)
#plt.axis([10,2e4,1e-9,1e-2])
plt.show()
