from plot1 import addPlot
from plotMC import addPlot as pltMC
import matplotlib.pyplot as plt
import numpy as np

r=2
MC = 'MC_h5_iso.samples'
HCs = ['HC_h5_iso.moments']
#HCs = ['HC_h1_iso.moments',
#       'HC_h3_iso.moments',
#       'HC_h5_iso.moments']

ref=[3001,1.0010574558596474e+00,7.1917726521375513e-05]

#addPlot(MC,'MC',ref=ref)
for h in HCs:
  addPlot(h,'HC',ref=ref,r=r)#+h.split('.')[0].split('_')[1],ref=ref)
pltMC(MC,'MC',ref=ref,r=r)

xs=np.linspace(2,32000,3)
ysMC = 1e-0/np.sqrt(xs)
ysHC1 = 5e1*xs**(-11)
#ysHC2 = 5e-5*xs**(-0.5)
plt.loglog(xs,ysMC,'k:',label=r'$c\ \eta^{-1/2}$')
plt.loglog(xs,ysHC1,'k-.',label=r'$c\ \eta^{-11}$')
#plt.loglog(xs,ysHC2,'k:')#,label=r'$C_2\eta^{-1/2}$')
plt.title(r'Error in var($k$); $N$=1, $h\sim$0.018')
plt.xlabel(r'PDE Solves $\eta$')
plt.ylabel('Rel. Error')
plt.legend(loc=4)
plt.axis([1,2e4,1e-9,1e0])
plt.show()
