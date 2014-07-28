from plot1 import addPlot
from plotMC import addPlot as pltMC
import matplotlib.pyplot as plt
import numpy as np

MC = 'MC_h1_simple10.samples'
HCs = ['HC_h1_simple10.moments']
TDs = ['TD_h1_simple10.moments']

ref=[0,4.345037968406584659766087558051122681529155352e-12,
       0] #TODO

for h in HCs:
  addPlot(h,'HC_iso',ref=ref)#+h.split('.')[0].split('_')[1],ref=ref)
for h in TDs:
  addPlot(h,'TD_iso',ref=ref)#+h.split('.')[0].split('_')[1],ref=ref)
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
plt.legend(loc=3)
#plt.axis([1,5e4,1e-11,1e-0])
plt.show()
