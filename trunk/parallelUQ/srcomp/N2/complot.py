from plot1 import addPlot
from plotMC import addPlot as pltMC
import matplotlib.pyplot as plt
import numpy as np

MC = 'MC_h47_src.samples'
HC = 'TD_h47_src.moments'

ref=[0,0.019354814806,3.85792294471e-07]

#addPlot(MC,'MC',ref=ref)
addPlot(HC,'HC',ref=ref)
pltMC(MC,'MC',ref=ref)

xs=np.linspace(2,32000,3)
ysMC = 1e-2/np.sqrt(xs)
ysHC = 1e-3*xs**(-0.3)
plt.loglog(xs,ysMC,'k:',label=r'$C_1\eta^{-1/2}$')
plt.loglog(xs,ysHC,'k-.',label=r'$C_2\eta^{-3/10}$')
plt.title(r'Error in $k$; $N$=2, $h$=3/55')
plt.xlabel(r'PDE Solves $\eta$')
plt.ylabel('Rel. Error')
plt.legend()
plt.show()
