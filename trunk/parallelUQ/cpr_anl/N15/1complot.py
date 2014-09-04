from plot1 import addPlot
from plotMC import addPlot as pltMC
from plotHDMR import addPlot as addHDMR
import matplotlib.pyplot as plt
import numpy as np

MC = 'MC_h1_thirty15.samples'
HCs = []#'HC_h5_N5_iso.moments']
TDs = ['TD_h1_thirty15.moments']
hdmr = ['hdmr_TD_N15_H1.out',
        'hdmr_TD_N15_H2.out',
        'hdmr_TD_N15_H3.out']
       #'hdmr_TD_N5_H2.out']
#HCs = ['HC_h1_iso.moments',
#       'HC_h3_iso.moments',
#       'HC_h5_iso.moments']

N=15
ref= [0,(0.2*N*(np.exp(-1./N)-np.exp(-6./N)))**N]

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

#xs=np.linspace(2,32000,3)
#ysMC = 1e-2/np.sqrt(xs)
#ysHC1 = 1e-0*xs**(-3)
#ysHC2 = 5e-5*xs**(-0.5)
#plt.loglog(xs,ysMC,'k:',label=r'$c\ \eta^{-1/2}$')
#plt.loglog(xs,ysHC1,'k-.',label=r'$c\ \eta^{-4}$')
#plt.loglog(xs,ysHC2,'k:')#,label=r'$C_2\eta^{-1/2}$')
plt.title(r'Error in $k$; $N$=15, $h\sim$0.018')
plt.xlabel(r'PDE Solves $\eta$ (~$L$)')
plt.ylabel('Rel. Error')
plt.legend(loc=3)
plt.axis([1,1e5,1e-8,1e-0])
plt.show()
