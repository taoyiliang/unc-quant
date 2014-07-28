from plot1 import addPlot
from plotMC import addPlot as pltMC
import matplotlib.pyplot as plt
import numpy as np
import cPickle as pk

def plot(case,ref,r):
  MC  = 'MC_'+case+'.samples'
  HCs = ['HC_'+case+'.moments']
  TDs = ['TD_'+case+'.moments']
  xs=[]
  ys=[]
  lbl=[]

  x,y=pltMC(MC,'MC',ref=ref,r=r)
  xs.append(x)
  ys.append(y)
  lbl.append('MC_'+case)

  for h in HCs:
    x,y=addPlot(h,'HC_iso',ref=ref,r=r)#+h.split('.')[0].split('_')[1],ref=ref)
    xs.append(x)
    ys.append(y)
    lbl.append('HC_iso_'+case)
  for h in TDs:
    x,y=addPlot(h,'TD_iso',ref=ref,r=r)#+h.split('.')[0].split('_')[1],ref=ref)
    xs.append(x)
    ys.append(y)
    lbl.append('TD_iso_'+case)

  pk.dump([xs,ys,lbl],file(case+'.pk','w'))

  plt.title(case)
  plt.xlabel(r'PDE Solves $\eta$')
  plt.ylabel('Rel. Error')
  plt.legend(loc=3)
