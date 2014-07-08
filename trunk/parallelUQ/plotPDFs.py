import matplotlib.pyplot as plt
import numpy as np

from tools import makePDF

def plotMC(inFile):
  avg=[]
  var=[]
  for line in file(inFile,'r'):
    if line.startswith('N'):continue
    entry = line.split(',')
    avg.append(float(entry[1]))
    var.append(float(entry[2]))
  abn,act = makePDF(avg)
  vbn,vct = makePDF(var)
  plt.plot(act,abn,'k:',label='MC')


def plotHC():
  pass

plotMC('MC_h5_N1_iso.samples')
plt.show()
