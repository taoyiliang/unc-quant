import matplotlib.pyplot as plt
import numpy as np
import cPickle as pk
from tools import makePDF

def avg(v):
  if len(v)%2!=0:
    return [v[0]]+avg(v[1:])
  else:
    v1=[]
    v2=[]
    for i in range(0,len(v)-1,2):
      v1.append(v[i])
      v2.append(v[i+1])
    return list(0.5*(np.array(v1)+np.array(v2)))

def plotMC(inFile):
  av=[]
  var=[]
  i=0
  tot=0
  for line in file(inFile,'r'):
    if line.startswith('N'):continue
    i+=1
    entry = line.split(',')
    if i>1:
      av.append(float(entry[1])-tot)
      var.append(float(entry[2])-tot)
    else:
      av.append(float(entry[1]))
      var.append(float(entry[2]))
    tot+=av[-1]
  abn,act = makePDF(av,bins=50)
  #vbn,vct = makePDF(var)
  plt.plot(avg(act),avg(abn),'k',label='MC')


def plotSC(inFile):
  data=pk.load(file(inFile,'r'))
  for d in data:
    plt.plot(avg(avg(d[1])),avg(avg(d[2])),label=r'$\eta$='+str(d[0]))




plotMC('MC_h5_N5_iso.samples')
plotSC('TD_h5_N10_iso.ROMpdf')
plt.legend()
plt.show()
