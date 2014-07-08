import numpy as np
import matplotlib.pyplot as plt
from tools import makePDF

def MCpdf(inName):
  inFile = file(inName,'r')
  u1=[]
  u2=[]
  a=0
  b=0
  for line in inFile:
    if line.startswith('N'):
      continue
    a1,b1=line.split(',')[1:]
    u1.append(float(a1)-a)
    u2.append(float(b1)-b)
    a=float(a1)
    b=float(b1)
  l1,c1=makePDF(u1)
  l2,c2=makePDF(u2)
  return [[l1,c1],[l2,c2]]



if __name__=='__main__':
  inFile = 'MC_h5_iso.samples'
  data = MCpdf(inFile)
  plt.plot(data[0][1],data[0][0],label='1')
  plt.plot(data[1][1],data[1][0],label='2')
  plt.legend()
  plt.show()
