import matplotlib.pyplot as plt
import numpy as np

inFile = file('conv.text','r')
n=[]
h=[]
k=[]
for line in inFile:
  if line.startswith('h'):continue
  n.append(int(line.split(',')[0]))
  h.append(1./float(n[-1]))
  k.append(1.0+float(line.split(',')[1]))

ks=np.array(k)
act = k[-1]
err = []
for k,val in enumerate(ks[1:-1]):
  err.append(abs(val-act)/abs(val-ks[k-1]))

#plt.plot(h[:-1],np.log(err))
plt.plot(h[:-2],err)
plt.show()
