import numpy as np
import matplotlib.pyplot as plt

M = int(1e4)
c = 1.0

def f(x):
  return x*x

samps = np.random.rand(10000)
vals = f(samps)

tot = 0.
means=[]
high=[]
low=[]
for m in range(1,M+1):
  tot += vals[m-1]
  mean = tot/float(m)
  means.append(mean)
  err = 2.*c/np.sqrt(m)
  high.append(mean*(1.0 + err))
  low.append(mean*(1.0 - err))

plt.semilogx(range(M),means,'b-',label="average")
plt.semilogx(range(M),high,'b-',alpha=0.5,label="95 %")
plt.semilogx(range(M),low,'b-',alpha=0.5,label="5 %")
plt.legend()

plt.savefig('conv_mc.pdf')
plt.show()

