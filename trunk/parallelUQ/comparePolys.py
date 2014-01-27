import numpy as np
import matplotlib.pyplot as plt

coeffs={}
coeffs[18]=-1.684e-6
coeffs[16]= 2.273e-5
coeffs[14]=-0.0001961
coeffs[12]= 0.001387
coeffs[10]=-0.008333
coeffs[8] = 0.04167
coeffs[6] =-0.1667
coeffs[4] = 0.5
coeffs[2] =-1.0
coeffs[0] = 1.0

def actual(x):
  return np.exp(-x*x)

def poly(cof,x):
  tot=0
  for key in cof.keys():
    tot+=cof[key] * (x**key)
  return tot

xs=np.linspace(-1,1,100)
plt.plot(xs,actual(xs),'k-')
plt.plot(xs,poly(coeffs,xs),'go')
plt.show()
