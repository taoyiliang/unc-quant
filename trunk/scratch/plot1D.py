import numpy as np
import matplotlib.pyplot as p

xs=np.arange(-1,1,0.001)

def func(x):
  return 4*np.exp(-x*x/2.)

def apx(x):
  rt2=np.sqrt(2.)
  tot=2.*rt2
  tot-=1./rt2*(x*x-1.)
  tot+=1./(8.*rt2) * (x**4-6.*x*x+3.)
  tot-=1./(96.*rt2) * (x**6-15.*x**4+45.*x*x-15.)
  return tot

act=func(xs)
apr=apx(xs)

p.plot(xs,act,'k--')
p.plot(xs,apr,'b-')

print np.sum(abs(apr**2-act**2))/len(xs)

p.show()
