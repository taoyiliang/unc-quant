import matplotlib.pyplot as plt
from matplotlib import rc
import numpy as np
import sys

leta = {}
leta[1]=5
leta[5]=45
leta[10]=91
leta[15]=361
leta[20]=597
leta[30]=1495
leta[50]=4483
leta[70]=9149
leta[90]=16207


def doError(h):
  err1=[]
  err2=[]
  ref = h[-1]
  for entry in h[:-1]:
    err1.append(abs(ref[1]-entry[1])/ref[1])
    err2.append(abs(ref[2]-entry[2])/ref[2])
  return err1,err2

def doPlot(etas,errs,h):
  size=1./(h*2.0)
  xs=[]
  for e in etas:
    xs.append(leta[e])
  x=np.log(xs)
  y=np.log(errs)
  #for i,mx in enumerate(x):
  #  print np.exp(mx),np.exp(y[i])
  plt.plot(x,y,'-o',label='h=%1.2e' %size)

def doStuff(h,num,j):
  h=np.array(h)
  e1,e2=doError(h)
  if j==0:
    doPlot(h[:-1,0],e1,float(num))
  else:
    doPlot(h[:-1,0],e2,float(num))

def plotMC(fac=100):
  xs=np.array(leta.values())
  hs=[]
  for x in xs:
    hs.append(x**(-0.5)/fac)
  hs=np.log(np.array(hs))
  xs=np.log(xs)
  plt.plot(xs,hs,'k:',label=r'$\eta^{-1/2}$')



# DATA
#   h      m1          m2

h2=[
  [1 ,,],
  [5 ,,],
  [10,,],
  [15,,],
  [20,,],
  [30,,],
  [50,,],
  [70,,],
  [90,,],
  ]

h5=[
  [1 ,,],
  [5 ,,],
  [10,,],
  [15,,],
  [20,,],
  [30,,],
  [50,,],
  [70,,],
  [90,,],
  ]

h30=[
  [1 ,,],
  [5 ,,],
  [10,,],
  [15,,],
  [20,,],
  [30,,],
  [50,,],
  [70,,],
  [90,,],
  ]

h50=[
  [1 ,,],
  [5 ,,],
  [10,,],
  [15,,],
  [20,,],
  [30,,],
  [50,,],
  [70,,],
  [90,,],
  ]

doStuff(h2,2,0)
doStuff(h5,5,0)
doStuff(h30,30,0)
doStuff(h50,50,0)
plotMC()
plt.title('First Moment Error')
plt.legend(loc=3)
plt.xlabel(r'log $\eta$')
plt.ylabel('log error')

plt.figure()
doStuff(h2,2,1)
doStuff(h5,5,1)
doStuff(h30,30,1)
doStuff(h50,50,1)
plotMC(30)
plt.title('Second Moment Error')
plt.legend(loc=3)
plt.xlabel(r'log $\eta$')
plt.ylabel('log error')


plt.show()
sys.exit()
#h50=[
#  [1 ,,],
#  [5 ,,],
#  [10,,],
#  [15,,],
#  [20,,],
#  [30,,]]

#toplot = [h2]
#for h in toplot:
#  plt.figure()
#  h=np.array(h)
