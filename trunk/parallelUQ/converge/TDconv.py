import matplotlib.pyplot as plt
from matplotlib import rc
import numpy as np
import sys

leta = {}
leta[1]=
leta[5]=
leta[10]=
leta[15]=
leta[20]=
leta[30]=
leta[50]=
#leta[70]=
#leta[90]=


def doError(h,ref=None):
  err1=[]
  err2=[]
  if ref==None:
    ref = h[-1]
  for entry in h[:-1]:
    err1.append(abs(ref[1]-entry[1])/ref[1])
    err2.append(abs(ref[2]-entry[2])/ref[2])
  return err1,err2

def doPlot(etas,errs,h):
  size=200/(h*2.0)
  xs=[]
  for e in etas:
    xs.append(leta[e])
  x=np.log(xs)
  y=np.log(errs)
  for i,mx in enumerate(x):
    print np.exp(mx),np.exp(y[i])
  plt.plot(x,y,'-o',label='h=%1.3e' %size)

def doStuff(h,num,ref=None):
  h=np.array(h)
  e1,e2=doError(h,ref=ref)
  doPlot(h[:-1,0],e1,float(num))

def plotMC():
  xs=np.array(leta.values())
  hs=[]
  for x in xs:
    hs.append(x**(-0.5)/100)
  hs=np.log(np.array(hs))
  xs=np.log(xs)
  plt.plot(xs,hs,'k:',label='MC')



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
#  [70,,],
#  [90,,],
  ]

doStuff(h2,2)#,ref=h50[-1])
#doStuff(h5,5)#,ref=h50[-1])
#doStuff(h30,30,ref=h50[-1])
#doStuff(h50,50,ref=h50[-1])
plotMC()

#plt.rc('text',usetex=True)
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
