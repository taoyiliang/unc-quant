import matplotlib.pyplot as plt
from matplotlib import rc
import numpy as np
import sys

leta = {}
leta[1]=2
leta[5]=12
leta[10]=31
leta[30]=241
leta[50]=651
leta[70]=1261
leta[90]=2071
leta[120]=3661
leta[150]=5701


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
  [1  ,1.00369460152,1.01059835867],
  [5  ,1.00370284043,1.01063657116],
  [10 ,1.00309569978,1.00889283548],
  [30 ,1.00347773739,1.00999461201],
  [50 ,1.00356589991,1.01024886814],
  [70 ,1.00360508326,1.01036187086],
  [90 ,1.00362723036,1.01042574196],
  [120,1.003646836,1.01048228359],
  [150,1.00365870257,1.01051650616],
  ]

h5=[
  [1 ,1.00367008259,1.01054897075],
  [5 ,1.00367832124,1.01058718108],
  [10,1.00307119779,1.00884353752],
  [30,1.00345322456,1.0099452558],
  [50,1.00354138459,1.01019949848],
  [70,1.00358056682,1.01031249523],
  [90,1.0036027133,1.01037636296],
  [120,1.00362231838,1.0104329016],
  [150,1.00363418461,1.01046712236],
  ]

h30=[
  [1 ,1.00366421416,1.01053714361],
  [5 ,1.00367245272,1.0105753533],
  [10,1.00306533434,1.0088317348],
  [30,1.00344735792,1.00993343723],
  [50,1.00353551721,1.01018767625],
  [70,1.00357469912,1.01030067137],
  [90,,],
  [120,,],
  [150,,],
  ]

#h50=[
#  [1 ,,],
#  [5 ,,],
#  [10,,],
#  [15,,],
#  [20,,],
#  [30,,],
#  [50,,],
#  [70,,],
#  [90,,],
#  ]

doStuff(h2,2,0)
doStuff(h5,5,0)
#doStuff(h30,30,0)
#doStuff(h50,50,0)
plotMC()
plt.title('First Moment Error')
plt.legend(loc=3)
plt.xlabel(r'log $\eta$')
plt.ylabel('log error')

plt.figure()
doStuff(h2,2,1)
doStuff(h5,5,1)
#doStuff(h30,30,1)
#doStuff(h50,50,1)
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
