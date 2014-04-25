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
  [1 ,0.995287906476,0.997929811519],
  [5 ,0.99661420503 ,1.00048862618],
  [10,0.997833392903,1.00144206642],
  [15,0.997478924161,1.00168506129],
  [20,0.997514098507,1.00155534068],
  [30,0.997457443635,1.00230109582],
  [50,0.997427253767,1.00182789109],
  [70,0.997372560839,1.00213650875],
  [90,0.997326377499,1.00218367713],
  ]

h5=[
  [1 ,0.995263602019,0.997881085725],
  [5 ,0.996589865396,1.00043976562],
  [10,0.99789023944,1.00139315698],
  [15,0.997454562516,1.00163613804],
  [20,0.997489736394,1.00150642451],
  [30,0.997433081203,1.00225213959],
  [50,0.997338530327,1.00177895988],
  [70,0.997348200574,1.00208756101],
  [90,0.997302018155,1.00213472676],
  ]

h30=[
  [1 ,0.995257790432,0.997869439503],
  [5 ,0.996584044263,1.00042808293],
  [10,0.997803195654,1.00138146096],
  [15,0.997427253767,1.00162443839],
  [20,0.99748390956,1.00149472675],
  [30,0.997427253767,1.00224043101],
  [50,0.997332703973,1.00176725823],
  [70,0.997342373822,1.00207585493],
  [90,0.997296191579,1.00212301999]
  ]

h50=[
  [1 ,0.995257367762,0.997868591872],
  [5 ,0.996583620968,1.00042723282],
  [10,0.997802771563,1.00138060948],
  [15,0.99744831182,1.0016235868],
  [20,0.997483485653,1.00149387525],
  [30,0.997426829963,1.00223957894],
  [50,0.997334311772,1.00177083454]]

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
