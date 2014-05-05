import matplotlib.pyplot as plt
from matplotlib import rc
import numpy as np
import sys

leta = {}
leta[1]=5
leta[3]=29
leta[5]=101
leta[7]=265
leta[10]=821
leta[12]=1513
leta[13]=1989
leta[14]=2577
leta[15]=3281
leta[16]=4129
leta[17]=5125
#leta[18]=
#leta[19]=
leta[20]=9241
leta[30]=41761
#leta[40]=
#leta[50]=


def doError(h):
  err1=[]
  err2=[]
  ref = h[-1]
  for entry in h[:-1]:
    err1.append(abs(ref[1]-entry[1])/ref[1])
    err2.append(abs(ref[2]-entry[2])/ref[2])
  print err1
  print err2
  return err1,err2

def doPlot(etas,errs,h):
  size=1./(h*2.0)
  xs=[]
  for e in etas:
    xs.append(leta[e])
  x=np.log(xs)
  y=np.log(errs)
  for i,mx in enumerate(x):
    print np.exp(mx),np.exp(y[i])
  plt.plot(x,y,'-o',label='h=%1.3e' %size)

def doStuff(h,num):
  h=np.array(h)
  e1,e2=doError(h)
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
  [1 ,0.995287906475,0.997929811517],
  [3 ,0.996057325093,1.00014946374],
  [5 ,0.996365889899,1.00104103429],
  [7 ,0.996532194141,1.00152182409],
  [10,0.997077840503,1.0017382807],
  [12,0.997070582473,1.00192164286],
  [13,0.996759128614,1.00217816917],
  [14,0.997067101163,1.00206261227],
  [15,0.996798481079,1.00229201514],
  [16,0.997065533907,1.00217425898],
  [17,0.996829404277,1.00238148131],
#  [18,,],
#  [19,,],
  [20,0.997065003787,1.00233972828],
  [30,0.997067691826,1.00257845249],
  #[50,,],
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
