import numpy as np
import classVar as var
import classQuadrature as q
import classUQ as uq
import os

class Solver:
  def __init__(self):
    pass
  def run(self,w,x,y,z):
    return w*x*x*(1.-y)*(1.-z*z)

def runMFS():
# set non-random variables, set None for randoms
  x=None
  y=None
  w=None
  z=None

  low=0.
  hiw=1.

  ux=3.
  sx=5.

  loy=0.
  hiy=1.
  ay=3.
  by=2.

  az=2.
  bz=3.
  loz=0
  p=[w,x,y,z]
# designate random variables
  randVars=[]
  randVars.append(var.uniVar('w',[low,hiw],paramIndex=0))
  randVars.append(var.normVar('x',[ux,sx],paramIndex=1))
  randVars.append(var.betaVar('y',[ay,by,loy,hiy],paramIndex=2))
  randVars.append(var.gammaVar('z',[az,bz,loz],paramIndex=3))

  solver=Solver()
  scsolve=uq.SC(p,randVars,solver)
  scsolve.propUQ(3)
  UQsoln=scsolve.UQsoln

# check against mfs
  xs=np.linspace(low,hiw,5)
  ys=xs.copy()
  zs=xs.copy()
  ws=xs.copy()
  #ys=np.arange(uy-2.*sy,uy+2.*sy,4.*sy/5.)[:-1]
  #ys=np.arange(0.,1.2,1./5.)
  us=np.zeros([len(ws),len(xs),len(ys),len(zs)])
  mfs=us.copy()
  for h in range(len(ws)):
    print 'on h=',h
    for i in range(len(xs)):
      for j in range(len(ys)):
        for k in range(len(zs)):
          us[h,i,j,k]=UQsoln([ws[h],xs[i],ys[j],zs[k]])
          mfs[h,i,j,k]=solver.run(ws[h],xs[i],ys[j],zs[k])
  l2norm=np.sqrt(np.sum((us-mfs)**2))
  print '\nL2norm of error is',l2norm

#visualize solution
#  import matplotlib.pyplot as plt
#  CS=plt.imshow(np.rot90(us),extent=[0,1,0,1])
#  CB=plt.colorbar(CS)
#  plt.title('Beta (%1.0f, %1.0f, %1.0f, %1.0f)\n& Gamma (%1.0f,%1.0f,%1.0f)'\
#      %(ax,bx,lox,hix,ay,by,loy))
#  plt.figure()
#  plt.subplot(2,2,4)
#  CS=plt.imshow(np.rot90(np.sqrt((us-mfs)**2)),extent=[0,1,0,1])
#  CB=plt.colorbar(CS)
#  plt.title('err, Gamma/Beta')
#  plt.xlabel('x')

#  plt.figure()
#  plt.plot(xs,us[:,-1])
#  plt.show()

if __name__=='__main__':
  runMFS()
