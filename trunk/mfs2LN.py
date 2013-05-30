import numpy as np
import classVar as var
import classQuadrature as q
import classUQ as uq
import os

class Solver:
  def __init__(self):
    pass
  def run(self,x,y):
    return (1.-x)*(1.-y)

def runMFS():
# set non-random variables, set None for randoms
  x=None
  y=None
  lox=0.
  hix=1.
  uy=0.5
  sy=0.3
  p=[x,y]
# designate random variables
  randVars=[]
  randVars.append(var.uniVar('x',[lox,hix],paramIndex=0))
  randVars.append(var.normVar('y',[uy,sy],paramIndex=1))

  solver=Solver()
  scsolve=uq.SC(p,randVars,solver)
  scsolve.propUQ(8)
  UQsoln=scsolve.UQsoln

# check against mfs
  xs=np.arange(lox,hix+.2,(hix-lox)/5.)
  #ys=np.arange(uy-2.*sy,uy+2.*sy,4.*sy/5.)[:-1]
  ys=np.arange(0.,1.2,1./5.)
  us=np.zeros([len(xs),len(ys)])
  mfs=us.copy()
  for i in range(len(xs)):
    for j in range(len(ys)):
      us[i,j]=UQsoln([xs[i],ys[j]])
      mfs[i,j]=solver.run(xs[i],ys[j])
      print xs[i],ys[j],us[i,j],mfs[i,j]
  l2norm=np.sqrt(np.sum((us-mfs)**2))
  print '\nl2norm of error is',l2norm

#visualize solution
  import matplotlib.pyplot as plt
  CS=plt.imshow(np.rot90(us),extent=[0,1,-.1,1.1])
  CB=plt.colorbar(CS)
  plt.title('Calculated')
  plt.figure()
  CS=plt.imshow(np.rot90(mfs),extent=[0,1,-.1,1.1])
  CB=plt.colorbar(CS)
  plt.title('MFS')

  plt.figure()
  plt.plot(xs,us[:,-1])
  plt.show()

if __name__=='__main__':
  runMFS()
