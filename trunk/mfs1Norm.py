import numpy as np
#import norm1var as solver
import classVar as var
import classQuadrature as q
import classUQ as uq
import os

class Solver:
  def __init__(self):
    pass
  def run(self,x,u=0,sig2=1):
    return 1./(np.sqrt(2.*np.pi*sig2))*\
        np.exp(-(x-u)**2/(2.*sig2))

def runMFS():
# set non-random variables, set None for randoms
  x=None
  u=3
  sig2=5
  p=[x,u,sig2]
# designate random variables
  randVars=[]
  randVars.append(var.normVar('x',[u,sig2],paramIndex=0,N=8))

  solver=Solver()
  scsolve=uq.SC(p,randVars,solver)
  scsolve.propUQ(2)
  UQsoln=scsolve.UQsoln

# check against mfs
  xs=np.arange(-10,10,0.5)
  mfs=solver.run(xs,u,sig2)
  us=np.zeros_like(xs)
  for i in range(len(xs)):
    us[i]=UQsoln([xs[i]])
  l2norm=np.sqrt(np.sum((us-mfs)**2))
  print '\nl2norm of error is',l2norm

#visualize solution
  import matplotlib.pyplot as plt
  plt.plot(xs,mfs,'k-')
  plt.plot(xs,us,'b.',markersize=9)
  plt.show()

if __name__=='__main__':
  runMFS()
