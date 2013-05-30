import numpy as np
import classVar as var
import classQuadrature as q
import classUQ as uq
import os

class Solver:
  def __init__(self):
    pass
  def run(self,x):
    return x**2+2.*x+1.

def runMFS():
# set non-random variables, set None for randoms
  x=None
  p=[x,u,sig2]
# designate random variables
  randVars=[]
  randVars.append(var.uniVar('x',[3,5],paramIndex=0,N=2))
  var=randVars[0]

  solver=Solver()
  scsolve=uq.SC(p,randVars,solver)
  scsolve.propUQ(3)
  UQsoln=scsolve.UQsoln

# check against mfs
  xs=np.arange(3,5,0.01)
  mfs=solver.run(xs)
#xs=np.arange(-1,1,0.01)
  us=np.zeros_like(xs)
  for i in range(len(xs)):
    us[i]=UQsoln([xs[i]])
  l2norm=np.sqrt(np.sum((us-mfs)**2))
  print '\nl2norm of error is',l2norm,'\n'


#visualize solution
  import matplotlib.pyplot as plt
  plt.plot(xs,mfs,'k-')
  plt.plot(xs,us,'b--')
  plt.show()

if __name__=='__main__':
  runMFS()
