import numpy as np
import classVar as var
import classQuadrature as q
import classUQ as uq

class Solver:
  def __init__(self):
    pass
  def run(self,x):
    return 1.-x
    #return x**2+2.*x+1.

def runMFS():
# set non-random variables, set None for randoms
  x=None
  p=[x]
# designate random variables
  randVars=[]
  low=0.
  hi=1.
  randVars.append(var.uniVar('x',[low,hi],paramIndex=0))

  solver=Solver()
  scsolve=uq.SC(p,randVars,solver)
  scsolve.propUQ(4)
  UQsoln=scsolve.UQsoln

# check against mfs
  xs=np.arange(low,hi,(hi-low)/50.)
  mfs=solver.run(xs)
  us=np.zeros_like(xs)
  for i in range(len(xs)):
    us[i]=UQsoln([xs[i]])
  l2norm=np.sqrt(np.sum((us-mfs)**2))
  print '\nav l2norm of error is',l2norm/len(xs),'\n'


#visualize solution
  import matplotlib.pyplot as plt
  plt.plot(xs,mfs,'k-',label='act')
  plt.plot(xs,us,'b.',label='comp',markersize=12)
  plt.legend()
  plt.show()

if __name__=='__main__':
  runMFS()
