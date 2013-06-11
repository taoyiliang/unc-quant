import numpy as np
import matplotlib.pyplot as plt
import classVar as var
import classQuadrature as q
import classUQ as uq

class Solver:
  def __init__(self):
    pass
  def run(self,x):
    return 1.-x
    #return x*x+2.*x+1.

def runMFS():
# set non-random variables, set None for randoms
  x=None
  p=[x]
# designate random variables
  randVars=[]
  low=0#-2.
  hi=1#5.
  randVars.append(var.uniVar([low,hi],paramIndex=0))

  solver=Solver()
  scsolve=uq.SC(p,randVars,solver)
  scsolve.propUQ(8)
  UQsoln=scsolve.UQsoln

# check against mfs
  xs=np.linspace(low,hi,20.)
  mfs=solver.run(xs)
  us=np.zeros_like(xs)
  for i in range(len(xs)):
    us[i]=UQsoln([xs[i]])
  l2norm=np.sqrt(np.sum((us-mfs)**2))
  print '\nav l2norm of error is',l2norm/len(xs),'\n'


#visualize solution
  plt.plot(xs,mfs,'k-',label='act')
  plt.plot(xs,us,'b.',label='comp',markersize=12)
  plt.legend(loc=2)
  plt.title('Uniform Test (%1.0f, %1.0f)' %(low,hi))
  plt.ylabel('f(x) = x*x + 2x + 1')
  #plt.show()

if __name__=='__main__':
  runMFS()
  plt.show()
