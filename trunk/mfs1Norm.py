import numpy as np
import classVar as var
import classQuadrature as q
import classUQ as uq

class Solver:
  def __init__(self):
    pass
  def run(self,x):
    return 1.-x#x*x+2.*x+1.
    #return 1./(np.sqrt(2.*np.pi*sig2))*\
    #    np.exp(-(x-u)**2/(2.*sig2))

def runMFS():
# set non-random variables, set None for randoms
  x=None
  p=[x]
# designate random variables
  randVars=[]
  u=0.5
  s=0.3
  randVars.append(var.normVar('x',[u,s],paramIndex=0))

  solver=Solver()
  scsolve=uq.SC(p,randVars,solver)
  scsolve.propUQ(4)
  UQsoln=scsolve.UQsoln

# check against mfs
  xs=np.arange(u-2.*s,u+2.*s,4.*s/50.)[:-1]
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
