import numpy as np
import classVar as var
import classQuadrature as q
import classUQ as uq

class Solver:
  def __init__(self):
    pass
  def run(self,x):
    #return 1.-x
    return x*x+2.*x+1.
    #return 1./(np.sqrt(2.*np.pi*sig2))*\
    #    np.exp(-(x-u)**2/(2.*sig2))

def runMFS():
# set parameters, set None for randoms
  x=None
  p=[x]
# designate random variables
  randVars=[]
  alpha=1.
  beta=2.
  a=3.
  randVars.append(var.gammaVar('x',[alpha,beta,a],paramIndex=0))

  solver=Solver()
  scsolve=uq.SC(p,randVars,solver)
  scsolve.propUQ(4)
  UQsoln=scsolve.UQsoln

# check against mfs
  xs=np.linspace(a,4.*(2.*alpha/beta-a),20.)
  #xs=np.linspace(0,10,20.)
  mfs=solver.run(xs)
  us=np.zeros_like(xs)
  for i in range(len(xs)):
    us[i]=UQsoln([xs[i]])
    print xs[i],mfs[i],us[i]
  l2norm=np.sqrt(np.sum((us-mfs)**2))
  print '\nav l2norm of error is',l2norm/len(xs),'\n'

#visualize solution
  import matplotlib.pyplot as plt
  plt.plot(xs,mfs,'k-',label='act')
  plt.plot(xs,us,'b.',label='comp',markersize=12)
  plt.title('Gamma Unc. (%1.0f, %1.0f, %1.0f)'%(alpha,beta,a))
  plt.ylabel('f(x)=x*x+2*x+1')
  plt.xlabel('x')
  plt.legend(loc=1)
  #plt.show()

if __name__=='__main__':
  runMFS()
