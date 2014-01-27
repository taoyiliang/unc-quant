import numpy as np
import scipy.stats as dists
import scipy.special as polys
import scipy.special.orthogonal as quads
from scipy.misc import factorial
import sys

def newVar(vartype,name,path='dudpath'):
  types=['uniform','normal']
  if vartype not in types:
    print 'ERROR: variable type',vartype,'not recognized for variable',path
    sys.exit()
  if vartype == 'uniform':
    newvar = Uniform(name,path)
  elif vartype == 'normal':
    newvar = Normal(name,path)
  return newvar

class Variable(object):
  def __init__(self,name,path):
    self.name=name
    self.path=path

  def setDist(self,args):
    print 'ERROR no setDist method for',self.distName,'var',self.path
    sys.exit()

  def setQuadrature(self,inputfile):
    print 'ERROR no setQuadrature method for',self.distName,'var',self.path
    sys.exit()

  def sample(self):
    return self.dist.rvs()

  def lowhi(self,pct=99.99):
    return list(self.dist.interval(pct))

  def prob(self,x):
    return self.dist.pdf(x)




class Uniform(Variable):
  def __init__(self,name,path):
    super(Uniform,self).__init__(name,path)
    self.distName='uniform'

  def convertToActual(self,x):
    return self.range*x+self.mean

  def convertToStandard(self,y):
    return (y-self.mean)/self.range

  def intvlShiftCoeff(self):
    return self.range

  def setDist(self,args,verbose=1):
    if len(args)==1:
      self.mean=float(args[0])
      frac=0.3
      low=(1.0-frac)*self.mean
      self.range = frac*self.mean
    elif len(args)==2:
      self.mean=float(args[0])
      if args[0]==0:
        self.range=float(args[1])
      else:
        self.range=float(args[1])
      low=-self.range
    else:
      print 'ERROR Unrecognized input arguments for',self.name,':',args
      sys.exit()
    if verbose:
      print 'set var',self.path,'type,mean,range:',self.distName,self.mean,self.range
    expr='self.dist=dists.uniform('+str(low)+','+str(2.0*self.range)+')'
    exec expr

  def setQuadrature(self,inputfile=None,order=2,verbose=False):
    if verbose:
      print 'set quadrature for',self.name
    #get order from input file
    if inputfile != None:
      self.order=inputfile('Variables/'+self.name+'/order',2)
    else:
      self.order=order
    #standard legendre quadrature
    self.pts,self.wts = quads.p_roots(self.order)
    #print 'set pts:',self.pts
    #print '    wts:',self.wts

  def evalNormPoly(self,x,n):
    norm = 1.0/np.sqrt(2.0/(2.0*float(n)+1.0))
    return norm*polys.eval_legendre(n,x)

  def oneDPoly(self,n):
    norm = 1.0/np.sqrt(2.0/(2.0*float(n)+1.0))
    return np.poly1d(norm*polys.legendre(n))







class Normal(Variable):
  def __init__(self,name,path):
    super(Normal,self).__init__(name,path)
    self.distName='normal'

  def convertToActual(self,x):
    return np.sqrt(2)*self.stdev*x+self.mean

  def convertToStandard(self,y):
    return (y-self.mean)/(np.sqrt(2)*self.stdev)

  def setDist(self,args):
    if len(args)==1:
      self.mean=float(args[0])
      self.stdev=0.3*self.mean
    elif len(args)==2:
      self.mean=float(args[0])
      self.stdev=float(args[1])
    else:
      print 'ERROR Unrecognized input arguments for',self.name,':',args
      sys.exit()
    self.variance = self.stdev*self.stdev
    expr='self.dist=dists.norm('+str(self.mean)+','+str(self.stdev)+')'
    exec expr
    print 'set var',self.path,'type,mean,stdev:',self.distName,self.mean,self.stdev

  def setQuadrature(self,inputfile=None,order=2,verbose=False):
    if verbose:
      print 'set quadrature for',self.name
    #get order from input file
    if inputfile != None:
      self.order=inputfile('Variables/'+self.name+'/order',2)
    else:
      self.order=order
    #standard legendre quadrature
    self.pts,self.wts = quads.h_roots(self.order)

  def evalNormPoly(self,x,n):
    norm = 1.0/(np.sqrt(np.sqrt(np.pi)*(2.0**n)*factorial(n)))
    return norm*polys.eval_hermite(n,x)

  def oneDPoly(self,n):
    norm = 1.0/(np.sqrt(np.sqrt(np.pi)*(2.0**n)*factorial(n)))
    return norm*polys.hermite(n)

