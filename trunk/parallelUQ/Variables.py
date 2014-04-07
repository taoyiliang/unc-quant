import numpy as np
import scipy.stats as dists
import scipy.special as polys
import scipy.special.orthogonal as quads
from scipy.misc import factorial
from math import ceil
import sys

def VariableFactory(vartype,name='x',path='dudpath'):
  types=['uniform','normal','lognormal']
  if vartype not in types:
    print 'ERROR: variable type',vartype,'not recognized for variable',path
    sys.exit()
  if vartype == 'uniform':
    newvar = Uniform(name,path)
  elif vartype == 'normal':
    newvar = Normal(name,path)
  elif vartype == 'lognormal':
    newvar = Lognormal(name,path)
  return newvar

class Variable(object):
  def __init__(self,name,path):
    self.name=name
    self.path=path

  def setDist(self,args):
    print 'ERROR no setDist method for',self.distName,'var',self.path
    sys.exit()

  def setQuadrature(self,maxOrder):
    self.quadOrd=maxOrder

  def sample(self,trunc=0):
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

  def probWeight(self,x):
    return 0.5/self.range

  def setDist(self,args,verbose=1):
    if len(args)==1:
      self.mean=float(args[0])
      frac=0.3
      low=(1.0-frac)*self.mean
      self.range = frac*self.mean
    elif len(args)==2:
      self.mean=float(args[0])
      self.range=float(args[1])
      low=self.mean-self.range
    else:
      print 'ERROR Unrecognized input arguments for',self.name,':',args
      sys.exit()
    if verbose:
      print 'set var',self.path,'type,mean,range:',
      print self.distName,self.mean,self.range
    expr='self.dist=dists.uniform(loc='+str(low)+',scale='+str(2.0*self.range)+')'
    exec expr
    print 'Value range (.9999 confidence):',self.dist.interval(0.9999)
    self.expval = self.dist.mean()
    self.secondmom = 0.5*(low*low + (self.mean+self.range)**2)

  def setQuadrature(self,maxOrder,verbose=False):
    super(Uniform,self).setQuadrature(maxOrder)
    pts,wts = quads.p_roots(self.quadOrd)

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

  def probWeight(self,x,scale='actual'):
    #if scale=='actual':
    #  coeff = 1#./(self.stdev*np.sqrt(2.*np.pi()))
    #  main = 1#np.exp(-(x-self.mean)**2/(2.*self.stdev**2))
    #  return 1#coeff*main
    #elif scale=='standard':
    return 1.0/np.sqrt(np.pi)#./np.sqrt(2.*np.pi())*np.exp(-x*x/2.)
    #else:
    #  msg = 'variable '+var.name+'.probWeight does not recognize '
    #  msg+= str(scale)+' (type '+str(type(scale))+')!\n'
    #  raise IOError(msg)

  def probdens(self,x):
    return np.exp(x*x)

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
    print 'set var',self.path,
    print 'type,mean,stdev:',self.distName,self.mean,self.stdev
    self.expval = self.dist.mean()
    self.secondmom = self.mean**2+self.stdev**2

  def setQuadrature(self,maxOrder):
    super(Normal,self).setQuadrature(maxOrder)
    #standard hermite quadrature
    pts,wts = quads.h_roots(self.quadOrd)
    self.pts,self.wts,self.quadOrds=super(Normal,self).checkPoints(pts,wts)
    self.quaddict = {}
    for o,order in enumerate(self.quadOrds):
      self.quaddict[order]=(self.pts[o],self.wts[o])

  def norm(self,n):
    return 1.0/(np.sqrt(np.sqrt(np.pi)*(2.0**n)*factorial(n)))

  def evalNormPoly(self,x,n):
    norm = self.norm(n)
    return norm*polys.eval_hermite(n,x)

  def oneDPoly(self,n):
    norm = self.norm(n)
    return norm*polys.hermite(n)

  def sample(self,trunc=0):
    #be a jerk and force it to be within x stddev
#TODO fix this! It's not reading from input file
    x = nDev #hand-fixed
    val = self.dist.rvs()
    if x > 0:
      if not (self.mean-x*self.stdev < val < self.mean+x*self.stdev):
        val = self.sample(trunc=trunc)
    return val


class Lognormal(Normal):
  def __init__(self,name,path):
    super(Lognormal,self).__init__(name,path)
    self.distName='lognormal'

  def convertToActual(self,x):
    #FIXME do I need to convert to N(0,1) first then do this?
    return np.exp(self.u+self.s*x)
  # TODO check on these values
  def convertToStandard(self,y):
    return (np.log(y)-self.u)/self.s

  def setDist(self,args):
    m = float(args[0])
    self.stdev = float(args[1])
    v=self.stdev**2
    self.u = np.log(m*m/np.sqrt(v+m*m))
    self.s = np.sqrt(np.log(1.+v/(m*m)))
    super(Lognormal,self).setDist([0.,1.])#[self.u,self.s])
    #expr='self.dist=dists.norm('+str(self.u)+','+str(self.stdev)+')'
    #exec expr
    self.mean = m
    self.stdev = np.sqrt(v)
    print '  lognorm mean, variance:',m,v
    print '    shape params mu,sigma',self.u,self.s
    #print 'm,v,u,s',m,v,self.u,self.s

  def sample(self,trunc=0,n=1):
    x = trunc #hand-fixed
    vals = self.dist.rvs(n)
    #vals=self.convertToStandard(vals)
    if n==1: return np.exp(self.u+self.s*vals)[0]
    return np.exp(self.u+self.s*vals)
    #return vals
    #cof = 1./(vals*self.stdev*np.sqrt(2.*np.pi))
    #exp = -(np.log(vals)-self.mean)**2/(2*self.stdev**2)
    if n==1: return (cof*np.exp(exp))[0]
    return cof*np.exp(exp)
    #if x > 0:
    #  if not (self.mean-x*self.stdev < val < self.mean+x*self.stdev):
    #    val = self.sample(trunc=trunc)
    #return val
