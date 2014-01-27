import numpy as np
import scipy.stats as dists
import scipy.special as polys
import scipy.special.orthogonal as quads
from scipy.misc import factorial
from math import ceil
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

  def setQuadrature(self,inputfile,order='dud',verbose=False):
    if verbose:
      print 'set quadrature for',self.name
    #get acceptable range from input file
    if inputfile != None:
      okRange = inputfile('Variables/'+self.name+'/okRange','-inf inf')
      low = okRange.split(' ')[0].strip()
      hi = okRange.split(' ')[1].strip()
      if low=='-inf':
        low=-1e30
      else:
        try: low=float(low)
        except ValueError:
          msg = 'okRange value not recognized for var "'
          msg+= self.name+'": "'+low+'"'
          raise IOError(msg)
      if hi=='inf':
        hi=1e30
      else:
        try: hi=float(hi)
        except ValueError:
          msg = 'okRange value not recognized for var "'
          msg+= self.name+'": "'+hi+'"'
          raise IOError(msg)
      self.okRange=(low,hi)
    else: #no input file?
      self.okRange=(-1e30,1e30)

  def checkPoints(self,pts,wts):
    retpts=[]
    retwts=[]
    orders=[]
    for p,pt in enumerate(pts):
      if self.convertToActual(pt) < self.okRange[0] or pt > self.okRange[1]:
        msg = 'WARNING: Sample point '+str(p)+' for variable "'
        msg+= self.name+'" discarded for being outside target range ('
        msg+= str(self.okRange[0])+','+str(self.okRange[1])+'): '
        msg+= str(pt)
        continue
      retpts.append(pt)
      retwts.append(wts[p])
      orders.append(p)
    return retpts,retwts,orders

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

#  def intvlShiftCoeff(self):
#    return self.range

  def probWeight(self,x,scale='actual'):
    if scale=='actual':
      return 0.5/self.range
    elif scale=='standard':
      return 0.5
    else:
      msg = 'variable '+var.name+'.probWeight does not recognize '
      msg+= str(scale)+' (type '+str(type(scale))+')!\n'
      raise IOError(msg)

  def probdens(self,x):
    return 1.0

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
    print low
    print 'Value range (.9999 confidence):',self.dist.interval(0.9999)
    self.expval = self.dist.mean()
    self.secondmom = 0.5*(low*low + (self.mean+self.range)**2)

  def setQuadrature(self,inputfile=None,expOrd=2,verbose=False):
    super(Uniform,self).setQuadrature(inputfile,order=expOrd,verbose=verbose)
    print 'okRange:',self.okRange
    #get order from input file
    if inputfile != None:
      self.expOrd=inputfile('Variables/'+self.name+'/order',2)
    else:
      self.expOrd=expOrd
    #quad of order N can do polys of order 2N-1
    self.quadOrd = int(ceil((self.expOrd+1)/2.0))
    #standard hermite quadrature
    pts,wts = quads.p_roots(self.quadOrd)
    self.pts,self.wts,self.quadOrds=super(Uniform,self).checkPoints(pts,wts)
    self.quaddict = {}
    for o,order in enumerate(self.quadOrds):
      self.quaddict[order]=(self.pts[o],self.wts[o])


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

  def setQuadrature(self,inputfile=None,expOrd=2,verbose=False):
    super(Normal,self).setQuadrature(inputfile,order=expOrd,verbose=verbose)
    print 'okRange:',self.okRange
    #get order from input file
    if inputfile != None:
      self.expOrd=inputfile('Variables/'+self.name+'/order',2)
    else:
      self.expOrd=expOrd
    #quad of order N can do polys of order 2N-1
    self.quadOrd = int(ceil((self.expOrd+1)/2.0))
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

  def sample(self,nDev=0):
    #be a jerk and force it to be within x stddev
    x = 2.0 #hand-fixed
    val = self.dist.rvs()
    if nDev > 0:
      if not (val>self.mean-x*self.stdev and val<self.mean+x*self.stdev):
        val = self.sample()
    return val