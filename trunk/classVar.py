import numpy as np
import sys
import random as rand
import classQuadrature as q
import scipy.special as sps
import scipy.stats as spst
from scipy.misc import factorial

#============================================================================\
class var():
  def __init__(self,name,unc,paramIndex=[],N=2):
    self.name=name
    self.uncType=''
    self.unc=unc #array of uncertainty info
    self.N=N
    self.paramIndex=paramIndex # so it can be changed in the list
    self.vals=[]
    self.weights=[]
    self.setUnc()
    self.setDistr()

  def setZeroToOneSample(self,toPPF,toPoly,order):
    self.toPPF=toPPF
    self.toPoly=toPoly
    self.fromPPF=self.dist.ppf
    self.shiftQuad=q.quadShiftLegendre(order)
    self.zeroToOneCoeffs=np.zeros(order)
    for o in range(order):
      self.zeroToOneCoeffs[o]=self.calcZeroOneCoeff(o)
    #self.zeroToOneCoeffs=self.zeroToOneCoeffs[::-1]
    #print np.poly1d(self.zeroToOneCoeffs[::-1])
    self.sampleZeroToOne=np.poly1d([])
    for c,cof in enumerate(self.zeroToOneCoeffs):
      new_coeffs=self.toPoly(c).c*cof*np.sqrt(2.*c+1.)
      self.sampleZeroToOne+=np.poly1d(new_coeffs)
    #now you can call sampleZeroToOne(x)
    self.z21quad=q.quadShiftLegendre(order)

  def calcZeroOneCoeff(self,o):
    #tot=0
    #for n in range(self.shiftQuad.order):
    #  tot+=self.shiftQuad.weights[n]*self.fromPPF(self.shiftQuad.ords[n])*\
    #      self.toPoly(o)(self.toPPF(self.shiftQuad.ords[n]))
    return np.sum(self.shiftQuad.weights*self.fromPPF(self.shiftQuad.ords)*\
          self.toPoly(o)(self.toPPF(self.shiftQuad.ords))*np.sqrt(2.*o+1.))

  def UQval(self,el):
    pass
  def UQwts(self,el):
    pass
  def UQpolySample(self,n,el):
    pass

  def setUnc(self):
    pass
  def setDistr(self):
    pass
  def sample(self,val=None):
    pass
  def samplePt(self):
    pass
  def sampleWt(self):
    pass
  def samplePoly(self):
    pass


#============================================================================\
class uniVar(var):
  def setUnc(self):
    self.uncType='uniform'
    if len(self.unc)!=2:
      print 'ERROR: Uniform-uncertainty variables require low and high values, got',self.unc
      sys.exit()
    if self.unc[0]>self.unc[1]:
      print 'WARNING: Low bound higher than high bound!  Swapping...'
      dummy=self.unc[1]
      self.unc[1]=self.unc[0]
      self.unc[0]=dummy
    self.low=self.unc[0]
    self.hi=self.unc[1]

  def setDistr(self):
    self.average=0.5*(self.low+self.hi)
    self.quad=q.quadLegendre(self.N)
    self.range=(self.average-self.low)
    self.dist=spst.uniform(self.low,2.*self.range)
    self.domain=[self.low,self.hi]
    #for i in range(self.N):
    #  self.vals.append(self.average+self.range*self.quad.ords[i])
    #  self.vals.append(self.quad.weights[i])

  def sample(self,val=None):
    if val==None:
      return self.low+rand.random()*(self.hi-self.low)
    return self.average+val*self.range

  def samplePt(self,x):
    '''samples equivalent point on [-1,1]'''
    return x*self.range+self.average

  def revertPt(self,x):
    '''returns from [-1,1] to [a,b]'''
    return (x-self.average)/self.range

  def sampleWt(self,o):
    '''samples equivalent weight for [-1,1]'''
    return self.quad.weights[o]/self.range

  def samplePoly(self,n,x,norm=True):
    '''samples Legendre polynomial of order n at x, default normalized.'''
    if norm: return sps.eval_legendre(n,x)*np.sqrt((2.*n+1.)/2.)
    else: return sps.eval_legendre(n,x)


#============================================================================\
class normVar(var):
  '''
  A variable with normally-distributed uncertainty.
  Establishes basic parameters and defines sampling for MonteCarlo and stochastic collocation.
  '''
  def setUnc(self):
    self.uncType='normal'
    if len(self.unc)!=2:
      print 'ERROR: Normal-uncertainty variables require mean and variance, got',self.unc
      sys.exit()
    self.mean=self.unc[0]
    self.var=self.unc[1]

  def setDistr(self):
    self.average=self.mean
    self.quad=q.quadStatHermite(self.N)
    self.dist=spst.norm(self.average,self.var)
    self.domain=[float('-inf'),float('inf')]

  def sample(self,val=None):
    if val==None:
      Z=np.sqrt(-2.0*np.log(rand.random()))*np.cos(2.0*np.pi*rand.random())
      return np.sqrt(self.var)*Z+self.mean
    return self.average+np.sqrt(2*self.var)*val

  def samplePt(self,x):
    return self.var*x+self.average

  def revertPt(self,x):
    return (x-self.average)/self.var

  def sampleWt(self,o):
    return self.quad.weights[o]

  def samplePoly(self,n,x,norm=True):
    if norm: return sps.eval_hermitenorm(n,x)*\
        (np.sqrt(np.sqrt(2.*np.pi)*factorial(n)))**(-1)
    else: return sps.eval_hermitenorm(n,x)
