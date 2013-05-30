import numpy as np
import sys
import random as rand
import classQuadrature as q
import scipy.special as sps
import scipy.stats as spst

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
  def sampleMC(self):
    pass
  def sampleSC(self):
    pass
  def sampleSCn(self):
    pass
  def samplePoly(self):
    pass


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
      return self.sampleMC()
    #return self.low+val*(self.hi-self.low)
    return self.average+val*self.range

  def sampleMC(self):
    return self.low+rand.random()*(self.hi-self.low)

  def sampleSC(self,n):
    return self.average+self.range*self.quad.ords[n]

  def samplePoly(self,n,x):
    return sps.eval_legendre(n,x)

  def samplePt(self,o):
    '''samples equivalent point on [-1,1]'''
    return self.quad.ords[o]*self.range+self.average

  def sampleWt(self,o):
    '''samples equivalent weight for [-1,1]'''
    return self.quad.weights[o]/self.range


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
    self.dist=spst.norm(self.average,np.sqrt(self.var))
    self.domain=[float('-inf'),float('inf')]

  def sampleMC(self):
    Z=np.sqrt(-2.0*np.log(rand.random()))*np.cos(2.0*np.pi*rand.random())
    return np.sqrt(self.var)*Z+self.mean

  def sampleSC(self,n):
    #return self.average+self.var*self.quad.ords[n]
    return self.average+np.sqrt(2.*self.var)*self.quad.ords[n]

  def sample(self,val=None):
    if val==None:
      return self.sampleMC()
    return self.average+np.sqrt(self.var)*val

  def samplePt(self,o):
    return self.quad.ords[o] #FIXME

  def sampleWt(self,o):
    return self.quad.weights[o] #FIXME
