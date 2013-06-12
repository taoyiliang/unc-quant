import numpy as np
import sys
import inspect
import random as rand
import classQuadrature as q
import scipy.special as sps
import scipy.stats as spst
from scipy.misc import factorial

#============================================================================\
class var():
  def __init__(self,unc,name='v',paramIndex=[],N=2):
    self.name=name
    self.uncType=''
    self.unc=unc #array of uncertainty info
    self.N=N
    self.paramIndex=paramIndex # so it can be changed in the list
    self.vals=[]
    self.weights=[]
    self.setUnc()
    self.setDistr()

  def setZtoSample(self,order,ps=[]):    
    quad=q.quadShiftLegendre(order) #quad for [0,1]
    toPoly=sps.sh_legendre
    def norm(o):
      return np.sqrt(2.*o+1)
    toDist=spst.uniform(0,1)
    toPPF=toDist.ppf
    coeffs=np.zeros(order)
    for o in range(order):
      coeffs[o]=np.sum(
            quad.weights*\
            self.dist.ppf(quad.ords)*\
            toPoly(o,*ps)(toPPF(quad.ords))*norm(o))
    sample=np.poly1d([])
    for c,cof in enumerate(coeffs):
      new_coeffs=toPoly(c,*ps).c*cof*norm(c)
      sample+=np.poly1d(new_coeffs)
    print 'sample in any\n',sample
    return sample,quad,toPoly,toDist,norm

  #TODO someday, translate onto other spaces besides uniform 0..1?

  def setStandardQuad(self,order):
    self.sampleOpt,self.quadOpt,self.polyOpt,self.distOpt,self.normOpt=self.setZtoSample(order)
    #For standard distros, overwrite this member in class

  def sample(self,val=None):
    if val==None:
      return self.dist.rvs()
    return self.dist.pdf(val)

  def sampleOpt(self,val=None):
    if val==None:
      return self.distOpt.rvs()
    return self.distOpt.pdf(val)
    
  def setUnc(self):
    pass
  def setDistr(self):
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
    self.low=float(self.unc[0])
    self.hi=float(self.unc[1])

  def setDistr(self):
    self.average=0.5*(self.low+self.hi)
    #self.quad=q.quadLegendre(self.N)
    self.range=(self.average-self.low)
    self.dist=spst.uniform(self.low,2.*self.range)
    self.domain=[self.low,self.hi]
    #for i in range(self.N):
    #  self.vals.append(self.average+self.range*self.quad.ords[i])
    #  self.vals.append(self.quad.weights[i])

#  def sample(self,val=None):
#    if val==None:
#      return self.low+rand.random()*(self.hi-self.low)
#    return self.average+val*self.range

#  def samplePt(self,x):
#    '''samples equivalent point on [-1,1]'''
#    return x*self.range+self.average

#FIXME do I really still need this?
  def revertPt(self,x):
    '''returns from [-1,1] to [a,b]'''
    #return (x-self.average)/(self.range)# [-1,1]
    return 0.5*(x-self.low)/self.range
    # this is analytic - numeric way to invert?

  def wtOpt(self,o):
    '''samples equivalent weight for [-1,1]'''
    return self.quadOpt.weights[o]/(2.*self.range)
    # -> self.range is analytic - way to get otherwise?
#
#  def samplePoly(self,n,x,norm=True):
#    '''samples Legendre polynomial of order n at x, default normalized.'''
#    if norm: return sps.eval_legendre(n,x)*np.sqrt((2.*n+1.)/2.)
#    else: return sps.eval_legendre(n,x)

  def sampleProbNorm(self,x):
    return self.quadOpt.probNorm(x)*self.range/0.5

#============================================================================\
class normVar(var):
  '''
  A variable with normally-distributed uncertainty.
  Establishes basic parameters and defines sampling.
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

#  def sample(self,val=None):
#    if val==None:
#      Z=np.sqrt(-2.0*np.log(rand.random()))*np.cos(2.0*np.pi*rand.random())
#      return np.sqrt(self.var)*Z+self.mean
#    return self.average+np.sqrt(2*self.var)*val

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

  def sampleProbNorm(self,x):
    return self.quad.probNorm(x)
    #return self.range*self.quad.probNorm(x)


#============================================================================\
class gammaVar(var):
  def setUnc(self):
    self.uncType='Gamma'
    if len(self.unc)!=3:
      print 'ERROR: Gamma-uncertainty variables require alpha,beta,left bound, got',self.unc
      sys.exit()
    self.alpha=self.unc[0]
    self.beta=self.unc[1]
    self.a=self.unc[2]

  def setDistr(self):
    self.average=self.alpha/self.beta
    self.quad=q.quadLaguerre(self.alpha,self.N)
    self.dist=spst.gamma(self.alpha,loc=self.a,scale=self.beta)
    self.domain=[self.a,float('inf')]

#  def sample(self,val=None):
#    if val==None:
#      return self.dist.rvs()
#    return self.dist.pdf(val)

  def samplePt(self,x):
    return x/self.beta+self.a

  def revertPt(self,x):
    return (x-self.a)*self.beta

  def sampleWt(self,o):
    return self.quad.weights[o] #*exp(a)?

  def samplePoly(self,n,x,norm=True):
    if norm: return sps.eval_genlaguerre(n,self.alpha,x)*\
        np.sqrt(factorial(n)/sps.gamma(n+self.alpha+1.0))
    else: return sps.eval_genlaguerre(n,self.alpha,x)

  def sampleProbNorm(self,x):
    return self.quad.probNorm(x)
    #return self.range*self.quad.probNorm(x)


#============================================================================\
class betaVar(var):
  def setUnc(self):
    self.uncType='Beta'
    if len(self.unc)!=4:
      print 'ERROR: Gamma-uncertainty variables require alpha,beta,left bound, got',self.unc
      sys.exit()
    self.alpha=self.unc[0]
    self.beta=self.unc[1]
    self.low=self.unc[2]
    self.hi=self.unc[3]

  def setDistr(self):
    self.average=(self.low+self.hi)*0.5
    self.range=self.average-self.low
    self.quad=q.quadJacobi(self.alpha,self.beta,self.N)
    self.dist=spst.beta(self.alpha,self.beta,scale=2.*self.range)
    self.domain=[self.low,self.hi]

  def samplePt(self,x):
    return x*self.range+self.average

  def revertPt(self,x):
    return (x-self.average)/self.range#(self.hi-self.low)

  def sampleWt(self,o):
    return self.quad.weights[o]/self.range

  def samplePoly(self,n,x,norm=True):
    a=self.alpha
    b=self.beta
    if norm: return sps.eval_jacobi(n,a,b,x)*\
        np.sqrt((2.*n+a+b+1.)/(2.**(a+b+1.)))*\
        np.sqrt(sps.gamma(n+a+b+1.)*factorial(n)/\
          (sps.gamma(n+a+1.)*sps.gamma(n+b+1.)))
    else: return sps.eval_jacobi(n,a,b,x)

  def sampleProbNorm(self,x):
    return self.quad.probNorm(x)*self.range
    #return self.range*self.quad.probNorm(x)





#============================================================================\
class triangVar(var):
  def setUnc(self):
    self.uncType='Triangular'
    if len(self.unc)!=3:
      print 'ERROR: Triangular-uncertainty variables require low,hi,frac, got',self.unc
      sys.exit()
    self.low=self.unc[0]
    self.hi=self.unc[1]
    self.frac=self.unc[2]

  def setDistr(self):
    self.center=(self.low+self.hi)*0.5
    self.range=self.center-self.low
    self.locHi=self.range*self.frac+self.low
    #FIXME - nonstandard quads
    # also, nonstandard polys - use quad's choice?
    # self.quad gets set only in zto sense
    
    self.dist=spst.triang(self.frac,loc=self.low,scale=2.*self.range)
    self.domain=[self.low,self.hi]

    # until zto, point to original functions
    self.samplePt=self.samplePtOrig
    self.revertPt=self.revertPtOrig

  def setZeroToOneDist(self,order):
    self.quad=q.quadShiftLegendre(order)
    self.toPPF=spst.uniform(0,1).ppf
    self.toPoly=sps.sh_legendre
    self.fromPPF=self.dist.ppf
    self.ztoCoeffs=np.zeros(order)
    for o in range(order):
      self.ztoCoeffs[o]=self.calcZtoCoeff(o)
    #redirect samples to new functions
    self.samplePt=np.poly1d([]) #new way to sample points
    for c,cof in enumerate(self.ztoCoeffs):
      new_coeffs=self.toPoly(c).c*cof*np.sqrt(2.*c+1)
      self.samplePt+=np.poly1d(new_coeffs)
 

  def calcZtoCoeff(self,o):
    return np.sum(self.quad.weights*self.fromPPF(self.quad.ords)*\
                  self.toPoly(o)(self.toPPF(self.quad.ords))*np.sqrt(2.*o+1.))
    
  def sampleOrig(self,val=None):
    if val==None:
      return self.dist.rvs()
    return self.dist.pdf(val)

  def samplePtOrig(self,x):
    return x*self.range+self.average

  def revertPtOrig(self,x):
    return (x-self.average)/self.range

  def sampleWtOrig(self,o):
    return self.quad.weights[o]/self.range

  def samplePoly(self,n,x,norm=True):
    #FIXME - nonstandard polys
    if norm: return sps.eval_jacobi(n,a,b,x)*\
        np.sqrt((2.*n+a+b+1.)/(2.**(a+b+1.)))*\
        np.sqrt(sps.gamma(n+a+b+1.)*factorial(n)/\
          (sps.gamma(n+a+1.)*sps.gamma(n+b+1.)))
    else: return sps.eval_jacobi(n,a,b,x)

  def sampleProbNorm(self,x):
    return self.quad.probNorm(x)*self.range
    #return self.range*self.quad.probNorm(x)
