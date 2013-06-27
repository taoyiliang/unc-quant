'''
Created Feb 10, 2013

@author: talbpaul
'''
import numpy as np
import sys
import scipy.special.orthogonal as orth
import scipy.stats as spst
import scipy.special as sps
from scipy.misc import factorial
from itertools import product

#============================================================================\
class Quadrature():
  '''
  Base class for developing quadrature of any dimension.
  '''
  def __init__(self,order):
    self.type = '' #string of quadrature type
    self.order = order #quadrature order
    self.qps = np.zeros(self.order) #points in quadrature
    self.wts = np.zeros(self.order) #wts in quadrature
    self.setQuad() #sets up points and wts
    self.setDist() #sets up associated distribution

  def setQuad(self):
    raise IOError('Quadrature.setQuad must be overwritten by individual quadrature!')
  def setDist(self):
    raise IOError('Quadrature.setDist must be overwritten by individual quadrature!')

  def resetQuad(self,order):
    self.order=order
    self.setQuad()

  def integrate(self,func,mult=1.0):
    #overwritten in the multiquad class, this is for single quads
    result=0
    for n in range(len(self.qps)):
      result+=self.wts[n]*func(self.qps[n])
    return result*mult

  def integrateArray(self,arr,mult=1.0):
    # same for individual quads or multiquad
    assert(len(arr) == len(self.wts))
    result=sum(arr*self.wts)
    return result*mult


#============================================================================\
class Legendre(Quadrature):
  def __init__(self,order=4,a=-1,b=1):
    self.a = a
    self.b = b
    Quadrature.__init__(self,order=order)

  def setQuad(self):
    self.type='Legendre'
    self.qps,self.wts = orth.p_roots(self.order)

  def setDist(self):
    self.dist=spst.uniform(-1,1)
    self.poly=sps.legendre

  def evNormPoly(self,o,x):
    return sps.eval_legendre(o,x)*np.sqrt((2.*o+1.)/2.)



class ShiftLegendre(Quadrature):
  def setQuad(self):
    self.type='ShiftedLegendre'
    self.qps,self.wts=orth.ps_roots(self.order)

  def setDist(self):
    self.dist=spst.uniform(0,1)
    self.poly=sps.sh_legendre

  def evNormPoly(self,o,x):
    return sps.eval_sh_legendre(o,x)*np.sqrt(2.*o+1.)


#============================================================================\
class Hermite(Quadrature):
  def setQuad(self):
    self.type='Hermite'
    self.qps,self.wts = orth.h_roots(self.order)

  def setDist(self):
    self.dist=spst.norm() #FIXME is this true?? exp(-x^2/<<2>>)
    self.poly=sps.hermite

  def evNormPoly(self,o,x):
    return sps.eval_hermitenorm(o,x)/np.sqrt(np.sqrt(np.pi)*2.**o*factorial(o))




class StatHermite(Quadrature):
  def setQuad(self):
    self.type='StatisticianHermite'
    self.qps,self.wts = orth.he_roots(self.order)

  def setDist(self):
    self.dist=spst.norm()
    self.poly=sps.hermitenorm

  def evNormPoly(self,o,x):
    return sps.eval_hermitenorm(o,x)/np.sqrt(np.sqrt(2.*np.pi)*factorial(o))




#============================================================================\
class Laguerre(Quadrature):
  def __init__(self,alpha,order):
    self.alpha=alpha
    Quadrature.__init__(self,order)

  def setQuad(self):
    self.quadType='GenLaguerre'
    self.qps,self.wts = orth.la_roots(self.order,self.alpha)

  def setDist(self):
    self.dist=spst.gamma(self.alpha) #shift from [a,inf] to [0,inf]?
    self.poly=sps.genlaguerre

  def evNormPoly(self,o,x):
    return sps.eval_genlaguerre(o,self.alpha,x)/\
        np.sqrt(sps.gamma(o+self.alpha+1.)/factorial(o))

#============================================================================\
class Jacobi(Quadrature):
  def __init__(self,alpha,beta,order=4):
    self.alpha=alpha
    self.beta=beta
    Quadrature.__init__(self,order)

  def setQuad(self):
    self.quadType='Jacobi'
    self.qps,self.wts = orth.j_roots(self.order,self.alpha,self.beta)

  def setDist(self):
    self.dist=spst.beta(self.alpha,self.beta)
    self.poly=sps.jacobi

  def evNormPoly(self,o,x):
    a=self.alpha
    b=self.beta
    gam=sps.gamma
    return sps.eval_jacobi(o,self.alpha,self.beta,x)/\
        np.sqrt(2**(a+b+1.)/(2.*o+a+b+1.))/\
        np.sqrt(gam(o+a+1.)*gam(o+b+1.)/(factorial(o)*gam(o+a+b+a+1.)))




#============================================================================\
class MultiQuad(Quadrature):
  '''
  Combines two or more quadratures to create ordinates, wts, integrate.
  '''
  def __init__(self,quads):
    self.quads=quads
    self.type='multi-'
    self.order=[]
    self.totOrder=np.product([q.order for q in self.quads])
    self.qps=np.zeros(self.totOrder,dtype=tuple)
    self.qps.copy() #tuple corresponding to qps
    self.wts=np.zeros_like(self.qps)
    #lookup dictionaries
    self.dict_quads={}
    self.indx_ord = {}
    self.indx_weight={}
    self.ord_indx={}
    self.ord_weight={}
    #Quadrature.__init__(self)

    for q in range(len(quads)):
      self.order.append(self.quads[q].order)
      self.dict_quads[quads[q]]=q
      self.type+=self.quads[q].type+'('+str(self.order[q])+')'
      if q<len(quads)-1:
        self.type+='-'

    #BETTER, use itertools.product to get all possible combos from lists
    self.indices=list(product(*[range(len(quad.qps)) for quad in self.quads]))
    self.qps=list(product(*[quad.qps for quad in self.quads]))
    wts=list(product(*[quad.wts for quad in self.quads]))
    self.wts = list(np.product(w) for w in wts)

    #make set of dictionaries
    self.indx_ord=dict(zip(self.indices,self.qps))
    self.indx_weight=dict(zip(self.indices,self.wts))
    self.ord_indx=dict(zip(self.qps,self.indices))
    self.ord_weight=dict(zip(self.qps,self.wts))

  def integrate(self,func,mult=1.0):
    result=0
    for n in range(len(self.qps)):
      result+=self.wts[n]*func(*self.qps[n])
    return result*mult




#utilities
def part_ndenum(arr,lvl):
  try:
    arr=np.array(arr)
    lvl=int(lvl)
  except ValueError: raise IOError('Inputs to Quadrature.part_ndenum must be array,integer!')
  assert lvl<=len(arr.shape),'Requested enumeration level > array dimension!'
  idx=np.ndindex(arr.shape[:lvl])
  for i in idx:
    yield i,arr[i]
