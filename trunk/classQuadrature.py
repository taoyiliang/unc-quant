import numpy as np
import scipy.special.orthogonal as orth
import orthopoly as op
import scipy.stats as spst
import scipy.special as sps

#============================================================================\
class quadrature():
  '''
  Base class for developing quadrature of any dimension.
  '''
  def __init__(self,order=4,a=-1,b=1):
    self.quadType=''
    self.a=a
    self.b=b
    self.order=order
    self.ords=np.zeros(self.order)
    self.weights=np.zeros(self.order)
    self.setQuad()
    self.setDist()

  def setQuad(self):
    pass #not used in base class, defined in subclasses
  def getProb(self,x):
    pass #overwritten in subclasses
  def setDist(self):
    pass #overwritten

  def integrate(self,func,mult=1.0):
    result=0
    typecomp=np.array([1.0,2.0])
    if type(self.ords[0]) in [type(np.pi),type(typecomp[0])]:
      for n in range(len(self.ords)):
        result+=self.weights[n]*func(self.ords[n])
    elif type(self.ords[0])in [type(typecomp)]:
      for n in range(len(self.ords)):
        result+=self.weights[n]*func(*self.ords[n])
    else:
      print 'ordinate type was not np64 or nparray:',type(self.ords[0])
    return result*mult

  def integrateArray(self,arr,mult=1.0):
    assert(len(arr)==len(self.ords))
    result=0
    for n in range(len(self.ords)):
      result+=arr[n]*self.weights[n]
    return result


#============================================================================\
class quadLegendre(quadrature):
  def setQuad(self):
    self.quadType='Legendre'
    self.ords,self.weights = orth.p_roots(self.order)

  def setDist(self):
    self.dist=spst.uniform(-1,1)
    self.poly=sps.legendre

  def evNormPoly(self,o,x):
    return sps.eval_legendre(o,x)*np.sqrt((2.*o+1.)/2.)

  def wtFunc(self,x):
    return 1.0

  def invWtFunc(self,x):
    return 1.0


class quadShiftLegendre(quadrature):
  def setQuad(self):
    self.quadType='ShiftedLegendre'
    self.ords,self.weights=orth.ps_roots(self.order)

  def setDist(self):
    self.dist=spst.uniform(0,1)
    self.poly=sps.sh_legendre

  def evNormPoly(self,o,x):
    return sps.eval_sh_legendre(o,x)*np.sqrt(2.*o+1.)

  def wtFunc(self,x):
    return 1.0

  def invWtFunc(self,x):
    return 1.0






#============================================================================\
class quadHermite(quadrature):
  '''
  Covers range [-infty,infty], weight function is exp(-x^2)
  '''
  def setQuad(self):
    self.quadType='Hermite'
    self.ords,self.weights = orth.h_roots(self.order)

  def setDist(self):
    self.dist=spst.norm() #FIXME is this true?? exp(-x^2/<<2>>)
    self.poly=sps.hermite

  def wtFunc(self,x):
    return np.exp(-x*x)

  def invWtFunc(self,x):
    return np.exp(x*x)





class quadStatHermite(quadrature):
  def setQuad(self):
    self.quadType='Statistician Hermite'
    self.ords,self.weights = orth.he_roots(self.order)

  def setDist(self):
    self.dist=spst.norm()
    self.poly=sps.hermitenorm

  def wtFunc(self,x):
    return np.exp(-x**2/2.)

  def invWtFunc(self,x):
    return np.exp(x**2/2.)





#============================================================================\
class quadLobatto(quadrature):
  def setQuad(self):
    self.quadType='Lobatto'
    alpha,beta=op.rec_jacobi(n,0,0)
    self.ords,self.weights=op.lobatto(alpha,beta,lo,hi)

  def getProb(self,x):
    pass #FIXME







#============================================================================\
class quadChebyshev(quadrature):
  def setQuad(self):
    self.quadType='Chebyshev'
    self.ords,self.weights = orth.t_roots(self.order)

  def getProb(self,x):
    return 1./np.sqrt(1.-x*x)






#============================================================================\
class quadMulti(quadrature):
  '''
  Combines two or more quadratures to create ordinates, weights, integrate.
  '''
  def __init__(self,quads):
    self.quads=quads
    self.quadType=''
    self.order=[]
    self.ords=[]
    self.weights=[]
    self.trackOrd=np.zeros(len(quads))
    self.trackWts=np.zeros(len(quads))
    self.ordsSet=0

    for q in range(len(quads)):
      self.order.append(self.quads[q].order)
      self.quadType+=self.quads[q].quadType+str(self.order[q])
      if q<len(quads)-1:
        self.quadType+='-'

    # use recursive function calls to establish ordinates
    counter=len(self.quads)
    self._quadLoop(counter)
    self.trackOrd=np.array([0,0])

    # get tricky with dimensions
    self.ordl=np.array(self.ords)
    self.weights=np.array(self.weights)

    self.ordm=self.ordl.copy()
    self.ordm.resize(self.order+[len(self.quads)])
    self.wtsm=self.weights.copy()
    self.wtsm.resize(self.order)

  def _quadLoop(self,counter):
    quad=self.quads[counter-1]
    for n in range(quad.order):
      self.trackOrd[counter-1]=quad.ords[n]
      self.trackWts[counter-1]=quad.weights[n]
      if counter==1:
        self.ords.append(np.array(self.trackOrd[:]))
        self.weights.append(np.product(self.trackWts))
        self.ordsSet+=1
      else: self._quadLoop(counter-1)

  def setQuad(self):
    pass #overwrote __init__ instead






#######################
# set Lobatto weights, quadrature

def lobatto(alpha,beta,xl1,xl2):
  # alpha, beta are recursion coeffs
  #   These come from whichever quadrature, like for Leg,
  # xl1 and xl2 are assigned node locations
  from scipy.linalg import solve_banded, solve
  n = len(alpha)-1
  en=np.zeros(n)
  en[-1]=1
  A1=np.vstack((np.sqrt(beta),alpha-xl1))
  J1=np.vstack((A1[:,0:-1],A1[0,1:]))
  A2=np.vstack((np.sqrt(beta),alpha-xl2))
  J2=np.vstack((A2[:,0:-1],A2[0,1:]))
  g1=solve_banded((1,1),J1,en)
  g1=solve_banded((1,1),J2,en)
  C=np.array(((1,-g1[-1]),(1,-g2[-1])))
  x1=np.array((xl1,xl2))
  ab=solve(C,x1)

  alpha1=alpha
  alpha1[-1]=ab[0]
  beta1=beta
  beta1[-1]=ab[1]
  x,w = gauss(alpha1,beta1)
  return x,w

def gauss(alpha,beta):
  from scipy.linalg import eig_banded
  A=np.vstack((np.sqrt(beta),alpha))
  x,V=eig_banded(A,lower=False)
  w=beta[0]*scipy.real(scipy.power(V[0,:],2))
  return x,w
