'''
Created June 17, 2013
@author: talbpw
'''

import unittest as ut
import numpy as np
from test import test_support
import scipy.special as polys
import scipy.stats.distributions as dist

#========================================================================\
# Quadrature.py tests
import Quadrature
class QuadratureTest():
  def testQps(self):
    for i in range(len(self.qps)):
      self.assertAlmostEqual(self.qps[i],self.quad.qps[i],places=12)
  def testWts(self):
    for i in range(len(self.wts)):
      self.assertAlmostEqual(self.wts[i],self.quad.wts[i],places=12)
  def testIntegrate(self):
    self.assertAlmostEqual(self.quad.integrate(self.func),self.intSoln,places=12)
  def testIntegrateAndScale(self):
    self.assertAlmostEqual(self.quad.integrate(self.func,mult=1./self.intSoln),1.,places=12)
  def testIntegrateArray(self):
    solns=[self.func(qp) for qp in self.quad.qps]
    self.assertAlmostEqual(self.quad.integrateArray(solns,mult=1.0),self.intSoln,places=12)
  def testIntegrateArrayAndScale(self):
    solns=[self.func(qp) for qp in self.quad.qps]
    self.assertAlmostEqual(self.quad.integrateArray(solns,mult=1./self.intSoln),1.,places=12)
  def testPolynomial(self):
    for o,order in enumerate(self.orders):
      for p,pt in enumerate(self.polypts):
        self.assertAlmostEqual(self.quad.evNormPoly(order,pt),self.poly[o][p],places=10)


      


class QuadratureLegendreTests(QuadratureTest,ut.TestCase):
  def setUp(self):
    self.quad=Quadrature.Legendre(4)
    self.qps=[-0.86113631159405246,
              -0.33998104358485637,
               0.33998104358485631,
               0.86113631159405224]
    self.wts=[ 0.34785484513745374,
               0.65214515486254632,
               0.65214515486254632,
               0.34785484513745374]
    self.orders=[0,2,5]
    self.polypts=[-5.,0,2.5]
    self.poly=[[     0.70710678118654757, 0.70710678118654757,   0.70710678118654757],
               [    58.502136713115021  ,-0.79056941504209488,  14.032607116997184],
               [-55171.015374923103     , 0.0                ,1493.924902408605]]
    self.intSoln=2./3.
    def func(x): return x*x
    self.func=func

class QuadratureShiftLegendreTests(QuadratureTest,ut.TestCase):
  def setUp(self):
    self.quad=Quadrature.ShiftLegendre(4)
    self.qps=[ 0.069431844202973769,
               0.33000947820757176,
               0.66999052179242846,
               0.93056815579702645]
    self.wts=[ 0.17392742256872706,
               0.32607257743127327,
               0.32607257743127299,
               0.17392742256872709]
    self.orders=[0,2,5]
    self.polypts=[-5.,0,2.5]
    self.poly=[[       1.0             , 1.0               ,    1.0],
               [     404.72830392746198, 2.2360679774997898,   52.547597471245062],
               [-4167839.8594249035    ,-3.3166247903553998,24912.827112754589]]
    self.intSoln=1./3.
    def func(x): return x*x
    self.func=func

class QuadratureStatHermiteTests(QuadratureTest,ut.TestCase):
  def setUp(self):
    self.quad=Quadrature.StatHermite(4)
    self.qps=[-2.334414218338976,
              -0.74196378430272669,
               0.74196378430272669,
               2.3344142183389769]
    self.wts=[ 0.11499371468450574,
               1.1383204226309944,
               1.1383204226309944,
               0.11499371468450574]
    self.orders=[0,2,5]
    self.polypts=[-5.,0,2.5]
    self.poly=[[   0.63161877774606467, 0.63161877774606467, 0.63161877774606467],
               [  10.718926100856025  ,-0.44662192086900115, 2.3447650845622556],
               [-112.43435200249249   , 0.0                ,-1.2162369807961937]]
    self.intSoln=np.sqrt(2.*np.pi)
    def func(x): return x*x#*np.exp(-x*x/2.),weight function
    self.func=func

class QuadratureLaguerre1Tests(QuadratureTest,ut.TestCase):
  def setUp(self):
    self.quad=Quadrature.Laguerre(1,4)
    self.qps=[ 0.74329192798143184,
               2.5716350076462811,
               5.7311787516891046,
              10.953894312683193]
              
    self.wts=[0.44687059321877698,
              0.47763577236386767,
              0.074177784731052118,
              0.0013158496863032462]
              
    self.orders=[0,2,5]
    self.polypts=[-5.,0,2.5]
    self.poly=[[  1.0              ,1.0               , 1.0],
               [ 17.609183210283589,1.7320508075688774,-0.79385662013573532],
               [337.12803619597094 ,2.4494897427831783, 0.36306456039950313]]
    self.intSoln=6.
    def func(x): return x*x#*np.exp(-x*x/2.),weight function
    self.func=func


class QuadratureLaguerre3Tests(QuadratureTest,ut.TestCase):
  def setUp(self):
    self.quad=Quadrature.Laguerre(3,4)
    self.qps=[ 1.7555216471854935,
               4.2656058656568261,
               8.0579406831380034,
              13.920931804019688]
              
    self.wts=[1.8603340741465004,
              3.3568910190289261,
              0.76445397284351757,
              0.01832093398106218]
              
    self.orders=[0,2,5]
    self.polypts=[-5.,0,2.5]
    self.poly=[[  0.40824829046386307,0.40824829046386307, 0.40824829046386307],
               [  6.1322236314950764 ,1.2909944487358056 , 0.08068715304598785],
               [104.94689348234934   ,3.0550504633038935 ,-0.2569316761014272]]
    self.intSoln=120.
    def func(x): return x*x#*np.exp(-x*x/2.),weight function
    self.func=func

class QuadratureJacobi22Tests(QuadratureTest,ut.TestCase):
  def setUp(self):
    self.quad=Quadrature.Jacobi(2,2,4)
    self.qps=[-0.69474659060686594,
              -0.25056280708573192,
               0.2505628070857317,
               0.69474659060686628]
              
    self.wts=[0.10170944689820199,
              0.4316238864351315,
              0.43162388643513139,
              0.10170944689820231]
              
    self.orders=[0,2,5]
    self.polypts=[-5.,0,2.5]
    self.poly=[[       5.3033008588991057, 5.3033008588991057,    5.3033008588991057],
               [    1091.8413346269688   ,-6.2749501990055663,  268.25412100748804],
               [-1781866.2059148131      , 0.0               ,50525.670799426152]]
    self.intSoln=16./105.
    def func(x): return x*x#*np.exp(-x*x/2.),weight function
    self.func=func

class QuadratureJacobi35Tests(QuadratureTest,ut.TestCase):
  def setUp(self):
    self.quad=Quadrature.Jacobi(3,5,4)
    self.qps=[-0.48258873560366899,
              -0.065607351868878838,
               0.34598814963599311,
               0.70220793783655422]

    self.wts=[0.048974422820318206,
              0.36408082731446229,
              0.48455458955691122,
              0.11826317618132405]
              
    self.orders=[0,2,5]
    self.polypts=[-5.,0,2.5]
    self.poly=[[31.217483082401113 ,31.217483082401113 ,31.217483082401113],
               [10038.076027992116 ,-22.865776829139218,2020.7630272751783],
               [-23042129.310517259,-44.160960470618157,473685.64233801217]]
    self.intSoln=64./495.
    def func(x): return x*x#*np.exp(-x*x/2.),weight function
    self.func=func

#========================================================================\
# Distributions.py tests
import Distributions
class DistributionTests(): #base class tests
  #parameter tests
  def testLowBound(self):
    self.assertAlmostEqual(self.dist.distribution.interval(1)[0],self.lowbound)
  def testHiBound(self):
    self.assertAlmostEqual(self.dist.distribution.interval(1)[1],self.hibound)
  def testRange(self):
    self.assertAlmostEqual(self.dist.range,self.range)

  #pdf tests
  def testPDFInRange(self):
    for i,x in enumerate(self.inpdf_range):
      self.assertAlmostEqual(self.dist.distribution.pdf(x),self.pdf[i])
  def testPDFUnderRange(self):
    for i,x in enumerate(self.ltpdf_range):
      self.assertAlmostEqual(self.dist.distribution.pdf(x),self.ltpdf[i])
  def testPDFAboveRange(self):
    for i,x in enumerate(self.gtpdf_range):
      self.assertAlmostEqual(self.dist.distribution.pdf(x),self.gtpdf[i])
  def testMean(self):
    self.assertAlmostEqual(self.dist.distribution.mean(),self.mean)
  def testVariance(self):
    self.assertAlmostEqual(self.dist.distribution.var(),self.variance)

  #polynomial tests
  def testNorms(self):
    for order in range(len(self.polycoeff)):
      self.assertAlmostEqual(self.dist.poly_norm(order),self.norms[order])

  # test shifting to and from standard and actual distribution
  def testShiftToActualDistribution(self):
    for x in range(len(self.std_pts)):
      self.assertAlmostEqual(self.dist.actual_point(self.std_pts[x]),self.act_pts[x])

  def testShiftToStandardDistribution(self):
    for x in range(len(self.act_pts)):
      self.assertAlmostEqual(self.dist.std_point(self.act_pts[x]),self.std_pts[x])

  def testShiftedWeights(self):
    for i,w in enumerate(self.wtrange):
      self.assertAlmostEqual(self.dist.actual_weight(w),self.adjwts[i])
  def testProbabilityNorm(self):
    for i,x in enumerate(self.wtrange):
      self.assertAlmostEqual(self.dist.probability_norm(x),self.probnorm[i])

class DistributionUniformTests(DistributionTests,ut.TestCase):
  def setUp(self):
    self.dist=Distributions.Uniform()
    self.dist.low=-4.
    self.dist.hi=8.
    self.dist.inDistr()
    #var specifics for tests
    self.std_pts=[-1.0, -0.5, 0.0, 0.5, 1.0]
    self.act_pts=[-4.0, -1.0, 2.0, 5.0,  8.0]
    self.inpdf_range=np.linspace(self.dist.low,self.dist.hi,5)
    self.ltpdf_range=np.linspace(-20,self.dist.low-0.0001,5)
    self.gtpdf_range=np.linspace(self.dist.hi+0.0001,20,5)
    self.mean=2.
    self.variance=12.
    self.wtrange=np.linspace(0,5,20)
    self.adjwts=self.wtrange/6.
    #benchmarks
    self.lowbound=-4.
    self.hibound=8.
    self.range=12.
    self.pdf=[1./12.]*5
    self.ltpdf=[0.0]*5
    self.gtpdf=self.ltpdf
    self.polycoeff=[[1.0],
                    [1.0, 0.0],
                    [1.5, 0.0, -0.5],
                    [2.5, 0.0, -1.5, 0.0]]
    self.norms = np.sqrt(np.array([0.5,1.5,2.5,3.5]))
    self.probnorm=[12.]*len(self.wtrange)

class DistributionContNormalTests(DistributionTests,ut.TestCase):
  def setUp(self):
    self.dist=Distributions.Normal()
    self.dist.mean=3.
    self.dist.sigma=np.sqrt(2.)
    self.dist.inDistr()
    #var specifics for tests
    self.std_pts=[-1.0,0.0,1.0] #TODO how to pick other points appropriately?
    self.act_pts=[ 2.0,3.0,4.0] #FIXME take care of it when I make sure it matches mfs
    self.inpdf_range=np.linspace(-1,7,5)
    self.ltpdf_range=np.linspace(-100,-1,5)
    self.gtpdf_range=np.linspace(7,100,5)
    self.mean=3.
    self.variance=2.
    self.wtrange=np.linspace(0,5,20)
    self.adjwts=self.wtrange
    #benchmarks
    self.lowbound=float('-inf')
    self.hibound=float('inf')
    self.range=float('inf')
    self.pdf=[0.0051667463385230176,
              0.10377687435514871,
              0.28209479177387814,
              0.10377687435514871,
              0.0051667463385230176]
    self.ltpdf=[0.0,0.0,0.0,0.0,0.0051667463385230176]
    self.gtpdf=[0.0051667463385230176,0.0,0.0,0.0,0.0]
    self.polycoeff=[[1.0],
                    [1.0, 0.0],
                    [1.0, 0.0, -1.],
                    [1.0, 0.0, -3, 0.0]]
    self.norms = [0.63161877774606467,
                  0.63161877774606467,
                  0.44662192086900115,
                  0.25785728623970555]
    self.probnorm=[1.]*len(self.wtrange)
  def testRange(self): #overwrite, since range doesn't make exist for hermite
    pass

class DistributionGamma2Tests(DistributionTests,ut.TestCase):
  def setUp(self):
    self.dist=Distributions.Gamma()
    self.dist.low=1.
    self.dist.alpha=2.
    self.dist.beta=1.
    self.dist.inDistr()
    #var specifics for tests
    self.std_pts=[0.0] #TODO how to pick other points appropriately?
    self.act_pts=[3.0] #FIXME take care of it when I make sure it matches mfs
    self.inpdf_range=np.linspace(1,20,5)
    self.ltpdf_range=np.linspace(-100,1,5)
    self.gtpdf_range=np.linspace(20,100,5)
    self.mean=3.
    self.variance=2.
    self.wtrange=np.linspace(0,5,20)
    self.adjwts=self.wtrange
    #benchmarks
    self.lowbound=1.0
    self.hibound=float('inf')
    self.range=float('inf')
    self.pdf=[0.0,
              0.041095552214823014,
              0.00071109238393315542,
              9.2282318505751397e-06,
              1.0645313231320826e-07]
    self.ltpdf=[0.0,0.0,0.0,0.0,0.0]
    self.gtpdf=[1.0645313231320826e-07,0.0,0.0,0.0,0.0]
    self.polycoeff=[[1.0],
                    [-1.0 , 3.0],
                    [0.5  ,-4.0, 6.],
                    [-1./6., 2.5,-10,10.0]]
    self.norms = [0.70710678118654757,
                  0.40824829046386302,
                  0.28867513459481287,
                  0.22360679774997896]
    self.probnorm=[1.]*len(self.wtrange)
  def testRange(self): #overwrite, since range doesn't make sense
    pass
  def testPolynomialCoeffs(self):
    for order in range(len(self.polycoeff)):
      for e,entry in enumerate(self.polycoeff[order]):
        self.assertAlmostEqual(self.dist.polynomial(order,self.dist.alpha).c[e],entry)
  def testPolynomialType(self):
    self.assertTrue(isinstance(self.dist.polynomial(1,self.dist.alpha),polys.orthogonal.orthopoly1d))

#========================================================================\
# Samplers.py tests
import Samplers
class StochCollSamplerTests(ut.TestCase):
  def setUp(self):
    self.sc=Samplers.StochasticCollocation()
    self.sc.min_poly_order=2
    self.sc.generateQuadrature()

    self.sc.distDict={}

    self.sc.distDict['a']=Distributions.Uniform()
    self.sc.distDict['a'].low=-4.
    self.sc.distDict['a'].hi=8.
    self.sc.distDict['a'].inDistr()

    self.sc.distDict['b']=Distributions.Normal()
    self.sc.distDict['b'].mean=3.
    self.sc.distDict['b'].sigma=np.sqrt(2)
    self.sc.distDict['b'].inDistr()

    self.sc.var_poly_order={}
    self.sc.var_poly_order['a']=2
    self.sc.var_poly_order['b']=4

    self.sc.generateQuadrature()

    self.qps=[(-0.57735026918962573,-2.334414218338976),
              (-0.57735026918962573,-0.74196378430272669),
              (-0.57735026918962573, 0.74196378430272669),
              (-0.57735026918962573, 2.3344142183389769),
              ( 0.57735026918962573,-2.334414218338976),
              ( 0.57735026918962573,-0.74196378430272669),
              ( 0.57735026918962573, 0.74196378430272669),
              ( 0.57735026918962573, 2.3344142183389769)]

    self.wts=[0.114993714685,
              1.13832042263,
              1.13832042263,
              0.114993714685,
              0.114993714685,
              1.13832042263,
              1.13832042263,
              0.114993714685]


  def testMultiQuadPoints(self):
    for q,qp in enumerate(self.sc.quad.qps):
      for i,n in enumerate(qp):
        self.assertAlmostEqual(self.qps[q][i],n)

  def testMultiQuadWeights(self):
    for w,wt in enumerate(self.sc.quad.wts):
      self.assertAlmostEqual(self.wts[w],wt)

  def testGenerateInputs(self):
    for i in range(len(self.wts)):
      pass




import Models
class SimpleUQ(ut.TestCase):
  def setUp(self):
    self.sampler
    self.distDict



def test_main():
  test_support.run_unittest(
                            QuadratureLegendreTests,
                            QuadratureShiftLegendreTests,
                            QuadratureStatHermiteTests,
                            QuadratureLaguerre1Tests,
                            QuadratureLaguerre3Tests,
                            QuadratureJacobi22Tests,
                            QuadratureJacobi35Tests,
                            DistributionUniformTests,
                            DistributionContNormalTests,
                            DistributionGamma2Tests,
                            #DistributionGamma5Tests,
                           )

if __name__== '__main__':
  #test_main()
  ut.main(verbosity=2)
