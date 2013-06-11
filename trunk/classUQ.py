import numpy as np
import classQuadrature as q
import classVar as cv
import scipy.special as sps
import scipy.stats as spst
import os
import random as rand

class UQ:
  pass

class SC(UQ):
  def __init__(self,params,rands,solver,stype=float):
    self.params=params
    self.randVars=rands
    self.solver=solver
    self.stype=stype #type for soln storage

    self.quad=None #multi-D quadrature for vars
    self.samples=np.array([]) #variable realizations
    self.coeffs=np.array([])

    self.trackPolys=np.zeros(len(self.randVars))
    self.trackProbNorm=np.zeros(len(self.randVars))
    self.trackIndx=np.zeros(len(self.randVars))
    self.trackVarVals=np.zeros(len(self.randVars))
    self.trackWeights=np.zeros(len(self.randVars))

  #FUNCTIONALITY
  def propUQ(self,order):
    # make multiquad
    quadOrder=int(round(0.5+(order+1)/2))
    quads=[]
    #stVars=(cv.uniVar,cv.normVar,
    #         cv.betaVar),cv.gammaVar)
    for var in self.randVars:
      var.setStandardQuad(order) # also maps non-optimal uncertainties
      #if not isinstance(var,stVars):
      #  var.setZeroToOneDist(order)
      quads.append(var.quadOpt)
    self.quad=q.quadMulti(quads)
    self.coeffs=np.zeros_like(self.quad.wtsm)
    # get variables on correct domain
    #ppf=spst.uniform(0,1).ppf
    #poly=sps.sh_legendre
    #for v,var in enumerate(self.randVars):
    #  self.randVars[v].setZeroToOneSample(ppf,poly,quadOrder)
      #sample them as var.sampleZeroToOne(x), x on [0,1]
    # calculate coefficients
    for indx,val in np.ndenumerate(self.coeffs):
      counter=len(self.randVars) #tracks which variable is being eval'd
      self.coeffs[indx]=self.calcCoeffLoop(counter,indx)
    #print self.coeffs

  def calcCoeffLoop(self,counter,indx):
    counter-=1
    n=indx[counter]
    #print 'coeff',n
    var=self.randVars[counter]
    toRet = 0
    for el,absc in enumerate(var.quadOpt.ords):
      self.trackIndx[counter]=el
      #self.trackVarVals[counter]=var.samplePt(absc)
      self.params[var.paramIndex]=var.sampleOpt(absc)
      self.trackWeights[counter]=var.wtOpt(el)
      self.trackPolys[counter]=var.polyOpt(n)(absc)
      self.trackProbNorm[counter]=var.sampleProbNorm(absc)
      
      if counter==0: # each variable has a set value
        tot=self.solver.run(*self.params)
        toRet+= tot*np.prod(self.trackWeights)*np.prod(self.trackPolys)*\
            np.prod(self.trackProbNorm)
      else: toRet+= self.calcCoeffLoop(counter,indx)
        
    return toRet


  def UQsoln(self,vals):
    tot=0
    for indx,cof in np.ndenumerate(self.coeffs):
      temp=cof
      #if np.isnan(cof):print 'cof is nan!'
      for v,val in enumerate(vals):
        var=self.randVars[v]
        val=var.revertPt(val) #this might be difficult
        temp*=var.polyOpt(indx[v])(val)
      print val,cof, tot
      tot+=temp
    return tot
