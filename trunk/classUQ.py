import numpy as np
import classQuadrature as q
import writeTransIn as w
import scipy.special as sps
import scipy.stats as spst
import os
import random as rand

class UQ:
  pass

class MC(UQ):
  def __init__(self,params,rands,solver,N=10000,stype=float):
    self.params=params
    self.case=str(params[0])
    self.randVars=rands
    self.solver=solver
    self.stype=stype
    self.N=N

  def run(self):
    self.soln=np.zeros(self.N)
    for n in range(self.N):
      if n in np.arange(0,self.N,self.N/10):
        print 'running',n
      rands=[]
      for r in range(len(self.randVars)):
        rands.append(self.randVars[r].sample())
      self.params[0]=self.case+'UQMC'
      for v,var in enumerate(self.randVars):
        if len(var.paramIndex)==1:
          self.params[var.paramIndex[0]]=rands[v]
        elif len(var.paramIndex)==3:
          self.params[var.paramIndex[0]][var.paramIndex[1]]\
              [var.paramIndex[2]]=rands[v]
      w.writeInput(self.params)
      phi,k=self.solver.run(self.params[0],doplot=False)
      self.soln[n]=k
      os.system('rm -f '+self.params[0]+'.*')
    import matplotlib.pyplot as plt
    plt.hist(self.soln)
    plt.show()


class SC(UQ):
  def __init__(self,params,rands,solver,stype=float):
    self.params=params
    self.case=str(params[0])
    self.randVars=rands
    self.solver=solver
    self.stype=stype #type for soln storage

    self.quad=None #multi-D quadrature for vars
    self.samples=np.array([]) #variable realizations
    self.coeffs=np.array([])

    self.trackPolys=np.zeros(len(self.randVars))
    #self.trackInvWtFunc=np.zeros(len(self.randVars))
    self.trackIndx=np.zeros(len(self.randVars))
    self.trackVarVals=np.zeros(len(self.randVars))
    self.trackWeights=np.zeros(len(self.randVars))
    #self.build_solution_matrix()
    #self.build_sample_matrix()
#    self.run()

  #UTILITIES
  #def scPoly(self,c,poly):
  #  print c,poly.c
  #  new_coeffs=poly.c*c
  #  return np.poly1d(new_coeffs)

  #FUNCTIONALITY
  def propUQ(self,order):
    # ASSUME all are uniform vars on -1,1 -> FIXME
    # make multiquad
    quadOrder=int(round(0.5+(order+1)/2))
    quads=[]
    for var in self.randVars:
      var.quad.order=quadOrder
      var.quad.setQuad()
      quads.append(var.quad)
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

  def calcCoeffLoop(self,counter,indx):
    counter-=1
    n=indx[counter]
    #print 'coeff',n
    var=self.randVars[counter]
    toRet = 0
    for el,absc in enumerate(var.ords):
      self.trackIndx[counter]=el
      #self.trackVarVals[counter]=var.samplePt(el)
      #self.trackWeights[counter]=var.sampleWt(el)
      #self.trackPolys[counter]=\
      #    self.evalNormLeg(n,var.quad.ords[el])
      #self.trackInvWtFunc[counter]=var.quad.invWtFunc(var.quad.ords[el])
      self.trackVarVals[counter]=var.UQval(el)
      self.trackWeights[counter]=var.UQwts(el)
      self.trackPolys[counter]=var.UQpolySample(n,el)
      #self.trackInvWtFunc[counter]=var.quad.invWtFunc(var.quad.ords[el])
      if counter==0: # each variable has a set value
        for v,val in enumerate(self.trackVarVals):
          var=self.randVars[v]
          self.params[var.paramIndex]=val
          #print 'x,ord',self.params,absc
        tot=self.solver.run(*self.params)
        #print tot
        toRet+= tot*np.prod(self.trackWeights)*np.prod(self.trackPolys)#*\
            #np.prod(self.trackInvWtFunc)
        #print self.trackWeights,self.trackPolys
      else: toRet+= self.calcCoeffLoop(counter,indx)
    #print 'toRet',toRet
    return toRet


  def UQsoln(self,vals):
    tot=0
    for indx,cof in np.ndenumerate(self.coeffs):
      temp=cof
      #if np.isnan(cof):print 'cof is nan!'
      for v,val in enumerate(vals):
        var=self.randVars[v]
        #temp*=sps.eval_legendre(indx[v],val)*self.norm(indx[v])
        temp*=self.evalNormShLeg(indx[v],val)
        #*self.randVars[v].range
        # FIXME range doesn't make sense for norm var
      tot+=temp
    return tot

  def evalNormLeg(self,o,x):
    return np.sqrt((2.*o+1.)/2.)*sps.eval_legendre(o,x)

  def evalNormShLeg(self,o,x):
    return np.sqrt((2.*o+1.))*sps.eval_sh_legendre(o,x)

