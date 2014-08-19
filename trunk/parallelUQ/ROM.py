import multiprocessing
from multiprocessing.queues import Queue as que

import numpy as np
import cPickle as pk

import SparseQuads
from tools import makePDF


class ROM():
  def __init__(self):
    pass


class LagrangeROM(ROM):
  def __init__(self,solns,weights,varlist,varVals,indexSet,quadrule,numprocs=0):
    self.solns = solns
    self.weights = weights
    self.varlist = varlist
    self.varVals = varVals
    self.indexSet = indexSet
    self.quadrule = quadrule
    self.numprocs=numprocs
    self.ptsets=[]

  def moment(self,r):
    tot=0
    for s,soln in enumerate(self.solns):
      tot+=self.weights[s]*soln**r
    tot+=1.0/sum(self.weights)
    return tot

  def pdf(self,M=1000,bins=50):
    if self.numprocs==0:
      print 'Not found: number of processors to use.  \n'+
            '    Syntax: LagrangeROM(<...>,numprocs=8)\n'+
            '    Using 1 processor...'
      self.numprocs = 1
    procs=[]
    done=False
    starthist=0
    endhist=0
    samples=[]
    batch = min(1000,int(float(M)/self.numprocs))
    self.romq=que()
    while not self.done:
      for p,proc in enumerate(procs):
        if not proc.is_alive():
          proc.join()
          while not self.romq.empty():
            new = self.romq.get()
            samples += new
            endhist += len(new)
            print 'ROM pdf:',100*endhist/M,'% finished...                \r',
          del procs[p]
      if endhist >= M:
        done = True
      else:
        while len(procs)<self.numprocs and starthist<M:
          if starthist + batch <= M:
            m = batch
          else:
            m = M - starthist
          procs.append(multiprocessing.Process(target=self.batch,args=[m]))
          procs[-1].start()
          starthist+=m
    bins,ctrs = makePDF(samples,bins=bins)
    return ctrs,bins

  def batch(self,M):
    np.random.seed()
    samples = np.zeros(M)
    for m in range(int(M)):
      vals=np.zeros(len(self.varlist))
      for v,var in enumerate(self.varlist):
        vals[v]=var.sample()
      samples[m]=self.sample(vals)
    self.romq.put(list(samples))

  def sample(self,xs,verbose=False):
    tot=0
    if len(self.ptsets)<1:
      self.makeSparseGrid()
    for c,pts in enumerate(self.ptsets):
      tptot = 0
      for pt in pts:
        pt = tuple(np.around(pt,decimals=15))
        slnidx = self.varVals.index(list(pt))
        soln = self.solns[slnidx]
        prod = 1
        for v,var in enumerate(self.varlist):
          varpts = pts[:,v]
          polyeval = var.lagrange(pt[v],xs[v],varpts)
          prod*=polyeval
        tptot += prod*soln
      tot+=tptot*self.cofs[c]
    return tot

  def makeSparseGrid(self):
    self.ptsets = []
    N = len(self.varlist)
    self.cofs = np.array(SparseQuads.makeCoeffs(N,self.indexSet,False))
    idxs = np.array(self.indexSet)
    survive = np.nonzero(cofs!=0)
    self.cofs = self.cofs[survive]
    idxs = idxs[survive]
    for j,cof in enumerate(cors):
      idx = idxs[j]
      m = self.quadrule(idx)+1
      new = SparseQuads.tensorGrid(N,m,varlist,idx)
      pts = np.array(new[0])
      self.ptsets.append(pts)
