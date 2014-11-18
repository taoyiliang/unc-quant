import multiprocessing
from multiprocessing.queues import Queue as que

import numpy as np
import cPickle as pk

import SparseQuads
from tools import makePDF


class ROM():
  def __init__(self):
    pass

  def pdf(self,M=1000,bins=50):
    if self.numprocs==0:
      print 'Not found: number of processors to use.  \n'+\
            '    Syntax: LagrangeROM(<...>,numprocs=8)\n'+\
            '    Using 1 processor...'
      self.numprocs = 1
    procs=[]
    self.done=False
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
        self.done = True
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
      vals=np.zeros(len(self.varDict))
      for v,var in enumerate(self.varDict.values()):
        vals[v]=var.sample()
      samples[m]=self.sample(vals)
    self.romq.put(list(samples))


class HDMR_ROM(ROM):
  def __init__(self,ROMs,varDict):
    self.ROMs=ROMs
    self.varDict=varDict

  def numRunsToCreate(self):
    return sum(r.numRunsToCreate() for r in self.ROMs.values())

  def serializable(self):
    store=[]
    storeROMs=[]
    for rom in self.ROMs.values():
      storeROMs.append(rom.serializable())
    store.append(storeROMs)
    keys=self.varDict.keys()[:]
    store.append(keys)
    vals=self.varDict.values()[:]
    storevals=[]
    for v in vals:
      storevals.append(v.serializable())
    store.append(storevals)
    return store

  @classmethod
  def unserialize(self,store):
    #roms
    roms=[]
    storeROMs=store[0]
#    for entry in storeROMs:
#TODO FIXME
#      newrom
    #vardict

  def sample(self,xs,lvl,verbose=False):
    samples={}
    ref=0
    for i in range(lvl+1):
      for romk,romv in self.ROMs.iteritems():
        if len(romv.varDict)!=i:continue
        if verbose: print 'i,romk:',i,romk
        samples[romk]=romv.sample(xs)
        if i==0:
          ref=samples[romk]
        else:
          samples[romk]-=ref
        for j in range(1,i):
          #if verbose:print '  trying j =',j
          for sampk,sampv in samples.iteritems():
            if len(sampk.split('_'))-1!=j:continue
            if all(s in romk.split('_') for s in sampk.split('_'))\
                and sampk.split('_')[-1]!='':
                #or sampk.split('_')[-1]=='':#len(sampk.split('_'))==1:
              if verbose: print '  subtracting',sampk
              samples[romk]-=sampv
    tot=sum(samples.values())
    if verbose: print 'total',tot,'\n'
    return samples,tot

#  def mean(self,lvl,verbose=False):
#    samples={}
#    ref=0
#    if verbose: print '\nMEAN'
#    for i in range(lvl+1):
#      for romk,romv in self.ROMs.iteritems():
#        if len(romv.varDict)!=i:continue
#        if verbose: print 'i,romk:',i,romk
#        samples[romk]=romv.moment(1)
#        if i==0:
#          ref=samples[romk]
#        else:
#          samples[romk]-=ref
#        for j in range(1,i):
#          #if verbose:print '  trying j =',j
#          for sampk,sampv in samples.iteritems():
#            if len(sampk.split('_'))-1!=j:continue
#            if all(s in romk.split('_') for s in sampk.split('_'))\
#                and sampk.split('_')[-1]!='':
#                #or sampk.split('_')[-1]=='':#len(sampk.split('_'))==1:
#              if verbose: print '  subtracting',sampk
#              samples[romk]-=sampv
#    tot=sum(samples.values())
#    if verbose: print 'mean',tot,'\n'
#    return tot

  def moment(self,lvl,r=1,verbose=False): #r is moment, defaults to mean
    samples={}
    ref=0
    if verbose: print '\nMEAN'
    for i in range(lvl+1):
      for romk,romv in self.ROMs.iteritems():
        if len(romv.varDict)!=i:continue
        if verbose: print 'i,romk:',i,romk
        samples[romk]=romv.moment(r)
        if i==0:
          ref=samples[romk]
        else:
          samples[romk]-=ref
        for j in range(1,i):
          #if verbose:print '  trying j =',j
          for sampk,sampv in samples.iteritems():
            if len(sampk.split('_'))-1!=j:continue
            if all(s in romk.split('_') for s in sampk.split('_'))\
                and sampk.split('_')[-1]!='':
                #or sampk.split('_')[-1]=='':#len(sampk.split('_'))==1:
              if verbose: print '  subtracting',sampk
              samples[romk]-=sampv
    tot=sum(samples.values())
    if verbose: print 'moment',r,':',tot,'\n'
    return tot



class LagrangeROM(ROM):
  def __init__(self,solns,weights,probs,varDict,varVals,indexSet,quadrule,numprocs=0):
    self.solns = solns
    self.weights = weights
    self.probs = probs
    self.varDict = varDict
    self.varVals = varVals
    self.indexSet = indexSet
    self.quadrule = quadrule
    self.numprocs=numprocs
    self.ptsets=[]

  def case(self):
    case=''
    for v in self.varDict.keys():
      case+=v+'_'
    return case[:-1]

  def numRunsToCreate(self):
    return len(self.solns)

  def moment(self,r):
    tot=0
    for s,soln in enumerate(self.solns):
      tot+=self.weights[s] * soln**r #* self.probs[s]
    tot*=1.0/sum(self.weights)
    return tot

  def sample(self,xs,verbose=False):
    tot=0
    if len(self.ptsets)<1:
      self.makeSparseGrid()
    for c,pts in enumerate(self.ptsets):
      tptot = 0
      #get distinct points
      pts_by_var = {}
      for v,(name,var) in enumerate(self.varDict.iteritems()):
        pts_by_var[name]=[]
        for pt in pts:
          if pt[v] not in pts_by_var[name]:
            pts_by_var[name].append(pt[v])
      for pt in pts:
        pt = tuple(np.around(pt,decimals=15))
        slnidx = self.varVals.index(list(pt))
        soln = self.solns[slnidx]
        prod = 1
        for v,(key,var) in enumerate(self.varDict.iteritems()):
          polyeval = var.lagrange(pt[v],xs[key],pts_by_var[key],verbose)
          prod*=polyeval
        tptot += prod*soln
      tot+=tptot*self.cofs[c]
    return tot

  def makeSparseGrid(self):
    self.ptsets = []
    N = len(self.varDict)
    self.cofs = np.array(SparseQuads.makeCoeffs(N,self.indexSet,False))
    idxs = np.array(self.indexSet)
    survive = np.nonzero(self.cofs!=0)
    self.cofs = self.cofs[survive]
    idxs = idxs[survive]
    for j,cof in enumerate(self.cofs):
      idx = idxs[j]
      m = self.quadrule(idx)+1
      new = SparseQuads.tensorGrid(N,m,self.varDict,idx)
      pts = np.array(new[0])
      self.ptsets.append(pts)

#  def sample(self,xs,verbose=False):
#    varlist = self.varDict
#    tot=0
#    N=len(varlist)
#    #get coefficients, index points
#    cofs = np.array(SparseQuads.makeCoeffs(N,self.indexSet,verbose))
#    idxs = np.array(self.indexSet)
#    survive = np.nonzero(cofs!=0)
#    cofs=cofs[survive]
#    idxs=idxs[survive]
#    #for each idx point j, get little tensor product
#    for j,cof in enumerate(cofs):
#      idx = idxs[j]
#      m = self.quadrule(idx)+1 #TODO quadrule
#      new = SparseQuads.tensorGrid(N,m,varlist,idx)
#      pts = np.array(new[0])
#      #get distinct points for each variable
#      pts_by_var = list([] for x in range(len(self.varDict)))
#      for pt in pts:
#        for v,var in enumerate(self.varDict.values()):
#          if pt[v] not in pts_by_var[v]:
#            pts_by_var[v].append(pt[v])
#      #evaulate Lagrange product for each ordinate
#      if verbose:
#        print 'idx:',idx
#        print 'cof:',cof
#        print 'pts:',pts
#      tptot = 0
#      for pt in pts:
#        pt = tuple(np.around(pt,decimals=15))
#        #print self.histories['varVals']
#        slnidx = self.varVals.index(list(pt))
#        soln = self.solns[slnidx]
#        if verbose:
#          print '  pt:',pt
#          print '  soln:',soln
#        prod = 1
#        for v,(key,var) in enumerate(self.varDict.iteritems()):
#          #varpts = pts[:,v] #TODO this is a problem spot
#          #polyeval = var.lagrange(pt[v],xs[key],varpts,verbose)
#          polyeval = var.lagrange(pt[v],xs[key],pts_by_var[v],verbose)
#          prod*=polyeval
#          if verbose:
#            print '    var',key,' poly:',polyeval
#        if verbose:
#          print '  prod:',prod
#          print '  prod*soln:',prod*soln
#        tptot+=prod*soln
#      if verbose:
#        print 'tptot*cof:',tptot*cof,'\n'
#      tot+=tptot*cof
#    if verbose:
#      print 'romsample',xs,'total:',tot
#    return tot


  def serializable(self):
    store=[]
    store.append(self.ptsets)
    store.append(self.cofs)
    store.append(self.solns)
    store.append(self.weights)
    keys=self.varDict.keys()[:]
    store.append(keys)
    vals=self.varDict.values()[:]
    storevals=[]
    for v in vals:
      storevals.append(v.serializable())
    store.append(storevals)
    store.append(self.indexSet)
    #store.append(self.quadrule)
    store.append(self.numprocs)
    return store

  @classmethod
  def unserialize(cls,store):
    varDict={}
    for v in range(len(store[5])):
      varDict[store[5][v]]=store[6][v]
    ret = cls(store[3],store[4],varDict,store[7],store[8],store[9],store[10])
    ret.ptsets=store[0]
    ret.cofs=store[1]
    ret.quadrule=store[2]
    return ret

