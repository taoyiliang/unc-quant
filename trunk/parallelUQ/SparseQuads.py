import numpy as np
from itertools import product as allcombos
from multiprocessing.queues import Queue as que
import multiprocessing
import scipy.weave as weave
#from memory_profiler import profile

def BasicSparse(N,L,indexset,quadrule,varlist):
  c=np.array(makeCoeffs(N,indexset))
  indexset=np.array(indexset)
  survive=np.nonzero(c!=0)
  #for i in range(len(indexset)):
  #  print indexset[i],c[i]
  #print 'pre-c',c#indexset
  c=c[survive]
  indexset=indexset[survive]
  print '  ...coefficient set established...'
  #print 'survive',indexset
  SG=[] #holds list (multi-quad pt, wt)
  for j,cof in enumerate(c):
    idx = indexset[j]
    #print '\nindex point:',idx
    m = quadrule(idx)+1
    #print '  m',m
    #print '  c',int(c[j])
    new = tensorGrid(N,m,varlist,idx)
    #print '  points and weights:'
    for i in range(len(new[0])):
      #print '    ',new[0][i],' | ',new[1][i]
      SG.append( [new[0][i],new[1][i]] )
      SG[-1][1]*=c[j]
  #TODO DEBUG check sum of weights
  #print '\n\n'
  return SG

def parBasicSparse(nump,N,L,indexset,quadrule,varlist):
  coeffMaker = CoeffMaker();
  c=np.array(coeffMaker.parMakeCoeffs(nump,N,indexset))
  indexset=np.array(indexset)
  survive=np.nonzero(c!=0)
  c=c[survive]
  indexset=indexset[survive]
  print '  ...coefficient set established...'
  SG=[]
  for j,cof in enumerate(c):
    idx = indexset[j]
    m = quadrule(idx)+1
    new = tensorGrid(N,m,varlist,idx)
    for i in range(len(new[0])):
      SG.append( [new[0][i],new[1][i]] )
      SG[-1][1]*=c[j]
  return SG

def oldmakeCoeffs(N,indexset,verbose=True):
  NI=len(indexset)
  c=np.zeros(NI)
  #set up index set as iset
  iset=indexset[:]
  #iset=np.array(indexset)
  #for i,entry in enumerate(iset):
  #  iset[i]=np.array(entry)
  #print 'index set:',iset
  #set up potential js as jset
  zerone=[0,1]
  sets=[]
  for n in range(N):
    sets.append(zerone)
  jset=list(allcombos(*sets))
  #jset=np.array(jset)
  for j,entry in enumerate(jset):
    jset[j]=np.array(entry)
  if verbose:
    print '\n\nmaking coeffs...'
  #TODO this is really slow.
  for i,ix in enumerate(iset):
    #print '  i:',i,ix
    ix=np.array(ix)
    for j,jx in enumerate(jset):
      #jx=np.array(jx)
      #print '    j:',j,jx
      comb = tuple(jx+ix)
      #print '  i,j,ix+jx',ix,jx,comb,comb in iset,
      if comb in iset:
        #print ' adding in',(-1)**sum(jx)
        c[i]+=(-1)**sum(jx)
      #else: print ' not adding.'
    #print 'final c[%i]:' %i,c[i]
    #for j in range(i+1,NI):
    #  d = indexset[j]-indexset[i]
    #  print '    d:',d
    #  if d.all()>=0 and d.all()<=1:
    #    c[i]+=(-1)**sum(d)
    #print '  c:',c[i]
  #print 'c:'
  #for i,cof in enumerate(c):
  #  print i,cof
  return c

def makeCoeffs(N,indexset,verbose=True):
  NI=len(indexset)
  c=np.zeros(NI)
  iset=indexset[:]
  zerone=[0,1]
  sets=[]
  for n in range(N):
    sets.append(zerone)
  #instead of a list, I want an iterator
  jiter=allcombos(*sets)
  #for j,entry in enumerate(jset):
    #jset[j]=np.array(entry)
  if verbose:
    print '\n\nmaking coeffs...'
  while True:
    try:
      jx = np.array(jiter.next())
      for i,ix in enumerate(iset):
        ix=np.array(ix)
        comb = tuple(jx+ix)
        if comb in iset:
          c[i]+=(-1)**sum(jx)
    except StopIteration:
      break
  return c

class CoeffMaker:
  def __init__(self):
    self.coefq = que() #TODO

  def parMakeCoeffs(self,nump,N,indexset,verbose=True):
    procs=[]
    done=False
    allStarted=False
    NI=len(indexset)
    c=np.zeros(NI)
    iset=indexset[:]
    jiter=allcombos([0,1],repeat=N)
    numdone=0
    numtodo = 2**N
    if verbose:
      print '\n\nmaking coeffs...'
    while not done:
      for p,proc in enumerate(procs):
        if not proc.is_alive():
          proc.join()
          while not self.coefq.empty():
            n,cofs = self.coefq.get()
            numdone+=n
            c+=cofs
          del procs[p]
          print '  ...finished jx %i/%i (%i pct)...\r'%(numdone,numtodo,100*float(numdone)/float(numtodo)),
      if allStarted and len(procs)==0:
        done=True
        break
      while len(procs)<nump and not allStarted:
        try:
          jxs=[]
          for i in range(500):
            jxs.append( np.array(jiter.next() ))
          procs.append(multiprocessing.Process(target=self.makeCoeffBatch,args=[jxs,NI,iset]))
          procs[-1].start()
        except StopIteration:
          allStarted=True
          if len(jxs)>0:
            procs.append(multiprocessing.Process(target=self.makeCoeffBatch,args=[jxs,NI,iset]))
            procs[-1].start()
          break
    return c

  def makeCoeffBatch(self,jxs,NI,iset):
    c=np.zeros(NI)
    for jx in jxs:
      for i,ix in enumerate(iset):
        ix=np.array(ix)
        comb = tuple(jx+ix)
        if comb in iset:
          c[i]+=(-1)**sum(jx)
    self.coefq.put([len(jxs),c])

  def makeCoeffChild(self,jx,NI,iset):
    c=np.zeros(NI)
    for i,ix in enumerate(iset):
      ix=np.array(ix)
      comb = tuple(jx+ix)
      if comb in iset:
        c[i]+=(-1)**sum(jx)
    self.coefq.put(c)

def tensorGrid(N,m,varlist,idx):
  #print '\n',idx
  quadptlists=[]
  quadwtlists=[]
  for n in range(N):
    mn = m[n]
    #print 'mn:',mn
    var=varlist[varlist.keys()[n]]
    var.setQuadrature(mn)
    quadptlists.append(var.pts)
    quadwtlists.append(var.wts)
  quadpts = list(allcombos(*quadptlists))
  quadwts = list(allcombos(*quadwtlists))
  for k,wtset in enumerate(quadwts):
    quadwts[k]=np.product(wtset)
  #for i in range(len(quadpts)):
  #  print quadpts[i],quadwts[i]
  #print 'wts,pts',quadwts,quadpts
  return quadpts,quadwts

def removeDuplicates(grid,tol=15):
  ptwt={}
  for e,entry in enumerate(grid):
    npt = tuple(np.around(entry[0],decimals=tol))
    if npt in ptwt.keys():
      ptwt[npt]+=entry[1]
    else:
      ptwt[npt]=entry[1]
  return ptwt

