import numpy as np
from itertools import product as allcombos

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

def makeCoeffs(N,indexset,verbose=True):
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

