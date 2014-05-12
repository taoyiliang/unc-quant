import numpy as np
from itertools import product as allcombos

def BasicSparse(N,L,indexset,quadrule,varlist):
  c=np.array(makeCoeffs(N,indexset))
  indexset=np.array(indexset)
  survive=np.nonzero(c!=0)
  #print 'pre-c',indexset
  c=c[survive]
  indexset=indexset[survive]
  print '  ...coefficient set established...'
  #print 'survive',indexset
  SG=[] #holds list (multi-quad pt, wt)
  for j,cof in enumerate(c):
    idx = indexset[j]
    m = quadrule(idx)+1
    #print 'm',m
    new = tensorGrid(N,m,varlist,idx)
    for i in range(len(new[0])):
      SG.append( [new[0][i],new[1][i]] )
      SG[-1][1]*=c[j]
  #TODO DEBUG check sum of weights
  return SG

def makeCoeffs(N,indexset):
  NI=len(indexset)
  c=np.ones(NI)
  indexset=np.array(indexset)
  for e,entry in enumerate(indexset):
    indexset[e]=np.array(entry)
  #print indexset
  for i in range(NI):
    #print 'i:',i,indexset[i]
    for j in range(i+1,NI):
      #print '    j:',indexset[j]
      d = indexset[j]-indexset[i]
      #print 'd:',d
      if d.all()>=0 and d.all()<=1:
        c[i]+=(-1)**sum(d)
      #bln = d<=1
      #bln*= d>=0
      #print '    bln:',bln,d*bln
      #c[i]+=(-1)**sum(d*bln)
    #print '  c[i]',c[i]
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

