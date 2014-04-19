import numpy as np
from itertools import product as allcombos

def BasicSparse(N,L,indexset,quadrule,varlist):
  c=np.array(makeCoeffs(N,indexset))
  print
  indexset=np.array(indexset)
  survive=np.nonzero(c!=0)
  print 'pre-c',indexset
  c=c[survive]
  indexset=indexset[survive]
  print 'survive',indexset
  SG=[] #holds list (multi-quad pt, wt)
  for j,cof in enumerate(c):
    idx = indexset[j]
    m = quadrule(idx)+1
    print 'm',m
    new = tensorGrid(N,m,varlist,idx)
    for i in range(len(new[0])):
      SG.append( [new[0][i],new[1][i]] )
      SG[-1][1]*=c[j]
  #TODO DEBUG check sum of weights
  return SG

def makeCoeffs(N,indexset):
  L=len(indexset)
  c=np.ones(L)
  indexset=np.array(indexset)
  for entry in indexset:
    indexset=np.array(indexset)
  for i in range(L):
    for j in range(i+1,L):
      d = indexset[j]-indexset[i]
      bln = d<=1
      bln*= d>=0
      c[i]+=(-1)**sum(d*bln)
  return c

def tensorGrid(N,m,varlist,idx):
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
  print 'quadwts:',quadwts
  print 'quadpts:',quadpts
  #print 'wts,pts',quadwts,quadpts
  return quadpts,quadwts


def Lagrange(j,x,pts):
  #TODO outdated.  In Variables.py now.
  prod=1
  for m,pt in enumerate(pts):
    if j!=m:
      newprod= (x-pts[m])/(pts[j]-pts[m])
      prod*=newprod
  return prod
