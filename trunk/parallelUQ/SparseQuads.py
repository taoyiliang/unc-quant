import numpy as np
from itertools import product as allcombos

def BasicSparse(N,L,indexset,quadrule,varlist):
  #re-type index set as numpy array for convenience
  indexset=np.array(indexset)
  #get coefficients
  c=np.array(makeCoeffs(N,indexset))
  #eliminate index pts with coefficients=zero
  survive=np.nonzero(c!=0)
  c=c[survive]
  indexset=indexset[survive]
  print '  ...coefficient set established...'
  #now make the sparse grid for remaining points
  SG=[] #holds list (multi-quad pt, wt) -> (array, scalar)
  for j,cof in enumerate(c):
    idx = indexset[j]
    #need p(i)+1 quad points for each variable
    m = quadrule(idx)+1
    #get list of new points, weights for index point i
    new = tensorGrid(N,m,varlist,idx)
    #append new points to existing points and weights, includes duplicates
    for i in range(len(new[0])):
      SG.append( [new[0][i],new[1][i]] )
      #multiply weights by coefficients, so I don't have to keep coefficients
      SG[-1][1]*=c[j]
  return SG

def makeCoeffs(N,indexset):
  NI=len(indexset) #cardinality of index set
  c=np.zeros(NI) #initialize c-values
  iset=indexset[:] #make a copy so I don't change it accidentally
  #make possible j vectors
  zerone=[0,1]
  sets=[] #a list of the possible values for each entry in j
  for n in range(N):
    sets.append(zerone)
  jset=list(allcombos(*sets))
  print '\n\nmaking coeffs...'
  #check all j to see if i+j in index set; if so, add into c[i]
  for i,ix in enumerate(iset):
    ix=np.array(ix)
    for j,jx in enumerate(jset):
      jx=np.array(jx)
      comb = tuple(jx+ix) #tuple makes it easier to check if i+j in iset
      if comb in iset:
        c[i]+=(-1)**sum(jx)
  return c

def tensorGrid(N,m,varlist,idx):
  quadptlists=[]
  quadwtlists=[]
  #get the quadrature points for each variable
  for n in range(N):
    #quadrature order for variable n
    mn = m[n]
    #var is the instance of the Variable class
    var=varlist[varlist.keys()[n]]
    var.setQuadrature(mn)
    quadptlists.append(var.pts)
    quadwtlists.append(var.wts)
  #now do a tensor product of pts and wts
  quadpts = list(allcombos(*quadptlists))
  quadwts = list(allcombos(*quadwtlists))
  #just need product of weights for each point, not each weight individually
  for k,wtset in enumerate(quadwts):
    quadwts[k]=np.product(wtset)
  return quadpts,quadwts

def removeDuplicates(grid,tol=15):
  #DANGER: this does round points to 15 decimal places.
  #This should be worth it for the speed, but maybe not always.
  ptwt={}
  for e,entry in enumerate(grid):
    npt = tuple(np.around(entry[0],decimals=tol))
    if npt in ptwt.keys():
      ptwt[npt]+=entry[1] #combine weights with existing point
    else:
      ptwt[npt]=entry[1] #add a new point-weight pair
  return ptwt

