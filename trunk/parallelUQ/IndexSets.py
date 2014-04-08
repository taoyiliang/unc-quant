import numpy as np
from itertools import product as allcombos

# each of these functions should take a list of order lists, as
# ordlist = [[a,b],[c,d]]

def IndexSetFactory(N,order,setType):
  okSetTypes = ['HC','TD','TP']#,'ETP'] TODO
  if setType not in okSetTypes:
    raise IOError('ERROR: index set type '+setType+' not recognized')

  orderlist=[]
  for n in range(N):
    orderlist.append(range(order))

  if setType=='HC':
    return HyperbolicCross(orderlist,order)
  if setType=='TD':
    return TotalDegree(orderlist,order)
  if setType=='TP':
    return TensorProduct(orderlist,order)


def TensorProduct(orderlist,maxorder=0):
  '''
  Given [[a,b],[c,d]], returns all possible cominations:
    [[a,c],[a,d],[b,c],[b,d]]
  '''
  return list(allcombos(*orderlist))


def TotalDegree(orderlist,maxorder):
  '''
  Given [[a,b],[c,d]], returns only possible cominations
    where for combination (p1,p2), p1+p2<=maxorder
  '''
  start = list(allcombos(*orderlist))
  end=[]
  tossed=0
  for entry in start:
    if np.sum(entry)<=maxorder:
      end.append(entry)
    else:
      tossed+=1
  print 'Discarded',tossed,'indices from TP'
  return end


def HyperbolicCross(orderlist,maxorder):
  '''
  Given [[a,b],[c,d]], returns only possible cominations
    where for combination (p1,p2), p1*p2<=maxorder
  '''
  start = list(allcombos(*orderlist))
  end=[]
  tossed=0
  for entry in start:
    tot=1
    for e in entry:
      tot*=e+1
    if tot<=maxorder+1:
      end.append(entry)
    else:
      tossed+=1
  print 'Discarded',tossed,'indices from TP'
  return end

def chooseSet(orderlist,name,maxorder=0):
  if name=='TP': return TensorProduct(orderlist,maxorder)
  elif name=='TD': return TotalDegree(orderlist,maxorder)
  elif name=='HC': return HyperbolicCross(orderlist,maxorder)
  else: raise IOError('Set not recognized: '+name)

#for debugging
def makeOrdlist(ords):
  sets=[]
  for o in ords:
    sets.append(range(o+1))
  return sets
