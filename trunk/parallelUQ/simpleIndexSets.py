import numpy as np
from itertools import product as allcombos

# each of these functions should take a list of order lists, as
# ordlist = [[a,b],[c,d]]

def IndexSetFactory(N,order,setType,impwts=None):
  okSetTypes = ['HC','TD','TP']#,'ETP'] TODO
  if setType not in okSetTypes:
    raise IOError('ERROR: index set type '+setType+' not recognized')

  orderlist=[]
  for n in range(N):
    orderlist.append(range(order+1))

  if impwts==None:
    impwts=list(np.ones(N))

  if setType=='HC':
    return HyperbolicCross(orderlist,order,impwts)
  if setType=='TD':
    return TotalDegree(orderlist,order,impwts)
  if setType=='TP':
    return TensorProduct(orderlist,order)


def TensorProduct(orderlist,maxorder=0):
  '''
  Given [[a,b],[c,d]], returns all possible cominations:
    [[a,c],[a,d],[b,c],[b,d]]
  '''
  return list(allcombos(*orderlist))


def TotalDegree(orderlist,maxorder,impwts):
  '''
  Given [[a,b],[c,d]], returns only possible cominations
    where for combination (p1,p2), p1+p2<=maxorder
  '''
  print '  ...level:',maxorder
  start = list(allcombos(*orderlist))
  end=[]
  tossed=0
  for entry in start:
    vals=np.array(entry)*np.array(impwts)
    if np.sum(vals)<=sum(impwts)/float(len(vals))*maxorder:
      end.append(entry)
    else:
      tossed+=1
  #    print 'Discarded index',entry
  print '  ...discarded',tossed,'indices from TP'
  return end


def HyperbolicCross(orderlist,maxorder,impwts):
  '''
  Given [[a,b],[c,d]], returns only possible cominations
    where for combination (p1,p2), p1*p2<=maxorder
  '''
  #TODO importance weighting
  #print 'starting HC index generation...'
  #print len(orderlist),len(orderlist[0])
  #print orderlist
  start = list(allcombos(*orderlist))
  target = (maxorder+1)**(sum(impwts)/float(len(impwts)))
  #print 'target:',target
  end=[]
  tossed=0
  for entry in start:
    tot=1
    for e,val in enumerate(entry):
      tot*=(val+1)**impwts[e]
    if tot<=target:
      end.append(entry)
    else:
      tossed+=1
  print '  ...discarded',tossed,'indices from TP'
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
