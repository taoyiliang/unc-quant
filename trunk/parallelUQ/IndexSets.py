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
  def rule(i):
    tot=0
    for j,p in enumerate(i):
      tot+=p*impwts[j]
    return tot<=sum(impwts)/float(len(i))*maxorder
  end = multiidxGenerator(len(orderlist),rule,maxorder)
  #print 'in total degree, returning'
  #for e in end:
  #  print e
  #import sys; sys.exit()
  return end


def HyperbolicCross(orderlist,maxorder,impwts):
  '''
  Given [[a,b],[c,d]], returns only possible cominations
    where for combination (p1,p2), p1*p2<=maxorder
  '''
  print '  ...level:',maxorder
  target = (maxorder+1)**(sum(impwts)/max(1,float(len(impwts))))
  def rule(i):
    tot=1;
    for e,val in enumerate(i):
      tot*=(val+1)**impwts[e]
    return tot<=target;
  end = multiidxGenerator(len(orderlist),rule,maxorder)
  print 'end:'
  for e in end:
    print e
  return end

def multiidxGenerator(N,rule,L,I=None,MI=None):
  if I==None: I=[]
  if MI==None: MI=[]
  if len(I)!=N:
    i=0
    while rule(I+[i]):
      MI = multiidxGenerator(N,rule,L,I+[i],MI)
      i+=1
  else:
    MI.append(tuple(I))
  return MI



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
