import math
import sys

def uniform(x):
  #integrated -1..1 should give 2/3
  return x*x

def testOrthogonal(vartype,polyords=[2,4],quadord=3,target=0,verbose=False):
  import Variables as v
  var = v.newVar(vartype,'x')
  var.setDist([0,1.0])
  var.setQuadrature(order=quadord)
  tot = 0
  for l,val in enumerate(var.pts):
    tot+=var.evalNormPoly(val,polyords[0])*\
         var.evalNormPoly(val,polyords[1])*\
         var.wts[l]
  if abs(tot-target)<1.0e-9:
    if verbose:
      print 'Passed orders',polyords,'with quad order',quadord
    return True
  else:
    if verbose:
      print 'ERROR!  Failed orders',polyords,'with quad order',quadord
      print '    Got',tot,'instead of',target
    return False

def testUniformOrthogonals(maxord):
  ok = True
  for i in range(maxord):
    for j in range(i+1):
      if i==j: target=1
      else: target=0
      #print 'i,j',i,j,i==j
      quadord=int(math.ceil((i+j+1)/2.0))
      ok*=testOrthogonal('uniform',polyords=[i,j],quadord=quadord,target=target)
      if not ok:
        print 'ERROR: Failed on Legendre poly orders',i,j
        sys.exit()
  if ok:
    print 'Passed all orthogonal Legendre tests up to order [',i,',',j,']'
  else:
    print 'ERROR!  There were failures!'

def testNormalOrthogonals(maxord):
  ok = True
  for i in range(maxord):
    for j in range(i+1):
      if i==j: target=1
      else: target=0
      #print 'i,j',i,j,i==j
      quadord=int(math.ceil((i+j+1)/2.0))
      ok*=testOrthogonal('normal',polyords=[i,j],quadord=quadord,target=target,verbose=True)
      if not ok:
        print 'ERROR: Failed on Hermite poly orders',i,j
        sys.exit()
  if ok:
    print 'Passed all orthogonal hermite tests up to order [',i,',',j,']'
  else:
    print 'ERROR!  There were failures!'


if __name__=='__main__':
  testUniformOrthogonals(int(sys.argv[1]))
  testNormalOrthogonals(int(sys.argv[1]))
