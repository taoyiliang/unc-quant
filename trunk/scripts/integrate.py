import scipy.special.orthogonal as quads

def func(x):
  return x*x

pts,wts = quads.p_roots(1000)
actual = 2./3.

tot=0
for l in range(len(pts)):
  tot+=func(pts[l])*wts[l]

print 'total:',tot
print 'err:',abs(tot-actual)/actual
