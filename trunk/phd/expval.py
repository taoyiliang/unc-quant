import numpy as np
import matplotlib.pyplot as plt
import scipy.stats as dists

def f(x):
  return x

#dist = dists.uniform(3,4)
dist = dists.norm(0,1)
print 'mean',dist.mean(),
print 'stdev',dist.std()
tot1=0
tot2=0
N = 0
for i in range(10):
  print '  ',i
  num=int(1e7)
  N+=num
  rands = dist.rvs(num)
  tot1+=sum(rands)
  tot2+=sum(rands*rands)

mom1 = tot1/N
mom2 = tot2/(N*(N-1))

print '1',mom1
print '2',mom2
