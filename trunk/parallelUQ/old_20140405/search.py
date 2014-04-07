import numpy as np
from bisect import bisect_left

bounds=np.linspace(0,10,11)
bins=np.zeros(len(bounds)-1)

#check
#for b in range(len(bins)):
#  print 'bounds:',bounds[b],bounds[b+1]

for r in np.random.rand(10)*10:
  print 'random is',r
  i=bisect_left(bounds,r)
  print '  found index',i
  bins[i-1]+=1
print bins
