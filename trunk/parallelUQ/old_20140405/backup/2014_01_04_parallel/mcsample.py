import numpy as np
import scipy.stats as dists
import matplotlib.pyplot as plt
from bisect import bisect_left

xs=dists.uniform(0,2).rvs(1e6)
solns=np.exp(-xs*xs)
bins=100

low=min(solns)
hi=max(solns)
bounds=np.linspace(low,hi,bins+1)
bins=np.zeros(int(bins))
ctrs=np.zeros(len(bins))
for c in range(len(ctrs)):
  ctrs[c]=0.5*(bounds[c]+bounds[c+1])
for s in solns:
  i=bisect_left(bounds,s)
  bins[i-1]+=1
for b in range(len(bins)):
  bins[b]=float(bins[b])/float(len(solns))

plt.plot(ctrs,bins)
plt.title('Actual')
plt.show()
