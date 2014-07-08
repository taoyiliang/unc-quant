import numpy as np
from bisect import bisect_left

def makePDF(vals,bins=100):
    #print vals
    low = min(vals)
    hi = max(vals)
    bounds = np.linspace(low,hi,bins+1)
    bins=np.zeros(int(bins))
    #make bin centers
    ctrs=np.zeros(len(bins))
    widths=np.zeros(len(bins))
    for c in range(len(ctrs)):
      ctrs[c]=0.5*(bounds[c]+bounds[c+1])
      widths[c]=bounds[c+1]-bounds[c]
    #TODO vectorize this binning
    for v in vals:
      i=bisect_left(bounds,v)
      bins[i-1]+=1
    #normalize
    #for b in range(len(bins)):
    #  bins[b]=float(bins[b])/float(len(vals))/widths[b]
    return bins,ctrs
