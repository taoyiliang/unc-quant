import matplotlib.pyplot as plt
import numpy as np
from bisect import bisect_left
import scipy.stats as dists

#do samples
inFile = file('temp.out','r')

samples=[]
for line in inFile:
  samples.append(float(line))

low = min(samples)
hi = max(samples)

nbins = 50
bounds = np.linspace(low,hi,nbins+1)
ctrs = 0.5*(bounds[:-1]+bounds[1:])
bins = np.zeros(nbins)

for s in samples:
  i = bisect_left(bounds,s)
  bins[i-1]+=1

dist = dists.norm(loc=0.1003,scale=0.001)
fitdist=dists.norm(*dists.norm.fit(samples))
print dists.norm.fit(samples)
fitvals=fitdist.pdf(ctrs)
fitvals*=max(bins)/max(fitvals)
expvals=dist.pdf(ctrs)
expvals*=max(bins)/max(expvals)

plt.figure()
plt.subplot(2,1,1)
plt.plot(ctrs,bins,'b.')
#plt.plot(ctrs,fitvals,'go')
plt.plot(ctrs,expvals,'k:')
plt.title('Sampling frequency, %i samples, %i bins' %(len(samples),nbins))
plt.xlabel('Variable Value')
plt.ylabel('Number Samples')

#store varvals
varvals=samples[:]
#do solns

inFile = file('../diffusion/solns.out','r')

samples=[]
for line in inFile:
  samples.append(float(line))

low = min(samples)
hi = max(samples)

print 'soln range:',low,hi

#nbins = 50
bounds = np.linspace(low,hi,nbins+1)
ctrs = 0.5*(bounds[:-1]+bounds[1:])
bins = np.zeros(nbins)

newsamples=[]
for s in samples:
  i = bisect_left(bounds,s)
  try:
    bins[i-1]+=1
    newsamples.append(s)
  except IndexError: pass

samples=newsamples[:]

#dist = dists.norm(loc=0.1003,scale=0.01)
fitdist=dists.norm(*dists.norm.fit(samples))
print dists.norm.fit(samples)
fitvals=fitdist.pdf(ctrs)
fitvals*=max(bins)/max(fitvals)
#expvals=dist.pdf(ctrs)
#expvals*=max(bins)/max(expvals)

plt.subplot(2,1,2)
plt.plot(ctrs,bins,'b')
#plt.plot(ctrs,fitvals,'k:')
plt.title('Solution frequency, %i samples, %i bins' %(len(samples),nbins))
plt.xlabel('Solution Value')
plt.ylabel('Number Samples')

#store solutions
solns = samples[:]

plt.show()
