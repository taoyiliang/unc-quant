import Variables as v
import matplotlib.pyplot as plt
import numpy as np

test=v.newVar('lognormal')
test.setDist([10.,4.])
vals=test.sample(n=int(1e6))
hist,bounds=np.histogram(vals,bins=100)
ctrs=np.zeros(len(bounds)-1)
for c in range(len(ctrs)):
  ctrs[c]=0.5*(bounds[c]+bounds[c+1])

plt.plot(ctrs,hist)
plt.show()
print np.average(vals)
