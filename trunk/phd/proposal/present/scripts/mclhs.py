import matplotlib.pyplot as plt
import numpy as np

samps = 10
invl = 1./float(samps)

sx = []
sy = []
for i in range(samps):
  sx.append(np.random.beta(10,10))
  sy.append(np.random.beta(10,10))

plt.figure()
plt.plot(sx,sy,'bo',markersize=10)
plt.title('Monte Carlo')
plt.axis([0,1,0,1])
plt.savefig('mc.pdf')

plt.figure()
grid = np.linspace(0,1,samps+1)
for l in grid:
  plt.plot( [l,l],[0,1],'k')
  plt.plot( [0,1],[l,l],'k')

sx = []
sy = []
for lx in grid[:-1]:
  sx.append( invl*np.random.random_sample()+lx)
  sy.append( invl*np.random.random_sample()+lx)

np.random.shuffle(sy)
plt.plot(sx,sy,'bo',markersize=10)
plt.title('LHS')
plt.savefig('lhs.pdf')

plt.show()
