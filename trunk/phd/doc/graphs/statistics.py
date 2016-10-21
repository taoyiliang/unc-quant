import numpy as np
import matplotlib.pyplot as plt
import scipy.stats as dists


xs=np.linspace(-4,4,5000)

#gaussian, different mean
plt.figure()
plt.subplot(2,2,1)
#standard
ys = dists.norm.pdf(xs,loc=0)
plt.plot(xs,ys,'b',label='$\mu$=0')
ax = plt.gca()
ax.fill_between(xs,0,ys,color='b',alpha=0.5)
#shifted
offset = 1
ys = dists.norm.pdf(xs,loc=offset)
plt.plot(xs,ys,'r',label='$\mu$='+str(offset))
ax = plt.gca()
ax.fill_between(xs,0,ys,color='r',alpha=0.5)
ax.xaxis.set_major_formatter(plt.NullFormatter())
ax.yaxis.set_major_formatter(plt.NullFormatter())
#plt.legend(loc=0)
plt.title('Mean')

#gaussian, different sigma
plt.subplot(2,2,2)
#plt.figure()
#standard
ys = dists.norm.pdf(xs,loc=0)
plt.plot(xs,ys,'b',label='$\sigma$=1')
ax = plt.gca()
ax.fill_between(xs,0,ys,color='b',alpha=0.5)
#shifted
offset = 1.5
ys = dists.norm.pdf(xs,scale=offset)
plt.plot(xs,ys,'r',label='$\sigma$='+str(offset))
ax = plt.gca()
ax.fill_between(xs,0,ys,color='r',alpha=0.5)
ax.xaxis.set_major_formatter(plt.NullFormatter())
ax.yaxis.set_major_formatter(plt.NullFormatter())
#plt.legend(loc=0)
plt.title('Variance')

xs=np.linspace(0,1,1000)
#beta, different skewness
plt.subplot(2,2,3)
#plt.figure()
a = 3
b = 3
ys = dists.beta.pdf(xs,a,b)
plt.plot(xs,ys,'b',label='$\gamma_1$=0')
ax = plt.gca()
ax.fill_between(xs,0,ys,color='b',alpha=0.5)
a = 6
b = 3
ys = dists.beta.pdf(xs,a,b)
plt.plot(xs,ys,'r',label='$\gamma_1$=0.41')
ax = plt.gca()
ax.fill_between(xs,0,ys,color='r',alpha=0.5)
ax.xaxis.set_major_formatter(plt.NullFormatter())
ax.yaxis.set_major_formatter(plt.NullFormatter())
#plt.legend(loc=0)
plt.title('Skewness')

from scipy.special import gamma
def p7(x,k):
  k = float(k)
  m = 5./2. + 3./k
  a = np.sqrt(2. + 6./k)
  return gamma(m)/(a*np.sqrt(np.pi)*gamma(m-0.5)) * (1. + (x/a)**2)**(-m)

xs=np.linspace(-4,4,5000)
#pearson VII, different kurtosis
plt.subplot(2,2,4)
#plt.figure()
k=5e-2
ys = p7(xs,k)
plt.plot(xs,ys,'b',label='$\gamma_2$=1e-5')
ax = plt.gca()
ax.fill_between(xs,0,ys,color='b',alpha=0.5)
ax.xaxis.set_major_formatter(plt.NullFormatter())
ax.yaxis.set_major_formatter(plt.NullFormatter())
#other
k=1e30
ys = p7(xs,k)
plt.plot(xs,ys,'r',label='$\gamma_2$=1e10')
ax = plt.gca()
ax.fill_between(xs,0,ys,color='r',alpha=0.5)
#plt.legend(loc=0)
plt.title('Kurtosis')
plt.savefig('change_stats.pdf')

plt.show()
