import scipy.stats as dists
import numpy as np
import matplotlib.pyplot as plt


def testplot(xs,a,b,c=0,d=0.5):
  temp=dists.beta(a,b,loc=c,scale=d)
  ys=[]
  for x in xs:
    ys.append(temp.pdf(x))
  plt.plot(xs,ys,label='%i,%i' %(a,b))


# do battery of tests
xs=np.linspace(-4,4,1000)
nums=range(1,7)
for i in nums:
  #for j in nums:
  testplot(xs,i,i,-3,6)

norm = dists.norm()
ys=[]
for x in xs:
  ys.append(norm.pdf(x))
plt.plot(xs,ys,'k:',label='norm')

plt.legend()
plt.show()


