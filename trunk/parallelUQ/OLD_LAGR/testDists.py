import scipy.stats as dists
import matplotlib.pyplot as plt
import numpy as np

f=5.0

norm = dists.norm(0,1)
beta1= dists.beta(12,12,loc=-1.*f,scale=2.*f)
beta2= dists.beta(1,2)
beta3= dists.beta(2,1)
beta4= dists.beta(1,1)

x=np.linspace(-10,10,1000)
plt.plot(x,norm.pdf(x),'k-')
plt.plot(x,beta1.pdf(x),'b:')
#plt.plot(x,beta2.pdf(x),'g:')
#plt.plot(x,beta3.pdf(x),'r:')
#plt.plot(x,beta4.pdf(x),'k:')
plt.show()
