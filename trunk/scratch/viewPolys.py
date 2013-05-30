import numpy as np
import scipy.special as sps
import matplotlib.pyplot as plt

#legendre
xs=np.arange(-1,1.05,0.05)
leg=sps.legendre(1)
for i in range(2,5):
  leg+=sps.legendre(i)
ys=leg(xs)
plt.plot(xs,ys,label="legendre")
plt.legend()
plt.axis([-1,1,-5,5])
plt.show()
