import numpy as np
from SparseQuads import Lagrange
import matplotlib.pyplot as plt

#pts = np.random.rand(10)*0.1
pts = np.linspace(0,1,11)
print pts
xs = np.linspace(0,1,100)
for j in range(4):
  ys = Lagrange(j,xs,pts)
  plt.plot(xs,ys,label=str(j))
plt.legend()
plt.axis([0,1,-10,10])
plt.show()
