import matplotlib.pyplot as plt
import numpy as np
from scipy.stats import norm
from mpl_toolkits.mplot3d import Axes3D


# MONTE CARLO

m = 125

xs1 = np.random.normal(0,1,m)
xs2 = np.random.normal(0,1,m)
xs3 = np.random.normal(0,1,m)

fig = plt.figure()
ax = fig.add_subplot(111,projection='3d')
ax.scatter(xs1,xs2,xs3,marker='x')
ax.set_xlim([-4,4])
ax.set_ylim([-4,4])
ax.set_zlim([-4,4])
plt.savefig('mc_3d.pdf')

plt.figure()
plt.plot(xs1,xs2,'x')
plt.axis([-4,4,-4,4])
plt.savefig('mc_2d.pdf')

#GRID
from itertools import product
x = np.linspace(-4,4,int(m**(1./3.)))
x1,x2,x3 = zip(*product(x,repeat=3))

fig = plt.figure()
ax = fig.add_subplot(111,projection='3d')
ax.scatter(x1,x2,x3,marker='x')
ax.set_xlim([-5,5])
ax.set_ylim([-5,5])
ax.set_zlim([-5,5])
plt.savefig('grid_3d.pdf')

plt.figure()
plt.plot(x1,x2,'x')
plt.axis([-5,5,-5,5])
plt.savefig('grid_2d.pdf')

# LHS
xlin = np.linspace(-4,4,m)
x = []
y = []
z = []
for i in range(len(xlin)-1):
  r = xlin[i+1]-xlin[i]
  l = xlin[i]
  x.append(np.random.rand()*r + l)
  y.append(np.random.rand()*r + l)
  z.append(np.random.rand()*r + l)

fig = plt.figure()
ax = fig.add_subplot(111,projection='3d')
ax.scatter(x,y,z,marker='x')
ax.set_xlim([-5,5])
ax.set_ylim([-5,5])
ax.set_zlim([-5,5])
plt.savefig('lhs_3d.pdf')

plt.figure()
plt.plot(x,y,'x')
plt.axis([-5,5,-5,5])
plt.savefig('lhs_2d.pdf')
plt.show()
