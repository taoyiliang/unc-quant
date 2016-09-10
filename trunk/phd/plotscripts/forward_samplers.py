import numpy as np
import matplotlib.pyplot as plt

#grid
h = 11

xs=range(h)
ys=range(h)
dx = 1.0/float(h-1)
for i in range(h-1):
  lh = [(i*dx,(i+1)*dx)]

plt.figure()
for i in range(h):
  for j in range(h):
    plt.plot(xs[i],ys[j],marker='o',color='b',markersize='10')
plt.title('Grid')
plt.savefig('grid.pdf')

plt.figure()
for i in range(11*11):
  x = np.random.rand()
  y = np.random.rand()
  plt.plot(x,y,marker='o',color='r',markersize='10')
plt.title('MonteCarlo')
plt.savefig('mc.pdf')

plt.figure()
dx = 1.0/float(h-1)
x2 = range(h)
y2 = range(h)
lowsx = list(x*dx for x in x2[:-1])
lowsy = list(y*dx for y in y2[:-1])
for low in lowsx:
  plt.plot([low,low],[0,1],'k-')
for low in lowsy:
  plt.plot([0,1],[low,low],'k-')
while len(lowsx)>0:
  lowx = np.random.choice(lowsx)
  lowy = np.random.choice(lowsy)
  lowsx.remove(lowx)
  lowsy.remove(lowy)
  x = lowx + np.random.rand()*dx
  y = lowy + np.random.rand()*dx
  plt.plot(x,y,marker='^',color='g',markersize='10')
plt.title('Stratified (LHS)')
plt.savefig('lhs.pdf')

plt.figure()
pts1 = [(0,8),(1,8),(2,8),(3,7),(3,6),(3,5),(4,4),(5,3),(6,3),(7,2),(8,1),(9,1),(10,1)]
pts2 = [(0,9),(1,9),(2,9),(3,8),(4,7),(4,6),(4,5),(5,4),(6,4),(7,3),(8,2),(9,2),(10,2)]
for p in pts1:
  plt.plot(0.1*p[0],0.1*p[1],'ro',markersize=10)
for p in pts2:
  plt.plot(0.1*p[0],0.1*p[1],'g^',markersize=10)
plt.title('Limit Surface')
plt.axis([0,1,0,1])
plt.savefig('limit.pdf')
plt.show()
