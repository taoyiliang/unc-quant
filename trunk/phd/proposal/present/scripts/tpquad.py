import matplotlib.pyplot as plt
from numpy import polynomial as p

xpts = p.legendre.leggauss(5)[0]
ypts = p.legendre.leggauss(5)[0]

for x in xpts:
  print 'x:',x
  for y in ypts:
    print 'y:',y
    plt.plot(x,y,'bo',markersize=10)

plt.title('Tensor Quadrature')
plt.savefig('tpquad.pdf')

plt.show()
