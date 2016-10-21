from scipy.special import orthogonal as quads
import matplotlib.pyplot as plt

ms = 10

#tensor
plt.figure()
x,w = quads.p_roots(4)
for i in x:
  for j in x:
    plt.plot([i],[j],'bo',markersize=ms)
plt.savefig('tensor_quad.pdf')

#smolyak
plt.figure()
x2,w2 = quads.p_roots(2)
plt.plot([0]*len(x),x,'bo',markersize=ms)
plt.plot(x,[0]*len(x),'bo',markersize=ms)
for i in x2:
  for j in x2:
    plt.plot([i],[j],'bo',markersize=ms)
plt.savefig('smolyak_quad.pdf')
