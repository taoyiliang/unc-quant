import matplotlib.pyplot as plt
import numpy as np
from scipy.stats import norm

#read in files
def readf(name):
  vals=[]
  for l,line in enumerate(file(name,'r')):
    if l==0: continue
    vals.append(float(line.split(',')[0]))
  return vals

equal = readf('grid_value.csv')
cdf = readf('grid_cdf.csv')

xl = [0.13333,0.26666]
x = [xl[0]]*11
plt.plot(equal,x,'o-',markersize=10,label='value')
x = [xl[1]]*11
plt.plot(cdf,x,'o-',markersize=10,label='cdf')

xs = np.linspace(-4,4,1000)
ys = norm.pdf(xs)
plt.plot(xs,ys,'k-')

plt.yticks(xl,['Value','CDF'])

plt.axis([-5,5,0,0.4])
plt.title('Example Grid Spacing')

plt.savefig('grid_spacing.pdf')

plt.show()
