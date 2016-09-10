import numpy as np
import matplotlib.pyplot as plt
from scipy.stats import norm


#plot the normal distribution

xs = np.linspace(-4.1,4.1,10000)
ys = norm.pdf(xs)

def plotBox(val,stdev,c,height):
  label = str(stdev)+r'$\sigma$'
  #fill
  ax = plt.gca()
  if stdev == 1:
    low = np.argmax(xs>-stdev)
    hi = np.argmin(xs<stdev)
    ax.fill_between(xs[low:hi],0,ys[low:hi],color=c,alpha=0.5)
    #left
    plt.plot([-stdev,-stdev],[0,height],c+'-')
    #right
    plt.plot([stdev,stdev],[0,height],c+'-')
    #top
    plt.plot([-stdev,stdev],[height,height],c+'-',label=label)
  else:
    low1 = np.argmax(xs>-stdev)
    low2 = np.argmax(xs>-(stdev-1))
    hi1 = np.argmin(xs<stdev)
    hi2 = np.argmin(xs<(stdev-1))
    ax.fill_between(xs[low1:low2],0,ys[low1:low2],color=c,alpha=0.5)
    ax.fill_between(xs[hi2:hi1],0,ys[hi2:hi1],color=c,alpha=0.5)
    height = ys[low2]
    #left
    plt.plot([-stdev,-stdev],[0,height],c+'-')
    #right
    plt.plot([stdev,stdev],[0,height],c+'-')
    #top
    plt.plot([-stdev,stdev],[height,height],c+'-',label=label)

pctDict = {1:68,2:95,3:99,4:99.99}
cDict = {1:'b',2:'g',3:'y',4:'r'}


plt.figure()
plt.plot(xs,ys,'k-',label='normal')

for stdev,val in pctDict.items():
  plotBox(val,stdev,cDict[stdev],(5-stdev)*0.4/4.0)

plt.legend(loc=0)
plt.axis([xs[0],xs[-1],0,0.41])
plt.savefig('stdev_pct.pdf')
plt.show()
