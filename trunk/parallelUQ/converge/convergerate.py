from plot1 import addPlot
import matplotlib.pyplot as plt
import numpy as np

varNs = ['1xc2','1xf2',
         '2xc2','2xf2',
         '3xc2','3xf2',
         '4xc2','4xf2',
         '5D2','5xc2']

keys = ['b-o','b:o',
        'r-o','r:o',
        'g-o','g:o',
        'y-o','y:o',
        'k-o','k:o']

for v,var in enumerate(varNs):
  inpN = 'HC_h5_'+var+'.moments'
  addPlot(inpN,var,key=keys[v])

plt.legend(loc=1)
plt.axis([1,20,1e-11,1e-2])
plt.show()
