import numpy as np
import norm1var as solver
#import uni1var as solver
import classVar as var
import classQuadrature as q
from writeTransIn import writeInput
import classUQ as uq
import os

# set non-random variables, set None for randoms
x=None
p=[x]
# designate random variables
randVars=[]
randVars.append(var.normVar('x',[2,4],paramIndex=0,N=8))
#randVars.append(var.uniVar('x',[3,7],paramIndex=0,N=2))
var=randVars[0]

scsolve=uq.SC(p,randVars,solver)
scsolve.propUQ(8)
UQsoln=scsolve.UQsoln
#visualize solution
import matplotlib.pyplot as plt
zsa=np.arange(0,1.001,0.001)
zs=np.array(randVars[0].z21quad.ords)
xs=np.zeros_like(zsa)
for i in range(len(xs)):
  xs[i]=randVars[0].sampleZeroToOne(zsa[i])
print xs
usa=np.exp(-xs*xs)
#xs=np.arange(var.low,var.hi+var.range/10.,var.range/10.)
us=np.zeros(len(zs))
np.set_printoptions(suppress=True)
for i,z in enumerate(zs):
  us[i]=UQsoln([z])
  print np.array([z,us[i]])
plt.plot(zs,us)
plt.plot(zsa,usa,'k--')
#plt.imshow(np.rot90(zs),extent=[-1,1,-1,1])

plt.figure()
xs=np.arange(2-4,2+4,0.001)
plt.plot(xs,np.exp(-xs*xs))


plt.show()
