import numpy as np
import triple as solver
import classVar as var
import classQuadrature as q
from writeTransIn import writeInput
import classUQ as uq
import os

# set non-random variables, set None for randoms
x=None
y=None
z=None
p=[x,y,z]
# designate random variables
randVars=[]
randVars.append(var.uniVar('x',[-1,1],paramIndex=0,N=2))
randVars.append(var.uniVar('y',[-1,1],paramIndex=1,N=2))
randVars.append(var.uniVar('z',[-1,1],paramIndex=2,N=2))

scsolve=uq.SC(p,randVars,solver,stype=list)
scsolve.propUQ(3)
UQsoln=scsolve.UQsoln
#visualize solution
import matplotlib.pyplot as plt
xs=np.arange(-1,1.5,0.5)
ys=xs.copy()
zs=xs.copy()
us=np.zeros([len(xs),len(ys),len(zs)])
np.set_printoptions(suppress=True)
for i,x in enumerate(xs):
  for j,y in enumerate(ys):
    for k,z in enumerate(zs):
      us[i,j,k]=UQsoln([x,y,z])
      print np.array([x,y,z,us[i,j,k]])
#plt.imshow(np.rot90(zs),extent=[-1,1,-1,1])
print us[3,3,3]
#plt.show()
