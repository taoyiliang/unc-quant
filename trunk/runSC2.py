import numpy as np
import double as solver
import classVar as var
import classQuadrature as q
from writeTransIn import writeInput
import classUQ as uq
import os

# set non-random variables, set None for randoms
x=None
y=None
p=[x,y]
# designate random variables
randVars=[]
randVars.append(var.uniVar('x',[0,7],paramIndex=0,N=2))
randVars.append(var.uniVar('y',[0,10],paramIndex=1,N=2))

scsolve=uq.SC(p,randVars,solver)
scsolve.propUQ(3)
UQsoln=scsolve.UQsoln
print scsolve.coeffs
#visualize solution
import matplotlib.pyplot as plt
xs=np.arange(0,1.2,0.2)
ys=xs.copy()
us=np.zeros([len(xs),len(ys)])
np.set_printoptions(suppress=True)
for i,x in enumerate(xs):
  for j,y in enumerate(ys):
    us[i,j]=UQsoln([x,y])
    print np.array([x,y,us[i,j]])
plt.imshow(np.rot90(us),extent=[-1,1,-1,1])
plt.show()
