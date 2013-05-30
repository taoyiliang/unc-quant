import numpy as np
import matplotlib.pyplot as p

# set non-random variables, set None for randoms

rt2=np.sqrt(2)

xs=np.arange(-1,1,0.01)
ys=xs.copy()
zs=np.zeros([len(xs),len(ys)])
zc=zs.copy()
for i,x in enumerate(xs):
  for j,y in enumerate(ys):
    zs[i,j]=(1.-x)*(np.exp(-y*y))
    zc[i,j]=1./rt2-x/rt2+1./(82.*rt2)*(8.*y*y-2.)-\
        x/(82.*rt2)*(8.*y*y-2.)




p.imshow(np.rot90(zs-zc),extent=[-1,1,-1,1])
print np.max(zs-zc)
p.show()
