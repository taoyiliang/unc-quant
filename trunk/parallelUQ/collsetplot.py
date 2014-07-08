import SparseQuads as sq
import Variables as v
import IndexSets as isets
import matplotlib.pyplot as plt
import numpy as np
from itertools import product as allcombos

def quadrule(x):
  return x


varlist={}
varlist[1]=v.VariableFactory('uniform')
varlist[1].setDist([0,1])
varlist[2]=v.VariableFactory('uniform')
varlist[2].setDist([0,1])

fx=[]
fy=[]
for i in range(5):
  pts,wts = sq.tensorGrid(2,[i+1,i+1],varlist,[i,i])
  for pt in pts:
    fx.append(pt[0])
    fy.append(pt[1])

plt.subplot(1,3,1)
plt.axis([-1,1,-1,1])
plt.plot(fx,fy,'kx',markersize=8)
plt.title('Tensor Product (%i points)' %53)

plt.subplot(1,3,2)

iset = isets.IndexSetFactory(2,4,'TD')
SG = sq.BasicSparse(2,4,iset,quadrule,varlist)
SG=np.array(SG)
pts=SG[:,0]
for p,pt in enumerate(pts):
  pts[p]=np.array(list(pt))
pts=np.array(pts)
xs=[]
ys=[]
for pt in pts:
  xs.append(pt[0])
  ys.append(pt[1])

plt.plot(xs,ys,'kx',mfc='None',markersize=8)
plt.title('Isotropic TD (%i points)' %51)
plt.axis([-1,1,-1,1])


plt.subplot(1,3,3)

iset = isets.IndexSetFactory(2,4,'HC')
SG = sq.BasicSparse(2,4,iset,quadrule,varlist)
SG=np.array(SG)
pts=SG[:,0]
for p,pt in enumerate(pts):
  pts[p]=np.array(list(pt))
pts=np.array(pts)
xs=[]
ys=[]
for pt in pts:
  xs.append(pt[0])
  ys.append(pt[1])

plt.plot(xs,ys,'kx',mfc='None',markersize=8)
plt.title('Isotropic HC (%i points)' %17)
plt.axis([-1,1,-1,1])

plt.show()


