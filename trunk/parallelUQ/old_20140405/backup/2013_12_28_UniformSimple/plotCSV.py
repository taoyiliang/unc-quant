import matplotlib.pyplot as plt
import numpy as np
import matplotlib.mlab as ml

inFile = file('test.out','r')
xs=[]
ys=[]
zs=[]
combo=[]
diffs=[]
for line in inFile:
  if line[0]=='R':
    xName,yName=line.split(',')[1:3]
    continue
  elif line[0]=='a':
    continue
  lar=line.split(',')
  tRun=int(lar[0])
  xs.append(float(lar[1]))
  ys.append(float(lar[2]))
  zs.append(float(lar[3]))
  combo.append((xs[-1],ys[-1],zs[-1]))
  diffs

#srt=zip(zs,xs,ys)
#srt.sort()
#for s in srt:
#  print s

xs=np.array(xs)
ys=np.array(ys)
zs=np.array(zs)

xi=np.linspace(np.min(xs),np.max(xs),100)
yi=np.linspace(np.min(ys),np.max(ys),100)
zi=ml.griddata(xs,ys,zs,xi,yi)

plt.contour(xi,yi,zi,15,linewidth=0.5,colors='k')
plt.pcolormesh(xi,yi,zi,cmap=plt.get_cmap('rainbow'))

plt.colorbar()
plt.scatter(xs,ys,marker='o',c='b',s=1,zorder=10)
plt.xlim(np.min(xs),np.max(xs))
plt.ylim(np.min(ys),np.max(ys))

#print 'lims',np.min(xs),np.max(xs),np.min(ys),np.max(ys),np.min(zs),np.max(zs)

plt.xlabel(yName)
plt.ylabel(xName)
plt.title('k_eff ('+str(tRun)+' runs)')
plt.show()
