import matplotlib.pyplot as plt
import numpy as np

#read in coarsest
inFiles = ['1.out','3.out','9.out','27.out','81.out','243.out']

inFile = file(inFiles[0],'r')
runs={}
for line in inFile:
  if line.startswith('DIM'):
    h=200./float(line.strip().split(',')[-1])
    crs = h
    runs[h]={}
  elif line.startswith('k'):
    runs[h]['k']=1.+float(line.strip().split(',')[1].split('+')[1])
  else:
    lnar=line.strip().split(' | ')
    spts = lnar[0].split(',')
    sflx = lnar[1].split('  , ')
    for s,pt in enumerate(spts):
      spts[s]=float(pt)
    for s,flx in enumerate(sflx):
      sflx[s]=float(flx)
      runs[h][tuple(spts)]={}
      runs[h][tuple(spts)][1]=sflx[0]
      runs[h][tuple(spts)][2]=sflx[1]
#for item in runs[h].keys():
#  print item,runs[h][item]

#now read in all others
for inFileName in inFiles[1:]:
  inFile = file(inFileName,'r')
  for line in inFile:
    if line.startswith('DIM'):
      h=200./float(line.strip().split(',')[-1])
      runs[h]={}
    elif line.startswith('k'):
      runs[h]['k']=1.+float(line.strip().split(',')[1].split('+')[1])
    else:
      lnar=line.strip().split(' | ')
      spts = lnar[0].split(',')
      for s,pt in enumerate(spts):
        spts[s]=float(pt)
      if tuple(spts) in runs[crs].keys():
        sflx = lnar[1].split('  , ')
        for s,flx in enumerate(sflx):
          sflx[s]=float(flx)
        runs[h][tuple(spts)]={}
        runs[h][tuple(spts)][1]=sflx[0]
        runs[h][tuple(spts)][2]=sflx[1]


#do g1 errs
err1={}
err2={}
errk={}
hs=runs.keys()[:]
hs.sort()
print 'hs: h, N/s'
for h in hs:
  print h,200./h
ref = hs[0]
for h in hs[1:]:
  err1[h]=1e-14
  err2[h]=1e-14
  for pt in runs[h].keys():
    if pt=='k':
      errk[h]=abs(runs[h]['k']-runs[ref]['k'])
    else:
      err1[h]=max(err1[h],abs(runs[h][pt][1]-runs[ref][pt][1])/runs[ref][pt][1])
      err2[h]=max(err2[h],abs(runs[h][pt][2]-runs[ref][pt][2])/runs[ref][pt][2])

p1 = zip(np.array(err1.keys()),np.array(err1.values()))
p2 = zip(np.array(err2.keys()),np.array(err2.values()))
pk = zip(np.array(errk.keys()),np.array(errk.values()))
p1.sort()
p2.sort()
pk.sort()
p1=np.array(p1)
p2=np.array(p2)
pk=np.array(pk)

def addSlopes(x,y):
  for i in range(len(x)-1):
    avx = 0.5*(x[i]+x[i+1])
    avy = 0.5*(y[i]+y[i+1])
    pty = y[i+1]
    slope = (y[i+1]-y[i])/(x[i+1]-x[i])
    plt.text(avx,pty,'%1.2f' %slope,ha='center')
    plt.plot([avx,avx],[avy,pty],'k--')



def doPlot(x,y,label):
  plt.figure()
  x=np.log(x)
  y=np.log(y)
  plt.plot(x,y,'-o',label=label)
  #plt.loglog(x,y,'-o',label=label)
  plt.title(label)
  plt.xlabel('log h')
  plt.ylabel('log err')
  addSlopes(x,y)
  #fit linear
  z=np.polyfit(x,y,1)
  f=np.poly1d(z)
  plt.plot(x,f(x),'k:',label='fit slope: %1.2f' %z[0])
  print '\n',label
  for i in range(len(x)):
    print x[i],y[i]
  plt.legend(loc=0)

doPlot(p1[:,0],p1[:,1],'g1 soln')
doPlot(p2[:,0],p2[:,1],'g2 soln')
doPlot(pk[:,0],pk[:,1],'k')

plt.show()
#plt.plot(p1[0],p1[1],label='g1')
#plt.plot(err2.keys(),err2.values(),label='g2')
#plt.show()
