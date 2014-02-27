import numpy as np
import matplotlib.pyplot as plt

d=dict()
xs=[]
ys=[]

inFile = file('flux.out','r')
for line in inFile:
  if line[:3]=='DIM':
    dud,nX,nY=line.strip().split(',')
    nX=int(nX)
    nY=int(nY)
    print 'Dims:',nX,nY
    continue
  if line[:1]=='k':
    ent=line[:-1].split(',')
    k=float(ent[1])
    print ent[0],'=',ent[1]
    continue
  temp=line.split('|')[0]
  newX=float(temp.split(',')[0])
  newY=float(temp.split(',')[1])
  #newY=float(temp.split('|')[0])
  if newX not in xs:
    xs.append(newX)
  if newY not in ys:
    ys.append(newY)
  ds = line.split('|')[1]
  d[(newX,newY)]=[]
  d[(newX,newY)].append(float(ds.split(',')[0].strip()))
  d[(newX,newY)].append(float(ds.split(',')[1].strip()))

#for key in d.keys():
#  print key,d[key]

xs.sort()
ys.sort()
#print xs
#print ys
s1=np.zeros([len(xs),len(ys)])
s2=np.zeros([len(xs),len(ys)])
for i,x in enumerate(xs):
  for j,y in enumerate(ys):
    s1[i,j]=d[(x,y)][0]
    s2[i,j]=d[(x,y)][1]
#s1=np.rot90(s1)
#s2=np.rot90(s2)

xs=np.array(xs)
ys=np.array(ys)

twxs=np.zeros(len(xs)-1)
for i,x in enumerate(twxs):
  twxs[i]=0.5*(xs[i]+xs[i+1])

twys=np.zeros(len(ys)-1)
for j,y in enumerate(twys):
  twys[j]=0.5*(ys[j]+ys[j+1])

X,Y=np.meshgrid(xs,ys)

from matplot import giveMap,addRegions

plt.figure()
#CS1 = plt.contourf(X,Y,s1,100)
CS1 = plt.imshow(np.rot90(s1),extent=[0,nX,0,nY],interpolation='none')
plt.colorbar(CS1)

if nX<100 and nY<100 and False:
  for x in twxs:
    tempx = np.ones(len(ys))*x
    plt.plot(tempx,ys,'k-')
  for y in twys:
    tempy = np.ones(len(xs))*y
    plt.plot(xs,tempy,'k-')

addRegions(nX/4)
plt.title('G1 flux, '+str(nX)+'x'+str(nY)+', k = '+str(k))

plt.figure()
#CS2 = plt.contourf(X,Y,s2,100)
CS2 = plt.imshow(np.rot90(s2),extent=[0,nX,0,nY],interpolation='none')
plt.colorbar(CS2)

addRegions(nX/4)

plt.title('G2 flux, '+str(nX)+'x'+str(nY)+', k = '+str(k))
#plt.axis([0,165,0,165])
#giveMap()
plt.show()
