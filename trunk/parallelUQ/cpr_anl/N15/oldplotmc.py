import matplotlib.pyplot as plt
import numpy as np

def getEntries(title):
  inFile = file(title,'r')
  entries=[]
  for line in inFile:
    if line=='\n' or line.startswith('N,'):
      continue
    entries.append([])
    vals=line.strip().split(',')
    entries[-1].append(int(vals[0]))
    entries[-1].append(float(vals[1])/float(entries[-1][0]))
    entries[-1].append(float(vals[2])/float(entries[-1][0])-entries[-1][1]**2)
  entries=np.array(entries)
  entries=entries[entries[:,0].argsort()]
  return entries

def getSparseEntries(title,vec):
  inFile = file(title,'r')
  entries=[]
  i=0
  for l,line in enumerate(inFile):
    if line=='\n' or line.startswith('N,'):
      continue
    #if l>100:break
    if l//10**i in vec and l%10**i==0:
      #print 'DEBUG',l,10**i,l//10**i
      entries.append([])
      vals=line.strip().split(',')
      entries[-1].append(int(vals[0]))
      entries[-1].append(float(vals[1])/float(entries[-1][0]))
      entries[-1].append(float(vals[2])/float(entries[-1][0])-entries[-1][1]**2)
    if l//10**i==vec[-1]:
      i+=1
      print 'FINISHED',i-1
  entries=np.array(entries)
  entries=entries[entries[:,0].argsort()]
  return entries

def linearize(xs,ys):
  m,b=np.polyfit(np.log(xs),np.log(ys),1)
  #nx = np.linspace(2,max(xs))
  ny = np.exp(m*np.log(xs)+b)
  plt.loglog(xs,ny,'k:',label='MC slope %1.2e' %m)


def addPlot(title,lbl,ref=None,r=1):
  entries=getEntries(title)
  #entries=getSparseEntries(title,range(10))
  if ref==None:
    errs=np.zeros([len(entries)-1,3])
    errs[:,0] = entries[:-1,0]
    errs[:,1] = abs(entries[:-1,1]-entries[-1,1])/entries[-1,1]
    errs[:,2] = abs(entries[:-1,2]-entries[-1,2])/entries[-1,2]
  else:
    errs=np.zeros([len(entries),3])
    errs[:,0] = entries[:,0]
    errs[:,1] = abs(entries[:,1]-ref[1])/ref[1]

  errs=errs[errs[:,0].argsort()]
  errs=zip(*errs)
  plt.loglog(errs[0],errs[r],'-',label=lbl,alpha=0.1)
  linearize(errs[0][:-1],errs[r][:-1])


if __name__=='__main__':
  title = 'MC_h5_iso.samples'
  r=2
  addPlot(title,'',r=r)
  plt.xlabel('Number of Solves')
  plt.ylabel('Error')
  plt.title('h1, Moment %i' %r)
  plt.show()

