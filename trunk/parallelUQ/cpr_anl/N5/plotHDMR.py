import matplotlib.pyplot as plt
import numpy as np

def addPlot(title,lbl,ref=None,r=1,slnfig=None):
  #if r>1:
  #  print 'HDMR cannot do r>1, tried to do r =',r
  #  return
  inFile = file(title,'r')

  entries=[]
  for line in inFile:
    if line=='\n' or line.startswith('Runs') or line.startswith('N,')\
        or line.startswith('#'):
      continue
    entries.append([])
    vals=line.strip().split(',')
    entries[-1].append(int(vals[0])) #comp solves
    entries[-1].append(int(vals[1])) #level
    entries[-1].append(float(vals[2])) #mean
    entries[-1].append(float(vals[3])) #variance

  entries=np.array(entries)
  entries=entries[entries[:,0].argsort()]

  if slnfig:
    errfig = plt.gcf()
    entr=entries.copy()
    #entr[:,1]=(1.-entr[:,1])/entr[:,1]
    #entr[:,2]=(1.-entr[:,2])/entr[:,2]
    entr=zip(*entr)
    plt.figure(slnfig.number)
    plt.plot(entr[0],entr[r+1],'-x',label=lbl)
    plt.figure(errfig.number)

  errs=np.zeros([len(entries),3])
  if ref==None:
    errs[:,0] = entries[:-1,0]
    errs[:,1] = abs(entries[:-1,r+1]-entries[-1,r+1])/entries[-1,r+1]
  else:
    errs[:,0] = entries[:,0]
    #errs[:,1] = abs(entries[:,1]-ref[1])/ref[1]
    errs[:,1] = abs(entries[:,r+1]-ref[r])/ref[r]

  #for e in errs:
  #  print e
  errs=errs[errs[:,0].argsort()]
#  print '\n\n'
  #for e in errs:
  #  print e
  errs=zip(*errs)
  plt.loglog(errs[0],errs[1],'-o',label=lbl)




if __name__=='__main__':
  title = 'HC_h3_iso.moments'
  r=2
  addPlot(title,'',r=r)
  plt.xlabel('Number of Solves')
  plt.ylabel('Error')
  plt.title('h1, Moment %i' %r)
  plt.show()

