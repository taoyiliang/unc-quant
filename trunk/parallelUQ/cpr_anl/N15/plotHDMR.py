import matplotlib.pyplot as plt
import numpy as np

def addPlot(title,lbl,ref=None,r=1):
  if r>1:
    print 'HDMR cannot do r>1, tried to do r =',r
    return
  inFile = file(title,'r')

  entries=[]
  for line in inFile:
    if line=='\n' or line.startswith('Runs') or line.startswith('N,')\
        or line.startswith('#'):
      continue
    entries.append([])
    vals=line.strip().split(',')
    entries[-1].append(int(vals[0]))
    entries[-1].append(int(vals[1]))
    entries[-1].append(float(vals[2]))

  entries=np.array(entries)
  entries=entries[entries[:,0].argsort()]
  errs=np.zeros([len(entries),3])
  if ref==None:
    errs[:,0] = entries[:-1,0]
    #errs[:,1] = abs(entries[:-1,1]-entries[-1,1])/entries[-1,1]
    errs[:,1] = abs(entries[:-1,2]-entries[-1,2])/entries[-1,2]
  else:
    errs[:,0] = entries[:,0]
    #errs[:,1] = abs(entries[:,1]-ref[1])/ref[1]
    errs[:,1] = abs(entries[:,2]-ref[1])/ref[1]

  #for e in errs:
  #  print e
  errs=errs[errs[:,0].argsort()]
#  print '\n\n'
  #for e in errs:
  #  print e
  errs=zip(*errs)
  plt.loglog(errs[0],errs[r],'-o',label=lbl)




if __name__=='__main__':
  title = 'HC_h3_iso.moments'
  r=2
  addPlot(title,'',r=r)
  plt.xlabel('Number of Solves')
  plt.ylabel('Error')
  plt.title('h1, Moment %i' %r)
  plt.show()

