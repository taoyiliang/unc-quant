import matplotlib.pyplot as plt
import numpy as np

def addPlot(title,lbl,ref=None):
  inFile = file(title,'r')

  entries=[]
  for line in inFile:
    if line=='\n' or line.startswith('Moments') or line.startswith('N,'):
      continue
    entries.append([])
    vals=line.strip().split(',')
    entries[-1].append(int(vals[0]))
    entries[-1].append(float(vals[1]))
    entries[-1].append(float(vals[2]))

  entries=np.array(entries)
  entries=entries[entries[:,0].argsort()]
  if ref==None:
    errs=np.zeros([len(entries)-1,3])
    errs[:,0] = entries[:-1,0]
    errs[:,1] = abs(entries[:-1,1]-entries[-1,1])/entries[-1,1]
    errs[:,2] = abs(entries[:-1,2]-entries[-1,2])/entries[-1,2]
  else:
    errs=np.zeros([len(entries),3])
    errs[:,0] = entries[:,0]
    errs[:,1] = abs(entries[:,1]-ref[1])/ref[1]
    errs[:,2] = abs(entries[:,2]-ref[2])/ref[2]

  #for e in errs:
  #  print e
  errs=errs[errs[:,0].argsort()]
#  print '\n\n'
#  for e in errs:
#    print e
  errs=zip(*errs)
  print errs[1]
  plt.loglog(errs[0],errs[1],'-o',label=lbl)




if __name__=='__main__':
  title = 'HC_h5_aniso.moments'
  fitx=np.linspace(1,10000,10)
  scale=1e-3
  fit2=scale*fitx**(-2.)
  fit3=scale*fitx**(-3.)
  fit4=scale*fitx**(-4.)
  fith=scale*fitx**(-0.5)
  plt.loglog(fitx,fit2,':',label='O(2)')
  plt.loglog(fitx,fit3,':',label='O(3)')
  plt.loglog(fitx,fit4,':',label='O(4)')
  plt.loglog(fitx,fith,':',label='O(1/2)')
  plt.legend(loc=3)
  addPlot(title,'')
  plt.xlabel('Number of Solves')
  plt.ylabel('Error')
  plt.title('h1, Mean')
  plt.axis([0,10000,1e-9,1e-4])
  plt.show()

