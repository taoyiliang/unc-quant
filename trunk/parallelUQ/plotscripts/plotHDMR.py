import matplotlib.pyplot as plt
import numpy as np

def addPlot(title,inFile,lbl,ref=None,r=1,slnfig=None):
  inFile.seek(0)

  if lbl.split('_')[1]=='TD': marker='-o'
  else: marker = '-^'

  level = int(lbl.split('_')[2][1:])
  if level==1: color='r'
  elif level==2: color='m'
  elif level==3: color='y'
  else: color='k'

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
    entr=zip(*entr)
    plt.figure(slnfig.number)
    plt.plot(entr[0],entr[r+1],color+marker,label=lbl)
    plt.figure(errfig.number)

  errs=np.zeros([len(entries),3])
  if ref==None:
    errs[:,0] = entries[:-1,0]
    errs[:,1] = abs(entries[:-1,r+1]-entries[-1,r+1])/entries[-1,r+1]
  else:
    errs[:,0] = entries[:,0]
    errs[:,1] = abs(entries[:,r+1]-ref[r])/ref[r]

  errs=errs[errs[:,0].argsort()]
  errs=zip(*errs)

  plt.loglog(errs[0],errs[1],color+marker,label=lbl)




if __name__=='__main__':
  title = 'HC_h3_iso.moments'
  r=2
  addPlot(title,'',r=r)
  plt.xlabel('Number of Solves')
  plt.ylabel('Error')
  plt.title('h1, Moment %i' %r)
  plt.show()

