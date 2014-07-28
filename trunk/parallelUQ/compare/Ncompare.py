import matplotlib.pyplot as plt
import cPickle as pk

def plotCompare(fnames):
  for name in fnames:
    xs,ys,lbls = pk.load(file(name,'r'))
    for i in range(len(xs)):
      plt.loglog(xs[i],ys[i],label=lbls[i])

if __name__=='__main__':
  files=['h1_simple3.pk',
         'h1_simple5.pk',
         'h1_simple10.pk']
  plotCompare(files)
  plt.show()
