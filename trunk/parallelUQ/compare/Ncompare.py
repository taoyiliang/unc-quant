import matplotlib.pyplot as plt
import cPickle as pk

def plotCompare(fnames):
  clrs=['b','g','r']
  styl=[':','-x','-o']
  size=[1,5,5]
  for n,name in enumerate(fnames):
    xs,ys,lbls = pk.load(file(name,'r'))
    for i in range(len(xs)):
      plt.loglog(xs[i],ys[i],clrs[n]+styl[i],
                 ms=size[i],mec=clrs[n],mfc=clrs[n],
                 label=lbls[i])
  plt.legend(loc=3)

if __name__=='__main__':
  files=['h1_simple3.pk',
         'h1_simple5.pk',
         'h1_simple10.pk']
  plotCompare(files)
  plt.show()
