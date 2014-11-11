import matplotlib.pyplot as plt
import cPickle as pk
from pdfMCplot import MCpdf

def ROMpdf(inName):
  inFile = file(inName,'r')
  u=[]
  x,y=pk.load(inFile)
  norm = float(max(y))
  for i,ys in enumerate(y):
    y[i]=float(ys)/norm
  return x,y

def normalize(x):
  norm = float(max(x))
  for i,xs in enumerate(x):
    x[i]=float(xs)/norm
  return x

[MCy,MCx],[dud,dd] = MCpdf('MC_h5_N5.samples')
TDx,TDy = ROMpdf('TD_h5_N5.ROMpdf.pk')
HCx,HCy = ROMpdf('HC_h5_N5.ROMpdf.pk')

MCy = normalize(MCy)
TDy = normalize(TDy)

plt.plot(MCx,MCy,label='MC')
plt.plot(TDx,TDy,label=r'TD, $L=4$')
plt.plot(HCx,HCy,label=r'HC, $L=10$')

plt.title(r'PDF of $k, N=5$')
plt.xlabel(r'Value of $k$')
plt.ylabel('Frequency')
plt.legend()

plt.show()
