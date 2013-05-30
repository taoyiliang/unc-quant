import classVar as v
import scipy.stats as spst
import scipy.special as sps
import numpy as np

test=v.normVar('x',[2,4])
test.setZeroToOneSample(spst.uniform(0,1).ppf,sps.sh_legendre,24)

for i in range(3,24,3):
  z=test.z21quad.ords[i]
  x=test.sampleZeroToOne(z)
  Ua=np.exp(-x*x)
  print x,z,Ua
