import matplotlib.pyplot as plt
import numpy as np

inputs=['1xc2.cof',
        '1xf2.cof',
        '4xc2.cof',
        '4xf2.cof',
        '5D2.cof']

for inp in inputs:
  inFile = file(inp,'r')
  ords=[]
  cofs=[]
  for line in inFile:
    terms=line.split(' | ')
    newCof=abs(float(terms[1]))
    if newCof > 1e-12:
      ords.append(int(terms[0].strip('(').strip(')').strip(',')))
      cofs.append(newCof)
  cofs=np.log10(np.array(cofs))
  plt.plot(ords,cofs,label=inp.split('.')[0])

plt.legend()
plt.xlabel('Expansion Moment')
plt.ylabel('log Coefficient')
plt.show()

