import matplotlib.pyplot as plt
import numpy as np

#inputs=['1.cof',
#        '5.cof',
#        '10.cof',
#        '20.cof']
inputs=[]
#rng=[1,2,7,8,9,10,20]
rng=range(1,6)
for i in rng:
  inputs.append(str(i)+'.cof')

for inp in inputs:
  try:
    inFile = file(inp,'r')
  except IOError:continue
  ords=[]
  cofs=[]
  for line in inFile:
    terms=line.split(' | ')
    newCof=abs(float(terms[1]))
    if newCof > 1e-12:
      ords.append(int(terms[0].strip('(').strip(')').strip(',')))
      cofs.append(newCof)
  cofs=np.log10(np.array(cofs))
  plt.plot(ords,cofs,'-',label=inp.split('.')[0])

plt.legend()
plt.xlabel('Expansion Moment')
plt.ylabel('log |Coeff|')
plt.title('1xc2')
plt.show()

