import matplotlib.pyplot as plt

title='Uniform Sampling by Order'
inNames = ['mcu']#,'sc2u','sc4u','sc8u','sc16u']
#inNames = ['mcn','sc2n','sc4n','sc8n','sc16n']
labels = ['MonteCarlo','SC2','SC4','SC8','SC16']

plt.figure()

for n,name in enumerate(inNames):
  inFile = file('../simple/'+name+'.out','r')
  ctrs=[]
  vals=[]
  for line in inFile:
    if line.startswith('Centers') or line.strip()=='':
      continue
    elif line.startswith('Average'):
      avg = float(line.strip().split(',')[1])
    elif line.startswith('2nd Moment'):
      mom2= float(line.strip().split(',')[1])
    elif line.startswith('Variance'):
      var = float(line.strip().split(',')[1])
    else:
      ctr,val=line.strip().split(',')
      ctrs.append(float(ctr))
      vals.append(float(val))
  plt.plot(ctrs,vals,label=labels[n])
plt.legend()
plt.title(title)
plt.xlabel('Solution Value')
plt.ylabel('Solution Frequency')
plt.show()
