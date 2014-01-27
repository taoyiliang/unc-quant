import cPickle as pk
import matplotlib.pyplot as plt


numlist = [2,4,8,16,64]

plt.figure()
for num in numlist:
  c,b = pk.load(open(str(num)+'.pk','rb'))
  plt.plot(c,b,label=str(num))

plt.title('Monte Carlo Sample of ROM by Expansion Order')
plt.ylabel('Probability')
plt.xlabel('Solution Value')
plt.axis([0.95,1.05,0,0.15])
plt.legend()
plt.show()
