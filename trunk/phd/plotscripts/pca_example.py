import numpy as np
import matplotlib.pyplot as plt
import cPickle as pk
import os

# linear system:
# x1 = a + 2b
# x2 = 3a - b
# x3 = 2a
# obviously one variable is redundant

def calc_x1(a,b):
  return 1.*a+2.*b

def calc_x2(a,b):
  return 3.*a-1.*b

def calc_x3(a,b):
  return 2.*a

#construct covariance matrix
if not os.path.isfile('./covar.pk'):
  x1=[]
  x2=[]
  x3=[]
  for i in range(100000):
    x = np.random.rand()
    y = np.random.rand()
    x1.append(calc_x1(x,y))
    x2.append(calc_x2(x,y))
    x3.append(calc_x3(x,y))
  covar = np.cov(np.array([x1,x2,x3]))
  pk.dump((covar,x1,x2,x3),file('covar.pk','w'))
else:
  covar,x1,x2,x3 = pk.load(file('covar.pk','r'))
print 'covar:',covar

#plot relationships
plt.figure()
plt.plot(x1[:1000],x2[:1000],'bo',label='x1,x2')
plt.legend()

plt.figure()
plt.plot(x1[:1000],x3[:1000],'bo',label='x1,x3')
plt.legend()

plt.figure()
plt.plot(x2[:1000],x3[:1000],'bo',label='x2,x3')
plt.legend()

#perform PCA
U,S,V = np.linalg.svd(covar,full_matrices=True)
# Q = U*sqrt(S)
ss = np.sqrt(S)
Q = U.copy()
Q[0:] *= ss
Q[1:] *= ss
Q[2:] *= ss
print ''
print 'S',S
print ''
print 'Q',Q
print ''

#sanity test
ys=[0,0,0]
print ys,':',np.matmul(Q,ys)
ys=[1,1,1]
print ys,':',np.matmul(Q,ys)
ys=[1,1,0]
print ys,':',np.matmul(Q,ys)

plt.show()
