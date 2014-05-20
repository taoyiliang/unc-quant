import scipy.stats as dists
import numpy as np
import matplotlib.pyplot as plt

mu=1.
var=1.

def getCofs(u,var,alpha):
  s=np.sqrt(var)
  term=np.sqrt(2.*alpha+1.)
  c=-s*term + u
  d=2.*s*term
  return c,d

xs=np.linspace(-5,5,1000)

normal=dists.norm(mu,var)
yN=[]
for x in xs:
  yN.append(normal.pdf(x))
plt.plot(xs,yN,'k:')

def pickA(a,b,mu,var):
  c,d=getCofs(mu,var,a)
  print 'a,L,R:',a,c,c+2.*d
  beta=dists.beta(a,b,loc=c,scale=d)
  yB=[]
  for x in xs:
    yB.append(beta.pdf(x))
  plt.plot(xs,yB,label='a='+str(a))

def pickB(c,d,var):
  a=d*d/8./var-0.5
  beta=dists.beta(a,a,loc=c,scale=d)
  yB=[]
  for x in xs:
    yB.append(beta.pdf(x))
  plt.plot(xs,yB,label='c='+str(c))


#hold alpha,beta and vary c,d
for i in range(1,100,20):
  pickA(i,i,mu,var)
plt.legend()
plt.title('Vary L,R')

plt.figure()
plt.plot(xs,yN,'k:')
#hold c,d and vary alpha
cs=np.array([-2,-2.5,-3,-3.5])
cs=cs+mu
for c in cs:
  pickB(c,2.*abs(mu-c),var)
plt.legend()
plt.title('Vary Alpha=Beta')
plt.show()


