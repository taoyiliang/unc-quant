import numpy as np
from scipy.special.orthogonal import he_roots
from scipy.misc import factorial
from scipy.special import eval_hermitenorm as eHe


def cw(i,n,xi):
  num=factorial(n)*np.sqrt(2.*np.pi)
  den=n*n*(eHe(n-1.,xi))**2
  return num/den




n=128.
os,wsa = he_roots(n)

wsc=np.zeros_like(wsa)
for i in range(len(wsa)):
  wsc[i]=cw(i,n,os[i])
  if abs(wsa[i]-wsc[i])>1e-13:
    print wsa[i]-wsc[i]
