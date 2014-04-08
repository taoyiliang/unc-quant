from multiprocessing import Pool
import numpy as np
import time

start = time.time()

def f(arg):
  x=arg[0]
  y=arg[1]
  return y,x*x

pool = Pool(processes=4)
arg=range(int(1e7)+1)
ins = zip(arg,arg)
one = time.time()
print 'one',one-start
res = pool.imap(f,ins,int(1e3))
two = time.time()
print 'two',two-one
for r in res:
  pass
#print res[-1]
three = time.time()
print three-two,'secs'
