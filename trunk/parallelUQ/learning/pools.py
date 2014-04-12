from multiprocessing import Pool
import numpy as np
import time

start = time.time()

class dummy(object):
  def __init__(self):
    self.counter=-1

  def __iter__(self):
    return self

  def next(self):
    self.counter+=1
    if self.counter>100:
      raise StopIteration
    else:
      return self.makeres()

  def makeres(self):
      res = {}
      res['val']=self.counter
      res['squ']=self.counter**2
      return res

class doPool(object):
  def __init__(self):
    pool = Pool(processes=4)
    arg=range(int(1e7)+1)
    ins = zip(arg,arg)
    one = time.time()
    print 'one',one-start
    fc = dummy()
    res = pool.imap(self.fun,fc,int(1e3))
    pool.close()
    pool.join()
    two = time.time()
    print 'two',two-one
    for r in res:
      print r
    three = time.time()
    print three-two,'secs'

  def fun(self,x):
    toprint = 'val: %i, squ: %i' %(x['val'],x['squ'])
    return toprint

test = doPool()

