import numpy as np

class AdaptiveSparse(object):
  def __init__(self,N):
    self.N=N
    self.O={}
    self.A={}
    self.r=0

  def AdaptiveSparseQuadStep(A,O,N,i):
   #find max gk index
   maxErr = max(self.A.values())
   maxErrIx = self.A.values().index(maxErr))
   maxIx = self.A.keys()[maxErrIx]
   self.O[maxIx]=maxErr
   del self.A.keys()[maxIx]
   self.r-=maxErr
   for range(N):
    j=
