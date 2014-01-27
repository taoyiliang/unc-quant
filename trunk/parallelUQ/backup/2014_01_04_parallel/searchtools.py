
#binomial search, I think this is
def findClosest_LowHi(vec,val):
  '''
  Given a sorted array vec, will find the index of the
  nearest value to val within vec.

  Inputs: array vec, float val
  Outputs: int index
  '''
  maxI=len(vec)-1
  loc=maxI
  lastloc=0
  curMax=maxI
  curMin=0
  lastDir='up'
  found = False
  track=0
  while not found:
    track+=1
    print 'At',loc,':',vec[loc]
    if vec[loc]==val:
      return loc
    elif vec[loc]<val: #search up
      if lastDir=='down':
        curMin=loc
      if loc==maxI:
        return loc
      adj = int((curMax-loc)*0.5)
      if adj==0:#or adj==1:
        if vec[loc+1]-val < val-vec[loc]:
          return loc+1
        else:
          return loc
      lastdir='up'
    else: #search down
      if lastDir=='up':
        curMax=loc
      if loc==0:
        return loc
      adj = -int(0.5*(curMin+loc))
      if adj==0:# or adj==-1:
        if val-vec[loc-1] < vec[loc]-val:
          return loc-1
        else:
          return loc
    lastloc=loc
    loc=loc+adj
