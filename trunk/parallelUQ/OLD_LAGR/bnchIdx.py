import newIndexSets as nis
import IndexSets as ois
import time

L=4
one = range(L)
lists=[]
for i in range(10):
  lists.append(one)

startold = time.time()
oldset = ois.IndexSetFactory(len(lists),L,'HC')
timeold = time.time()-startold
print 'len old:',len(oldset)
print 'old time:',timeold

startnew = time.time()
newset = nis.IndexSetFactory(len(lists),L,'HC')
timenew = time.time()-startnew
print 'len new:',len(newset)
print 'new time:',timenew

print 'check:',oldset==newset
