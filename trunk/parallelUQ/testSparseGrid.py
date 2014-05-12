import Variables as v
import SparseQuads as sq
import IndexSets as isets

N=3
L=3
isetname='TD'

def qrule(x):
  return x

varlist={}
varlist['x']=v.VariableFactory('uniform',name='x',path='/x')
varlist['y']=v.VariableFactory('uniform',name='y',path='/y')
varlist['z']=v.VariableFactory('uniform',name='z',path='/z')

iset=isets.IndexSetFactory(N,L,isetname)

for var in varlist.values():
  var.setDist([0,1])



SG = sq.BasicSparse(N,L,iset,qrule,varlist)
ptwt = sq.removeDuplicates(SG)
for pt in ptwt.keys():
  print pt,ptwt[pt]
print 'TOTAL'
print sum(ptwt.values())
