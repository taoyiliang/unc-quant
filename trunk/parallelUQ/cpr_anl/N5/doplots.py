from matplotlib.pyplot import show
from collections import OrderedDict
import os,sys

cursys=sys.path
sys.path.append('/'.join(os.getcwd().split('/')[:-2])+'/plotscripts')
from slnerrplot import slnerrplot

cases=OrderedDict({})
for line in file('list.inp','r'):
  fname=line.strip()
  cases[fname]=file(fname,'r')

reffile = file('ref.inp','r')
for line in reffile:
  if line.startswith('N'):N=int(line.strip().split('=')[1])
  elif line.startswith('title'):title=line.strip().split('=')[1]
  elif line.startswith('xlim'):exec(line.strip())
  elif line.startswith('s1ylim'):exec(line.strip())
  elif line.startswith('e1ylim'):exec(line.strip())
  elif line.startswith('s2ylim'):exec(line.strip())
  elif line.startswith('e2ylim'):exec(line.strip())
  elif line.startswith('ref'):
    ref = line.strip().split('=')[1].split(',')
    ref[0]=int(ref[0]) #dud
    ref[1]=float(ref[1]) #mean
    ref[2]=float(ref[2]) #variance

print 'starting first moment...'
slnerrplot(cases,title,N,xlim,s1ylim,e1ylim,ref,1)
print 'starting second moment...'
slnerrplot(cases,title,N,xlim,s2ylim,e2ylim,ref,2)

for cfile in cases.values():
    cfile.close()

show()

sys.path=cursys
