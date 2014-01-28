from GetPot import GetPot
import numpy as np
import sys
import time

def g(x):
  #return x*x - 2.*x + 0.5
  return 1.+2.*x
  #return x*x 

def f(x,y):
  #return x+y
  return x*x + x*y + y*y

cl = GetPot(sys.argv)
if cl.search('-i'):
  inpFileName = ''
  inpFileName = cl.next(inpFileName)
  input_file=GetPot(Filename=inpFileName)
else: raise IOError('Requires an input file using -i.')

x = input_file('Problem/x',1.0)
y = input_file('Problem/y',1.0)
outFileName = input_file('Output/file','test.out')

result = g(x)
#result = f(x,y)

writeFile = file(outFileName,'w')
writeFile.writelines('res,'+str(result)+'\n')
writeFile.close()
