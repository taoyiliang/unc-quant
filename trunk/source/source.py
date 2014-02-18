from GetPot import GetPot
import numpy as np
import sys
import time

def g(xa=0.75,loc=2.0,D=0.5,S=1.0):
  #S = 1.0 #source strength
  #D = 0.5 #diffusion coefficient
  #loc = 2.0 #depth in cm
  cof = S/xa
  #print '...in source.py...'
  #print 'D,xa:',D,xa
  L = np.sqrt(D/xa)
  flux = cof*(1.0-np.exp(-loc/L))
  return flux

if __name__=='__main__':
  cl = GetPot(sys.argv)
  if cl.search('-i'):
    inpFileName = ''
    inpFileName = cl.next(inpFileName)
    input_file=GetPot(Filename=inpFileName)
  else: raise IOError('Requires an input file using -i.')

  xa = input_file('Problem/xa' ,1.0)
  D  = input_file('Problem/D'  ,1.0)
  loc= input_file('Problem/loc',1.0)
  S  = input_file('Problem/S'  ,1.0)

  outFileName = input_file('Output/file','test.out')
  print 'outFileName:',outFileName

  result = g(xa,D,loc,S)

  writeFile = file(outFileName,'w')
  writeFile.writelines('res,'+str(result)+'\n')
  writeFile.close()
