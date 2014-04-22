from matplotlib.pyplot import show
import time
import datetime as dt
import os
import sys
from GetPot import GetPot
import Executor as Exec
import multiprocessing
import numpy as np

class Driver(object):
  '''Since I want a different executor for MC, MLMC, PCESC,
  adding a driver to manage executors makes sense to me.'''
  def __init__(self,argv):
    self.starttime=time.time()
    print '\nStarting UQ at',dt.datetime.fromtimestamp(self.starttime)
    self.loadInput(argv)
    self.loadExec()
    self.finishUp()

  def loadInput(self,argv):
    cl = GetPot(argv)
    print argv
    if cl.search('-i'):
      self.unc_inp_file=cl.next('')
      print 'Selected uncertainty input',self.unc_inp_file,'...'
      #check desired file exists
      if not os.path.isfile(self.unc_inp_file):
        raise IOError('Uncertainty input not found: '+self.unc_inp_file)
      print '  ...found uncertainty input file.'
    else:
      raise IOError('Requires and input file using -i.')
    self.input_file = GetPot(Filename=self.unc_inp_file)

  def loadExec(self):
    self.ex_type = self.input_file('Problem/executor','not found')
    print 'Loading executor',self.ex_type
    self.ex=Exec.ExecutorFactory(self.ex_type,self.input_file)

  def finishUp(self):
    elapsed=time.time()-self.starttime
    print 'Driver run time:',elapsed,'sec'
    print '\nStarting postprocessing...'
    makePDF = self.input_file('Backends/makePDF',0)
    if makePDF:
      print '...sampling ROM...'
      numSamples = self.input_file('Backends/PDFsamples',-1)
      if numSamples==-1:
        print '...Backends/PDFsamples not found; using 1e4...'
        numSamples = int(1e4)
      self.makePDF(numSamples)

    #TODO DEBUG
    samp = self.ex.ROM([1.,1.])
    print 'sparseU(1,1) =',samp

    show()
    print '\nDriver complete.\n'

  def makePDF(self,M):
    wantprocs = self.input_file('Problem/numprocs',1)
    numprocs = min(wantprocs,multiprocessing.cpu_count())
    nBins = self.input_file('Backends/PDF/bins',10)
    binmin = self.input_file('Backends/PDF/min',-10)
    binmax = self.input_file('Backends/PDF/max',10)
    bins=np.linspace(binmin,binmax,nBins+1)
    self.ex.makePDF(numprocs,M,bins)


if __name__=='__main__':
  print sys.argv,type(sys.argv)
  drv = Driver(sys.argv)
