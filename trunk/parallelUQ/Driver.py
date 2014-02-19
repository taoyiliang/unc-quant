from matplotlib.pyplot import show
import time
import datetime as dt
import os
import sys
from GetPot import GetPot
import parExecutor as Exec

class Driver(object):
  '''Since I want a different executor for MC, MLMC, PCESC,
  adding a driver to manage executors makes sense to me.'''
  def __init__(self,argv):
    self.starttime=time.time()
    self.restart=False
    print '\nStarting UQ at',dt.datetime.fromtimestamp(self.starttime)
    self.loadInput(argv)
    self.loadExec()
    self.finishUp()

  def loadInput(self,argv):
    cl = GetPot(argv)
    if cl.search('-r'):
      self.restart=True
      print 'Restart option detected...'
    else:
      print 'No restart option detected...'
    if cl.search('-i'):
      self.unc_inp_file=cl.next('')
      print 'Selected uncertainty input',self.unc_inp_file,'...'
      #check desired file exists
      if not os.path.isfile(self.unc_inp_file):
        raise IOError('Uncertainty input not found: '+self.unc_inp_file)
      print '  ...found uncertainty input file!'
    else:
      raise IOError('Requires and input file using -i.')
    self.input_file = GetPot(Filename=self.unc_inp_file)

  def loadExec(self):
    self.ex_type = self.input_file('Problem/executor','not found')
    print 'Loading executor',self.ex_type
    self.ex=Exec.newExecutor(self.ex_type,self.input_file,self.restart)

  def finishUp(self):
    elapsed=time.time()-self.starttime
    print 'Driver run time:',elapsed,'sec'
    show()
    print '\nDriver complete.\n'


if __name__=='__main__':
  drv = Driver(sys.argv)
