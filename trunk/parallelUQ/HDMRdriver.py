import os
import sys
import time
import datetime as dt
import multiprocessing
from itertools import combinations as combos

import numpy as np
from GetPot import GetPot
from matplotlib.pyplot import show

import Executor
import Variables
import InputEditor

class Driver(object):
  def __init__(self,argv):
    self.starttime=time.time()
    print '\nStarting HDMR UQ at',dt.datetime.fromtimestamp(self.starttime)
    self.loadInput(argv)
    self.loadVars()
    self.setRuns()
    self.createROMs()
    #self.finishUp()

  def loadInput(self,argv):
    cl = GetPot(argv)
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
    self.hdmr_level = self.input_file('HDMR/level',0)
    if self.hdmr_level==0:
      print 'HDMR level not specified.  Using 2...'
      self.hdmr_level=2

  def loadVars(self):
    print '\nLoading uncertain variables...'
    uVars = self.input_file('Variables/names','').split(' ')
    self.varDict={}
    for var in uVars:
      path = self.input_file('Variables/'+var+'/path',' ')
      dist = self.input_file('Variables/'+var+'/dist',' ')
      args = self.input_file('Variables/'+var+'/args',' ').split(' ')
      for a,arg in enumerate(args):
        args[a]=float(arg)
      impwt = self.input_file('Variables/'+var+'/weight',1.0)
      self.varDict[var]=Variables.VariableFactory(dist,var,path,impwt)
      self.varDict[var].setDist(args)

  def setRuns(self):
    #set up the runs that we want to do
    #singles
    self.todo = {}
    varnames = self.varDict.keys()
    for i in range(1,self.hdmr_level+1):
      new = self.addHDMRLevel(varnames,i)
      for key,value in new.iteritems():
        self.todo[key]=value
    #for i in range(1,5):
    #  print '\nSet',i
    #  print list(x for x in self.getRunSet(i).keys())

  def getRunSet(self,num):
    ret={}
    for key,value in self.todo.iteritems():
      if len(key)==num:
        ret[key]=value
    return ret

  def addHDMRLevel(self,varnames,lvl):
    todo={}
    nameset = combos(varnames,lvl)
    for entry in nameset:
      todo[entry]=[]
      for e,ent in enumerate(entry):
        todo[entry].append(self.varDict[ent])
    return todo

  def createROMs(self):
    self.ROMs={}
    ie = InputEditor.HDMR_IO()
    #do reference problem TODO this assumes they left it mean in the first place
    chlist,ident = self.makeCase({})
    runfile = ie.writeInput(self.unc_inp_file,chlist,ident)
    inp_file = GetPot(Filename=runfile)
    ex = Executor.ExecutorFactory('SC',{},inp_file)
    ex.run()
    #END reference case
    self.ROMs[ident]=ex.ROM
    numruns={}
    for i in range(1,self.hdmr_level+1):
      nr=0
      print '\n=================================='
      print '      STARTING ORDER %i RUNS' %i +'      '
      print '==================================\n'
      new=self.createROMLevel(i,ie)
      for key,value in new.iteritems():
        nr+=1
        self.ROMs[key]=value
      print '\nnumber of order %i runs: %i' %(i,nr)
      numruns[i]=nr
    xs={}
    for key,value in self.varDict.iteritems():
      xs[key]=1
    for key,value in self.ROMs.iteritems():
      print 'sampled',key,':',self.ROMs[key].sample(xs)
    print '\nROMs per level:'
    for key,value in numruns.iteritems():
      print ' ',key,value

  def createROMLevel(self,lvl,ie):
    ROMs={}
    runs = self.getRunSet(lvl)
    for run,vrs in runs.iteritems():
      print '\nStarting run:',run
      chlist,ident = self.makeCase(run)
      runfile = ie.writeInput(self.unc_inp_file,chlist,ident)
      inp_file = GetPot(Filename=runfile)
      exdict={}
      for i in range(len(run)):
        exdict[run[i]]=vrs[i]
      ex = Executor.ExecutorFactory('SC',exdict,inp_file)
      ex.run()
      ROMs[ident]=ex.ROM
    return ROMs

  def makeCase(self,chvars):
    changelist={}
    changelist['Variables/names']=' '.join(chvars)
    ident = 'hdmr_'+'_'.join(chvars)
    changelist['Backend/outLabel']=ident
    return changelist,ident

  def finishUp(self):
    elapsed=time.time()-self.starttime
    print 'Driver run time:',elapsed,'sec'
    print '\nStarting postprocessing...'
    makePDF = self.input_file('Backend/makePDF',0)
    if makePDF:
      print '...sampling ROM...'
      numSamples = self.input_file('Backends/PDFsamples',-1)
      if numSamples==-1:
        print '...Backends/PDFsamples not found; using 1e4...'
        numSamples = int(1e4)
      self.makePDF(numSamples)

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
  #print sys.argv,type(sys.argv)
  drv = Driver(sys.argv)
