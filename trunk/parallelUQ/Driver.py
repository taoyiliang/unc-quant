import os
import sys
import time
import datetime as dt
import multiprocessing
import cPickle as pk
from itertools import combinations as combos

import numpy as np
from GetPot import GetPot
from matplotlib.pyplot import show

import Executor
import Variables
import InputEditor
import ROM

class Driver(object):
  def __init__(self,argv,verbose={})
    self.verbose={'Backends':False,
                  'Driver':False,
                  'warning':True,
                  'Executor':False,
                  'HDMR':False,
                  'IndexSets':False,
                  'InputEditor':False,
                  'ROM':False,
                  'run':True,
                  'Samplers':False,
                  'SparseQuads':False,
                  'Variables':False}):
    for key,val in verbose.iteritems():
      self.verbose[key]=val
    self.starttime=time.time()
    print '\nStarting UQ at',dt.datetime.fromtimestamp(self.starttime)
    self.loadInput(argv)
    self.loadVars()
    self.setRuns()
    self.createROMs()
    self.finishUp()

  def loadInput(self,argv):
    cl = GetPot(argv)
    if cl.search('-i'):
      self.unc_inp_file=cl.next('')
      if verbose['Driver']:
        print 'Selected uncertainty input',self.unc_inp_file,'...'
      #check desired file exists
      if not os.path.isfile(self.unc_inp_file):
        raise IOError('Uncertainty input not found: '+self.unc_inp_file)
      if verbose['Driver']:print '  ...found uncertainty input file.'
    else:
      raise IOError('Requires and input file using -i.')
    self.input_file = GetPot(Filename=self.unc_inp_file)

  def loadVars(self):
    if verbose['Driver']:print '\nLoading uncertain variables...'
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
    self.todo={}
    self.todo[tuple(self.varDict.keys())]=self.varDict.values()

  def createROMs(self):
    self.ROMs={}
    ident = self.todo.keys()[0]
    if self.verbose['run']:print '\nStarting run:',ident
    inp_file = GetPot(Filename=self.unc_inp_file)
    ex_type = inp_file('Problem/executor','')
    if self.verbose['Driver']:print 'ex type:',ex_type
    ex = Executor.ExecutorFactory(ex_type,self.varDict,inp_file)
    ex.run(verbose=self.verbose)
    try:
      self.ROMs[ident]=ex.ROM
      if verbose['Driver']:print 'sampled mean:',ex.ROM.moment(1)
    except AttributeError:
      pass #MC doesn't store a rom at this point
    self.ex = ex
    #xs={}
    #for key,value in self.varDict.iteritems():
    #  xs[key]=1
    #print 'sampled:',ex.ROM.sample(xs,verbose=False)
    return ex

  def finishUp(self):
    elapsed=time.time()-self.starttime
    if self.verbose['Driver']:
      print 'Driver run time:',elapsed,'sec'
      print '\nStarting postprocessing...'
    #makePDF = self.input_file('Backend/makePDF',0)
    #if makePDF:
    #  print '...sampling ROM...'
    #  numSamples = self.input_file('Backends/PDFsamples',-1)
    #  if numSamples==-1:
    #    print '...Backends/PDFsamples not found; using 1e4...'
    #    numSamples = int(1e4)
    #  self.makePDF(numSamples)

    needWrite = bool(self.input_file('Backend/writeOut',0))
    if needWrite:
      self.ex.writeOut()

    show()
    if self.verbose['Driver']:print '\nDriver complete.\n'

  def makePDF(self,M):
    if self.verbose['Driver']:print '\n...constructing pdf...\n'
    wantprocs = self.input_file('Problem/numprocs',1)
    numprocs = min(wantprocs,multiprocessing.cpu_count())
    nBins = self.input_file('Backends/PDF/bins',10)
    binmin = self.input_file('Backends/PDF/min',-10)
    binmax = self.input_file('Backends/PDF/max',10)
    bins=np.linspace(binmin,binmax,nBins+1)
    self.ex.makePDF(numprocs,M,bins)




class HDMR_Driver(Driver):
  def loadInput(self,argv):
    super(HDMR_Driver,self).loadInput(argv)
    self.hdmr_level = self.input_file('HDMR/level',0)
    if self.hdmr_level==0:
      if verbose['warning']:print 'HDMR level not specified.  Using 2...'
      self.hdmr_level=2

  def setRuns(self):
    self.todo = {}
    varnames = self.varDict.keys()
    for i in range(1,self.hdmr_level+1):
      new = self.addHDMRLevel(varnames,i)
      for key,value in new.iteritems():
        self.todo[key]=value

  def addHDMRLevel(self,varnames,lvl):
    todo={}
    nameset = combos(varnames,lvl)
    for entry in nameset:
      todo[entry]=[]
      for e,ent in enumerate(entry):
        todo[entry].append(self.varDict[ent])
    return todo

  def getRunSet(self,num):
    ret={}
    for key,value in self.todo.iteritems():
      if len(key)==num:
        ret[key]=value
    return ret

  def createROMs(self):
    self.ROMs={}
    ie = InputEditor.HDMR_IO()
    #reference run
    chlist,ident = self.makeCase({})
    runfile = ie.writeInput(self.unc_inp_file,chlist,ident)
    inp_file = GetPot(Filename=runfile)
    ex = Executor.ExecutorFactory('SC',{},inp_file)
    ex.run()
    os.system('rm '+runfile)
    self.ROMs[ident]=ex.ROM
    # rest of runs
    numruns={}
    for i in range(1,self.hdmr_level+1):
      nr=0
      if self.verbose['HDMR']:
        print '\n=================================='
        print '      STARTING ORDER %i RUNS' %i +'      '
        print '==================================\n'
      new=self.createROMLevel(i,ie)
      for key,value in new.iteritems():
        nr+=1
        self.ROMs[key]=value
      if self.verbose['HDMR']:print '\nnumber of order %i runs: %i' %(i,nr)
      numruns[i]=nr
    xs={}
    for key,value in self.varDict.iteritems():
      xs[key]=1
    if self.verbose['HDMR']:
      for key,value in self.ROMs.iteritems():
        print 'mean',key,':',value.moment(1)
      print '\nROMs per level:'
      for key,value in numruns.iteritems():
        print ' ',key,value
    #for rom in self.ROMs.values():
    #  pk.dump(rom.serializable(),file('hdmr_'+rom.case()+'.pk','w'))
    self.HDMR_ROM=ROM.HDMR_ROM(self.ROMs,self.varDict)
    if self.verbose['HDMR']:
      print 'Total Det. Runs:',self.HDMR_ROM.numRunsToCreate()
      print 'HDMR sampled',self.HDMR_ROM.sample(xs,self.hdmr_level)[1]
    #store output
    case = 'hdmr'
    case+= '_'+self.input_file('Sampler/SC/indexSet','')
    case+= '_N'+str(len(self.varDict))
    case+= '_H'+str(self.hdmr_level)
    #case+= '_L'+self.input_file('Sampler/SC/expOrd','')
    mean = self.HDMR_ROM.mean(self.hdmr_level,verbose=False)
    outFile = file(case+'.out','a')
    outFile.writelines('\nRuns,SC Level,Mean\n')
    outFile.writelines(str(self.HDMR_ROM.numRunsToCreate())+',')
    outFile.writelines(self.input_file('Sampler/SC/expOrd','')+',')
    outFile.writelines('%1.15e \n' %mean)
    outFile.close()
#TODO FIXME how to calculate moments?


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
      os.system('rm '+runfile)
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
    if self.verbose['Driver']:
      print 'Driver run time:',elapsed,'sec'
      print '\nDriver complete.\n'

if __name__=='__main__':
  #print sys.argv,type(sys.argv)
  if '-hdmr' in sys.argv:
    drv = HDMR_Driver(sys.argv)
  else:
    drv = Driver(sys.argv,verbose={'Driver':True})
