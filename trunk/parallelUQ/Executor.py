import os
import time
import multiprocessing
from multiprocessing.queues import Queue as que
import cPickle as pk
import datetime as dt
from sys import *
from math import ceil
from itertools import product as allcombos
from bisect import bisect_left

import numpy as np
import scipy.stats as stt
import matplotlib.pyplot as plt
from GetPot import GetPot
#from memory_profiler import profile

import InputEditor
import Variables
import Samplers as spr
import Backends as be
import IndexSets
import SparseQuads
import ROM
from tools import makePDF


def ExecutorFactory(exec_type,varlist,inp_file):
  oktypes=['PCESC','SC','MC','MLMC']
  if exec_type not in oktypes:
    msg = 'Desired exec type <'+exec_type+'> not recognized.'
    msg+= '\n    Options include '+str(oktypes)
    raise IOError(msg)
  todo = 'ex = '+exec_type+'Exec(varlist,inp_file)'
  exec todo
  return ex

class Executor(object): #base class
  def __init__(self,varDict,inp_file):
    '''Constructor'''
    #print '\nInitializing Executor...'
    self.outq = que()
    self.input_file = inp_file
    self.varDict = varDict

  def run(self,verbose=False):
    self.loadInput(verbose)
    if verbose: print '...input loaded...'
    self.setDirs()
    if verbose: print '...directories set...'
    self.loadSampler()
    if verbose: print '...sampler loaded...'
    self.loadBackends()
    if verbose: print '...backends loaded...'
    self.clearInputs(verbose)
    if verbose: print '...old inputs cleared...'
    try:
      self.runParallel(verbose)
    except KeyboardInterrupt:
      self.checkConverge(True,verbose)
    if verbose: print '...parallel run complete...'
    return self.finish()

  def loadInput(self,verbose=False):
    self.templateFile=self.input_file('Problem/templateFile','')
    self.execDir    = self.input_file('Problem/execDir'  ,'')
    self.inputDir   = self.input_file('Problem/inputDir' ,'')
    self.input_editor = self.input_file('Problem/input_editor','')
    torun = 'self.ie = InputEditor.'+self.input_editor+'(\''+self.execDir+'\')'
    exec torun
    mesh_size = self.input_file('Problem/mesh_factor',-1)
    self.meshFactor = mesh_size
    self.setSetType()
    self.setCase()
    print '...case set <'+self.case+'>...'

  def setSetType(self):
    pass

  def setDirs(self):
    #set directories
    self.uncDir  = os.getcwd()
    os.chdir(self.execDir)
    self.execDir = os.getcwd()
    os.chdir(self.uncDir)
    #print 'in dir:',inputDir
    os.chdir(self.inputDir)
    self.inputDir = os.getcwd()
    os.chdir(self.uncDir)

  def setCase(self):
    pass #overwritten

  def loadSampler(self):
    pass

  def loadBackends(self):
    pass

  def clearInputs(self,verbose=False):
    if verbose: print '\nAttempting to clear old inputs...'
    os.chdir(self.inputDir)
    fail=os.system('rm '+self.templateFile+'.unc*')
    if (not fail) and verbose:
      print '...successfully cleared old input files.'

  #@profile
  def runParallel(self,verbose=False):
    wantprocs = self.input_file('Problem/numprocs',1)
    self.numprocs = min(wantprocs,multiprocessing.cpu_count())
    if verbose: print '\n...using %i processers...' %self.numprocs
    try: self.total_runs
    except AttributeError: self.total_runs=-1
    procs=[]
    self.done=False
    try: self.histories
    except AttributeError: self.histories={}
    self.histories['varNames'] = self.varDict.keys()
    self.histories['vars']=self.varDict.values()
    self.histories['varVals']=[]
    self.histories['nRun']=[]
    self.histories['soln']=[]
    self.histories['solwt']=[]
    self.histories['varPaths']=[]
    if verbose: print '...uncertain variables:',self.histories['varNames'],'...'
    for var in self.varDict.values():
      self.histories['varPaths'].append(var.path)
    #preset some file change stuff
    self.totalProcs=0
    try: finished=self.histories['entries']
    except KeyError: finished = 0
    while not self.done:
      for p,proc in enumerate(procs):
        if not proc.is_alive():
          proc.join()
          while not self.outq.empty():
            #TODO taking values needs to be specific to executor type
            self.readSolutions(self.outq.get())
            #n,wt,soln = self.outq.get()
            #self.histories['soln'][n]=soln
            #self.histories['solwt'][n]=wt
            finished+=1
          print 'Runs finished:',finished,
          del procs[p]
          self.checkConverge(verbose)
      if not self.sampler.converged:
        while len(procs)<self.numprocs and not self.sampler.converged:
        #while len(procs)<self.numprocs:
          try:
            runDict = self.sampler.next()
            #print 'added runvals:',runDict['varVals']
          except StopIteration: break
          self.totalProcs+=1
          self.total_runs+=1
          print '...runs started:',self.total_runs,'...\r',
          #self.histories['nRun'].append(self.total_runs)
          runDict['nRun'] = self.total_runs
          self.setupHistories(runDict)
          runDict['fileChange']={}
          runDict['outFileName']='run.out'+str(self.total_runs)
          runDict['fileChange']['Output/file']=runDict['outFileName']
          if self.meshFactor > 0:
            runDict['fileChange']['Mesh/nx_per_reg']=self.meshFactor
            runDict['fileChange']['Mesh/ny_per_reg']=self.meshFactor
          runDict['inp_file'] = self.ie.writeInput(self.templateFile,
                                self.inputDir,
                                self.histories['varPaths'],
                                runDict['varVals'],
                                runDict['fileChange'],
                                self.total_runs)
          procs.append(multiprocessing.Process(target=self.runSample,
                           args=[runDict]))
          procs[-1].start()
      self.done = len(procs)==0 and self.sampler.converged
      time.sleep(0.1)

  def runSample(self,runDict):
    self.ie.runSolve(runDict['inp_file'])
    soln = self.ie.storeOutput(runDict['outFileName'])
    return soln

  def runPDFSample(self,M,bounds):
    np.random.seed()
    bins = np.zeros(len(bounds)-1)
    runs = 0
    hi=-1e14
    low=1e14
    numhi=0
    numlow=0
    while runs < M:
      uvars = self.varDict.values()
      sampVals = np.zeros(len(uvars))
      for i in range(len(sampVals)):
        sampVals[i]=uvars[i].sample()
      soln = self.ROM(sampVals)
      low=min(low,soln)
      hi=max(hi,soln)
      if soln > bounds[-1]:
        numhi+=1
      elif soln < bounds[-1]:
        numlow+=1
      else:
        indx = bisect_left(bounds,soln)
        bins[i-1]+=1
      runs+=1
    self.pdfque.put([bins,numlow,numhi,low,hi])

  def finish(self):
    os.chdir(self.uncDir)
    print 'Executor complete.'

class SC(Executor):
  '''Intermediary before choosing SC type'''

  def setCase(self):
    inp = self.input_file('Backend/outLabel','')
    self.case = self.settype+'_h'+str(self.meshFactor)+'_'+inp

  def setSetType(self):
    self.settype=self.input_file('Sampler/SC/indexSet','dud')
    if self.settype=='dud':
      print 'Index set not specified; using tensor product.'
      self.settype='TP'

  def loadSampler(self):
    self.loadIndexSet()
    self.loadQuadSet()

  def loadIndexSet(self):
    self.expOrder = self.input_file('Sampler/SC/expOrd',-1)
    if self.expOrder==-1:
      print '...expansion order not set in Sampler/SC.  Using 2...'
      self.expOrder=2
    impwts=[]
    for var in self.varDict.values():
      impwts.append(var.impwt)
    self.indexSet = IndexSets.IndexSetFactory(len(self.varDict.keys()),
                                              self.expOrder,
                                              self.settype,
                                              impwts)
    print '...%i expansion indices used...' %len(self.indexSet)

  def loadQuadSet(self,verbose=False):
    def single(x):
      return x
    def double(x):
      return 2**x
    rule = self.input_file('Sampler/SC/quadrule','<None>')
    okRules=['single','double']
    if rule not in okRules:
      print '...quadrule not recognized in Sampler/SC. Using single...'
      rule='single'
    todo = 'quadrule='+rule
    exec todo
    self.quadrule = quadrule
    #for i in self.indexSet:
    #  print 'index point:',i
    if verbose: print '...constructing sparse grid...'
    wantprocs = self.input_file('Problem/numprocs',1)
    self.numprocs = min(wantprocs,multiprocessing.cpu_count())
    grid = SparseQuads.parBasicSparse(self.numprocs,
                                   len(self.varDict.keys()),
                                   self.expOrder,
                                   self.indexSet,
                                   quadrule,
                                   self.varDict,
                                   verbose)
    run_samples = {}
    run_samples['varNames']=self.varDict.keys()
    run_samples['variables']=self.varDict.values()
    run_samples['quadpts']=[]
    run_samples['weights']={}
    if verbose:
      print '  ...removing duplicate quadrature points...'
      print '    ...number of pts including duplicates: %i' %len(grid)
  # NEW ROUNDED POINT METHOD
    for e,entry in enumerate(grid):
      npt = tuple(np.around(entry[0],decimals=15))
      if npt in run_samples['quadpts']:
        run_samples['weights'][npt]+=entry[1]
      else:
        run_samples['quadpts'].append(npt)
        run_samples['weights'][npt]=entry[1]
        numsamp = len(run_samples['quadpts'])
        if verbose:
          print '    ...new number of quad pts collected:',
          print numsamp,
          print '(%i pct complete)' %(int(100.*(e+1)/len(grid))),
          print '(%i duplicates)' %(e+1-numsamp),'\r',
    #exit()
    #TODO DEBUG
    #for pt in run_samples['quadpts']:
      #print 'quad pts,wts: (%1.4e), %1.4e' %(pt[0],run_samples['weights'][pt])
    #print 'sum wts:',sum(run_samples['weights'].values())
    if verbose:print '...constructing sampler...'
    self.sampler = spr.StochasticPoly(self.varDict,
                                      run_samples)
    self.numquadpts = len(run_samples['quadpts'])
    if verbose:print '...%i quadrature points used...' %self.numquadpts

  def setupHistories(self,runDict):
    self.histories['soln'].append(0)  #placeholders
    self.histories['solwt'].append(0)
    for key in runDict.keys():
      try:self.histories[key].append(runDict[key])
      except KeyError: self.histories[key]=[runDict[key]]

  def runSample(self,runDict):
    soln = super(SC,self).runSample(runDict)
    self.outq.put([runDict['nRun'],runDict['quadWts'],soln])

  def readSolutions(self,out):
    n,wt,soln = out
    self.histories['soln'][n]=soln
    self.histories['solwt'][n]=wt

  def checkConverge(self,force=False,verbose=False):
    print '                        \r',

  def writeOut(self):
    mean = self.ROMmoment(1)
    r2 = self.ROMmoment(2)
    var = r2 - mean*mean
    name = self.case+'.moments'
    print '...writing to',name,'...'
    outFile = file(name,'a')
    outFile.writelines('\nMoments\nN,mean,var\n')
    expv='%1.16e' %mean
    svar='%1.16e' %var
    outFile.writelines(','.join([str(self.numquadpts),
      expv,svar])+'\n')
    outFile.close()
    if self.done:
      print '\n'
      print 'N, mean, var:'
      print self.numquadpts,mean,var
      print '\n'


class SCExec(SC):
  def setCase(self):
    inp = self.input_file('Backend/outLabel','')
    self.case = self.settype+'_h'+str(self.meshFactor)+'_'+inp

  def finish(self):
    self.ROM = ROM.LagrangeROM(self.histories['soln'],
                               self.histories['solwt'],
                               self.varDict,
                               self.histories['varVals'],
                               self.indexSet,
                               self.quadrule,
                               self.numprocs)
    super(SC,self).finish()



#####################################################
#####################################################
#####################################################


class MC(Executor):
  def setCase(self):
    inp = self.input_file('Backend/outLabel','')
    self.case = 'MC_h'+str(self.meshFactor)+'_'+inp

  def loadSampler(self):
    self.storeSolns=bool(self.input_file('Sampler/MC/writeSamples',0))
    loadSolns=bool(self.input_file('Sampler/MC/loadSamples',0))

    if loadSolns:
      try:
        inFile = file(self.uncDir+'/'+self.case+'.samples','r')
        for line in inFile:
          pass
        n,one,two=line.strip().split(',')
        n=int(n)
        one=float(one)
        two=float(two)
        self.histories={}
        self.histories['entries']=n
        self.histories['first']=one
        self.histories['second']=two
        self.total_runs=n
        print '...loaded %i previous samples...' %n
      except IOError: pass


    self.maxM = self.input_file('Sampler/MC/maxSamples','0')
    self.maxM = int(float(self.maxM))
    self.targetTol = self.input_file('Sampler/MC/convergence',0.0)
    if self.targetTol==0.0:
      if self.maxM==0: #neither tol nor M is set
        raise IOError('In Sampler/MC/ neither maxSamples nor convergence '+\
            'criteria are set.  Please set one.')
      else: #M is set, but tol is not -> use M criteria only
        self.targetTol = 1e-40
    else:
      if self.maxM==0: #tol is set, but not M -> use tol criteria only
        self.maxM = int(1e20)
    print '...max samples to run: %1.0e...' %self.maxM
    print '...mean tol to converge: %1.1e...' %self.targetTol
    self.timesConverged = 0
    self.sampler = spr.MonteCarlo(self.varDict)
    self.N = 0
    self.mean = 1e14
    self.var = 1e14

  def setupHistories(self,runDict):
    pass #nothing needed

  def runSample(self,runDict):
    #TODO make this run a batch of samples
    soln = super(MC,self).runSample(runDict)
    self.outq.put([runDict['nRun'],soln])

  def readSolutions(self,out):
    n,soln = out
    try: self.histories['entries']+=1
    except KeyError: self.histories['entries']=1

    try: self.histories['first']+=soln
    except KeyError: self.histories['first']=soln

    try: self.histories['second']+=soln*soln
    except KeyError: self.histories['second']=soln*soln

    if self.storeSolns:
      name = self.case+'.samples'
      #if os.path.isfile(name) and self.histories['entries']==1:
      try: self.slnFile
      except AttributeError: self.slnFile = file(self.uncDir+'/'+name,'a')
      if self.histories['entries']==1:
        #msg='  WARNING: sample file exists!  Append solutions?'
        #ctu = raw_input(msg)
        #if ctu.strip().lower()[-1] in ['y','t']:
        #else:
        #  self.slnFile = file(name,'w')
        self.slnFile.writelines('N,sum,sum^2\n')
      self.slnFile.writelines('%i,%1.16e,%1.16e' %(
                              self.histories['entries'],
                              self.histories['first'],
                              self.histories['second'])+'\n')
      self.slnFile.flush()

  def writeOut(self):
    name = self.case+'.moments'
    print '...writing to',name,'...'
    outFile = file(name,'a')
    outFile.writelines('\nMoments\nN,mean,var\n')
    outFile.writelines(','.join(
      [str(self.N),str(self.mean),str(self.var)])+'\n')
    outFile.close()
    if self.done:
      print '\n'
      print 'N, mean, 2nd moment:'
      print self.N,self.mean,self.var+self.mean*self.mean
      print '\n'


class MCExec(MC):
  #def setCase(self):
  #  case = ''
  #  case+= self.templateFile.split('.')[0]+'_AMC_'
  #  case+= str(len(self.varDict.keys()))+'var'
  #  self.case=case

  def moments(self):
    N = self.histories['entries']
    mean = self.histories['first']/float(N)
    second = self.histories['second']/float(N)
    var = second - mean*mean
    #print '%i %1.4e %1.4e' %(N,mean,var),
    return N,mean,var
    #print '\n\nStatistics:'
    #print 'N Runs:',N
    #print 'Mean  :',mean
    #print 'Var   :',var

  def checkConverge(self,force=False):
    if force:
      print '\n\nEXECUTOR TERMINATED\n'
      print 'N,mean,var | convergence mean, var'
      print self.N,self.mean,self.var
      print ''
      raise KeyboardInterrupt
    oN = self.N
    oMean = self.mean
    oVar = self.var
    self.N,self.mean,self.var = self.moments()
    if self.N<5 or self.sampler.converged:
      return
    #if self.N<5:
    #  return
    #print 'N,mean,var',self.N,self.mean,self.var
    convMean = abs(self.mean - oMean)/self.mean
    convVar = abs(self.var - oVar)/self.var
    print '| %1.4e %1.4e' %(convMean,convVar),'\r',
    #TODO what about variance convergence?
    if convMean <= self.targetTol:# and convMean!=0.0:
      self.timesConverged+=1
      #print '\nconverged:',self.timesConverged,convMean
      #print ''
    else: self.timesConverged = 0
    if self.N>=self.maxM or self.timesConverged >= 10:# and self.done:
      self.sampler.converged = True

class MLMCExec(MC):
  def setCase(self):
    case = ''
    case+= self.templateFile.split('.')[0]+'_MLMC_'
    case+= str(len(self.varDict.keys()))+'var'
    self.case=case

  #TODO overwrite the runParallel process
