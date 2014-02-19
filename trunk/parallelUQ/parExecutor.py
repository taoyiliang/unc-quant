import numpy as np
from math import ceil
import scipy.stats as stt
import time
import datetime as dt
import os
from GetPot import GetPot
from sys import *
import InputEditor
import Variables
import Samplers as spr
import Backends as be
import matplotlib.pyplot as plt
import cPickle as pk
from itertools import product as allcombos

import multiprocessing
from multiprocessing.queues import SimpleQueue as sque

def newExecutor(exec_type,inp_file,restart):
  oktypes=['PCESC','MC','MLMC']
  if exec_type not in oktypes:
    msg = 'Desired exec type <'+exec_type+'> not recognized.'
    msg+= '\n    Options include '+str(oktypes)
    raise IOError(msg)
  todo = 'ex = '+exec_type+'Exec(inp_file,restart)'
  exec todo
  return ex

class Executor(object): #base class
  def __init__(self,inp_file,restart):
    '''Constructor'''
    print '\nInitializing Executor...'
    self.input_file = inp_file
    self.restart=restart
    #run
    self.loadInput()
    self.setDirs()
    self.loadVars()
    self.setCase()
    self.loadSampler()
    self.loadBackends()
    self.clearInputs()
    self.parallelRun(restart=self.restart)
    self.runBackends()

  def savestate(self,savedict):
    print 'Executor has no implemented savestate function!'

  def loadRestart(self):
    print 'Executor has no implemented loadRestart function!'

  def loadInput(self):
    self.templateFile=self.input_file('Problem/templateFile','')
    self.execDir    = self.input_file('Problem/execDir'  ,'')
    self.inputDir   = self.input_file('Problem/inputDir' ,'')
    self.input_editor = self.input_file('Problem/input_editor','')
    torun = 'self.ie = InputEditor.'+self.input_editor+'(\''+self.execDir+'\')'
    exec torun

  def setDirs(self):
    '''
    Defines working directories, for UQ files, run file, IO files
    Input: none
    Output: none
    '''
    #set directories
    self.uncDir  = os.getcwd()
    os.chdir(self.execDir)
    self.execDir = os.getcwd()
    os.chdir(self.uncDir)
    #print 'in dir:',inputDir
    os.chdir(self.inputDir)
    self.inputDir = os.getcwd()
    os.chdir(self.uncDir)

  def loadVars(self):
    '''
    Fills dictionary with variables to have their U Q'd
    Input: none
    Output: none
    '''
    print '\nLoading uncertain variables...'
    uVars = self.input_file('Variables/names','').split(' ')
    self.varDict={}
    for var in uVars:
      path=self.input_file('Variables/'+var+'/path',' ')
      dist=self.input_file('Variables/'+var+'/dist',' ')
      args=self.input_file('Variables/'+var+'/args',' ').split(' ')
      for a,arg in enumerate(args):
        args[a]=float(arg)
      self.varDict[var]=Variables.newVar(dist,var,path)
      self.varDict[var].setDist(args)

  def clearInputs(self):
    '''
    Deletes old run input files from UQ runs
    Input: none
    Output: none
    '''
    print '\nAttempting to clear old inputs...'
    os.chdir(self.inputDir)
    os.system('rm '+self.templateFile+'.unc*')

  def cprToActual(self,x):
    '''
    Optional "actual" solution to compare to
    Input: variable value
    Output: function evaluated at value
    '''
    #return np.exp(-x*x)
    #return x*x - 2.*x + 0.5
    return 1.+2.*x

  def runBackends(self):
    '''
    Runs each post-processing step requested
    Input: none
    Output: none
    '''
    trunc = self.input_file('Backend/ROM/truncate',0)
    for btype in self.backends.keys():
      backend = self.backends[btype]
      if btype in ['PDF']:
        backend.MCcalc(self.histories)
      elif btype in ['plot2D']:
        backend.compareVars(self.histories)
      elif btype in ['ROM','rom']:
        backend.makeROM(self.histories)
        samples=self.input_file('Backend/ROM/samples',100.)
        samples=int(samples)
        checkSoln=self.input_file('Backend/ROM/checkSolution',0)
        makePDF=self.input_file('Backend/ROM/createPDF',0)
        numbins=self.input_file('Backend/ROM/bins',100)
        doStats=self.input_file('Backend/ROM/stats',0)
        MCROMs=None
        #sample ROM
        if checkSoln:
          print 'Sampling ROM...'
          MCROMs,correct = backend.MCsampleDict(self.histories,\
                                                trials=samples,\
                                                cprFunc=self.cprToActual)
          print 'Checking convergence to actual soln...'
          numbad=0
          for i,key in enumerate(MCROMs.keys()):
            if not np.abs(MCROMs[key]-correct[key])<1e-9:
              print 'Mismatch!',key,MCROMs[key],correct[key]
              numbad+=1
          if numbad==0: print '  ...Converged!'
          else:
            print  '  ...ERROR: not converged!',numbad,'errors found'

        if doStats:
          backend.setStats(self.histories)
        #pdf rom

        if makePDF:
          print '\nConstructing discrete pdf...'
          if MCROMs == None:
            bins,ctrs=backend.makePDF(self.histories,
                                      numbins,
                                      samples,
                                      trunc=trunc)
          else:
            bins,ctrs=backend.makePDF(self.histories,
                                      numbins,
                                      samples,
                                      MCROMs.values())
          plt.plot(ctrs,bins)
          plt.title('PDF solution, 1e'+str(int(np.log10(samples)))+' samples')
          plt.xlabel('Solutions')
          plt.ylabel('Probability')
          #plt.axis([0.9,1.2,0,0.05])
          #print 'dumping to',os.getcwd(),'4.pk'
          #pk.dump([ctrs,bins],open('4.pk','wb'))

        #write stats out to file
        #if doStats and makePDF:
        #  outFile=file(self.BEoutFile,'w')
        #  for i,ctr in enumerate(ctrs):
        #    msg=','.join([str(ctr),str(bins[i])])
        #    outFile.writelines(','.join([str(ctr),str(bins[i]),])+'\n')
        #  outFile.close()

    print 'Executioner complete...'







class PCESCExec(Executor):
  def setCase(self):
    case = ''
    case+= self.templateFile.split('.')[0]+'_SC_'
    case+= str(len(self.varDict.keys()))+'var'
    self.case=case

  def loadSampler(self):
    self.sampler = spr.StochasticPoly(self.varDict,self.input_file)

  def loadBackends(self):
    '''
    Sets up post-processing
    Input: none
    Output: none
    '''
    backendTypes = ['ROM']#self.input_file('Backend/active','').split(' ')
    self.BEoutFile = self.case+'.SC'
    #self.input_file('Backend/outFile','BE.out')
    self.backends={}
    for beType in backendTypes:
      if beType == 'ROM':
        #TODO bad full tensor product, so fix this
        orderlist=[]
        for var in self.varDict.values():
          orderlist.append(range(var.expOrd))
        tensorprod=list(allcombos(*orderlist))
        backend = be.ROM(tensorprod,self.BEoutFile)
        self.backends['ROM']=backend

  def parallelRun(self,restart=False):
    '''
    Runs Sampler in parallel and collects solns
    Input: none
    Output: none
    '''
    self.outq=sque() #to dump (nRun,soln) in
    ps=[]
    wantprocs = self.input_file('Problem/numprocs',1)
    self.numprocs = min(wantprocs,multiprocessing.cpu_count())
    if not self.restart:
      print '\nRunning samples in parallel...'
      self.total_runs = -1
      ps=[]
      self.done=False
      self.histories={}
      self.histories['varNames']=self.varDict.keys()
      self.histories['vars']=self.varDict.values()
      self.histories['varPaths']=[]
      for var in self.varDict.values():
        #need to preset these so they can be filled with nRun?
        self.histories['varPaths'].append(var.path)
        self.histories['varVals']=[]
        self.histories['nRun']=[]
        self.histories['soln']=[]
      print '  ...uncertain variables:',self.histories['varNames']
      tempOutFile = file('solns.out','w')
    else: #start from restart
      print '\nStarting from restart...'
      print '  Starting after run',self.total_runs
    print '  ...using',self.numprocs,'processors...'
    self.numPs=0
    finished=0
    while not self.done:
      #remove dead processses
      for p,proc in enumerate(ps):
        if not proc.is_alive():
          proc.join()
          while not self.outq.empty():
            n,sln = self.outq.get()
            self.histories['soln'][n]=sln
            finished+=1
          print 'Runs finished:',finished,
          del ps[p]
          #save state
      if not self.sampler.converged:
        while len(ps)<self.numprocs and not self.sampler.converged:
          self.numPs+=1
          self.total_runs+=1
          #TODO print single line, started/finished runs
          print 'Runs started:',self.total_runs,'\r',
          try:
            tot1 = self.backends['PDF'].tot1
            tot2 = self.backends['PDF'].tot2
            N = self.backends['PDF'].N
            runDict = self.sampler.giveSample([N,tot1,tot2])
          except KeyError:
            runDict = self.sampler.giveSample()
          self.histories['nRun'].append(self.total_runs)
          self.histories['soln'].append(0) #placeholder
          #print self.sampler.type
          for key in runDict.keys():
            try: self.histories[key].append(runDict[key])
            except KeyError:
              self.histories[key]=[runDict[key]]
          #  print '  ',key,self.histories[key][-1]
          #add flux output identifier to change dict
          #TODO this expects the same output block for everything!
          runDict['fileChange']={}
          outFileName='run.out'+str(self.total_runs)
          runDict['fileChange']['Output/file']=outFileName
          mesh_size = self.input_file('Problem/mesh_factor',1)
          runDict['fileChange']['Mesh/nx_per_reg']=mesh_size
          runDict['fileChange']['Mesh/ny_per_reg']=mesh_size
          inp_file = self.ie.writeInput(self.templateFile,
                               self.inputDir,
                               self.histories['varPaths'],
                               runDict['varVals'],
                               runDict['fileChange'],
                               self.total_runs)
          ps.append(multiprocessing.Process(\
                      target = self.runSample,\
                      args=(self.total_runs,inp_file,outFileName))\
                   )
          ps[-1].start()
      self.done = (len(ps)==0) and self.sampler.converged

  def runSample(self,nRun,infile,outFile):
    '''
    Child process to run a single sample
    Inputs:
      - nRun: identifier for sorting run into histories later
      - inp_file: the UQ-made version of the input file to run
      - outFile: the provided name for the executable to put the soln in
    Output: none
    '''
    self.ie.runSolve(infile)
    soln = self.ie.storeOutput(outFile)
    self.outq.put([nRun,soln])
    #print 'Finished',str(nRun)+':',soln





class MCExec(Executor):
  def setCase(self):
    case = ''
    case+= self.templateFile.split('.')[0]+'_MC_'
    case+= str(len(self.varDict.keys()))+'var'
    self.case=case

  def savestate(self,savedict):
    curcwd = os.getcwd()
    os.chdir(self.uncDir)
    outfile=file('MC.backup.pk','w')
    pk.dump(savedict,outfile,2)
    outfile.close()
    os.chdir(curcwd)

  def loadRestart(self):
    print '\nLoading MC restart from',os.getcwd()+'/MC.backup.pk'
    outfile=file('MC.backup.pk','r')
    toret = pk.load(outfile)
    outfile.close()
    return toret

  def loadSampler(self):
    self.solnRange=self.input_file('Sampler/MC/solnRange','-1e30 1e30').split()
    for i,r in enumerate(self.solnRange):
      self.solnRange[i]=float(r)
    self.sampler = spr.MonteCarlo(self.varDict,self.input_file)

  def loadBackends(self):
    '''
    Sets up post-processing, including building ROM for SC
    Input: none
    Output: none
    '''
    backendTypes = ['PDF']#self.input_file('Backend/active','').split(' ')
    self.BEoutFile = self.case+'.MC.stats'
    #self.input_file('Backend/outFile','BE.out')
    self.backends={}
    for beType in backendTypes:
      if beType == 'PDF':
        bins = self.input_file('Backend/PDF/bins',100)
        low  = self.input_file('Backend/PDF/low',0.0)
        hi   = self.input_file('Backend/PDF/hi',10.0)
        backend = be.Stats(low,hi,self.BEoutFile,bins)
        if self.restart:
          loaddict=self.loadRestart()
          runs=backend.loaddict(loaddict)
          self.total_runs=runs
          self.sampler.restart(runs)
        self.backends['PDF']=backend

  def parallelRun(self,restart=False):
    '''
    Runs Sampler in parallel and collects solns
    Input: none
    Output: none
    '''
    self.outq=sque() #to dump (nRun,soln) in
    ps=[]
    wantprocs = self.input_file('Problem/numprocs',1)
    self.numprocs = min(wantprocs,multiprocessing.cpu_count())
    trackThousand = 0
    trials = int(self.sampler.totalSamples)
    print '\nRunning %1.0e samples in parallel...' %trials
    ps=[]
    self.numPs=0
    self.done=False
    self.histories={}
    self.histories['varNames']=self.varDict.keys()
    self.histories['vars']=self.varDict.values()
    self.histories['varPaths']=[]
    for var in self.varDict.values():
      #need to preset these so they can be filled with nRun?
      self.histories['varPaths'].append(var.path)
      self.histories['varVals']=[]
    print '  ...uncertain variables:',self.histories['varNames'],'...'
    if not restart:
      self.total_runs = -1
      trialsLeft = trials
      #tempOutFile = file('solns.out','w')
    else: #start from restart
      print '\nStarting from restart...'
      print '  Starting after run',self.total_runs
      trialsAtRst = self.total_runs
      trialsLeft = trials - self.total_runs
    print '  ...using',self.numprocs,'processors...'

    printFreq = self.input_file('Sampler/MC/printfreq',1000)
    print '  ...print frequency is',printFreq,'...'

    trialsPerProc = int(trials/float(self.numprocs))
    trialsPerProc = min(trialsPerProc,int(ceil(printFreq/4)))
    if trialsPerProc > 1000:
      trialsPerProc = 1000
    mesh_size = self.input_file('Problem/mesh_factor',1)

    self.done=False
    runDict={}
    starttime=time.time()
    rcvProc=0
    lastPrint=0
    doAPrint = False
    thrown=0
    print '\nFinished Run | Time Elapsed | Est. Remaining',
    print '| Number Discarded Solutions'
    while not self.done:
      #remove dead processses
      for p,proc in enumerate(ps):
        if not proc.is_alive():
          proc.join()
          del ps[p]
          rcvProc+=1
          while not self.outq.empty():
            slns,newthrown = list(self.outq.get())
            thrown+=newthrown
            lastPrint+=len(slns)
            if lastPrint >= printFreq:
              self.backends['PDF'].addToBins(slns,True)
              doAPrint=True
              lastPrint=0
            else:
              self.backends['PDF'].addToBins(slns)
          self.savestate(self.backends['PDF'].savedict())
          if rcvProc==self.numprocs:
            rcvProc=0
            if doAPrint:
              doAPrint = False
              lastPrint=0
              #print progress
              #FIXME fix for restart case
              finished = trials-trialsLeft
              if self.restart:
                finished -= trialsAtRestart
              elapTime = time.time()-starttime
              dpdt = float(finished)/float(elapTime)
              toGo = dt.timedelta(seconds=(int(trialsLeft/dpdt)))
              elapTime = dt.timedelta(seconds=int(elapTime))
              print '%12i | %12s | %12s | %9i' %(finished,elapTime,toGo,thrown),
              print '                      \r',
      if trialsLeft > 0:
        while len(ps)<self.numprocs and not self.done:
          if trialsLeft > trialsPerProc: newtrials = trialsPerProc
          else: newtrials = trialsLeft
          trialsLeft -= newtrials
          self.numPs+=1
          runDict['fileChange']={}
          runDict['fileChange']['Mesh/nx_per_reg']=mesh_size
          runDict['fileChange']['Mesh/ny_per_reg']=mesh_size
          ps.append(multiprocessing.Process(\
                      target = self.runSample,\
                      args=(trackThousand,newtrials,runDict)\
                   ))
          ps[-1].start()
          trackThousand+=1
      self.done = self.outq.empty() and len(ps)==0
    print '\n'

  def runSample(self,prefix,batch,runDict):
    np.random.seed()
    solns=[]
    nThrown = 0
    for i in range(batch):
      ident = '_'+str(prefix)+'_'+str(i)
      outFileName='run.out'+ident
      runDict['fileChange']['Output/file']=outFileName
      runDict['varVals']=self.sampler.giveSample()
      inp_file = self.ie.writeInput(self.templateFile,
                               self.inputDir,
                               self.histories['varPaths'],
                               runDict['varVals'],
                               runDict['fileChange'],
                               ident)
      self.ie.runSolve(inp_file)
      soln = self.ie.storeOutput(outFileName)
      if self.solnRange[0] < soln < self.solnRange[1]:
        solns.append(soln)
      else:
        nThrown+=1
        #print 'WARNING! Solution outside of range.  Tossed',soln,'\n'
    self.outq.put([solns,nThrown])


#if __name__=='__main__':
  #run syntax: python parExecutor -i <UQ input file>
#  if '-r' in argv:
#    restart = True
#    print 'setting up restart...'
#  else: restart = False
#  if 'mc' in argv:
#    ex = MCExec(argv,restart=restart)
#  elif 'sc' in argv:
#    ex = PCESCExec(argv,restart)
#  else:
#    print 'Need "mc" or "sc" in argument list!'
