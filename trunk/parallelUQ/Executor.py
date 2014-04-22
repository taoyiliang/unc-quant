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

import InputEditor
import Variables
import Samplers as spr
import Backends as be
import IndexSets
import SparseQuads


def ExecutorFactory(exec_type,inp_file):
  oktypes=['PCESC','SC','MC','MLMC']
  if exec_type not in oktypes:
    msg = 'Desired exec type <'+exec_type+'> not recognized.'
    msg+= '\n    Options include '+str(oktypes)
    raise IOError(msg)
  todo = 'ex = '+exec_type+'Exec(inp_file)'
  exec todo
  return ex

class Executor(object): #base class
  def __init__(self,inp_file):
    '''Constructor'''
    print '\nInitializing Executor...'
    self.outq = que()
    self.input_file = inp_file
    #run
    self.loadInput()
    print '...input loaded...'
    self.setDirs()
    print '...directories set...'
    self.loadVars()
    print '...variables loaded...'
    self.setCase()
    print '...case set <'+self.case+'>...'
    self.loadSampler()
    print '...sampler loaded...'
    self.loadBackends()
    print '...backends loaded...'
    self.clearInputs()
    print '...old inputs cleared...'
    self.runParallel()
    print '...parallel run complete...'
    self.finish()

  def loadInput(self):
    self.templateFile=self.input_file('Problem/templateFile','')
    self.execDir    = self.input_file('Problem/execDir'  ,'')
    self.inputDir   = self.input_file('Problem/inputDir' ,'')
    self.input_editor = self.input_file('Problem/input_editor','')
    torun = 'self.ie = InputEditor.'+self.input_editor+'(\''+self.execDir+'\')'
    exec torun

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


  def loadVars(self):
    print '\nLoading uncertain variables...'
    uVars = self.input_file('Variables/names','').split(' ')
    self.varDict={}
    for var in uVars:
      path=self.input_file('Variables/'+var+'/path',' ')
      dist=self.input_file('Variables/'+var+'/dist',' ')
      args=self.input_file('Variables/'+var+'/args',' ').split(' ')
      for a,arg in enumerate(args):
        args[a]=float(arg)
      self.varDict[var]=Variables.VariableFactory(dist,var,path)
      self.varDict[var].setDist(args)

  def setCase(self):
    pass #overwritten

  def loadSampler(self):
    pass

  def loadBackends(self):
    pass

  def clearInputs(self):
    print '\nAttempting to clear old inputs...'
    os.chdir(self.inputDir)
    fail=os.system('rm '+self.templateFile+'.unc*')
    if not fail:
      print '...successfully cleared old input files.'

  def runParallel(self):
    wantprocs = self.input_file('Problem/numprocs',1)
    self.numprocs = min(wantprocs,multiprocessing.cpu_count())
    print '\n...using %i processers...' %self.numprocs
    #
    # set up processor pool TODO get pickling error - can't use instance method
    #
    #pool = multiprocessing.Pool(processes=self.numprocs)
    #res = pool.imap(self.runSample,self.sampler,int(np.sqrt(self.numquadpts)))
    #pool.close()
    #pool.join()
    #
    # use the old queue method
    #
    self.total_runs = -1
    procs=[]
    self.done=False
    self.histories={}
    self.histories['varNames'] = self.varDict.keys()
    self.histories['vars']=self.varDict.values()
    self.histories['varVals']=[]
    self.histories['nRun']=[]
    self.histories['soln']=[]
    self.histories['solwt']=[]
    self.histories['varPaths']=[]
    print '...uncertain variables:',self.histories['varNames'],'...'
    for var in self.varDict.values():
      self.histories['varPaths'].append(var.path)
    #preset some file change stuff
    mesh_size = self.input_file('Problem/mesh_factor',-1)
    self.totalProcs=0
    finished=0
    while not self.done:
      for p,proc in enumerate(procs):
        if not proc.is_alive():
          proc.join()
          while not self.outq.empty():
            n,wt,soln = self.outq.get()
            self.histories['soln'][n]=soln
            self.histories['solwt'][n]=wt
            finished+=1
          print 'Runs finished:',finished,
          del procs[p]
      if not self.done:
        while len(procs)<self.numprocs and not self.sampler.converged:
          try:
            runDict = self.sampler.next()
            print 'added runvals:',runDict['varVals']
          except StopIteration: break
          self.totalProcs+=1
          self.total_runs+=1
          print '...runs started:',self.total_runs,'...\r',
          #self.histories['nRun'].append(self.total_runs)
          runDict['nRun'] = self.total_runs
          self.histories['soln'].append(0)  #placeholders
          self.histories['solwt'].append(0)
          for key in runDict.keys():
            try:self.histories[key].append(runDict[key])
            except KeyError: self.histories[key]=[runDict[key]]
          runDict['fileChange']={}
          runDict['outFileName']='run.out'+str(self.total_runs)
          runDict['fileChange']['Output/file']=runDict['outFileName']
          if mesh_size > 0:
            runDict['fileChange']['Mesh/nx_per_reg']=mesh_size
            runDict['fileChange']['Mesh/ny_per_reg']=mesh_size
          runDict['inp_file'] = self.ie.writeInput(self.templateFile,
                                self.inputDir,
                                self.histories['varPaths'],
                                runDict['varVals'],
                                runDict['fileChange'],
                                self.total_runs)
          procs.append(multiprocessing.Process(target=self.runSample,
                           args=[runDict]))
          procs[-1].start()
        self.done = (len(procs)==0) and self.sampler.converged
        if self.done: print ''


  def makePDF(self,P,M,bounds):
    for n,sln in enumerate(self.histories['soln']):
      print self.histories['varVals'][n],'|',sln
    print 'Creating PDF by MC sample of ROM...'
    print '...using %i processors...' %P
    print '...using %i bins from %1.3e to %1.3e...' \
                             %(len(bounds)-1,bounds[0],bounds[-1])
    total_runs_finished = 0
    total_runs_started = 0
    self.pdfque = que()
    procs=[]
    bad=[0,0] #discarded solns
    rge=[1e14,-1e14]
    bins=np.zeros(len(bounds)-1)
    print 'Runs Started / Finished'
    while total_runs_finished < M:
      #collect finished solutions
      for p,proc in enumerate(procs):
        if not proc.is_alive():
          proc.join()
          del procs[p]
          while not self.pdfque.empty():
            newbins,newlow,newhi,newmin,newmax = list(self.pdfque.get())
            total_runs_finished+=len(newbins)
            bins+=newbins
            bad[0]+=newlow
            bad[1]+=newhi
            rge[0]=min(rge[0],newmin)
            rge[1]=max(rge[1],newmax)
            print '%i / %i' %(total_runs_started,total_runs_finished),'\r',
      #queue new runs
      if total_runs_started < M:
        runs_left = M - total_runs_started
        while len(procs)<P and runs_left > 0:
          if runs_left > 10: #TODO make this an input
            new_runs = 10
          else: new_runs = runs_left
          runs_left -= new_runs
          total_runs_started+=new_runs
          procs.append(multiprocessing.Process(
            target = self.runPDFSample, args=(new_runs,bounds)))
          procs[-1].start()
    print '\n'
    #normalize results
    Mgood = M - bad[0] - bad[1]
    for b,bn in enumerate(bins):
      bins[b] = bn/Mgood
    #printout
    print 'Range of solutions: %1.3e -> %1.3e' %(rge[0],rge[1])
    #plot it
    #centers = 0.5*(bounds[:-1]+bounds[1:])
    #plt.figure()
    #plt.plot(centers,bins)
    #plt.title('ROM PDF by MC, %i bins' %len(bins))
    #plt.xlabel('Solution Value')
    #plt.ylabel('Frequency')

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

  def loadSampler(self):
    self.loadIndexSet()
    self.loadQuadSet()

  def loadIndexSet(self):
    self.expOrder = self.input_file('Sampler/SC/expOrd',-1)
    if self.expOrder==-1:
      print '...expansion order not set in Sampler/SC.  Using 2...'
      self.expOrder=2
    settype=self.input_file('Sampler/SC/indexSet','dud')
    if settype=='dud':
      print 'Index set not specified; using tensor product.'
      settype='TP'
    self.indexSet = IndexSets.IndexSetFactory(len(self.varDict.keys()),
                                              self.expOrder,
                                              settype)
    print '...%i expansion indices used...' %len(self.indexSet)

  def loadQuadSet(self):
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
    for i in self.indexSet:
      print 'index point:',i
    grid = SparseQuads.BasicSparse(len(self.varDict.keys()),
                                   self.expOrder,
                                   self.indexSet,
                                   quadrule,
                                   self.varDict)
    run_samples = {}
    run_samples['varNames']=self.varDict.keys()
    run_samples['variables']=self.varDict.values()
    run_samples['quadpts']=[]
    run_samples['weights']={}
    for entry in grid:
      npt = entry[0]
      nwt = entry[1]
      alreadyThere = len(run_samples['quadpts'])>0
      for pt in run_samples['quadpts']:
        alreadyThere = True and len(run_samples['quadpts'])>0
        #if abs(npt[0])<1e-12 or abs(npt[1])<1e-12:
          #print 'checking',npt,pt
        for i,dud in enumerate(pt):
          #print npt[i]-pt[i]
          alreadyThere*= abs(npt[i]-pt[i])<1e-13
        if alreadyThere:
          npt = pt
          break
      if not alreadyThere:
        run_samples['quadpts'].append(npt)
        run_samples['weights'][npt]=nwt
        print 'SG new entry:',npt,nwt
      else:
        #print '...duplicate point',npt,'- combining weights.'
        run_samples['weights'][npt]+=nwt
    #exit()
    self.sampler = spr.StochasticPoly(self.varDict,
                                      run_samples)
    self.numquadpts = len(run_samples['quadpts'])
    print '...%i quadrature points used...' %self.numquadpts

  def runSample(self,runDict):
    self.ie.runSolve(runDict['inp_file'])
    soln = self.ie.storeOutput(runDict['outFileName'])
    self.outq.put([runDict['nRun'],runDict['quadWts'],soln])
    #print 'Finished',str(runDict['nRun'])+':',soln

class SCExec(SC):
  def setCase(self):
    case =''
    case+= self.templateFile.split('.')[0]+'_SC_'
    case+= str(len(self.varDict.keys()))+'var'
    self.case=case
    return

  def ROM(self,xs):
    #TODO this is a strange place for this, but idk where to do
    #gather all the pts at which a variable has been evaluated
    varvals = []
    for i in self.histories['varVals'][0]:
      varvals.append([])
    for run in self.histories['varVals']:
      for v,val in enumerate(run):
        if val not in varvals[v]:
          varvals[v].append(val)
    for e,entry in enumerate(varvals):
      print 'Y_%i quad pts:' %e,entry

    tot=0
    for j in range(self.numquadpts):
      evalpts = self.histories['varVals'][j]
      print 'k: %i, Y^(k):' %j,evalpts
      prod=1
      for v,var in enumerate(self.histories['vars']):
        print 'n: %1i, varname:' %(v+1),var.name
        polyeval = var.lagrange(evalpts[v],xs[v],varvals[v])
        prod*=polyeval
        print '  L_k%1i:' %v,polyeval
      print 'barL:',prod
      sln = self.histories['soln'][j]
      print 'u:',sln
      print 'u*barL:',sln*prod
      tot+=sln*prod
      print ''
    return tot




class PCESCExec(SC):
  pass





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
        self.backends['PDF']=backend

  def parallelRun(self):
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
    if not self.restart:
      self.total_runs = -1
      trialsLeft = trials
      trialsAtRestart = 0
      #tempOutFile = file('solns.out','w')
    else: #start from restart
      print '\nStarting from restart...'
      print '  Starting after run',self.total_runs
      trialsAtRestart = self.total_runs
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
              totDone = finished + trialsAtRestart
              if self.restart:
                finished -= trialsAtRestart
              elapTime = time.time()-starttime
              dpdt = float(finished)/float(elapTime)
              toGo = dt.timedelta(seconds=(int(trialsLeft/dpdt)))
              elapTime = dt.timedelta(seconds=int(elapTime))
              print '%12i | %12s | %12s | %9i' %(totDone,elapTime,toGo,thrown),
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
