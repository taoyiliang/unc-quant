import numpy as np
import numpy.random as rd
import scipy.stats as stats
import time
from GetPot import GetPot
from itertools import product as allcombos
from sys import *
import InputEditor as ie
import IndexSets

def convertTime(timetogo):
  try: hrs = int(timetogo/3600)
  except ZeroDivisionError: hrs=0
  try: mnt = int((timetogo%(hrs*3600))/60)
  except ZeroDivisionError: mnt = int(timetogo/60)
  try: sec = int((timetogo%(hrs*3600+mnt*60)))
  except ZeroDivisionError: sec = int(timetogo)
  return hrs,mnt,sec

class Sampler:
  def __init__(self,varDict):
    self.varDict=varDict
    self.counter=0
    self.converged=False

  def giveSample(self):
    raise TypeError (str(self)+' is missing a giveSample method...')

class MonteCarlo(Sampler):
  def __init__(self,varDict,input_file):
    Sampler.__init__(self,varDict)
    self.type = 'MonteCarlo sampler'
    self.varDict=varDict
    self.varlist = varDict.values()
    self.totalSamples = input_file('Sampler/MC/totalSamples','')
    self.totalSamples = float(self.totalSamples)
    #self.convVar = input_file('Sampler/MC/convergeVar','1e0')
    #self.convAvg = input_file('Sampler/MC/convergeAvg','1e0')
    #self.convVar = float(self.convVar)
    #self.convAvg = float(self.convAvg)
    #self.batchSize = input_file('Sampler/MC/batch',1)
    self.print_freq = input_file('Sampler/MC/printfreq',1)
    self.truncDevs = input_file('Sampler/MC/truncate',0)
    self.oldAvg = 0
    self.oldVar = 0
    #self.lastone = False
    #self.sincebatch = 0
    self.sincePrint = 0
    self.tempFile = file('temp.out','w')

  def __del__(self):
    self.tempFile.close()

  def restart(self,N):
    self.counter=N
    self.starttime=time.time()

  def giveSample(self,moms=None):
    self.counter+=1
    #if self.counter==1:
    #  self.starttime=time.time()
    #  self.curTime = self.starttime
    if self.counter >= min(self.totalSamples,1e14):
      self.converged = True
    #update statistics
    #if (self.counter >= self.batchSize and moms != None) or self.lastone:
    #  if self.sincebatch < 2:
    #    self.sincebatch = self.batchSize
    #self.sincePrint+=1
    #if ((self.sincePrint >= self.print_freq) and moms!=None) or self.converged:
    #  self.sincePrint = 0
    #  N = moms[0]
    #  newM1 = moms[1]/N
    #  newM2 = moms[2]/N
    #  newAvg = newM1
    #  newVar = (newM2-newM1*newM1)#/(N-1) pop var
    #  checkVar = abs(newVar - self.oldVar)/newVar
    #  checkAvg = abs(newAvg - self.oldAvg)/newAvg
    #  self.oldVar = newVar
    #  self.oldAvg = newAvg

      #self.lastpct = self.pct
    #  pct = self.counter/float(self.totalSamples)
    #  toGopct = 1.0-pct
      #estimate completion
      #self.lastTime=self.curTime
    #  curTime = time.time()
      #timeforbatch = self.curTime - self.lastTime
      #dpct = (self.pct - self.lastpct)/timeforbatch
    #  elapsed = curTime - self.starttime
    #  dpct = pct/elapsed #pct per sec
    #  toGoT = toGopct/dpct
    #  hrs,mnt,sec = convertTime(toGoT)
    #  elapsed = convertTime(elapsed)

    #  pct = float(100*pct)
    #  print '\nMC run: %2.2f%%' %pct
    #  print '  Run %i/%1.0e' %(self.counter,self.totalSamples)
    #  print '  Average : %1.4e' %self.oldAvg
    #  print '  Moment2 : %1.4e' %newM2
    #  print '  Variance: %1.4e' %self.oldVar
    #  print '    Converge average : %1.2e' %checkAvg
    #  print '    Converge variance: %1.2e' %checkVar
    #  print '  time elapsed (hms): %2i:%2i:%2i ' %(elapsed)
    #  print '  est. time remain  : %2i:%2i:%2i ' %(hrs,mnt,sec)
      #if checkVar < self.convVar and checkAvg < self.convAvg:
      #  if self.lastone: self.converged = True
      #  else: self.lastone = True
      #else:
      #  self.converged = False
      #  self.lastone = False
    vals=[]
    for v,var in enumerate(self.varlist):
      vals.append(var.sample())
    #self.tempFile.writelines(str(runDict['varVals'])[1:-1]+'\n')
    #self.tempFile.flush()
    return vals

class StochasticPoly(Sampler):
  def __init__(self,varDict,input_file):
    Sampler.__init__(self,varDict)
    self.type = 'StochasticPoly sampler'
    for var in varDict.values():
      var.setQuadrature(input_file)
    self.varlist = varDict.values()
    #choose the index set
    #FIXME this has nothing to do with sampling, this is backend!
    # ask prinja about this!
    indexSet=input_file('Sampler/SC/indexSet','dud')
    if indexSet=='dud':
      print 'Index set not specified; using tensor product.'
      indexSet='TP'
    maxorder=input_file('Sampler/SC/maxorder',0)
    self.pickGaussPoints(indexSet,maxorder)

  def pickGaussPoints(self,iset,maxorder):
    #TODO change this to use a particular indexing system
    varlist = self.varlist
    orderlist = []
    for var in varlist:
      orderlist.append(var.quadOrds)
    #get the desired index set
    if maxorder==0:
      maxorder=np.max(np.max(orderlist))
      print 'Max order for index set not specified;',
      print 'using max from vars:',maxorder
    self.runords=IndexSets.chooseSet(orderlist,iset,maxorder)
    print '...size of index set:',len(self.runords),'...'
    #self.runords = list(allcombos(*orderlist))

  def giveSample(self):
    ords = self.runords[self.counter]
    runDict={}
    runDict['varVals']=[]
    runDict['quadWts']=[]
    runDict['runOrds']=ords
    for v,var in enumerate(self.varlist):
      #print 'Trying to access:'
      #print '  var num :',v
      #print '  var name:',var.name
      #print '  counter :',self.counter
      #print '  run[ctr]:',self.runords[self.counter]
      #print '  point   :',len(var.pts)
      #std_pt = var.pts[self.runords[self.counter][v]]
      std_pt,wt = var.quaddict[ords[v]]
      act_pt = var.convertToActual(std_pt)
      runDict['varVals'].append(act_pt)
      runDict['quadWts'].append(wt)
    self.counter+=1
    if self.counter>=len(self.runords):
      self.converged = True
    return runDict
