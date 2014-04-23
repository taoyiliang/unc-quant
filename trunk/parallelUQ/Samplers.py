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
    self.counter=-1
    self.converged=False

  def __iter__(self):
    return self

  def next(self):
    self.counter+=1
    res = self.giveSample()
    if res == 'END':
      self.converged=True
      raise StopIteration
    else:
      return self.giveSample()

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
    self.print_freq = input_file('Sampler/MC/printfreq',1)
    self.truncDevs = input_file('Sampler/MC/truncate',0)
    self.oldAvg = 0
    self.oldVar = 0
    self.sincePrint = 0
    self.tempFile = file('temp.out','w')


  def giveSample(self,moms=None):
    if self.counter >= min(self.totalSamples,1e14):
      self.converged = True
    vals=[]
    for v,var in enumerate(self.varlist):
      vals.append(var.sample())
    return vals

class StochasticPoly(Sampler):
  def __init__(self,varDict,run_set):
    Sampler.__init__(self,varDict)
    self.type = 'StochasticPoly sampler'
    self.varlist = varDict.values()
    self.run_set = run_set
    #print 'run set vals:',run_set['quadpts']

  def giveSample(self):
    try:
      quadpts = self.run_set['quadpts'][self.counter]
    except IndexError: return 'END'
    weight = self.run_set['weights'][quadpts]
    runDict={}
    runDict['varVals']=[]
    runDict['quadWts']=[]
    #runDict['runOrds']=ords
    for v,var in enumerate(self.varlist):
      #std_pt = quadpts[v]
      #wt = weights[v]
      # TODO Integrated now act_pt = var.convertToActual(std_pt)
      runDict['varVals'].append(quadpts[v])
    runDict['quadWts']=weight
    #if self.counter>=len(self.runords):
    #  self.converged = True
    return runDict
