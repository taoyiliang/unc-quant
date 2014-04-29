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
  def __init__(self,varDict):#,input_file):
    Sampler.__init__(self,varDict)
    self.type = 'MonteCarlo sampler'
    self.varDict=varDict
    self.varlist = varDict.values()
    #either tell me how many samples or target error
    #self.totalSamples = input_file('Sampler/MC/totalSamples',0)
    #self.targetError = input_file('Sampler/MC/targetError',1.0)
    #if self.totalSamples==0 and self.targetError==1.0:
    #  raise IOError('Neither totalSamples nor targetError were '+\
    #      'set in Sampler/MC.  Please set one.')

  def giveSample(self):
    if self.converged:
      return 'END'
    runDict={}
    runDict['varVals']=[]
    for v,var in enumerate(self.varlist):
      runDict['varVals'].append(var.sample())
    return runDict

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
    except IndexError:
      return 'END'
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
