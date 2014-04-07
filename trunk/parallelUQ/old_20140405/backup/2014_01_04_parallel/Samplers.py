import numpy as np
import scipy.stats as stats
import time
from GetPot import GetPot
from itertools import product as allcombos
from sys import *


import InputEditor as ie

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
    self.varDict=varDict
    self.totalSamples   = input_file('Sampler/MC/totalSamples',10)
    self.samplesPerIter = input_file('Sampler/MC/samplesPerIter',1)
    self.targetVariance = input_file('Sampler/MC/targetVar',1.0)

  def giveSample(self):
    runDict={}
    runDict['varVals']={}
    for key in self.varDict.keys():
      runDict['varVals'][self.varDict[key].path]=self.varDict[key].dist.rvs()
      #print self.counter,'Sampled',key,'as',runDict['varVals'][self.varDict[key].path]
    self.counter+=1
    return runDict

class StochasticPoly(Sampler):
  def __init__(self,varDict,input_file):
    Sampler.__init__(self,varDict)
    for var in varDict.values():
      var.setQuadrature(input_file)
    self.varlist = varDict.values()
    self.pickGaussPoints()

  def pickGaussPoints(self):
    varlist = self.varlist
    orderlist = []
    for var in varlist:
      orderlist.append(range(var.order))
    #tensor product of runs - not the best solve
    self.runords = list(allcombos(*orderlist))

  def giveSample(self):
    runDict={}
    #runDict['varvals']=[]
    runDict['varVals']=[]
    runDict['quadWts']=[]
    for v,var in enumerate(self.varlist):
      std_pt = var.pts[self.runords[self.counter][v]]
      act_pt = var.convertToActual(std_pt)
      #print '  std | act:',std_pt,act_pt
      #old runDict['varVals'].append(std_pt)
      runDict['varVals'].append(act_pt)
      runDict['quadWts'].append(var.wts[self.runords[self.counter][v]])
    self.counter+=1
    if self.counter>=len(self.runords):
      self.converged = True
    #print 'sampling',list((runDict['varVals'][i] for i in range(len(runDict['varVals']))))
    return runDict;
    #TODO
