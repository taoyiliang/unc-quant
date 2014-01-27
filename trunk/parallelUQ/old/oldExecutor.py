import numpy as np
import scipy.stats as stt
import time
import os
from GetPot import GetPot
from sys import *
import InputEditor as ie
from Variables import Variable
import Samplers as spr
import Backends as be
import matplotlib.pyplot as plt

starttime=time.time()

#load uncertainty input
cl = GetPot(argv)
if cl.search('-i'):
  inpFileName=''
  inpFileName=cl.next(inpFileName)
  input_file=GetPot(Filename=inpFileName)
else: raise IOError('Requires an input file using -i.')

#load problem-specific information
templateFile=input_file('Problem/templateFile','')

#load variables
uVars = input_file('Variables/names','').split(' ')
varDict={}
for var in uVars:
  path=input_file('Variables/'+var+'/path',' ')
  dist=input_file('Variables/'+var+'/dist',' ')
  args=input_file('Variables/'+var+'/args',' ').split(' ')
  for a,arg in enumerate(args):
    args[a]=float(arg)
  varDict[var]=Variable(var,path)
  varDict[var].setDist(dist,args)

#build sampler
samplerType = input_file('Sampler/active','')
if samplerType in ['mc','MC','mC','MonteCarlo','Monte Carlo']:
  sampler = spr.MonteCarlo(varDict,input_file)
elif samplerType in ['SC','sc','StochPoly']:
  sampler = spr.StochasticPoly(varDict,input_file)
else:
  raise IOError ('Sampler type not recognized: |'+samplerType+'|')


#"converged" is a property of the sampler type.

# clear old unc files, if any
os.system('rm '+templateFile+'.unc*')

# run samples
total_runs = -1
curTime=time.time()
histories={}
histories['varVals']=[]
histories['nRun']=[]
histories['soln']=[]
histories['vars']=varDict.keys()
while not sampler.converged and total_runs<20:
  total_runs+=1
  print 'run:',total_runs
  runDict = sampler.giveSample()
  histories['nRun'].append(total_runs)
  for key in runDict.keys():
    try: histories[key].append(runDict[key])
    except IndexError:
      histories[key]=[runDict[key]]
  histories['varVals'].append(runDict['varVals'])

  #get input file name, write input file
  inp_file = ie.writeInput(templateFile,runDict['varVals'],total_runs)

  #run file
  os.system('./TwoDProblem -i '+inp_file+' > /dev/null')
  os.system('rm '+inp_file)
  histories['soln'].append(ie.storeOutput())

#build backend
backendTypes = input_file('Backend/active','').split(' ')
for beType in backendTypes:
  if beType in ['plot2D']:
    backend = be.Plotter2D()
    cprvarsX=input_file('Output/Plot2D/xvars','').split(' ')
    cprvarsY=input_file('Output/Plot2D/yvars','').split(' ')
    for i in range(len(xvars)):
      backend.compareVars(cprvarsX[i],cprvarsY[i],histories)
    plt.show()
    #TODO
  elif beType in ['ROM','rom']:
    pass
  elif backendType in ['CSV','csv']:
    backend = be.CSVMaker()
    outFileName = input_file('Backend/outFile','dud')
    if outFileName == 'dud': outFileName='out.csv'
    backend.makeCSV(histories,outFileName)
