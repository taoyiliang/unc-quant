import numpy as np
import scipy.stats as stt
import time
import os
from GetPot import GetPot
from sys import *
import InputEditor
import Variables
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
#maxSamples = input_file('Problem/totalSamples',0)
storeFlux  = input_file('Problem/storeFlux',0)
execDir    = input_file('Problem/execDir'  ,'')
inputDir   = input_file('Problem/inputDir' ,'')
input_editor = input_file('Problem/input_editor','')

torun = 'ie = InputEditor.'+input_editor+'()'
exec torun

#set directories
uncDir  = os.getcwd()
os.chdir(execDir)
execDir = os.getcwd()
os.chdir(uncDir)
#print 'in dir:',inputDir
os.chdir(inputDir)
inputDir = os.getcwd()
os.chdir(uncDir)

#load variables
uVars = input_file('Variables/names','').split(' ')
varDict={}
for var in uVars:
  path=input_file('Variables/'+var+'/path',' ')
  dist=input_file('Variables/'+var+'/dist',' ')
  args=input_file('Variables/'+var+'/args',' ').split(' ')
  for a,arg in enumerate(args):
    args[a]=float(arg)
  varDict[var]=Variables.newVar(dist,var,path)
  varDict[var].setDist(args)

#build sampler
samplerType = input_file('Sampler/active','')
if samplerType in ['mc','MC','mC','MonteCarlo','Monte Carlo']:
  sampler = spr.MonteCarlo(varDict,input_file)
elif samplerType in ['SC','sc','StochPoly']:
  sampler = spr.StochasticPoly(varDict,input_file)
else:
  raise IOError ('Sampler type not recognized: |'+samplerType+'|')

#for var in varDict.values():
#  print var.name,var.pts,var.wts
#print sampler.runords

#"converged" is a property of the sampler type.

# clear old unc files, if any
os.chdir(inputDir)
os.system('rm '+templateFile+'.unc*')

# run samples
total_runs = -1
curTime=time.time()
histories={}
histories['varNames']=varDict.keys()
histories['vars']=varDict.values()
histories['varPaths']=[]
for var in varDict.values():
  histories['varPaths'].append(var.path)
histories['varVals']=[]
histories['nRun']=[]
histories['soln']=[]
while not sampler.converged:# and total_runs<maxSamples:
  total_runs+=1
  print 'Starting run',total_runs
  runDict = sampler.giveSample()
  #if total_runs==0:
  #  histories['varPaths']=runDict['varVals']
  histories['nRun'].append(total_runs)
  for key in runDict.keys():
    try: histories[key].append(runDict[key])
    except KeyError:
      histories[key]=[runDict[key]]
  #histories['varVals'].append(runDict['varVals'])
  #print '  histories:'
  #for key in histories.keys():
  #  print '    ',key,histories[key]

  #get input file name, write input file
  inp_file = ie.writeInput(templateFile,
                           inputDir,
                           histories['varPaths'],
                           runDict['varVals'],
                           total_runs)

  #run file
  ie.runSolve(inp_file)
  #osstat = os.system('./TwoDProblem -i '+inp_file+' > /dev/null')
  #if osstat != 0:
  #  print 'Run attempt failed with error code',osstat
  #  exit()
  #os.system('rm '+inp_file)
  histories['soln'].append(ie.storeOutput())
  #print 'solution:',histories['soln'][-1]

#build backend
def cprToActual(y,x):
  return x*x + x*y + y*y
  #return x+y
#def cprToActual(x):
#  return 1.+2.*x

backendTypes = input_file('Backend/active','').split(' ')
for beType in backendTypes:
  if beType in ['plot2D']:
    backend = be.Plotter2D()
    #cprvarsX=input_file('Output/Plot2D/xvars','').split(' ')
    #cprvarsY=input_file('Output/Plot2D/yvars','').split(' ')
    #for i in range(len(xvars)):
    backend.compareVars(histories)
    plt.show()

  elif beType in ['ROM','rom']:
    backend = be.ROM(sampler.runords)
    backend.makeROM(histories)
    samples=input_file('Backend/ROM/samples',100.)
    samples=int(samples)
    checkSoln=input_file('Backend/ROM/checkSolution',0)
    makePDF=input_file('Backend/ROM/createPDF',0)
    numbins=input_file('Backend/ROM/bins',100)

    MCROMs=None
    #sample ROM
    if checkSoln:
      print 'sampling ROM...'
      MCROMs,correct = backend.MCsampleDict(histories,trials=samples,cprFunc=cprToActual)
      print 'Checking convergence to actual soln...'
      ok = np.product(np.abs(np.array(MCROMs.values())-np.array(correct.values()))<1e-9)==1
      numbad=0
      for i,key in enumerate(MCROMs.keys()):
        if not np.abs(MCROMs[key]-correct[key])<1e-9:
          numbad+=1
          print key,(MCROMs[key]-correct[key])
      if ok: print '  ...Converged!'
      else:
        print  '  ...ERROR: not converged!',numbad,'errors found'

  #pdf rom
    if makePDF:
      print 'constructing discrete pdf'
      if MCROMs == None:
        bins,ctrs=backend.makePDF(histories,numbins,samples)
      else:
        bins,ctrs=backend.makePDF(histories,numbins,samples,MCROMs.values())
      plt.plot(ctrs,bins)
      plt.title('PDF of solution, 1e'+str(int(np.log10(samples)))+' samples')
      plt.xlabel('Solutions')
      plt.ylabel('Probability')
      plt.show()


  elif beType in ['CSV','csv']:
    backend = be.CSVMaker()
    outFileName = input_file('Backend/outFile','dud')
    if outFileName == 'dud': outFileName='out.csv'
    backend.makeCSV(histories,outFileName)

print 'Executioner complete.'
