import numpy as np
import scipy.stats as stt
import time
from GetPot import GetPot
from sys import *
import os

import InputEditor as ie

import matplotlib.pyplot as plt
import matplotlib.mlab as ml
plt.ion()

def doplots(x,y,datapts,kx,ky,show=False):
  plt.close(1)
  plt.figure(1)
  datapts=np.array(datapts)
#create dictionary to get back soln values

#x=datapts[:,0]
#y=datapts[:,1]
  x=np.array(x)
  y=np.array(y)
  z=np.array(datapts)

  xi=np.linspace(np.min(x),np.max(x),100)
  yi=np.linspace(np.min(y),np.max(y),100)
  zi=ml.griddata(x,y,z,xi,yi)

  plt.contourf(xi,yi,zi,15)#,linewidth=0.0,colors='k')
  plt.pcolormesh(xi,yi,zi,cmap=plt.get_cmap('rainbow'))

  plt.colorbar()
  plt.scatter(x,y,marker='o', c='b',s=1,zorder=10)
  plt.xlim(np.min(x),np.max(x))
  plt.ylim(np.min(y),np.max(y))

  plt.xlabel(':'.join(kx.split('/')))
  plt.ylabel(':'.join(ky.split('/')))
  plt.title('Total Runs:'+str(totruns))
  if show: plt.show()
  else: plt.draw()

class Variable:
  def __init__(self,name,path):
    self.name=name
    self.path=path
    print 'Added variable',self.path

  def setDist(self,distName,args):
    if distName == 'uniform':
      if len(args)==1: #only given mean
        self.mean = float(args[0])
        self.range = 0.3*self.mean
        args=[self.mean-self.range,2.*self.range]
      elif len(args)==2:
        self.mean = float(args[0])
        if float(args[1])<0:
          self.mean *= 0.5
          self.range = self.mean
        else:
          self.range = float(args[1])*self.mean
        args=[self.mean-self.range,2.*self.range]
    expr='self.dist=stt.'+distName+'('+str(args)[1:-1]+')'
    exec expr

  def sample(self,n=1):
    if n==1:
      return self.dist.rvs()
    else:
      return self.dist.rvs(n)

  def sampleLow(self):
    return self.mean - self.range
  def sampleHi(self):
    return self.mean + self.range


starttime=time.time()

# load uncertainty input file
cl = GetPot(argv)
if cl.search('-i'):
  inpFileName=''
  inpFileName=cl.next(inpFileName)
  input_file=GetPot(Filename=inpFileName)
else: raise IOError('Require an input file using -i.')
print 'Loading uncertainty input:',inpFileName,'...'

uVars = input_file('Variables/names','').split(' ')
totalSamples   = input_file('Problem/totalSamples',0)
samplesPerIter = input_file('Problem/samplesPerIter',0)
targetVariance = input_file('Problem/targetVariance',0.0)
templateFile   = input_file('Problem/templateFile',' ')
storeFlux      = input_file('Problem/storeFlux',0)

print 'Using template input:',templateFile,'...'

#clear old flux files
if storeFlux:
  print 'Clearing old flux files...'
  os.system('rm -f store_flux/flux*.out')

# set up dictionary for random variables
varDict={}
for var in uVars:
  path=input_file('Variables/'+var+'/path',' ')
  dist=input_file('Variables/'+var+'/dist',' ')
  args=input_file('Variables/'+var+'/args',' ').split(' ')
  for a,arg in enumerate(args):
    args[a]=float(arg)
  varDict[var]=Variable(var,path)
  varDict[var].setDist(dist,args)

converged = False
runDict={}
datapts=[]
totruns=0
curTime=time.time()
x=[]
y=[]
print 'Running',totalSamples,'histories...'

#sample extremum first
runFileName,cx,cy,kx,ky = ie.writeInput(templateFile,runDict,totruns)

runDict[varDict[varDict.keys()[1]].path]=varDict[varDict.keys()[0]].sampleLow()
runDict[varDict[varDict.keys()[0]].path]=varDict[varDict.keys()[0]].sampleLow()

while not converged:
  # build sample values
  #if totruns < 4: #do an extreme case
  #  if totruns == 0: #low,low
  #    runDict[varDict[varDict.keys()[0]].path]=varDict[varDict.keys()[0]].sampleLow()
  #    runDict[varDict[varDict.keys()[1]].path]=varDict[varDict.keys()[1]].sampleLow()
  #  elif totruns == 1:
  #    runDict[varDict[varDict.keys()[0]].path]=varDict[varDict.keys()[0]].sampleHi()
  #    runDict[varDict[varDict.keys()[1]].path]=varDict[varDict.keys()[1]].sampleLow()
  #  elif totruns == 2:
  #    runDict[varDict[varDict.keys()[0]].path]=varDict[varDict.keys()[0]].sampleLow()
  #    runDict[varDict[varDict.keys()[1]].path]=varDict[varDict.keys()[1]].sampleHi()
  #  elif totruns == 3:
  #    runDict[varDict[varDict.keys()[0]].path]=varDict[varDict.keys()[0]].sampleHi()
  #    runDict[varDict[varDict.keys()[1]].path]=varDict[varDict.keys()[1]].sampleHi()
  #else:
  for key in varDict.keys():
    runDict[varDict[key].path]=varDict[key].sample()
  # run samples
  runFileName,cx,cy,kx,ky = ie.writeInput(templateFile,runDict,totruns)
  exitst = os.system('./TwoDProblem -i '+runFileName+' > /dev/null')
  #exitst = os.system('./TwoDProblem -i '+runFileName+' | grep "Final keff"')
  if exitst != 0:
    print 'TwoDProblem run failed with exit status',exitst
    exit()
  os.system('rm '+runFileName)

  #store flux
  if storeFlux:
    curPath = os.getcwd()
    if not os.path.isdir(curPath + '/store_flux'):
      os.system('mkdir store_flux')
    os.system('mv flux.out store_flux/')

  #collect k value
  newdata = ie.storeOutput()
  if newdata < 0:
    print 'WARNING: negative k found!'
    print ' ',kx,'=',cx,' |',ky,'=',cy,' |  k =',newdata
    print '  Excluding data point.'
    continue
  x.append(cx)
  y.append(cy)
  datapts+=[newdata]
  totruns+=1
  if totruns % (totalSamples/10.) == 0:
    print 'Completed',totruns,'samples,',totruns/(time.time()-starttime),'samples/sec'
    curTime=time.time()
    doplots(x,y,datapts,kx,ky)
  if totruns>=totalSamples:
    converged = True

print 'Total time:',curTime-starttime
print 'Avg run time:',(curTime-starttime)/float(totruns)
plt.ioff()
doplots(x,y,datapts,kx,ky,show=True)
