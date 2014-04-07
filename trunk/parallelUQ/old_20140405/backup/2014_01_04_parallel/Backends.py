import matplotlib.pyplot as plt
import matplotlib.mlab as mlab
import numpy as np
from bisect import bisect_left
from sys import *
import os

import multiprocessing
from multiprocessing.queues import SimpleQueue as sque

class Backend:
  def __init__(self):
    pass

class Plotter2D(Backend):
  def __init__(self):
    #super().__init__()
    pass

  def compareVars(self,histories):
    #find x,y in histories
    nameX=histories['vars'][0]
    nameY=histories['vars'][1]
    pathX=histories['varPaths'][0]
    pathY=histories['varPaths'][1]
    xs=np.zeros(histories['nRun'][-1]+1)
    ys=np.zeros(histories['nRun'][-1]+1)
    for e,edict in enumerate(histories['varVals']):
      xs[e]=edict[pathX]
      ys[e]=edict[pathY]
    zs=np.array(histories['soln'])

    xi=np.linspace(np.min(xs),np.max(xs),100)
    yi=np.linspace(np.min(ys),np.max(ys),100)

    plt.figure()
    zi=mlab.griddata(xs,ys,zs,xi,yi)
    plt.contour(xi,yi,zi,15,linewidth=0.5,colors='k')
    plt.pcolormesh(xi,yi,zi,cmap=plt.get_cmap('rainbow'))

    plt.colorbar()
    plt.scatter(xs,ys,marker='o',c='b',s=1,zorder=10)
    plt.xlim(np.min(xs),np.max(xs))
    plt.ylim(np.min(ys),np.max(ys))

    plt.xlabel(nameX)
    plt.ylabel(nameY)
    plt.title('Total Runs: '+str(len(xs)))

class CSVMaker(Backend):
  def __init__(self):
    #super().__init__()
    pass

  def makeCSV(self,histories,outFileName):
    outFile=file(outFileName,'w')
    #place header
    header='Run Number'
    for val in histories['vars']:
      header+=','+val
    header+=',soln\n'
    outFile.writelines(header)
    #place individual runs
    avVals=np.zeros(len(histories['vars']))
    avSoln=0
    for i in range(len(histories['nRun'])):
      line=str(histories['nRun'][i])+','
      for v in range(len(histories['vars'])):
        avVals[v]+=histories['varVals'][i][v]
        line+=str(histories['varVals'][i][v])+','
      avSoln+=histories['soln'][i]
      line+=str(histories['soln'][i])+'\n'
      outFile.writelines(line)
    #place averages
    avVals/=len(histories['nRun'])
    avSoln/=len(histories['nRun'])
    avg='avg,'
    for v in range(len(histories['vars'])):
      avg+=str(avVals[v])+','
    avg+=str(avSoln)+'\n'
    outFile.writelines(avg)
    outFile.close()
    os.system('grep Run '+outFileName)
    os.system('grep avg '+outFileName)

class ROM(Backend):
  def __init__(self,ords):
    self.coeffs={}
    self.ords = ords
    self.outq=sque()

  def makeROM(self,histories,verbose=False):
    print 'Building ROM'
    c={}
    #for key in histories.keys():
    #  print '  ',key,histories[key]
    #print 'in histories: varvals'
    #print '  ',histories['varVals']
    #TODO need to store vars for poly evals later
    for i,ordset in enumerate(self.ords): #poly orders
      c[tuple(ordset)]=0
      for m,soln in enumerate(histories['soln']): #histories (quadpts)
        #valList = histories['varVals'][m]
        temp_order = 1
        for v,var in enumerate(histories['vars']): #variables
          temp_var = histories['quadWts'][m][v]
          std_val = var.convertToStandard(histories['varVals'][m][v])
          temp_var *=var.evalNormPoly(std_val,ordset[v])
          #temp_var *=var.intvlShiftCoeff()
          if verbose:
            print 'var',var.name,'order',ordset[v]
            print '  varval:',histories['varVals'][m][v]
            print '  polord:',ordset[v]
            print '  stdval:',std_val
            print '  weight:',histories['quadWts'][m][v]
            print '  poly:',var.evalNormPoly(std_val,ordset[v])
            print '  soln:',soln
            #print '  eval:',temp_var
          temp_order*=temp_var
        c[tuple(ordset)]+=temp_order*soln
        #FIXME this is hacky, but eliminates some roundoff errors
        if abs(c[tuple(ordset)])<1e-12:
          c[tuple(ordset)]=0
    #remove zero coefficients
    noelim=True
    for key in c.keys():
      if abs(c[key])<1e-12:
        print '  Coeff',(key),'is less than 1e-12 and was removed'
        del c[key]
        noelim*=False
    if noelim: print '  ...No coefficients were removed.'
    orderedKeys = c.keys()[:]
    orderedKeys.sort()
    varnames=''
    for var in histories['vars']:
      varnames+=var.name+' '
    print 'Coefficient solutions: ( ',varnames,')'
    for key in orderedKeys:
      print '  ',key,c[key]
    self.coeffs=c

  def sampleROM(self,vals,histories):
    tot = 0
    for c in self.coeffs.keys():
      coef = self.coeffs[c]
      temp = 1
      for v,var in enumerate(histories['vars']):
        std_val = var.convertToStandard(vals[v])
        temp *= var.evalNormPoly(std_val,c[v])
      tot+=coef*temp
    return tot

  def MCsampleDict(self,histories,trials=1000,cprFunc=None):
    solns={}
    actual={}
    trials = int(trials)
    #TODO this can be vectored? 
    for i in range(trials):
      vals=[]
      actvals=[]
      if trials%(0.1*trials)==0:
        print '  completed',i,'/',trials,'trials...'
      for var in histories['vars']:
        #actvals.append(var.sample())
        #vals.append(var.convertToStandard(actvals[-1]))
        vals.append(var.sample())
      solns[tuple(vals)]=self.sampleROM(vals,histories)
      if cprFunc!=None:
        actual[tuple(vals)]=cprFunc(*vals)
        #WARNING: this depends strongly on what order vals is stored!
    if cprFunc==None:
      return solns
    else:
      return solns,actual

  def MCsample(self,histories,trials=1000):
    trials = int(trials)
    solns=np.zeros(trials)
    #TODO this can be vectored? 
    for i in range(trials):
      if trials%(0.1*trials)==0:
        print '  completed',i,'/',trials,'trials...'
      vals=[]
      for var in histories['vars']:
        #actvals.append(var.sample())
        #vals.append(var.convertToStandard(actvals[-1]))
        vals.append(var.sample())
      solns[i]=self.sampleROM(vals,histories)
    return solns

  def MCsampleChild(self,num,histories,trials=1000):
    #print 'MCchild,',trials,'trials'
    np.random.seed()
    trials = int(trials)
    solns=[]
    for i in range(trials):
      vals=[]
      for var in histories['vars']:
        vals.append(var.sample())
      solns.append(self.sampleROM(vals,histories))
    self.outq.put(solns)
    #print 'Child',num,'entered',trials,'solns in queue'

  def MCsampleParallel(self,histories,trials=1000):
    self.numprocs = multiprocessing.cpu_count()
    trials = int(trials)
    #turns out can only handle up to 1400 trials per go
    trialsPerProc = int(float(trials)/float(self.numprocs))
    if trialsPerProc >1000:
      trialsPerProc=1000

    trialsLeft = trials
    self.ps=[]
    solns=[]
    numPs=0
    done=False
    thqu=half=qu = False
    while not done:
      #printout flags
      if trialsLeft <= 0.75*trials and thqu==False:
        thqu = True
        print 'Trials complete:',trials-trialsLeft,'/',trials
      elif trialsLeft <= 0.5*trials and half==False:
        half = True
        print 'Trials complete:',trials-trialsLeft,'/',trials
      elif trialsLeft <= 0.25*trials and qu==False:
        qu = True
        print 'Trials complete:',trials-trialsLeft,'/',trials
      #remove dead processes
      for p,proc in enumerate(self.ps):
        if not proc.is_alive():
          proc.join()
          while not self.outq.empty():
            solns+=list(self.outq.get())
          #gather solns
          del self.ps[p]
      #check to see if more processes need to be spawned
      if trialsLeft > 0:
        while len(self.ps)<self.numprocs:
          #add new process
          if trialsLeft > trialsPerProc:
            newtrials = trialsPerProc
          else:
            newtrials = trialsLeft
          trialsLeft -= newtrials
          numPs+=1
          self.ps.append(multiprocessing.Process(\
                            target=self.MCsampleChild,\
                            args=(numPs,histories,newtrials))\
                        )
          self.ps[-1].start()
      #to be finished, 
      #  - outq must be empty,
      #  - all processes must be dead
      #  - no trials must be left to run
      done=self.outq.empty()
      for p in self.ps:
        done*=not p.is_alive()
      done *= (trialsLeft == 0)

    print 'ran',numPs,'child processes on',self.numprocs,'processers'
    #for p in self.ps:
    #  p.start()
    #i=0
    #for p in self.ps:
    #  p.join()
    #  print 'Child',i,'Done'
    #  i+=1

    #for i,p in enumerate(self.ps):
    #  print '  check',i,'is done:',not p.is_alive()
    if len(solns)==trials:
      print 'All solutions successfully recorded.'
    else:
      print 'Only recovered',len(solns),'solutions out of',trials
    return solns

  def makePDF(self,histories,bins,samples=1e6,solns=None):
    #TODO better if they don't get stored in memory!
    if solns==None:
      solns=self.MCsampleParallel(histories,int(samples))
    low = min(solns)
    hi = max(solns)
    bounds = np.linspace(low,hi,bins+1)
    bins=np.zeros(int(bins))
    #make bin centers
    ctrs=np.zeros(len(bins))
    for c in range(len(ctrs)):
      ctrs[c]=0.5*(bounds[c]+bounds[c+1])
    #TODO vectorize this binning
    for s in solns:
      i=bisect_left(bounds,s)
      bins[i-1]+=1
    #normalize
    for b in range(len(bins)):
      bins[b]=float(bins[b])/float(len(solns))
    return bins,ctrs

  def OneDPrintPoly(self,histories):
    poly = np.poly1d(0)
    for c in self.coeffs.keys():
      coef = self.coeffs[c]
      #temp = 1
      #for v,var in enumerate(histories['vars']): #there's only one
      temp = histories['vars'][0].oneDPoly(c[0])
      poly+=np.poly1d(coef*temp)
    #fix zeros
    cofs=poly.c
    for c,cof in enumerate(cofs):
      if abs(cof)<1e-10: cofs[c]=0.0
    return np.poly1d(cofs)
