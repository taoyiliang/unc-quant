import matplotlib.pyplot as plt
import matplotlib.mlab as mlab
import numpy as np
import cPickle as pk
from bisect import bisect_left
from sys import *
import os
import time
import datetime

import multiprocessing
from multiprocessing.queues import SimpleQueue as sque

class Backend(object):
  '''
  Meta-class for post-processors
  '''
  def __init__(self):
    pass

class Plotter2D(Backend):
  '''
  Used for producing 2D temperature plot of soln as a func of 2 uncertain vars
  '''
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
  '''
  Used to write data to a file in CSV format
  '''
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

class Stats(Backend):
  '''
  Used to calculate mean, variance, usually for MC runs
  '''
  def __init__(self,low,hi,outFileName,nbins=100):
    self.type='PDF'
    self.nbins = nbins
    self.bounds = np.linspace(low,hi,self.nbins+1)
    self.outFileName = outFileName
    self.bins=np.zeros(int(self.nbins))
    self.ctrs=np.zeros(len(self.bins))
    self.mom1 = 0
    self.mom2 = 0
    self.tot1 = 0
    self.tot2 = 0
    self.N = 0
    self.below=0
    self.above=0
    self.range=[1e10,-1e10]
    for c in range(len(self.ctrs)):
      self.ctrs[c]=0.5*(self.bounds[c]+self.bounds[c+1])

  def addToBins(self,vals,toprint=False):
    for val in vals:
      self.tot1+=val
      self.tot2+=val*val
      self.N+=1
      self.range[0]=min(val,self.range[0])
      self.range[1]=max(val,self.range[1])
      if val <= self.bounds[0]:
        self.below+=1
      elif val > self.bounds[-1]:
        self.above+=1
      else:
        i=bisect_left(self.bounds,val)
        self.bins[i-1]+=1
    #if toprint:
    #  print 'Finished run',self.N
    #print 'at',datetime.datetime.fromtimestamp(time.time())

  def savedict(self):
    #toSave = {}
    #items = ['tot1','tot2','N','range','below','above','bins',\
    #        'bounds','ctrs']
    #for mbr in items:
    #  exec('toSave[mbr]=self.'+mbr)
    #return toSave
    return self.__dict__

  def loaddict(self,loaddict):
    self.__dict__.update(loaddict)
    print '\nRestarted MC:'
    print '  runs:',self.N
    print '  sum1:',self.tot1
    print '  sum2:',self.tot2
    return self.N

  def MCcalc(self,histories):
    self.mom1=self.tot1/self.N
    self.mom2=self.tot2/self.N
    var = self.mom2-self.mom1*self.mom1

    for v,val in enumerate(self.bins):
      norm = (self.N-self.above-self.below)*(self.bounds[v+1]-self.bounds[v])
      self.bins[v]=val/norm

    height = max(self.bins)
    line = np.linspace(0,height,3)

    plt.figure()
    plt.plot(self.ctrs,self.bins,'b-')
    plt.plot(np.ones(len(line))*self.mom1,line,'g-')
    #plt.plot(np.ones(len(line))*var,      line,'k-')
    plt.text(self.mom1,0.95*height,'Avg='+str(self.mom1))
    plt.text(self.mom1,0.9*height,'Var='+str(var))

    #write it out to file, too
    outFile = file(self.outFileName,'w')
    outFile.writelines('Centers,Bin Values,Bin Width\n')
    outFile.writelines('Average,'+str(self.mom1)+'\n')
    outFile.writelines('2nd Moment,'+str(self.mom2)+'\n')
    outFile.writelines('Variance,'+str(var)+'\n')
    for i,ctr in enumerate(self.ctrs):
      msg=','.join([str(ctr),str(self.bins[i])])
      outFile.writelines(','.join([str(ctr),str(self.bins[i]),])+'\n')
    outFile.close()
    print 'Data was saved to',self.outFileName

    print '==================================================================='
    print 'Statistics:'
    print 'Mean:',self.mom1
    print 'Mom2:',self.mom2
    print 'Var :',var
    print 'Trys:',self.N
    print ''
    print 'The extreme bounds given for bins was',[self.bounds[0],self.bounds[-1]]
    print 'The extreme values obtained to bin is',self.range
    print 'Number of samples excluded from bins (hi/low):',self.above,self.below
    print '==================================================================='

  def makePDF(self,histories,bins):
    #TODO better if they don't get stored in memory!
    solns = np.array(histories['soln'])
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

class ROM(Backend):
  '''
  Uses SC runs to generate Reduced Order Model polynomial representation of
    soln as a function of the uncertain variables
  '''
  def __init__(self,ords,case):
    self.expOrds = ords
    self.case = case
    self.coeffs={}
    self.outq=sque()

  def makeROM(self,histories,verbose=False):
    print '\nBuilding ROM...'
    c={}
    #for key in histories.keys():
    #  print '  ',key,histories[key]
    #print 'in histories: varvals'
    #print '  ',histories['varVals']
    #TODO need to store vars for poly evals later
    cofFile = file(self.case+'.cof','w')
    #TODO do each ordset in parallel!
    for i,ordset in enumerate(self.expOrds): #poly orders
      ords = tuple(ordset)
      c[ords]=0
      print 'Constructing coefficient',ords,'\r',
      for m,soln in enumerate(histories['soln']): #histories (quadpts)
        #valList = histories['varVals'][m]
        temp_order = 1
        for v,var in enumerate(histories['vars']): #variables
          temp_var = histories['quadWts'][m][v]
          std_val = var.convertToStandard(histories['varVals'][m][v])
          temp_var *=var.evalNormPoly(std_val,ordset[v])
          #temp_var *=var.probdens(std_val)
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
        c[ords]+=temp_order*soln
        #this is hacky, but eliminates some roundoff errors
        if abs(c[ords])<1e-12:
          c[ords]=0
      cofFile.writelines(str(ords)+' | %1.8e' %c[ords] +'\n')
    cofFile.close()
    print 'Ceofficients written to',os.getcwd()+'/'+self.case+'.cof'
    #remove zero coefficients
    print 'Eliminating near-zero coefficients...'
    noelim=True
    for key in c.keys():
      if abs(c[key])<1e-11:
        print '  Coeff',(key),'is less than 1e-11 and was removed'
        del c[key]
        noelim*=False
    if noelim: print '  ...No coefficients were removed.'
    orderedKeys = c.keys()[:]
    orderedKeys.sort()
    varnames=''
    for var in histories['vars']:
      varnames+=var.name+' '
    if len(c)<=16:
      print '\nCoefficient solutions: ( ',varnames,')'
      for key in orderedKeys:
        print '  ',key,c[key]
    self.coeffs=c

  def setStats(self,histories):
    outFile = file(self.case+'.stats','w')
    #get mean
    base=[0]*len(histories['vars'])
    self.mean=float(self.evalSingleTerm(base,histories))
    outFile.writelines('Mean,'+str(self.mean)+'\n')

    #get second moment
    mom2=0
    for cof in self.coeffs.values():
      mom2+=float(cof)**2
    for var in histories['vars']:
      mom2*=var.probWeight(0,scale='standard')
    self.var = mom2-self.mean**2
    outFile.writelines('Variance,'+str(self.var)+'\n')
    outFile.writelines('2nd Moment,'+str(mom2)+'\n')
    print 'stats | mean,mom2,var'
    print self.mean,mom2,self.var
    outFile.close()

  def getcoeff(self,base):
    return self.coeffs[tuple(base)]

  def evalSingleTerm(self,base,histories):
    base=tuple(base)
    cof=self.coeffs[base]
    temp=1
    for v,var in enumerate(histories['vars']):
      val = 1.0
      add = var.evalNormPoly(val,base[v])
      temp *= add
    temp*=cof
    return temp

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

  #def MCsampleDict(self,histories,trials=1000,cprFunc=None):
  #  solns={}
  #  actual={}
  #  trials = int(trials)
    #TODO this can be vectored? 
  #  for i in range(trials):
  #    vals=[]
  #    actvals=[]
      #if trials%(0.1*trials)==0:
      #  print '  completed',i,'/',trials,'trials...'
  #    for var in histories['vars']:
        #actvals.append(var.sample())
        #vals.append(var.convertToStandard(actvals[-1]))
  #      vals.append(var.sample())
  #    solns[tuple(vals)]=self.sampleROM(vals,histories)
  #    if cprFunc!=None:
  #      actual[tuple(vals)]=cprFunc(*vals)
        #WARNING: this depends strongly on what order vals is stored!
  #  if cprFunc==None:
  #    return solns
  #  else:
  #    return solns,actual

  #def MCsample(self,histories,trials=1000):
  #  trials = int(trials)
  #  solns=np.zeros(trials)
    #TODO this can be vectored? 
  #  for i in range(trials):
      #if trials%(0.1*trials)==0:
      #  print '  completed',i,'/',trials,'trials...'
  #    vals=[]
  #    for var in histories['vars']:
        #actvals.append(var.sample())
        #vals.append(var.convertToStandard(actvals[-1]))
  #      vals.append(var.sample())
  #    solns[i]=self.sampleROM(vals,histories)
  #  return solns

  def MCsampleChild(self,num,histories,trials=1000,trunc=0):
    #print 'MCchild,',trials,'trials'
    np.random.seed()
    trials = int(trials)
    solns=[]
    for i in range(trials):
      vals=[]
      for var in histories['vars']:
        vals.append(var.sample(trunc=trunc))
        #print 'sample',vals[-1]
      solns.append(self.sampleROM(vals,histories))
    self.outq.put(solns)
    #print 'Child',num,'entered',trials,'solns in queue'

  def MCsampleParallel(self,histories,trials=1000,trunc=0):
    self.numprocs = multiprocessing.cpu_count()
    print '  ...using',self.numprocs,'processors...'
    trials = int(trials)
    print '  ...sampling 1e'+str(int(np.log10(trials))),'trials...'
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
    curTime = time.time()
    pkFile = self.case+'.samples.pk'
    while not done:
      #printout flags
      if trialsLeft <= 0.75*trials and thqu==False:
        thqu = True
        elapsed=time.time()-curTime
        curTime=time.time()
        print 'Trials complete:',trials-trialsLeft,'/',trials,',',elapsed,'secs'
        pk.dump(solns,open(pkFile,'wb'))
      elif trialsLeft <= 0.5*trials and half==False:
        half = True
        elapsed=time.time()-curTime
        curTime=time.time()
        print 'Trials complete:',trials-trialsLeft,'/',trials,',',elapsed,'secs'
        pk.dump(solns,open(pkFile,'wb'))
      elif trialsLeft <= 0.25*trials and qu==False:
        qu = True
        elapsed=time.time()-curTime
        curTime=time.time()
        print 'Trials complete:',trials-trialsLeft,'/',trials,',',elapsed,'secs'
        pk.dump(solns,open(pkFile,'wb'))
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
                            args=(numPs,histories,newtrials,trunc))\
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

    pk.dump(solns,open(pkFile,'wb'))
    print 'ran',numPs,'child processes on',self.numprocs,'processers'
    if len(solns)==trials:
      print 'All solutions successfully recorded.'
    else:
      print 'Only recovered',len(solns),'solutions out of',trials
    return solns

  def makePDF(self,histories,bins,samples=1e6,solns=None,trunc=0):
    #TODO better if they don't get stored in memory!
    outFile = file(self.case+'.stats','a')
    if solns==None:
      solns=self.MCsampleParallel(histories,int(samples),trunc)
    low = min(solns)
    hi = max(solns)
    bounds = np.linspace(low,hi,bins+1)
    bins=np.zeros(int(bins))
    #make bin centers
    ctrs=np.zeros(len(bins))
    widths=np.zeros(len(bins))
    for c in range(len(ctrs)):
      ctrs[c]=0.5*(bounds[c]+bounds[c+1])
      widths[c]=bounds[c+1]-bounds[c]
    #TODO vectorize this binning
    for s in solns:
      i=bisect_left(bounds,s)
      bins[i-1]+=1
    #normalize
    for b in range(len(bins)):
      bins[b]=float(bins[b])/float(len(solns))/widths[b]
    for c,ctr in enumerate(ctrs):
      outFile.writelines(str(ctr)+','+str(bins[i])+'\n')
    outFile.close()
    print '  ...PDF histogram stats added to',self.case+'.stats'
    #pk.dump([ctrs,bins],open(outFile,'wb'))
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
