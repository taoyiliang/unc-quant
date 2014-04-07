import matplotlib.pyplot as plt
import matplotlib.mlab as mlab
import numpy as np
from bisect import bisect_left
from sys import *
import os

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
          temp_var *=var.evalNormPoly(histories['varVals'][m][v],ordset[v])
          if verbose:
            print 'var',var.name,'order',ordset[v]
            print '  weight:',temp_var
            print '  value: ',histories['varVals'][m][v]
            print '  poly:',var.evalNormPoly(histories['varVals'][m][v],ordset[v])
            print '  eval:',temp_var
          temp_order*=temp_var
        c[tuple(ordset)]+=temp_order*soln
        #FIXME this is hacky, but eliminates some roundoff errors
        if abs(c[tuple(ordset)])<1e-12:
          c[tuple(ordset)]=0
    #remove zero coefficients
    noelim=True
    for key in c.keys():
      if c[key]<1e-12:
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
    #TODO only set for one variable; fix this!
    tot = 0
    for c in self.coeffs.keys():
      coef = self.coeffs[c]
      temp = 1
      for v,var in enumerate(histories['vars']):
        temp *= var.evalNormPoly(vals[v],c[v])
      #print 'key,coef[key]',c,coef
      tot+=coef*temp
    #act = 5. + 4.*vals[0] + 3.*vals[0]*vals[0] + 2.*vals[0]**3 + vals[0]**4
    #return tot,act
    return tot

  def MCsampleDict(self,histories,trials=1000,cprFunc=None):
    #should I make a separate non-dictionary version of this?
    solns={}
    actual={}
    trials = int(trials)
    #TODO this can be vectored? 
    for i in range(trials):
      vals=[]
      actvals=[]
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
      vals=[]
      for var in histories['vars']:
        #actvals.append(var.sample())
        #vals.append(var.convertToStandard(actvals[-1]))
        vals.append(var.sample())
      solns[i]=self.sampleROM(vals,histories)
    return solns

  def MCsampleVec(self,histories,trials=1000):
    trials = int(trials)
    solns=np.zeros(trials)
    #build value data TODO


  def makePDF(self,histories,bins,samples=1e6,solns=None):
    if solns==None:
      solns=self.MCsample(histories,int(samples))
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
