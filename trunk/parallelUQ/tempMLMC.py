

  def MLMCRun(self):
    print 'Starting MLMC runs...'
    #read in MLMC inputs
    self.tol = self.input_file('Sampler/MLMC/tol',0.0)
    if self.tol==0.0:
      print '  ...tol not specified, using 1e-4...'
      self.tol = 1e-4

    self.beta = self.input_file('Sampler/MLMC/beta',0)
    if self.beta==0:
      print '  ...beta not specified, using 2...'
      self.beta=2

    h0 = self.input_file('Sampler/MLMC/h0',0)
    if h0==0:
      print '  ...h0 not specified, using 1...'
      h0=1

    m0 = self.input_file('Sampler/MLMC/m0',0)
    if m0==0:
      print '  ...m0 not specified, using 1e3...'
      m0=1e3
    m0=int(m0)

    #TODO
    self.q=0
    #algorithm
    self.L=-1
    self.h=[] # l x 1, mesh size
    self.M=[] # l x 1, num solns
    self.G=[] # l x 1
    self.barG=[] #l x m
    self.v=[] #l x 1, var
    self.W = [] #l x 1, work
    self.runs=[] #m x {#vars}, input dicts
    self.solns=[] #l x m, soln vals
    converged=False
    while not converged and L<=10:
      L+=1
      self.h.append(h0*self.beta**(-self.L))
      self.M.append(m0)
      #TODO if m0 > len(inputs): generate inputs
      #TODO run samples
      #TODO G[l],barG[l]
      self.v.append(self.get_vL(L))
      #TODO loop update all M[l]
      #          if M[l] > inputs size, generate more inputs
      #          if M[l] got bigger, run more samples
      #          recalculate v[L]
      #     loop until sufficient M[L] done, and v[L] set
      #     calculate A
      #     calculate total error estimate
      #     check tolerance



  def get_vL(self,l):
    tot=0.0
    for m,soln in self.solns[l]:
      tot+=(self.G[l][m]-self.barG[l])**2
    return tot*1./self.M[l]

  def get_ML(self,l):
    #TODO c_alph?
    first = (self.tol/c_alph)**(-2)
    second = np.sqrt(self.v[l]/self.M[l])
    third=0.0
    for l in range(L+1):
      third+=np.sqrt(self.v[l]*self.W[l])
    return first*second*third

  def get_err(self):
    til_v = 0.0
    for l in range(L+1):
      til_v+=self.v[l]/self.M[l]
    c_w=max(get_cw(-1),get_cw(-2))
    err1 = c_w*self.h[-1]**self.q
    err2 = c_alph*np.sqrt(til_v)

  def get_cw(self,l):
    den=self.h[l]**self.q*(self.beta**self.q-1.0)
    return abs(barG[l])/den

  def generateInputs(self,num,prefix,start):
    new=[]
    totIndx = int(start)
    for m in range(len(num)):
      new.append({})
      totIndx+=1
      ident = '_'+str(prefix)+'_'+str(totIndx)
      new[-1]['id']=ident
      new[-1]['varVals']=self.sampler.giveSample() #inserts 'varVals'
    return new

  def runParallel(self,Ml,l,h):
    torun = self.inputs[:]
    #TODO is this a waste?
    self.solns[l]=np.zeros(len(torun))
    #set up multiprocessing queues
    self.inq=sque(1000) #inputs to be run
    self.outq=sque(1000) #finished solns
    #set up processors
    procs=[]
    wantprocs = self.input_file('Problem/numprocs',1)
    numprocs = min(wantprocs,multiprocessing.cpu_count())
    print '  ...using %i processors...' %numprocs
    #set up histories
    self.histories={}
    self.histories['varNames']=self.varDict.keys()
    self.histories['vars']=self.varDict.values()
    self.histories['varPaths']=[]
    for var in self.varDict.values()
      self.histories['varPaths'].append(var.path)
      self.histories['varVals']=[]
      self.histories['nRun']=[]
      self.histories['soln']=[]
    #start loop
    m=-1
    finished=False
    numFinished=0
    while not finished:
      #fill new inputs
      while True:
        m+=1
        try: self.inq.put([m,torun[m]])
        except Full:
          m-=1
          break
      #collect outputs
      for p,proc in enumerate(procs):
        if not proc.is_alive():
          proc.join()
          while not self.outq.empty():
            m,soln = self.outq.get()
            self.solns[l][m]=soln
            numFinished+=1
          del procs[p]
      if len(torun)>0:
        while len(procs)<numprocs and len(torun)>0:
          newrun=torun.pop()
          newrun['fileChange']['Output/file']='run.out'+str(newrun['id'])
          meshsize = int(1.0/self.h[l])
          newrun['fileChange']['Mesh/nx_per_reg']=meshsize
          newrun['fileChange']['Mesh/ny_per_reg']=meshsize
          inp_file = self.ie.writeInput(self.templateFile,
                                        self.inputDir,
                                        self.histories['varPaths'],
                                        newrun['varVals'],
                                        newrun['fileChange'],
                                        m)
          procs.append(multiprocessing.Process(\
                         target = self.runSample,\
                         args=(m,inp_file,outFileName)))
          procs[-1].start()
          self.runInput(newrun) #TODO run a new sample

  def runInput(self,run):
    pass #TODO

