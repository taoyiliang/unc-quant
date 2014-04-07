
def newExecutor(exec_type,inp_file,restart):
  oktypes=['PCESC','MC','MLMC']
  if exec_type not in oktypes:
    msg = 'Desired exec type <'+exec_type+'> not recognized.'
    msg+= '\n    Options include '+str(oktypes)
    raise IOError(msg)
  todo = 'ex = '+exec_type+'Exec(inp_file,restart)'
  exec todo
  return ex

class Executor(object):
  def __init__(self,inp_file,restart):
    print 'Initializing executor...'
    self.input_file = inp_file
    self.restart = restart

    self.loadUncInput()
    self.setDirs()
    self.loadDetSolver()
    self.loadVars()
    self.setup() #virtual
    self.loadSampler() #virtual
    self.runSets()

  def loadUncInput(self):
    print 'Reading uncertainty input...'
    self.templateFile=self.input_file('Problem/templateFile','')
    self.execDir    = self.input_file('Problem/execDir'  ,'')
    self.inputDir   = self.input_file('Problem/inputDir' ,'')
    self.input_editor = self.input_file('Problem/input_editor','')

  def setDirs(self):
    print 'Setting directories...'
    self.uncDir  = os.getcwd()
    os.chdir(self.execDir)
    self.execDir = os.getcwd()
    os.chdir(self.uncDir) #using relative paths
    os.chdir(self.inputDir)
    self.inputDir = os.getcwd()
    os.chdir(self.uncDir) #back to start dir
    print '  ...clearing old input files...'
    os.chdir(self.inputDir)
    fail=os.system('rm '+self.templateFile+'.unc*')
    if fail: print '    ...no old inputs found...'
    else: print    '    ...successfully cleared old inputs...'

  def loadDetSolver(self):
    print 'Loading deterministic solver...'
    torun = 'self.ie = InputEditor.'+self.input_editor+'(\''+self.execDir+'\')'
    exec torun

  def loadVars(self):
    print 'Loading uncertain variables...'
    uVars = self.input_file('Variables/names','').split(' ')
    self.varDict={}
    for var in uVars:
      path = self.input_file('Variables/'+var+'/path',' ')
      dist = self.input_file('Variables/'+var+'/dist',' ')
      args = self.input_file('Variables/'+var+'/args',' ')
      for a,arg in enumerate(args):
        args[a]=float(arg)
      self.varDict[var]=Variables.newVar(dist,var,path)
      self.varDict[var].setDist(args)

  def runSets(self):
    pass

class SCESC(Executor): #base class
  def setup(self):
    #set problem case
    self.case=self.templateFile.split('.')[0]+'_SC_'
    self.case+=str(len(self.varDict.keys()))+'var'
    print 'In setup for executor',self.case
    self.makeIndexSet()
    self.makeQuadSet()
    self.makeRunSet()

  def makeIndexSet(self):
    print '  ...generating index set...'
    #generate univariate expansion order lists
    indexranges=[]
    for var in self.varDict.values():
      var.setExpansion(self.input_file)
      indexranges.append(range(var.expOrd))
    maxorder = np.max(np.max(indexranges))
    #create index set
    iset = self.input_file('Sampler/SC/indexSet','dud')
    if iset == 'dud':
      print '  ...index set not specified; using tensor product...'
      iset = 'TP'
    self.indexSet = IndexSets.chooseSet(indexranges,iset,maxorder)
    print '  ...using %i expansion moments...' %len(self.indexSet)

  def makeQuadSet(self):
    multfac = self.input_file('Sampler/SC/quadFactor',0)
    if multfac==0:
      print '  ...quad mult fact not specified; using 2...'
      multfac = 2
    #find largest quad needed for each var, set quad
    self.quadSet=[]
    listSets=zip(*self.indexSet)
    maxOrds=[]
    for v,var in enumerate(self.varDict.values()):
      maxOrds.append(max(listedOrds[v])+1)
      var.setQuadrature(maxOrd[-1])
    #make run list
    #TODO this is where sparse quadratures should be made
    ranges=[]
    for mo in maxOrds:
      ranges.append(range(mo))
    self.runQuadOrds=list(allcombos(*ranges))
    print '  ...using %i quadrature points...' %len(self.quadSet)

  def loadSampler(self):
    self.sampler = spr.StochasticPoly(self.varDict,
                                      self.input_file,
                                      self.runQuadOrds)
