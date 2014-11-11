import os
import sys
import numpy as np

class InputEditor:
  def __init__(self,runpath):
    '''constructor'''
    pass

  def storeOutput(self,outFile):
    '''
    obtains desired value from one output file
    output: float, the desired value
    '''
    pass

  def runSolve(self,input_file,outFile):
    '''
    runs the problem
    '''
    pass

  def writeInput(self,templateName,inputDir,varList,valList,otherChange,ident):
    '''
    generates an input file to run.
    Inputs: templateName - the base input file name to work from
            inputDir     - the directory the base input file is located in
            varList      - list of variable namepaths to look for
            valList      - list of values to set vars in varList to
            ident        - the run number for the history
    '''
    curPath=''
    newSectionFlag=True

    startDir = os.getcwd()
    os.chdir(inputDir)
    #print os.getcwd()
    readFile=file(templateName,'r')
    writeFileName = templateName+'.unc'+str(ident)
    writeFile=file(writeFileName,'w')

    for line in readFile:
      if line.strip().startswith('#'):continue
      if newSectionFlag and line[0]=='[': #new section
        curPath=line.strip()[1:-1]
        newLine=line[:]
      elif line.strip()[0:5]=='[../]': #go up a level
        curPath='/'.join(curPath.split('/')[:-1])
        newLine=line[:]
      elif line.strip()[0:3]=='[./': #go down a level
        curPath=curPath+'/'+line.strip()[3:-1]
        newLine=line[:]
      elif line.strip()!='': #normal entry
        keyword=curPath+'/'+line.split('=')[0].strip()
        if keyword in varList: #uncertain parameter
          indx = varList.index(keyword)
          newLine=line.split('=')[0]+'= '+str(valList[indx])+'\n'
        elif keyword in otherChange.keys(): #not uncertain, but still change
          newLine=line.split('=')[0]+'= '+str(otherChange[keyword])+'\n'
        else: #no change to entry
          newLine = line[:]
      else: #not an interesting line
        newLine=line[:]
      writeFile.writelines(newLine)
    writeFile.close()
    readFile.close()
    os.chdir(startDir)
    return writeFileName

class IE_Simple(InputEditor):
  def __init__(self,runpath=''):
    self.type = 'InputOutput simple'
    self.inp=0
    self.out=0

  def writeInput(self,templateName,inputDir,varList,valList,otherChange,ident):
    self.inp = valList[0]
    return 'dud'

  def storeOutput(self,outFile):
    return self.out

  def runSolve(self,input_file):
    sys.path.insert(0,os.getcwd())
    from simple import g
    self.out=g(self.inp)

class IE_Double(InputEditor):
  def __init__(self,runpath=''):
    self.type = 'InputOutput simple'
    self.x=0
    self.y=0
    self.out=0

  def writeInput(self,templateName,inputDir,varList,valList,otherChange,ident):
    self.x = valList[0]
    self.y = valList[1]
    return 'dud'

  def storeOutput(self,outFile):
    return self.out

  def runSolve(self,input_file):
    sys.path.insert(0,os.getcwd())
    from simple import f
    self.out=f(self.y,self.x)

class IE_projectile(InputEditor):
  def __init__(self,runpath=''):
    self.type = 'InputOutput projectile'
    self.a=[]
    self.out=0

  def writeInput(self,templateName,inputDir,varList,valList,otherChange,ident):
    self.varList = varList
    self.valList = valList
    self.dt = otherChange['dt']
    return 'dud'

  def storeOutput(self,outFile):
    return self.out

  def runSolve(self,input_file):
    sys.path.insert(0,os.getcwd())
    from proj import R
    #default values
    m=0.145
    r=0.0336
    C=0.5
    rho=1.2
    v=50
    ang=35
    g=9.7988
    sy=0
    for p,path in enumerate(self.varList):
      var = path.split('/')[-1]
      expr = var+' = '+str(self.valList[p])
      exec(expr)
    self.out=R(m,r,C,rho,v,ang,g,sy=sy,dt=self.dt,retpts=False)

class IE_Thirty(InputEditor):
  def __init__(self,runpath=''):
    self.type = 'InputOutput simple'
    self.a=[]
    self.out=0

  def writeInput(self,templateName,inputDir,varList,valList,otherChange,ident):
    self.a = valList
    return 'dud'

  def storeOutput(self,outFile):
    return self.out

  def runSolve(self,input_file):
    sys.path.insert(0,os.getcwd())
    from simple import h5 as h
    self.out=h(self.a)

class IE_Thirty5(IE_Thirty):
  def runSolve(self,input_file):
    sys.path.insert(0,os.getcwd())
    from simple import h5 as h
    self.out=h(self.a)
class IE_Thirty10(IE_Thirty):
  def runSolve(self,input_file):
    sys.path.insert(0,os.getcwd())
    from simple import h10 as h
    self.out=h(self.a)
class IE_Thirty15(IE_Thirty):
  def runSolve(self,input_file):
    sys.path.insert(0,os.getcwd())
    from simple import h15 as h
    self.out=h(self.a)
class IE_Thirty30(IE_Thirty):
  def runSolve(self,input_file):
    sys.path.insert(0,os.getcwd())
    from simple import h30 as h
    self.out=h(self.a)

class IE_Source(InputEditor):
  def __init__(self,runpath=''):
    self.type = 'InputOutput source'
    self.varVals=[]
    self.varList=[]
    self.out=0

  def writeInput(self,templateName,inputDir,varList,valList,otherChange,ident):
    self.varList=varList[:]
    self.varVals=valList[:]
    try:
      self.d=otherChange['Mesh/nx_per_reg']
    except KeyError:
      pass
    return 'dud'

  def storeOutput(self,outFile):
    return self.out

  def runSolve(self,input_file):
    sys.path.insert(0,os.getcwd())
    from solver import solve
    expr='soln=solve('
    for v,var in enumerate(self.varList):
      vname=var.split('/')[-1]
      expr+=vname+'='+str(self.varVals[v])+','
    expr+='N='+str(self.d)+')'
    exec expr
    #norm=np.sqrt(1./len(soln)*sum(soln**2))
    #print 'soln',norm
    i = (len(soln)-1)/4
    self.out=soln[i]

class IE_OldSource(InputEditor):
  def __init__(self,runpath=''):
    self.type = 'InputOutput soure'
    self.varVals=[]
    self.varList=[]
    self.out=0

  def writeInput(self,templateName,inputDir,varList,valList,otherChange,ident):
    self.varList=varList[:]
    self.varVals=valList[:]
    return 'dud'

  def storeOutput(self,outFile):
    return self.out

  def runSolve(self,input_file):
    sys.path.insert(0,os.getcwd())
    from source import g
    expr='self.out=g('
    for v,var in enumerate(self.varList):
      vname=var.split('/')[-1]
      expr+=vname+'='+str(self.varVals[v])+','
    expr=expr[:-1]+')'
    exec expr
    #self.out=g(self.inp)

class IE_Diffusion(InputEditor):
  def __init__(self,runpath=''):
    self.type = 'InputOutput diffusion'

  def storeOutput(self,outFile):
    readFile=file(outFile,'r')
    for line in readFile:
      if line.split(',')[0]=='k':
        readFile.close()
        val = float(line.split(',')[1].strip())
        #val = 1.0+float(line.split(',')[1].split('1+')[1].strip())
        #print 'line:',line[:-1].strip(),
        #print '  grabbed',val
        if val > 0:
          os.system('rm '+outFile+'\n')
        return val#float(line.split(',')[1].strip())

  def runSolve(self,input_file):
    osstat = os.system('./TwoDProblem -i '+input_file+' > /dev/null 2>&1')
    if osstat != 0:
      print 'Run attempt failed with error code',osstat
      sys.exit()
    os.system('rm '+input_file)

class HDMR_IO(InputEditor):
  def __init__(self,path=''):
    self.type = 'HDMR_IO_editor'

  def writeInput(self,templateName,changelist,ident):
    curPath=''
    newSectionFlag=True

    readFile = file(templateName,'r')
    writeFileName = templateName+'.'+str(ident)
    writeFile = file(writeFileName,'w')

    for line in readFile:
      if line.strip().startswith('#'):continue
      if newSectionFlag and line[0]=='[': # new section
        curPath = line.strip()[1:-1]
        newLine=line[:]
      elif line.strip()[:5]=='[../]': # go up a level
        curPath='/'.join(curPath.split('/')[:-1])
        newLine=line[:]
      elif line.strip()[:3]=='[./': # add a level
        curPath = curPath + '/' + line.strip()[3:-1]
        newLine=line[:]
      elif line.strip()!='': # value entry
        keyword = curPath+'/'+line.split('=')[0].strip()
        if keyword in changelist.keys():
          newLine = line.split('=')[0] + '= ' + str(changelist[keyword])+'\n'
        else:
          newLine=line[:]
      else:
        newLine=line[:]
      writeFile.writelines(newLine)
    writeFile.close()
    readFile.close()
    return writeFileName

