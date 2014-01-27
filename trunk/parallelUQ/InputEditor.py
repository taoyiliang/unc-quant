import os
import sys

class InputEditor:
  def __init__(self):
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
  def __init__(self):
    self.type = 'InputOutput simple'

  def storeOutput(self,outFile):
    readFile=file(outFile,'r')
    for line in readFile:
      if line.split(',')[0]=='res':
        val=float(line.split(',')[1].strip())
        #print 'grabbed',val
        readFile.close()
        os.system('rm '+outFile+'\n')
        return val

  def runSolve(self,input_file):
    #print os.getcwd(),input_file
    osstat = os.system('python simple.py -i '+input_file+' > /dev/null')
    if osstat != 0:
      print 'Run attempt failed with error code',osstat
      sys.exit()
    os.system('rm '+input_file)

class IE_Source(InputEditor):
  def __init__(self):
    self.type = 'InputOutput soure'

  def storeOutput(self,outFile):
    readFile=file(outFile,'r')
    for line in readFile:
      if line.split(',')[0]=='res':
        val=float(line.split(',')[1].strip())
        #print 'grabbed',val
        readFile.close()
        os.system('rm '+outFile+'\n')
        return val

  def runSolve(self,input_file):
    #print os.getcwd(),input_file
    osstat = os.system('python source.py -i '+input_file+' > /dev/null')
    if osstat != 0:
      print 'Run attempt failed with error code',osstat
      sys.exit()
    os.system('rm '+input_file)

class IE_Diffusion(InputEditor):
  def __init__(self):
    self.type = 'InputOutput diffusion'

  def storeOutput(self,outFile):
    readFile=file(outFile,'r')
    for line in readFile:
      if line.split(',')[0]=='k':
        #print 'grabbed',line.split(',')[1].strip()
        readFile.close()
        val = float(line.split(',')[1].strip())
        if val > 0:
          os.system('rm '+outFile+'\n')
        return val#float(line.split(',')[1].strip())

  def runSolve(self,input_file):
    osstat = os.system('./TwoDProblem -i '+input_file+' > /dev/null')
    if osstat != 0:
      print 'Run attempt failed with error code',osstat
      sys.exit()
    os.system('rm '+input_file)
