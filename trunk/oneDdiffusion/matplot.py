import numpy as np
import matplotlib.pyplot as plt

nR = 4

def giveMap():
  matmap = [[1,1,1,1],
            [1,1,1,1],
            [1,1,1,1],
            [1,1,1,1]]

  matmap=np.array(matmap)
  plt.figure()
  CS=plt.imshow(np.rot90(matmap),extent=[0,nR,0,nR],interpolation='none')
  plt.colorbar(CS)
  plt.title('Material Layout by Number')
  addRegions()

def addRegions(cPr=1):
  myRange = np.array(range(nR+1))*cPr
  for i in myRange:
    plt.plot(np.ones(nR+1)*i,myRange,'k-')
    plt.plot(myRange,np.ones(nR+1)*i,'k-')
  plt.axis([0,nR*cPr,0,nR*cPr])
