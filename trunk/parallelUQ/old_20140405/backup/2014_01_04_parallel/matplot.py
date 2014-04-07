import numpy as np
import matplotlib.pyplot as plt

def giveMap():
  matmap = [[2,1,1,1,1,2,2,3,3,5,5],
            [1,1,1,1,1,1,1,3,3,5,5],
            [1,1,1,1,1,1,1,3,3,5,5],
            [1,1,1,1,1,1,1,3,3,5,5],
            [1,1,1,1,1,1,1,3,3,5,5],
            [2,1,1,1,1,2,2,3,3,5,5],
            [2,1,1,1,1,2,2,3,3,5,5],
            [3,3,3,3,3,3,3,4,5,5,5],
            [3,3,3,3,3,3,3,5,5,5,5],
            [5,5,5,5,5,5,5,5,5,5,5],
            [5,5,5,5,5,5,5,5,5,5,5]]

  matmap=np.array(matmap)
  plt.figure()
  CS=plt.imshow(np.rot90(matmap),extent=[0,11,0,11],interpolation='none')
  plt.colorbar(CS)
  plt.title('Material Layout by Number')
  addRegions()

def addRegions(cPr=1):
  myRange = np.array(range(12))*cPr
  for i in myRange:
    plt.plot(np.ones(12)*i,myRange,'k-')
    plt.plot(myRange,np.ones(12)*i,'k-')
  plt.axis([0,11*cPr,0,11*cPr])
