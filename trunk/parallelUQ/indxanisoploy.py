import IndexSets as isets
import matplotlib.pyplot as plt
import numpy as np

wts = [2,1]

TP = np.array(isets.IndexSetFactory(2,4,'TP',impwts=wts))
TD = np.array(isets.IndexSetFactory(2,4,'TD',impwts=wts))
HC = np.array(isets.IndexSetFactory(2,4,'HC',impwts=wts))


plt.figure(figsize=(13.6,3.7))
plt.subplot(1,3,1)
plt.plot(TP[:,0],TP[:,1],'kx',markersize=10)
plt.axis([-0.2,4.2,-0.2,4.2])
plt.subplot(1,3,2)
plt.plot(TD[:,0],TD[:,1],'kx',markersize=10)
plt.axis([-0.2,6.2,-0.2,6.2])
plt.subplot(1,3,3)
plt.plot(HC[:,0],HC[:,1],'kx',markersize=10)
plt.axis([-0.2,10.2,-0.2,10.2])
#plt.show()
plt.savefig('SGaniso.pdf',format='pdf')

