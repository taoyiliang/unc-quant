import IndexSets as isets
import matplotlib.pyplot as plt
import numpy as np

TP = np.array(isets.IndexSetFactory(2,4,'TP'))
TD = np.array(isets.IndexSetFactory(2,4,'TD'))
HC = np.array(isets.IndexSetFactory(2,4,'HC'))

plt.plot(TP[:,0],TP[:,1],'k.',markersize=5)
plt.plot(TD[:,0],TD[:,1],'k+',markersize=10)
plt.plot(HC[:,0],HC[:,1],'ko',markersize=10)
plt.axis([-0.2,4.2,-0.2,4.2])
plt.show()


