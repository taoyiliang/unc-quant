import matplotlib.pyplot as plt
import numpy as np

#uniform
xs = [4,8,16]
m = [#1.2547221522,
      1.25569029702,
      1.25569096924,
      1.25569096924]
v = [#0,
     0.049198975952,
     0.0492316191443,
     0.0492316191611]

ys = [0,16]
MCm =[1.25554670335]*2
MCv =[0.0492854851117]*2

plt.figure()
plt.subplot(2,1,1)
plt.plot(xs,m,'g-o')
#plt.plot(ys,MCm,'b:')
plt.title('Mean')

plt.subplot(2,1,2)
plt.plot(xs,v,'g-o')
#plt.plot(ys,MCv,'b:')
plt.title('Variance')

#plt.show()


#normal
#time
# 2 
# 4 
# 8 4.74
#16 6.88
xs=[2,4,8,16]
m=[1.2547221522,
   1.25569029702,
   1.25569096924,
   1.25569096924]

v=[0.0,
   0.049198975952,
   0.0492316191443,
   0.0492316191611]

ys = [0,16]
MCm =[1.24922240195]*2
MCv =[0.0488719424418]*2

plt.figure()
plt.subplot(2,1,1)
plt.plot(xs,m,'g-o')
#plt.plot(ys,MCm,'b:')
plt.title('Mean')

plt.subplot(2,1,2)
plt.plot(xs,v,'g-o')
#plt.plot(ys,MCv,'b:')
plt.title('Variance')

plt.show()
