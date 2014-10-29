import numpy as np
from proj import R

from matplotlib import pyplot as plt
from matplotlib import animation

m=0.145
r=0.0336
C=0.5
rho=1.2
v=50.
ang=35
g=9.81

sets=[]
for j in range(1000):#np.linspace(0.1,1,10):
  c = np.random.rand()
  pos,dist = R(m,r,c,rho,v,ang,g,dt=0.01,verbose=False)
  sets.append([pos[:,1],pos[:,2]])

fig = plt.figure()
ax = plt.axes(xlim=(0,350),ylim=(0,120))
line,=ax.plot([],[],lw=2)
#line=plt.plot(x,y)

def update_plot(i,sets,dud):
  ax.plot(sets[i][0],sets[i][1])
  #line.set_data(
  return line,

def init():
  line.set_data([],[])
  return line,

ani = animation.FuncAnimation(fig,update_plot,
                   init_func=init,
                   frames=len(sets),
                   interval=1,
                   fargs=(sets,1),blit=True,repeat=False)

ani.save('mc_C.mp4',fps=30)#,extra_args=['-vcodec','libx264'])

plt.show()

