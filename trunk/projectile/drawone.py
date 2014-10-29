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


pos,dist = R(m,r,C,rho,v,ang,g)
xx=pos[:,1]
yy=pos[:,2]


fig = plt.figure()
ax = plt.axes(xlim=(0,350),ylim=(0,120))
line,=ax.plot([],[],lw=2)
#line=plt.plot(x,y)

def update_plot(i,xx,yy):
  p=i*10
  x=xx[:p]
  y=yy[:p]
  line.set_data(x,y)
  return line,

def init():
  line.set_data([],[])
  return line,

ani = animation.FuncAnimation(fig,update_plot,
                   init_func=init,
                   frames=len(xx)/10,
                   interval=1,
                   fargs=(xx,yy),blit=True,repeat=False)

#ani.save('one_shot.mp4',fps=30)#,extra_args=['-vcodec','libx264'])

plt.show()

