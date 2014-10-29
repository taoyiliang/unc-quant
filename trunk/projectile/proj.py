import numpy as np

def R(m,r,C,rho,v,ang,g,sy=0,dt=0.001,tol=1e-12,verbose=False):
  D = rho*C*0.5*(np.pi*r*r)
  if C==0:D=0
  ang*=np.pi/180.
  t = 0
  sx = 0.
  #sy taken as input parameter
  vx = v*np.cos(ang)
  vy = v*np.cos(ang)
  nores = vx/g*(vy+np.sqrt(vx*vx+2*g*sy))
  if verbose: print 'Expected no-resistance value:',nores
  converged = False
  pos=[]
  while not converged:
    old,new = takeAStep(t,dt,vx,vy,sx,sy,D,m,g,verbose=False)
    if new[2]>0:
      pos.append(old)
      t=new[0]
      sx=new[1]
      sy=new[2]
      vx=new[3]
      vy=new[4]
    else:
      if dt>=tol:
        dt*=0.1
      else:
        converged=True
        pos.append(old)
        pos.append(np.array(new))
  pos=np.array(pos)
  if verbose:
    print 'END STATS'
    print '  Range:',pos[-1][1]
    print '  Max height:',np.max(pos[:,2])
    print '  Time of Flight:',pos[-1][0]
    print pos[-1]
  return pos,new[1]

def takeAStep(t,dt,vx,vy,sx,sy,D,m,g,verbose=True):
  v = np.sqrt(vx*vx + vy*vy)
  ax = -D/m*v*vx
  ay = -g-D/m*v*vy
  old = [t,sx,sy,vx,vy]
  vx += ax*dt
  vy += ay*dt
  sx += vx*dt + 0.5*ax*dt*dt
  sy += vy*dt + 0.5*ay*dt*dt
  if verbose:
    print '  At %1.2e, %1.2e, t,dt=%1.2e, %1.2e' %(sx,sy,t,dt)
  t += dt
  new = [t,sx,sy,vx,vy]
  return old,new
