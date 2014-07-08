
cofFile = '2G_2D_SC_5var.SC.cof'
rgs = [0.0055,0.011,0.0041,0.01,0.016]
factor = 0.5**5
for r in rgs:
  factor*= 1./r

mean=1.00014901416
#read in coefficients
tot=0
for line in file(cofFile,'r'):
  cof = float(line.split(' | ')[1])
  tot+=cof*cof * factor

var = tot - mean*mean
print 'mean is',mean
print 'var is',var
