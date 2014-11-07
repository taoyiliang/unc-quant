import os
import sys

runs=range(10)

#runs=[0   ,2   ,4   ,8  ,
#      10  ,20  ,40  ,80 ,
#      100 ,200 ,400 ,800,
#      1000,2000,4000,8000,
#      10000]

templateName = 'diffusion.unc'
for run in runs:
  #make new input
  procfile = file('processors.txt','r')
  for line in procfile:
    procs=int(line)
  template = file(templateName,'r')
  inname = templateName+'_'+str(run)
  infile = file(inname,'w')
  for line in template:
    if line.startswith('    expOrd'):
      infile.writelines(line.split('=')[0]+' = '+str(run)+'\n')
    elif line.startswith('  numprocs'):
      infile.writelines(line.split('=')[0]+'= '+str(procs)+'\n')
    else:
      infile.writelines(line)
  infile.close()
  template.close()
  osstat=os.system('python Driver.py -hdmr -i '+inname+' | grep STARTING')#-e \'quadrature points used\' -e \'level:\' -e processers')
  if osstat != 0:
    print 'ERROR driver failed with code',osstat
    sys.exit()
  print '\n\nRuns completed: ',runs[0:runs.index(run)+1]
  print 'Runs remaining: ',runs[runs.index(run)+1:],'\n\n'
  os.system('rm '+inname)


