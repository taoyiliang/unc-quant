import parExecutor as ex
import cPickle as pk

#load executor
print 'Loading executor for restart...'
rstDict=pk.load('executor.backup.pk','rb')
for key in rstDict.keys():
  print key,rstDict[key]
#newexec = Executor(restart=True,rstDict=rstDict)
