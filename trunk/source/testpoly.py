
class Sampler:
  def __init__(self,avg,rg):
    self.range=rg
    self.average=avg

  def convertToActual(self,x):
    return self.range*x+self.average

  def convertToStandard(self,x):
    return (x-self.average)/self.range

  def sample(self,y):
    x=self.convertToStandard(y)
    print 'sampling at',x
    return -0.06922*x*x*x + 0.1241*x*x - 0.0897*x + 1.003
