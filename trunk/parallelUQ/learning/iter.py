

class test(object):
  def __init__(self):
    self.counter=-1

  def __iter__(self):
    self.counter+=1
    return self.counter
