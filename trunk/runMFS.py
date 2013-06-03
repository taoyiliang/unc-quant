from mfs1Uni import runMFS as u1
from mfs1Norm import runMFS as n1
from mfs1Gamma import runMFS as g1
from mfs1Beta import runMFS as b1
from mfs2un import runMFS as un2
from mfs2gb import runMFS as gb2

import numpy as np
import matplotlib.pyplot as p

p.subplot(2,2,1)
u1()
p.subplot(2,2,2)
n1()
p.subplot(2,2,3)
g1()
p.subplot(2,2,4)
b1()

p.figure()
p.subplot(2,2,1)
un2()
p.subplot(2,2,2)
gb2()

p.show()
