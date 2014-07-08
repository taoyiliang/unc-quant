from itertools import product

one = range(6)
test=[]
for i in range(10):
  test.append(one)

#bm = list(product(*test))


for i in range(len(test)):
  if i==0:
    trial = [test[i]]
  else:
    trial = list(product(test[i],*trial))
