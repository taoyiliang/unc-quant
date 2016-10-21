import os

base_dir = '~/projects/raven/inputs/scgpc_stress/analytic_cases'

cases = {
    'linear_poly':['three','five','ten'],
    'sudret':['three','five'],
    'attenuate':['two','four','six'],
    'sfu_gauss_peak':['three','five'],
    'ishigami' :['three'],
    'sobolG' :['three','five']
         }

for case,nums in cases.items():
  for num in nums:
    dir = '/'.join([base_dir,case,num])
    os.system('cp %s/*pdf .' %dir)
