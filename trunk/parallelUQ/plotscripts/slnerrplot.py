from plot1 import addPlot
from plotMC import addPlot as pltMC
from plotHDMR import addPlot as addHDMR
import matplotlib.pyplot as plt
import numpy as np


def slnerrplot(cases,title,N,xlim,sylim,eylim,ref,mom):
    slnplot=plt.figure()
    errplot=plt.figure()
    plt.figure(errplot.number)

    for cname,cfile in cases.iteritems():
      ary=cname.split('.')[0].split('_')
      if ary[0]=='MC':
        pltMC(cname,cfile,'MC',ref=ref,r=mom,slnfig=slnplot)
      elif ary[0] in ['HC','TD']:
        ary=ary[0]#.remove(ary[1])
        addPlot(cname,cfile,ary,ref=ref,r=mom,slnfig=slnplot)
      elif ary[0]=='hdmr':
        namelist = cname.split('.')[0].split('_')
        name = '_'.join([namelist[0],namelist[1],namelist[-1]])
        addHDMR(cname,cfile,name,ref=ref,r=mom,slnfig=slnplot)


    plt.title(r'Error in $\mathbb{E}[u^%i]$; %s, $N$=%i'\
                       %(mom,title,N))
    plt.xlabel(r'PDE Solves $\eta$')
    plt.ylabel('Rel. Error')
    plt.xlim(xlim)
    plt.ylim(eylim)
    plt.legend(loc=3)

    plt.figure(slnplot.number)
    plt.title(r'Solution for $\mathbb{E}[u^%i]$; %s, $N$=%i'\
                       %(mom,title,N))
    plt.xlabel(r'PDE Solves $\eta$')
    plt.ylabel(r'$\mathbb{E}[u^%i]$' %mom)
    plt.legend(loc=4)
    plt.xlim(xlim)
    plt.ylim(sylim)
    plt.gca().set_xscale('log')

