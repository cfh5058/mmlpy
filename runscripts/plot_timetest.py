#!/usr/bin/python
import pysim
import os,copy,pprint
import matplotlib.pyplot as plt
import numpy as np
from mmlutils import *
from pysim.files import gadget

LIST_GDPTYPS=gadget.LIST_PTYPS

Npart_list=[long(1e4),long(1e5),long(1e6),long(1e7)]
Nproc_list=[1,10,50]#100]
compid='stampede'
methLIST=['gadget']
perpart=True
perstep=False
pertime=True

times={'4on1' :{'gadget':18.}, # done
       '4on10':{'gadget':10.}, # done
       '4on50':{'gadget':31.}, # done
       '5on1' :{'gadget':1015.}, # done
       '5on10':{'gadget':300.}, # done
       '5on50':{'gadget':139.}, # done
       '6on1' :{'gadget':3128.+3136.+73473.+18541.},
       '6on10':{'gadget':3082.+3115.+38425.}, # done
       '6on50':{'gadget':3097.+15798.}, # done
       '7on1' :{'gadget':0.}, # won't run
       '7on10':{'gadget':3467.+76313.},
       '7on50':{'gadget':3277.+74024.}}

print len(methLIST)

fig,axs=plt.subplots(1,len(methLIST))#,subplot_kw=dict(aspect='equal'))
for method in methLIST:
    iaxs=plt.subplot(1,1,1)
    # if len(methLIST)==1: iaxs=axs
    # else               : iaxs=axs[methLIST.index(method)]
    plotlist={}
    # Loop over particle number
    for iNpart in Npart_list:
        for iNproc in Nproc_list:
            ilogN=int(np.log10(iNpart))
            istr='{}on{}'.format(ilogN,iNproc)
            if ilogN not in plotlist: plotlist[ilogN]=[]
            if method=='gadget':
                sim=
            plotlist[ilogN].append(times[istr][method])
    # Plot                             
    for iNpart in Npart_list:
        ilogN=int(np.log10(iNpart))
        x=np.array(Nproc_list)
        y=np.array(plotlist[ilogN],dtype=float)
        if perpart: y/=float(iNpart)
        iaxs.plot(x,y,label='log(Npart)={}'.format(ilogN))
    # Legend
    plt.xlabel('Nproc')
    ylab='Time'
    if perpart: ylab+='/particle'
    if perstep: ylab+='/step'
    if pertime: ylab+='/time'
    plt.ylabel(ylab+' (s)')
    iaxs.set_xscale('log')
    iaxs.set_yscale('log')
    iaxs.legend(loc=1)
    fname='/home/langmm/code/python/mine/runscripts/timetest_'+method+'.png'
    plt.savefig(fname)
    print '    '+fname
