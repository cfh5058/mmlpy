#!/usr/bin/python 
import Image
import matplotlib
from matplotlib import pyplot
import pysim
import os,pprint
from mmlutils import mmlplot,mmlparam,mmlio,mmlmath
import numpy as np
from pysim import paperplots

# 115 3 armed spiral
# 135 peri
# 205 strongest bar
# 375 bar destruction
runtags=['MW0B7_MW0B7_Q1RP0p1E1I0',
         'MW0B7_MW0B6_Q10RP0p1E1I0',
         'MW0B7_MW0B7_Q1RP0p1E1I180']
#runtags=['MW0B7_MW0B7_Q1RP0p1E1I0']
#runtags=['MW0B7_MW0B6_Q10RP0p1E1I0']
#runtags=['MW0B7_MW0B7_Q1RP0p1E1I180']

#method='print_param'
#method='orbit'
method='barevol'
#method='testellipse'
#method='testbarfit'

#method='plotsnap'
plotlist=['image','ellipse']
plotlist_anno=['ellipse']

#extlist=[115,135,205,375]
extlist=[135,205,375]
#extlist=[135]
#extlist=[205]
#extlist=[375]
#extlist=[845]
#extlist=[585]
#extlist=[175]
extlist_cbar=[375]

overwrite=False
owdata=False
owhist=True

for runtag in runtags:
    if runtag=='MW0B7_MW0B6_Q10RP0p1E1I0': galList=['galaxy1','galaxy2']
    else                                 : galList=['galaxy1']
    typList=['disk']

    # Get simboj
    x=pysim.loadsim(runtag)

    # Print info on interaction & galaxies
    if method=='print_param': x.show_galmod()

    # Orbit
    elif method=='orbit': out=x.orbit_param(cmeth='pot',ftype='gadget',askuser=False,verbose=True)

    # Test ellipse fitting routines
    elif method=='testellipse':
        from mmlastro import ellipse
        ellipse.test('meth_derivt',atest=4.07,verbose=True,flipellip=True)
        break

    elif method=='testbarfit':
        from mmlastro import bar
        #bar.test('ellipse',Na=100,warn=False,rmin=0.5)
        bar.test('ellipse',owplot=True,meth_derivt='rings',flipellip=True,cutmeth='err')
        break

    # Create snapshot plots
    elif method=='plotsnap':
        for plot in plotlist:
            # Loop over galaxies
            for gal in galList:
                # Loop over types
                for typ in typList:
                    # Loop over extensions
                    for ext in extlist:
                        print runtag,gal,typ,ext
                        paperplots.snapshot(x,gal,typ,ext,plot,owhist=owhist,plotext='eps',
                                            colorbar=(ext in extlist_cbar),
                                            annotate=(plot in plotlist_anno))
    # Create evolution plots
    elif method=='barevol':
        print '\n\n'+runtag
        tlim=[0.,5.]
        #methods=['ellipse']
        methods=['scfAm2','ellipse']
        # Loop over galaxies
        for gal in galList:
            if gal=='galaxy2': rmin=0.3
            else             : rmin=1.0
            # Loop over types
            for typ in typList:
                paperplots.barevolution(x,gal,typ,methods=methods,tlim=tlim,plotext='eps',
                                        overwrite=overwrite,owhist=owhist,owdata=owdata,
                                        rmin=rmin)


