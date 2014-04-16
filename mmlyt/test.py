#!/usr/bin/env python
import glob
#import mmlyt
#import pynbody
import pysim
from mmlyt.main import runtag2simobj

runtag='MW_6M11LH11FG0'

simstr=runtag2simobj(runtag)

fns=sorted(glob.glob(simstr.fdict['gadget']['output']['snapbase']+'*'))
ifn=fns[0]
ifn=simstr.fdict['gadget']['input']['ic']

pysim.load(ifn)

#print mmlyt.mmlgadget.test_snap(ifn)
#print mmlyt.mmlgadget.get_blocklist(ifn)

#inb=mmlyt.mmlgadget.rw_snap('R',ifn)
#print inb
#inb=mmlyt.load(ifn)
