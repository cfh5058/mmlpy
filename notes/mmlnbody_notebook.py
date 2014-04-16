#!/usr/bin/python                                                                                                                                                                  
import sys,os,shutil,glob,copy,pprint,scipy,math,pylab
import numpy as np
#from mmlutils import *
#from mmlastro import *
import mmlnbody
#from mmlnbody import *

from yt.imods import *

# Select the run
runtag='MW_6M11LH11FG0'

# Get the simstr object
simstr=mmlnbody.main.runtag2simobj(runtag)

# Get the snapshot names
fns=sorted(glob.glob(simstr.fdict['gadget']['output']['snapbase']+'*'))
ifn=fns[0]

# Load snap
ipf = load(ifn)

# Print stats
ipf.h.print_stats()
ipf.h.field_list
ipf.h.derived_field_list
