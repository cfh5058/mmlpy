#!/usr/bin/python 
import sys
import pysim
import os
from mmlutils import mmlplot,mmlio,mmlmath,mmlplot
from matplotlib import pyplot as plt
import time
tic = time.clock()

runtag=None
#runtag='MW7_MW7_Q1RP0p05E1I0_SFTNEW'
#runtag='MW_6p52M11p52LH11FG0'
#runtag='MW7_MW7_Q1RP0p1E1I0_SFTBLG'
#runtag='MW7_MW7_Q1RP0p1E1I0_SFTNEW'
#runtag='MW0B_6M11LH11FG0'
runtag='MW0B_7M12LH11FG0'
#runtag='MW0B_8M12LH11FG0'
#runtag='MW0B_7M11LH11FG0'
#runtag='MW0B7_MW0B7_Q1RP0p1E1I0'
#runtag='MW0B7_MW0B7_Q1RP0p1E1I180'
#runtag='MW0B7_MW0B6_Q10RP0p1E1I0'
#runtag='MW0B8_MW0B8_Q1RP0p1E1I0_TEST'

method='status'
#method='testtype'
#method='fithalo'
#method='setlimits'
#method='status'
#method='subbuildgal'
#method='endbuildgal'
#method='subgadget'
#method='endgadget'
#method='delgadget'
#method='subscf'
method='endscf'
#method='plotimg'
#method='plotctr'
#method='plotell'
#method='plotscf'
#method='plotpnl'
#method='plotbar_Am2'
#method='plotbar_derAm2'
#method='plotbar_rot'
#method='plotbar_ellipse'
#method='plotorbit'
#method='orbitparam'
#method='plotbarprop_scfAm2'

# Pars commandline input
if len(sys.argv)>1:
    runtag=sys.argv[1]
if len(sys.argv)>2:
    method=sys.argv[2]
print runtag,method

# Aliases
if   method=='plotsnap'   : method='plotimg'
elif method=='plotellipse': method='plotell'
elif method=='plotpanel'  : method='plotpnl'
elif method.startswith('bar_'): method='plot'+method


plotmeth='contourf' ; histmeth='numpy'


#fext='ic'
#fext=205
#fext='0*5'
fext='*5'
family='disk'
galaxy='galaxy1'
#ow=True
ow=False
#owhist=True
owhist=False
# if method.endswith('scfAm2'):
#     owhist=True
#     ow=True

# import multiprocessing
# print multiprocessing.cpu_count()

pmethsing=None
plotkw={'annotate_arrow':False}
plotfext=None

# Print some info
if method!='status':
    pysim.statsim(runtag,fname0='stdout')
    x=pysim.loadsim(runtag)

# Pars method
if   method=='status'         : pysim.statsim(fname0='stdout')
elif method=='subgadget'      : pysim.files.subrun(x,'gadget',askuser=True)
elif method=='endgadget'      : pysim.files.endrun(x,'gadget')
elif method=='delgadget'      : pysim.files.endrun(x,'gadget',empty=True)
elif method=='subscf'         : pysim.files.subrun(x,'scf')
elif method=='endscf'         : pysim.files.endrun(x,'scf')
elif method=='delscf'         : pysim.files.endrun(x,'scf',empty=True)
elif method=='subbuildgal'    : pysim.files.subrun(x,'buildgal')
elif method=='endbuildgal'    : pysim.files.endrun(x,'buildgal')
elif method=='delbuildgal'    : pysim.files.endrun(x,'buildgal',empty=True)
elif method=='plotsnap'       : x.allsnaps(method='plot',ftype='gadget',fext=fext,
                                           tagstr=plttag,overwrite=ow,askuser=True)
elif method=='plotorbit'      : x.orbit_track(owplot=True,ftype='gadget',cmeth='pot',overwrite=False)
elif method=='orbitparam'     : x.orbit_param(ftype='gadget',cmeth='pot')
elif method.startswith('plot'): x.stdplots(method.split('plot')[-1],family,galaxy,
                                           fext=fext,overwrite=ow,owhist=owhist)
elif method=='print':
    if x['runtyp']=='galsim':
        from pysim.files import buildgal
        bgpar=buildgal.model2param(x.infodict)
        print bgpar.keys()
        for par in ['n','mp','mtot','mscl','rscl','zscl']:
            print par
            for ptyp in ['bulge','disk','halo']:
                print '    {:5s} = {}'.format(ptyp,bgpar[ptyp+'_'+par])
        print 'Mh/Md={}'.format(bgpar['halo_mtot']/bgpar['disk_mtot'])
    y=x.loadsnap('ic',ftype='gadget')
    print y.galaxies()
    print y.families()
    print 'npart = {}'.format(y.header.npart)
    print 'disk  = {}'.format(y.disk['mass'].sum())
    print 'halo  = {}'.format(y.dm['mass'].sum())
elif method=='testtype':
    y=x.loadsnap(fext)
    pysim.nbody.gadget.test_typeidx(y)
elif method=='fithalo':
    y=x.loadsnap('ic')
    fam=pysim.nbody.family.get_family('halo')
    gal=pysim.nbody.galaxy.get_galaxy('galaxy1')
    print pysim.nbody.analysis.langmm.fitprofile(y[fam][gal],'HERNQUIST',plotfile='fitpynbody.png')
elif method=='setlimits'  :
    pro=x.get_profpar(overwrite=True)
    lim=x.get_limits(overwrite=True) ; x.limits=lim
else: raise Exception('Invalid method: {}'.format(method))

toc = time.clock()
print 'Time ellapsed: {}'.format(toc-tic)

# if __name__ == '__main__':
#     main()
