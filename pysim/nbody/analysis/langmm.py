
"""
snapshot
========

This module implements the  :class:`~pynbody.snapshot.SimSnap` class which manages and stores snapshot data.
It also implements the :class:`~pynbody.snapshot.SubSnap` class (and relatives) which 
represent different views of an existing :class:`~pynbody.snapshot.SimSnap`.

"""

from .. import plot as nbplot
from .. import array
from .. import family, galaxy
from .. import units
from . import mapping,angmom

import numpy as np
import copy
import exceptions
import sys
import warnings
import threading
import re
import pprint

####################################################################################################################################
# BAR PROP
def bar(sim,method,verbose=False,
        filename=None,overwrite=False,histfile=None,owhist=False,
        galaxy=None,family=None,cmeth=None,
        rlim=(0.01,30.0),rscl='log',nbins=100,**kws):
    """Determine bar properties"""
    import mmlastro.bar as fbar
    from mmlutils import mmlio
    # Properties from ellipse fits
    if method=='ellipse':
        # Get histogram of mass distribution
        hist=mapping.hist2d(sim,'x','y',nbins=nbins,
                            xlim=(-rlim[1],rlim[1]),xscl='lin',
                            ylim=(-rlim[1],rlim[1]),yscl='lin',
                            wvar='mass',wlim=(1.0e4,1.0e10),wscl='log',
                            galaxy=galaxy,family=family,cmeth=cmeth,
                            units=['km s**-1','kpc','Msol'],
                            filename=histfile,overwrite=owhist)
        # Convert histogram to data used by ellipse
        kws.update(**fbar.parshist(method,hist))
        kws.update(Na=kws.pop('nellip',100),amin=rlim[0],amax=rlim[1],ascl=rscl,avgwin=4.*rlim[1]/nbins,
                   title=nbplot.generic.val2str(sim.properties['time'],label='t',units='Gyr',latex=False))
        # Fit bar with ellipses
        out=fbar.ellipse(**kws)
    # Properties from m=2 part of scf decomposition
    elif method=='scfAm2':
        # Get values
        ampthresh=kws.pop('ampthresh',0.001)
        # Get histogram of m=2 mode in phi vs rad
        hist=mapping.hist2d(sim,'R','modaz',nbins=nbins,
                            xlim=rlim,xscl=rscl,
                            ylim=(0.,np.pi),yscl='lin',
                            wvar='Am2',wlim=(ampthresh/10.,1.),wscl='symlog',
                            avar='N',
                            galaxy=galaxy,family=family,cmeth=cmeth,
                            units=['km s**-1','kpc','Msol'],
                            filename=histfile,overwrite=owhist)
        # Convert histogram to data used by ellipse
        kws.update(**fbar.parshist(method,hist))
        kws.update(title=nbplot.generic.val2str(sim.properties['time'],label='t',units='Gyr',latex=False))
        # Fit bar from m=2 mode
        out=fbar.scfAm2(**kws)
    # Error
    else: raise Exception('Invalid method: {}'.format(method))
    # Add time and fext
    out['fext']=sim.fext
    out['time']=float(sim.properties['time'].in_units('Gyr'))
    # Add entry to table
    if isinstance(filename,str):
        mmlio.rwtable('A',filename,out,overwrite=overwrite,uniquekey='fext')
    # Return output
    if verbose: pprint.pprint(out)
    return out

####################################################################################################################################
# METHOD FOR ALIGNING SELECT PART OF SIMULATION
def align(sim0,galaxy=None,family=None,cmeth=None,pNbody=False):
    """Select particles by galaxy/family and align angular momentum with z"""
    # Select galaxy/family to center
    if cmeth:
        cgal=galaxy if galaxy else 'galaxy1'
        cfam=family if family else None
    # Use pNbody to select particles and center
    if pNbody:
        if cmeth: sim0.center(cmeth,centertyp=cfam,centergal=cgal)
        sim1=sim0.select(family) if family else sim0
        sim2=sim1.select(galaxy) if galaxy else sim1
        sim=sim2
    # Use pynbody to select particles and center
    else:
        from ..family import get_family
        from ..galaxy import get_galaxy
        sim1=sim0[get_family(family)] if family else sim0
        sim2=sim1[get_galaxy(galaxy)] if galaxy else sim1
        sim=sim2
        
        if cmeth: angmom.faceon(sim1[get_galaxy(cgal)],mode=cmeth)
    # Return adjusted view
    return sim
    
####################################################################################################################################
# METHOD FOR MERGING TWO SIM OBJECTS
def merge(sim1,sim2,N1=None,N2=None,N1min=0,N2min=0,repcut=0.,typList=None,
          keepTypeRatio=False,pNbody=False,warn=False):
    """
    Merge two simulations removing repeats
        N1   : # of particles that you would like to select (adjusts if warn set)
        N1min: minimum # of particles that should be selected
    """
    from mmlutils import mmlio,mmlmath
    import string
    # Pars input
    if typList is None:
        if keepTypeRatio: typList=sim1.families()
        else            : typList=set(sim1.families()+sim2.families())
    # Get list of unique positions
    idxuni2,idxuni1=index_norepeats(sim2,sim1,cut=repcut,pNbody=pNbody)
    N1uni=len(idxuni1)
    N2uni=len(idxuni2)
    mmlio.verbose('    N (unique) = {: 8d} ({: 8d} + {: 8d})'.format(N1uni+N2uni,N1uni,N2uni))
    # Define warning string
    def warnstr(simnum,Nuni,Nreq,tstr=None,w=False,xstr=None):
        if tstr is None: tstr='overall'
        else           : tstr=string.ljust(tstr,5)
        wstr='    Sim {} {} has {: 8d} unique positions, but {: 8d} required.'.format(simnum,tstr,Nuni,Nreq)
        if xstr: wstr+=' '+xstr
        if w: mmlio.verbose(wstr)
        else: raise Exception(wstr)
    # Determine if there are enough unique positions overall
    if N1 is None: N1=N1uni
    if N2 is None: N2=N2uni
    if   N1uni<N1   : warnstr(1,N1uni,N1   ,w=warn)
    elif N1uni<N1min or N1uni==0: 
        warnstr(1,N1uni,N1min,w=warn,xstr='Use Sim 2')
        return sim2
    if   N2uni<N2   : warnstr(2,N2uni,N2   ,w=warn)
    elif N2uni<N2min or N2uni==0: 
        warnstr(2,N2uni,N2min,w=warn,xstr='Use Sim 1')
        return sim1
    # Loop over families
    rng1=np.arange(len(sim1),dtype=int)
    rng2=np.arange(len(sim2),dtype=int)
    tidx1={} ; tidxuni1={}
    tidx2={} ; tidxuni2={}
    valid=False
    while not valid:
        # Stop if invalid N1/N2
        if N1<N1min or N1==0:
            warnstr(1,N1,N1min,w=warn,xstr='Use Sim 2')
            return sim2
        if N2<N2min or N2==0:
            warnstr(2,N2,N2min,w=warn,xstr='Use Sim 1')
            return sim1
        # Reset validity and indicies for new N1/N2
        valid=True
        rep1=np.zeros(len(sim1),dtype=bool)
        rep2=np.zeros(len(sim2),dtype=bool)
        for t in typList:
            # Get indicies for each type
            if t not in tidx1 or t not in tidx2:
                if pNbody:
                    tslc1=sim1.get_index(ptyp=t)
                    tslc2=sim2.get_index(ptyp=t)
                else:
                    tslc1=sim1._get_family_slice(t)
                    tslc2=sim2._get_family_slice(t)
                tidx1[t]=rng1[tslc1]
                tidx2[t]=rng2[tslc2]
                # Combine to total type index
                tidxuni1[t]=mmlmath.intersect1d(idxuni1,tidx1[t],sort=False)
                tidxuni2[t]=mmlmath.intersect1d(idxuni2,tidx2[t],sort=False)
            # Get number of particles in each type
            Ntyp1uni=len(tidxuni1[t])
            Ntyp2uni=len(tidxuni2[t])
            Ntyp1   =N1   *len(tidx1[t])/len(sim1)
            Ntyp1min=N1min*len(tidx1[t])/len(sim1)
            if keepTypeRatio: 
                Ntyp2   =N2   *len(tidx1[t])/len(sim1)
                Ntyp2min=N2min*len(tidx1[t])/len(sim1)
            else            : 
                Ntyp2   =N2   *len(tidx2[t])/len(sim2)
                Ntyp2min=N2min*len(tidx2[t])/len(sim2)
            if Ntyp1==0 and Ntyp2==0: continue
            # Get random selection of unique indices if there are enough particles
            if Ntyp1uni>=Ntyp1 and Ntyp2uni>=Ntyp2:
                if Ntyp1uni>Ntyp1:
                    np.random.shuffle(tidxuni1[t])
                    tidxuni1[t]=tidxuni1[t][:Ntyp1]
                if Ntyp2uni>Ntyp2:
                    np.random.shuffle(tidxuni2[t])
                    tidxuni2[t]=tidxuni2[t][:Ntyp2]
                # Add to total sim index
                if Ntyp1!=0: rep1[tidxuni1[t]]=True
                if Ntyp2!=0: rep2[tidxuni2[t]]=True
            # Otherwise warn and adjust N1 & N2
            else:
                valid=False
                if Ntyp1uni<Ntyp1:
                    warnstr(1,Ntyp1uni,Ntyp1,tstr=t,w=warn)
                    N1=Ntyp1uni*len(sim1)/len(tidx1[t])
                    break
                if Ntyp2uni<Ntyp2:
                    warnstr(2,Ntyp2uni,Ntyp2,tstr=t,w=warn)
                    if keepTypeRatio: N2=Ntyp2uni*len(sim1)/len(tidx1[t])
                    else            : N2=Ntyp2uni*len(sim2)/len(tidx2[t])
                    break
    # Isolate sims
    if pNbody: 
        sim1_fin=sim1.selectc(rep1)
        sim2_fin=sim2.selectc(rep2)
    else     : 
        sim1_fin=sim1[rep1]
        sim2_fin=sim2[rep2]
    # Add and return
    sim=sim1_fin+sim2_fin
    mmlio.verbose(    'N (no rep) = {: 8d}'.format(len(sim)))
    return sim

####################################################################################################################################
# METHODS RETURNING INFOMATION ON REPEATED POSITIONS
def index_norepeats(sim,sim0=None,cut=0.,pNbody=False,printrep=False):
    """
    Returns the indicies of unique particle positions 
    """
    from mmlutils import mmlmath,mmlio
    # Recover positions
    if pNbody:
        pos=sim.pos
        if sim0: pos0=sim0.pos
    else:
        pos=sim['pos']
        if sim0: pos0=sim0['pos']
    if sim0: pos=np.append(pos0,pos,axis=0)
    # Get indices without repeats
    idx=mmlmath.index_norepeats(pos,cut=cut,axis=0)
    # Only include indices from sim
    if sim0:
        idx0=idx[idx< len(sim0)]
        idx1=idx[idx>=len(sim0)]-len(sim0)

        # Print some repeated postions
        if printrep and len(idx1)<len(sim):
            boolRep=np.ones(len(sim),dtype=bool)
            boolRep[idx1]=False
            if pNbody: print sim.pos[boolRep,:]
            else     : print sim['pos'][boolRep,:]
            if not mmlio.yorn('Continue?'):
                raise Exception('User stopped run')

        return idx1,idx0
    else:
        return idx

    # a=np.lexsort(sim['pos'].T)
    # diff=np.ones(len(sim))
    # diff[1:]=np.sqrt((np.diff(sim['pos'][a,:],axis=0)**2).sum(axis=1))
    # idx=np.sort(a[diff>cut])
    # return idx

def count_repeats(sim,**kws):
    """
    Count the number of repeated particle positions
    """
    idx_ntot=index_norepeats(sim,**kws)
    N=len(sim)-len(idx_ntot)
    return N

####################################################################################################################################
# METHODS FOR FITTING PROFILES
def fitprofile(sim,profid,coord=None,nbin=100,cmeth=None,rmax=None,**profkw):
    """
    Fit a 1D profile to the particles
    """
    from mmlastro import mmlprofiles
    # Pars input
    profid=profid.upper()
    if profid not in mmlprofiles.LIST_PROFILES:
        raise Exception('Invalid profid: {}'.format(profid))
    if   profid in mmlprofiles.LIST_PROFILES_SPH: coordlist=['r']
    elif profid in mmlprofiles.LIST_PROFILES_CYL: coordlist=['R']
    if profid=='HIYASHI': coordlist.append('z')
    if coord not in coordlist:
        if coord==None: coord=coordlist[0]
        else: raise Exception('Invalid coordinate: {}'.format(coord))
    # Center
    from . import angmom
    angmom.faceon(sim,mode=cmeth)
    # Get histogram
    if not rmax: rmax=sim[coord].max()
    x=sim[coord]
    m=sim['mass']
    mhist,xbins=np.histogram(x[x<=rmax],bins=nbin,weights=m[x<=rmax])
    # Select data to fit
    menc=np.cumsum(mhist)
    xenc=xbins[1:]
    # Fit data
    profkw[coord.lower()]=xenc
    profpar=mmlprofiles.fit(profid,menc,**profkw)
    profpar['xunits']=x.units
    profpar['munits']=m.units
    # Return parameters
    return profpar

####################################################################################################################################
# METHOD FOR PUTING GALAXY ON AN ORBIT
def add_orbit(sim,M1_msol,M2_msol,R1_pc=None,R2_pc=None,rinit=1.,rperi=0.5,ecc=1.5,incl=0.,**exkw):
    from mmlastro import mmlprofiles,mmlconst
    from mmlutils import mmlmath
    # Set constants
    #gconst=4.302e-3 # pc (km/s)**2 / Msol
    exkw['units']='galaxy' # pc/Msol
    gconst=mmlconst.main('phys',subtag=exkw['units'])['G']
    # Get virial radii
    if not R1_pc: R1_pc=mmlprofiles.mvir2rvir(M1_msol,**exkw)
    if not R2_pc: R2_pc=mmlprofiles.mvir2rvir(M2_msol,**exkw)
    # Give input units based on virial radii
    rtot_pc=R1_pc+R2_pc
    rini_pc=rinit*rtot_pc
    rper_pc=rperi*rtot_pc
    # Determine intial position/velocity of orbit
    xorb_pc=np.array([rini_pc,0.,0.])
    vorb_kms,phiorb=mmlmath.orbit(M1_msol,M2_msol,rini_pc,rper_pc,ecc,gconst)
    # Shift particle positions and velocities to be on orbit
    sim['pos']+=array.SimArray(xorb_pc ,'pc'      ).in_units(sim['pos'].units)
    sim['vel']+=array.SimArray(vorb_kms,'km s**-1').in_units(sim['vel'].units)
    # Incline the orbit
    sim.rotate_y(-180*incl/np.pi)
    # print sim.rotate(180*incl/np.pi,axis=[0,-1,0],ret_mat=True)
    # print sim.rotate(-180*incl/np.pi,axis='y',ret_mat=True)
    # print sim.rotate_y(-180*incl/np.pi,ret_mat=True)
    # Force pericenter to negative x-axis
    sim.rotate_z(180*phiorb/np.pi)
    #sim.rotate(180*phiorb/np.pi,axis=[0,0,1])
    # Print some info
    xorb=np.array(array.SimArray(xorb_pc ,'pc'      ).in_units('kpc'     ))
    vorb=np.array(array.SimArray(vorb_kms,'km s**-1').in_units('km s**-1'))
    verbose=True
    if verbose:
        if ecc == 0:              print 'CIRCULAR orbit'
        elif ecc > 0 and ecc < 1: print 'ELLIPTICAL orbit'
        elif ecc == 1:            print 'PARABOLIC orbit'
        elif ecc > 1:             print 'HYPERBOLIC orbit'
        print '    M1    = {:4.2e} Msol'.format(M1_msol)
        print '    M2    = {:4.2e} Msol'.format(M2_msol)
        print '    R1    = {:5f} kpc '.format(float(array.SimArray(R1_pc,'pc').in_units('kpc')))
        print '    R2    = {:5f} kpc '.format(float(array.SimArray(R2_pc,'pc').in_units('kpc')))
        print '    rinit = {:5f} kpc '.format(float(array.SimArray(rini_pc,'pc').in_units('kpc')))
        print '    rperi = {:5f} kpc '.format(float(array.SimArray(rper_pc,'pc').in_units('kpc')))
        print '    ecc   = {:2.2f}     '.format(ecc)
        print '    pos   = [{:5f},{:5f},{:5f}] kpc '.format(xorb[0],xorb[1],xorb[2])
        print '    vel   = [{:5f},{:5f},{:5f}] km/s'.format(vorb[0],vorb[1],vorb[2])
        print '    phi   = {:4.1f} deg '.format(phiorb*(180./np.pi))
        print '    theta = {:4.1f} deg '.format(incl*(180./np.pi))
    # Return control
    return

####################################################################################################################################
# KINEMATICALLY HEAT THE DISK BY SOME FACTOR
def heatdisk(sim,heatfact):
    print '[nbody.analysis.langmm.heatdisk] Check that the disk is being heated.'
    print 'before: {}'.format(sim['disk']['vel'][0])
    sim['disk']['vel']*=heatfact
    print 'after:  {}'.format(sim['disk']['vel'][0])

####################################################################################################################################
# DETERMINE THE STABILITY OF THE SIMULATION
def get_stability(sim,method):
    from mmlastro import mmldisk
    if   method=='op73'    : crit=mmldisk.stability_op73(sim['Erot'],sim['Ep'])
    elif method=='efstat82':
      G=sim.get_gconst()
      Vpeak=sim.get_Vpeak(ptyp=ptyp,pgal=pgal,idxtot=idxtot)
      Md=sim['mass'].sum()
      Rd=4. # CHANGE THIS!!!                                                                                                                                                         
      crit=mmldisk.stability_efstat82(Vpeak,Md,Rd,G)
    elif method=='toomre81':
      G=sim.get_gconst()
      R=4.
      sigr=sim.get_sigmaR(ptyp=ptyp,pgal=pgal,idxtot=idxtot)
      crit=None
#      kappa                                                                                                                                                                         
#      crit=mmldisk.stability_toomre81(R,sigr,kappa,rhoS,G)                                                                                                                          
    else: raise Exception('Invalid stability criterion: {}'.format(method))
    return crit


####################################################################################################################################

