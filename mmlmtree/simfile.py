#!/usr/bin/python
####################################################################################################################################
#
# MEAGAN LANG'S SIMFILE METHODS
#
####################################################################################################################################
import sys,os,shutil,glob,copy,pprint,scipy,math,itertools,re,struct,array
import numpy as np
import matplotlib as mplib
LIST_METHODS=['getfdict']
from mmlutils import *
import main as mmlmtree
import simlist,simstat
LIST_FILEMETHODS=copy.deepcopy(simlist.LIST_METHTYP)
for imtyp in simlist.LIST_METHTYP_NOF: LIST_FILEMETHODS.remove(imtyp)

def main(): return mmlnbody.walk(mtype='file')

####################################################################################################################################
####################################################################################################################################
# HIGH LEVEL METHODS REQUIRING MMLSIM

####################################################################################################################################
# METHOD FOR RUNNING DIFFERENT FILE METHODS
def run(simstr,method=None,**method_kw):
    """
    Provides interface for running different FILE methods
    """
    # Pars input
    if method not in LIST_METHODS: method=mmlio.askselect('Select a valid simfile method:',LIST_METHODS)
    # Initialize default output
    out=None
    # Proceed based on method
    if   method == 'getfdict' : out=files(simstr=simstr,**method_kw)
    else: raise Exception('Invalid method: {}'.format(method))
    # Return output
    return out

####################################################################################################################################
# METHOD TO RETURN DATA TYPES
def get_dtypes():
    """
    Returns a dictionary of data types used
    """
    dtypes={'long' :np.int32,
            'float':np.float32,
            'char' :np.int8}#np.dtype((str,1))}
    return dtypes

####################################################################################################################################
# METHOD TO READ/WRITE FOFCAT FILE
def readgrpcat(fname,grptyp=None,**exkw):
    """
    Reads data from group catelogue files files
    """
    out={} ; dts=get_dtypes()
    grptyp=mmlpars.mml_pars(grptyp,list=['fof','sub'],default='fof')
    varlist=[dict(name='len'   ,dim=1,dtype=dts['long']),
             dict(name='offset',dim=1,dtype=dts['long'])]
    if   grptyp=='fof': varlist+=[dict(name='subs'  ,dim=1,dtype=dts['long'])]
    elif grptyp=='sub': varlist+=[dict(name='parent',dim=1,dtype=dts['long'])]
    with open(fname,mode='rb') as file: 
        out['Ngrp']=np.fromfile(file,dtype=dts['long'],count=1)[0]
        if out['Ngrp']==0:
            for ivar in varlist: out[ivar['name']]=np.zeros((out['Ngrp']*ivar['dim']),dtype=ivar['dtype'])
        else:
            for ivar in varlist: out[ivar['name']]=np.fromfile(file,dtype=ivar['dtype'],count=out['Ngrp']*ivar['dim'])
        for ivar in varlist: 
            if ivar['dim']>1: out[ivar['name']]=out[ivar['name']].reshape((-1,ivar['dim']))
    return out

####################################################################################################################################
# METHOD TO READ GROUP PROPERTY FILES
def readgrpprop(fname,grptyp=None,**exkw):
    """
    Reads data from group property files
    """
    out={} ; dts=get_dtypes()
    grptyp=mmlpars.mml_pars(grptyp,list=['fof','sub'],default='fof')
    varlist=[dict(name='cm'  ,dim=3,dtype=dts['float']),
             dict(name='cmv' ,dim=3,dtype=dts['float']),
             dict(name='mtot',dim=1,dtype=dts['float']),
             dict(name='mgas',dim=1,dtype=dts['float'])]
    with open(fname,mode='rb') as file: 
        out['Ngrp']=np.fromfile(file,dtype=dts['long'],count=1)[0]
        if out['Ngrp']==0:
            for ivar in varlist: out[ivar['name']]=np.zeros((out['Ngrp']*ivar['dim']),dtype=ivar['dtype'])
        else:
            for ivar in varlist: out[ivar['name']]=np.fromfile(file,dtype=ivar['dtype'],count=out['Ngrp']*ivar['dim'])
        for ivar in varlist: 
            if ivar['dim']>1: out[ivar['name']]=out[ivar['name']].reshape((-1,ivar['dim']))
    return out

####################################################################################################################################
# METHOD TO READ GROUP TYPE FILES
def readgrptype(fname,grptyp=None,**exkw):
    """
    Reads data from group type files
    """
    out={} ; dts=get_dtypes()
    grptyp=mmlpars.mml_pars(grptyp,list=['fof','sub'],default='fof')
    varlist=[dict(name='type',dim=1,dtype=dts['char'])]
    with open(fname,mode='rb') as file: 
        out['Npart']=np.fromfile(file,dtype=dts['long'],count=1)[0]
        if out['Npart']==0:
            for ivar in varlist: out[ivar['name']]=np.zeros((out['Npart']*ivar['dim']),dtype=ivar['dtype'])
        else:
            for ivar in varlist: out[ivar['name']]=np.fromfile(file,dtype=ivar['dtype'],count=out['Npart']*ivar['dim'])
        for ivar in varlist: 
            if ivar['dim']>1: out[ivar['name']]=out[ivar['name']].reshape((-1,ivar['dim']))
    return out

####################################################################################################################################
# METHOD TO READ GROUP PARTICLE ID FILES
def readgrpids(fname,grptyp=None,**exkw):
    """
    Reads data from group particle ID files
    """
    out={} ; dts=get_dtypes()
    grptyp=mmlpars.mml_pars(grptyp,list=['fof','sub'],default='fof')
    varlist=[dict(name='ids',dim=1,dtype=dts['long'])]
    with open(fname,mode='rb') as file: 
        out['Npart']=np.fromfile(file,dtype=dts['long'],count=1)[0]
        if out['Npart']==0:
            for ivar in varlist: out[ivar['name']]=np.zeros((out['Npart']*ivar['dim']),dtype=ivar['dtype'])
        else:
            for ivar in varlist: out[ivar['name']]=np.fromfile(file,dtype=ivar['dtype'],count=out['Npart']*ivar['dim'])
        for ivar in varlist: 
            if ivar['dim']>1: out[ivar['name']]=out[ivar['name']].reshape((-1,ivar['dim']))
    return out

####################################################################################################################################
# METHOD TO READ GROUP PARTICLE POSITION FILES
def readgrppos(fname,grptyp=None,**exkw):
    """
    Reads data from group particle position files
    """
    out={} ; dts=get_dtypes()
    grptyp=mmlpars.mml_pars(grptyp,list=['fof','sub'],default='fof')
    varlist=[dict(name='pos',dim=3,dtype=dts['float'])]
    with open(fname,mode='rb') as file: 
        out['Npart']=np.fromfile(file,dtype=dts['long'],count=1)[0]
        if out['Npart']==0:
            for ivar in varlist: out[ivar['name']]=np.zeros((out['Npart']*ivar['dim']),dtype=ivar['dtype'])
        else:
            for ivar in varlist: out[ivar['name']]=np.fromfile(file,dtype=ivar['dtype'],count=out['Npart']*ivar['dim'])
        for ivar in varlist: 
            if ivar['dim']>1: out[ivar['name']]=out[ivar['name']].reshape((-1,ivar['dim']))
    return out

####################################################################################################################################
# METHOD TO READ GROUP PARTICLE VELOCITY FILES
def readgrpvel(fname,grptyp=None,**exkw):
    """
    Reads data from group particle velocity files
    """
    out={} ; dts=get_dtypes()
    grptyp=mmlpars.mml_pars(grptyp,list=['fof','sub'],default='fof')
    varlist=[dict(name='vel',dim=3,dtype=dts['float'])]
    with open(fname,mode='rb') as file: 
        out['Npart']=np.fromfile(file,dtype=dts['long'],count=1)[0]
        if out['Npart']==0:
            for ivar in varlist: out[ivar['name']]=np.zeros((out['Npart']*ivar['dim']),dtype=ivar['dtype'])
        else:
            for ivar in varlist: out[ivar['name']]=np.fromfile(file,dtype=ivar['dtype'],count=out['Npart']*ivar['dim'])
        for ivar in varlist: 
            if ivar['dim']>1: out[ivar['name']]=out[ivar['name']].reshape((-1,ivar['dim']))
    return out

####################################################################################################################################
# METHOD TO READ/WRITE MTREE TABLES
def rwtable(rwid,fname,method=None,**exkw):
    """
    Read/writes interaction lists from/to files
    """
    if   method=='intlist':
        exkw['keylineno']=0
        exkw['typlineno']=1
    elif method=='inttab' :
        exkw['keylineno']=1
        exkw['typlineno']=2
    elif method=='snaptab':
        exkw['keylineno']=2
        exkw['typlineno']=3
    elif method==None     : pass
    else: raise Exception('Invalid method: {}'.format(method))
    out=mmlio.rwtable(rwid,fname,**exkw)
    return out

####################################################################################################################################
# METHOD TO RETURN SIMULATION FILES
def files(simstr=None,verbose=False,**input_kw):
    """
    Generates a dictionary of simulation file names
    """
    # Pars input
    fdict=pars_fileinput(simstr=simstr,**input_kw)
    # Add associated files
    fdict['gdpm']='/home/langmm/code/Gadget/default.param'
    fdict['reds']=os.path.join(fdict['grpdir'],'redshift')
    fdict['ints']=os.path.join(fdict['grpdir'],'fof_to_subhalo_mergers.txt')
    # Add base files
    simdir=os.path.dirname(fdict['grpdir'])
    fdict['snapbase']=os.path.join(simdir,'snapshots','{}mpc{}_'.format(mmlstring.dec2str(simstr['boxsize']),simstr['npart3'])+'{:03d}')
    fdict['tabbase']=os.path.join(fdict['grpdir'],'table_{:03d}.txt')
    fdict['grpbase']=os.path.join(fdict['grpdir'],'groups_{:03d}')
    grplist=['fofcat','subcat','subprop','fofprop','types','ids','pos']
    for igrp in grplist: fdict['grp'+igrp]=fdict['grpbase']+'.'+igrp
    fdict['grpvel']=os.path.join(fdict['rundir'],'vels','groups_{:03d}.vel')
    # Add files that require loading files
    fdict['redlist']=mmlio.rwlist('R',fdict['reds'])
    fdict['snap0'  ]=fdict['redlist'].index(0)-1
    if verbose: mmlio.verbose('z=0 @ snapshot {}'.format(fdict['snap0']))
    fdict['tab0'   ]=fdict['tabbase'].format(fdict['snap0'])
    # Create directories
    for ifmeth in LIST_FILEMETHODS:
        fdict[ifmeth]={'dir' :os.path.join(fdict['rundir'],ifmeth),'pfix':fdict['pfix']}
    # Return output
    return fdict

####################################################################################################################################
####################################################################################################################################
# SUPPORTING METHOD FOR GENERATING FILE NAMES FROM MMLSIM OBJECTS

####################################################################################################################################
# METHOD TO RETURN INPUT FOR FILE METHODS
def pars_fileinput(simstr=None,rundir=None,runtag=None,
                   topdir=None,simdir=None,grpdir=None,
                   memcomp=None,checkaccess=None,**extra_kw):
    """
    Returns a dictionary of default file input variables
    """
    # Get memcomp from rundir
    if isinstance(rundir,str):
        memcompDEF=mmlinfo.dir2compid(rundir)
    else:
        memcompDEF=mmlinfo.hostname2compid()
    # Get runtyp from runtag
    runtyp='mtree'
    # Pars simstr
    if simstr==None:
        simstr={'memcomp':memcompDEF,'runtag':'','grpdir':''}
    memcomp=mmlpars.mml_pars(memcomp,default=simstr['memcomp'],type=str)
    runtag=mmlpars.mml_pars(runtag,default=simstr['runtag'],type=str)
    grpdir=mmlpars.mml_pars(grpdir,default=simstr['grpdir'],type=str)
    if len(runtag)==0: pfix=''
    else             : pfix=runtag+'.'
    # Get directories
    topdirDEF=get_topdir(memcomp=memcomp,checkaccess=checkaccess)
    topdir=mmlpars.mml_pars(topdir,default=topdirDEF,type=str)
    simdirDEF=get_simdir(topdir=topdir)
    simdir=mmlpars.mml_pars(simdir,default=simdirDEF,type=str)
    rundirDEF=get_rundir(runtag=runtag,simdir=simdir)
    rundir=mmlpars.mml_pars(rundir,default=rundirDEF,type=str)
    # Create and return dictionary
    input_kw=dict(rundir=rundir,runtag=runtag,
                  topdir=topdir,simdir=simdir,grpdir=grpdir,
                  memcomp=memcomp,pfix=pfix)
    return input_kw

####################################################################################################################################
# METHOD TO RETURN THE DEFAULT DIRECTORY FOR SIMULATIONS
def get_topdir(memcomp=None,checkaccess=None):
    """
    Returns the default directory for run data given computer info
    """
    # Pars input
    memcomp=mmlpars.mml_pars(memcomp,default=mmlinfo.hostname2compid(),type=str)
    checkaccess=mmlpars.mml_pars(checkaccess,default=False,type=bool)
    # Get top directory
    topdir=mmlinfo.computers(memcomp=memcomp,checkaccess=checkaccess)['simdir']
    # Return directory
    return topdir

####################################################################################################################################
# METHOD TO RETURN THE DEFAULT SIMULATION DIRECTORY
def get_simdir(topdir=None,**topdir_kw):
    """
    Returns the properly assigned simulation type directory
    """
    # Pars input
    topdirDEF=get_topdir(**topdir_kw)
    topdir=mmlpars.mml_pars(topdir,default=topdirDEF,type=str)
    runtyp='mtree'
    # Create simulation type directory
    if len(runtyp)==0:
        simdir=topdir
    else:
        simdir=os.path.join(topdir,runtyp.lower())
    # Return directory
    return simdir

####################################################################################################################################
# METHOD TO RETURN THE DEFAULT RUN DIRECTORY
def get_rundir(runtag=None,simdir=None,**simdir_kw):
    """
    Returns the propertly assigned run directory for a simulation
    """
    # Pars input
    simdirDEF=get_simdir(**simdir_kw)
    simdir=mmlpars.mml_pars(simdir,default=simdirDEF,type=str)
    runtag=mmlpars.mml_pars(runtag,default='',type=str)
    # Create rundir based on runtag
    if len(runtag)==0:
        rundir=simdir
    else:
        rundir=os.path.join(simdir,runtag)
    # Return directory
    return rundir

