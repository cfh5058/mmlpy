#!/usr/bin/python
####################################################################################################################################
#
# MEAGAN LANG'S SIMFILE METHODS
#
####################################################################################################################################
import sys,os,shutil,glob,copy,pprint,scipy,math
import numpy as np
import matplotlib as mplib
import pNbody
LIST_SIMKEYS=['runsheet','runtag']
LIST_RUNKEYS=['runsheet','runtag','simtyp','addruntag','subruntag','mkinfo',
              'inclbgtag','inclgdtag','inclscftag','newstyletag',
              'memcomp','nprocbg','nprocgd','nprocscf','ntot']
SIMLIST_DIR='/home/langmm/analysis'
from mmlutils import *
import main as mmlnbody
import simlist,simfile,simcalc,simplot,mmlbuildgal,mmlgadget,mmlscf,mmlgalfit,icgen

####################################################################################################################################
####################################################################################################################################
# METHODS FOR RETURNING LISTS
def list_runs(runtyp):
    """
    Returns a list of valid runsheets for a specified simulation type
    """
    runlist=simfile.fpar_taglist('icgen',runtyp)
    return runlist

####################################################################################################################################
####################################################################################################################################
# METHODS FOR ACCESSING RUNSHEET DATA

####################################################################################################################################
# METHOD FOR READ/WRITING RUNSHEET DATA
def runsheet_rw(rwid,runsheet,simdict=None,overwrite=None):
    """
    Reads/writes runsheet info from/to a runsheet file
    """
    # Pars input
    rwid=mmlpars.mml_pars(rwid,list=['R','W'])
    runsheet=mmlpars.mml_pars(runsheet,type=str)
    simdict=mmlpars.mml_pars(simdict,default={},type=dict)
    overwrite=mmlpars.mml_pars(overwrite,default=False,type=bool)
    # Get some info
    keylist=LIST_RUNKEYS
    fname=runsheet2fname(runsheet)
    # Read
    if   rwid=='R':
        if overwrite: return simdict
        # Read simdict
        if os.path.isfile(fname):
            simdict_in=mmlio.rwdict('R',fname)
        else:
            simdict_in={}
        # Populated input structure according to overwrite
        for ikey in keylist:
            if ikey not in simdict: simdict[ikey]=None
            if ikey in simdict_in and not overwrite:
                simdict[ikey]=simdict_in[ikey]
        # Return dictionary
        return simdict
    # Write
    elif rwid=='W':
        # Write simdict
        simdict['keylist']=keylist
        for ikey in keylist:
            if simdict[ikey]==None: print '[runsheet_rw: {}] Key {} is none.'.format(simdict['runsheet'],ikey)
        mmlio.rwdict('W',fname,simdict,overwrite=overwrite)
        return
        
####################################################################################################################################
# METHOD FOR RETURNING SIMULATION INFO
def get_runsheet(runtag=None):
    """
    Returns runsheet info for a given runsheet/runtag
    """
    # Get dictionary
    simdict=askuser_siminfo(runtag=runtag)
    # Return runsheet
    return simdict

####################################################################################################################################
# METHOD FOR RETURN MMLSIM OBJECT
def get_mmlsim(runtag=None):
    """
    Returns an mmlsim object for a given runsheet/runtag
    """
    # Get runsheet dictionary
    simdict=get_runsheet(runtag=runtag)
    # Get mmlsim object
    simobj=runsheet2mmlsim(**simdict)
    return simobj

####################################################################################################################################
####################################################################################################################################
# METHODS FOR CREATING HIGH-LEVEL OBJECTS FROM A RUNSHEET

####################################################################################################################################
# METHOD TO CREATE MMLSIM OBJECT FROM RUNSHEET
def runsheet2mmlsim(dontask=False,noinfodict=False,**simdict):
    """
    Creates a mmlsim object from a runsheet
    """
    # Fill in missingkeys & generate mmlsim object
    if not dontask: simdict=askuser_siminfo(**simdict)
    simobj=dict(runtag      = simdict['runtag'],
                runtyp      = simdict['simtyp'],
                addtag      = simdict['subruntag'],
                memcomp     = simdict['memcomp'],
                nprocbg     = simdict['nprocbg'],
                nprocgd     = simdict['nprocgd'],
                nprocscf    = simdict['nprocscf'],
                inclbgtag   = simdict['inclbgtag'],
                inclgdtag   = simdict['inclgdtag'],
                inclscftag  = simdict['inclscftag'],
                newstyletag = simdict['newstyletag'],
                ntot        = simdict['ntot'])
    simobj=mmlnbody.mmlsim(**simobj)
    # Get infodict
    if not noinfodict:
        # Get simulation type object & add it to the mmlsim object
        if   simdict['simtyp']=='galaxy'  : simtypobj=runsheet2galsim(dontask=dontask,**simdict)
        elif simdict['simtyp']=='interact': simtypobj=runsheet2intsim(dontask=dontask,**simdict)
        else: raise Exception('Invalid simulation type: {}'.format(simdict['simtyp']))
        simobj['infodict']=simtypobj
    # Return object
    return simobj

####################################################################################################################################
# METHOD TO CREATE GALSIM OBJECT FROM RUNSHEET
def runsheet2galsim(dontask=False,**simdict):
    """
    Creates a galsim object from a runsheet
    """
    # Fill in missing keys
    if 'mkinfo' not in simdict:
        if dontask: simdict['mkinfo']={}
        else      : simdict=askuser_siminfo(**simdict)
    if not dontask: simdict['mkinfo']=askuser_galsim(runsheet=simdict['runsheet'],**simdict['mkinfo'])
    # Create object
    galobj=dict(model     = simdict['mkinfo']['model'   ],
                haloprof  = simdict['mkinfo']['haloprof'],
                concen    = simdict['mkinfo']['concen'  ],
                ntot      = simdict['mkinfo']['ntot'    ],
                mvir      = simdict['mkinfo']['mvir'    ],
                fgas      = simdict['mkinfo']['fgas'    ])
    galobj=mmlnbody.galsim(**galobj)
    # Return object
    return galobj

####################################################################################################################################
# METHOD TO CREATE INTSIM OBJECT FROM RUNSHEET
def runsheet2intsim(dontask=False,**simdict):
    """
    Creates an intsim object from a runsheet
    """
    # Fill in missing keys
    if 'mkinfo' not in simdict:
        if dontask: simdict['mkinfo']={}
        else      : simdict=askuser_siminfo(**simdict)
    if not dontask: simdict['mkinfo']=askuser_intsim(runsheet=simdict['runsheet'],**simdict['mkinfo'])
    # Create object
    intobj=dict(galsim1 = runsheet2galsim(runsheet=simdict['mkinfo']['run1']),
                galsim2 = runsheet2galsim(runsheet=simdict['mkinfo']['run2']),
                rperi   = simdict['mkinfo']['rperi'],
                ecc     = simdict['mkinfo']['ecc'  ],
                incl    = simdict['mkinfo']['incl' ],
                ntot    = simdict['mkinfo']['ntot' ])
    intobj=mmlnbody.intsim(**intobj)
    # Return object
    return intobj

####################################################################################################################################
# METHOD TO ALLOW USER TO CREATE MMLSIM OBJECT FROM COMMAND LINE
def askuser_siminfo(owflag=None,**simdict):
    """
    Generates a runsheet by asking the user for information not supplied
    """
    # Set constants
    keylist=LIST_RUNKEYS
    slist_gal=list_runs('galaxy')
    slist_int=list_runs('interact')
    slist_tot=slist_gal+slist_int
    simtypLIST=simlist.LIST_RUNTYPS
    # Pars input
    simdict=mmlpars.mml_pars(simdict,default={},type=dict)
    owflag=mmlpars.mml_pars(owflag,default=False,type=bool)
    # Fill missing values with Nones
    for ikey in keylist:
        if ikey not in simdict: simdict[ikey]=None
    # Add simtyp key if it's missing
    if simdict['runtyp'] not in simtypLIST: simdict['runtyp']=runtag2runtyp(simdict['runtag'])
    slist=list_runs(simdict['runtyp'])
    # Get fname key if it's missing
    if not isinstance(simdict['runsheet'],str):
        runsheet=mmlio.askselect('Select a simulation runsheet:',slist+['newfile'])
        while runsheet=='newfile':
            runsheet=mmlio.askquest('Enter name for new simulation info file:',default='newfile',dtype=str)
            if runsheet in slist:
                owflag=mmlio.yorn('That file already exists. Overwrite?')
                if not owflag: runsheet='newfile'
        simdict['runsheet']=runsheet
    # Load run info from file
    if os.path.isfile(runsheet2fname(simdict['runsheet'])) and not owflag:
        simdict=runsheet_rw('R',simdict['runsheet'],simdict=simdict,overwrite=owflag)
    # Determine if anything should be written
    if any([simdict[ikey]==None for ikey in keylist]): writeflag=True
    else: writeflag=False
    # Get simulation type dependent info
    if simdict['mkinfo'      ]==None: simdict['mkinfo'     ]={}
    simdict['mkinfo']=askuser_mkinfo(simdict['simtyp'],runsheet=simdict['runsheet'],**simdict['mkinfo'])
    # Get general info
    simdict['fname']=runsheet2fname(simdict['runsheet'])
#    if simdict['fname'       ]==None: simdict['fname'      ]=runsheet2fname(simdict['runsheet'])
    if simdict['addruntag'   ]==None: simdict['addruntag'  ]=mmlio.askquest('[{}] String to append to end of runtag?'.format(simdict['runsheet']),default='',dtype='str')
    if simdict['subruntag'   ]==None: simdict['subruntag'  ]=mmlio.askquest('[{}] String identifying subset of run files?'.format(simdict['runsheet']),default='',dtype='str')
    if simdict['inclbgtag'   ]==None: simdict['inclbgtag'  ]=mmlio.yorn('[{}] Include Buildgal # of processors tag?'.format(simdict['runsheet']))
    if simdict['inclgdtag'   ]==None: simdict['inclgdtag'  ]=mmlio.yorn('[{}] Include Gadget # of processors tag?'.format(simdict['runsheet']))
    if simdict['inclscftag'  ]==None: simdict['inclscftag' ]=mmlio.yorn('[{}] Include SCF # of processors tag?'.format(simdict['runsheet']))
    if simdict['newstyletag' ]==None: simdict['newstyletag']=mmlio.yorn('[{}] Use new style for tagging files?'.format(simdict['runsheet']))
    if simdict['memcomp'     ]==None: simdict['memcomp'    ]=mmlio.askquest('[{}] VPAC machine to use for storage?'.format(simdict['runsheet']),default='vpac38',dtype='str')
    if simdict['nprocbg'     ]==None: simdict['nprocbg'    ]=mmlio.askquest('[{}] Number of BUILDGAL processors?'.format(simdict['runsheet']),default=1,dtype='int')
    if simdict['nprocgd'     ]==None: simdict['nprocgd'    ]=mmlio.askquest('[{}] Number of GADGET processors?'.format(simdict['runsheet']),default=32,dtype='int')
    if simdict['nprocscf'    ]==None: simdict['nprocscf'   ]=mmlio.askquest('[{}] Number of SCF processors?'.format(simdict['runsheet']),default=50,dtype='int')
    if simdict['runtag'      ]==None: simdict['runtag'     ]=runsheet2runtag(**simdict)
    if simdict['ntot'        ]==None: simdict['ntot'       ]=simdict['mkinfo']['ntot']
    # Save runsheet
    if writeflag:
        mmlio.verbose('Saving runsheet to file: '+simdict['fname'])
        runsheet_rw('W',simdict['runsheet'],simdict=simdict,overwrite=True)
    # Add new file to the list
    simlist_rw('W',simdict['simtyp'],**simdict)
    # Return dictionary
    return simdict

####################################################################################################################################
# METHOD TO ALLOWER USER TO CREATE INTSIM/GALSIM OBJECTS FROM COMMAND LINE
def askuser_mkinfo(simtyp1,runsheet=None,cleankeys=None,**mkinfo):
    """
    Provides command line access to creation of intsim/galsim objects
    """
    # Set constants
    simtypLIST=simlist.LIST_RUNTYPS
    # Pars input
    simtyp1=mmlpars.mml_pars(simtyp1,list=simtypLIST)
    runsheet=mmlpars.mml_pars(runsheet,default='',type=str)
    # Proceed based on simulation type
    if   simtyp1=='galaxy'  : mkinfo=askuser_galsim(runsheet=runsheet,cleankeys=cleankeys,**mkinfo)
    elif simtyp1=='interact': mkinfo=askuser_intsim(runsheet=runsheet,cleankeys=cleankeys,**mkinfo)
    else: raise Exception('Invalid simulation type: {}'.format(simtyp1))
    # Return make dictionary
    return mkinfo

####################################################################################################################################
# METHOD TO ALLOW USER TO CREATE GALSIM OBJECT FROM COMMAND LINE
def askuser_galsim(runsheet=None,cleankeys=None,**mkinfo):
    """
    Provides command line acces to creation of GALSIM objects
    """
    # Pars input
    runsheet=mmlpars.mml_pars(runsheet,default='',type=str)
    cleankeys=mmlpars.mml_pars(cleankeys,default=False,type=bool)
    inclextra=not cleankeys
    # Get info on galaxy from user
    mkinfo=simfile.askfpar('icgen','galsim',fpar=mkinfo,inclextra=inclextra,noload=True)
    # Return make info
    return mkinfo
    
####################################################################################################################################
# METHOD TO ALLOW USER TO CREATE INTSIM OBJECT FORM THE COMMAND LINE
def askuser_intsim(runsheet=None,cleankeys=None,**mkinfo):
    """
    Provides command line acces to creation of INTSIM objects
    """
    # Pars input
    runsheet=mmlpars.mml_pars(runsheet,default='',type=str)
    cleankeys=mmlpars.mml_pars(cleankeys,default=False,type=bool)
    inclextra=not cleankeys
    # Get info on galaxy from user
    mkinfo=simfile.askfpar('icgen','intsim',fpar=mkinfo,inclextra=inclextra,noload=True)
    # Get number of particles
    if mkinfo['ntot'       ]==None:
        simdict1=askuser_siminfo(runsheet=mkinfo['run1'])
        simdict2=askuser_siminfo(runsheet=mkinfo['run2'])
        mkinfo['ntot'       ]=simdict1['ntot']+simdict2['ntot']
    # Return dictionary
    return mkinfo

####################################################################################################################################
####################################################################################################################################
# SIMLIST METHODS
def simlist_typefile(simtyp):
    """
    Returns the name of the simulation list file for a given simulation type
    """
    simdir=SIMLIST_DIR
    if   simtyp=='galaxy'  : fname=os.path.join(simdir,'simlist_galaxy')
    elif simtyp=='interact': fname=os.path.join(simdir,'simlist_interact')
    else: raise Exception('Invalid simtyp: {}'.format(simtyp))
    return fname

def simlist_rw(rwid,simtyp1,**listkw):
    """
    Edits the list of simulations for a given type
    """
    # Set constants
    keyLIST=LIST_SIMKEYS
    # Pars input
    rwid=mmlpars.mml_pars(rwid,list=['R','W','D'])
    simtyp1=mmlpars.mml_pars(simtyp1,list=simlist.LIST_RUNTYPS)
    # Get info
    fname=simlist_typefile(simtyp1)
    # Read in existing data
    if os.path.isfile(fname):
        simdict=mmlio.rwtable('R',fname)
    else:
        simdict={ikey:[] for ikey in keyLIST}
    # Read
    if   rwid=='R': return simdict
    # Write
    elif rwid=='W':
        flagdict={}
        flaglist=[]
        indxdict={}
        for ikey in keyLIST:
            flagdict[ikey]=(listkw[ikey] in simdict[ikey])
            flaglist.append(flagdict[ikey])
            if flagdict[ikey]: indxdict[ikey]=simdict[ikey].index(listkw[ikey])
            else:              indxdict[ikey]=-1
        # Old
        if all(flaglist):
            key0=keyLIST[0]
            for ikey in keyLIST[1:]:
                if not indxdict[ikey]==indxdict[key0]:
                    raise Exception('{} {} is already listed, but not with {} {}'.format(ikey,listkw[ikey],key0,listkw[key0]))
        # Mixed
        elif any(flaglist):
            for ikey in keyLIST:
                if flagdict[ikey]: print '{} {} is already listed.'.format(ikey,listkw[ikey])
            raise Exception('Above keywords already listed.')
        # New
        else:
            mmlio.verbose('Adding runsheet to the list: {}'.format(listkw['runsheet']))
            for ikey in keyLIST:
                simdict[ikey].append(listkw[ikey])
            mmlio.rwtable('W',fname,tabdict=simdict,overwrite=True)
        return
    # Remove entry
    elif rwid=='D':
        for ikey in keyLIST:
            if ikey in listkw:
                key0=ikey
                break
        if listkw[key0] in simdict[key0]:
            rmidx=simdict['key0'].index(listkw[key0])
            for ikey in keyLIST:
                simdict[ikey].remove(simdict[ikey][rmidx])
        return
    # Error handling
    else:
        raise Exception('There is not an option for an rwid of {} in runsheet.simlist_rw'.format(rwid))
    return
            
####################################################################################################################################
####################################################################################################################################
# METHOD FOR CHANGING BETWEEN DIFFERENT VARIABLES
def runtag2fname(runtag):
    """
    Converts runsheets into file names
    """
    return simfile.ftag2fname('icgen','mmlsim',runtag)
def fname2runtag(fname):
    """
    Converts runsheets file names into runsheets
    """
    return simfile.fname2ftag('icgen','mmlsim',fname)
def runtag2runtyp(runtag):
    """
    Determines the simulation type from the runsheet
    """
    # Get lists
    runlist_gal=simfile.fpar_taglist('icgen','galsim')
    runlist_int=simfile.fpar_taglist('icgen','intsim')
    runtypLIST=simlist.LIST_RUNTYPS
    # Determine what type of simulation it is
    if   runtag in runlist_gal: runtyp='galsim'
    elif runtag in runlist_int: runtyp='intsim'
    else: runtyp=mmlio.askselect('[{}] Select a simulation type:'.format(runsheet),runtypLIST)
    # Return runtype
    return runtyp
def runtag2fdict(runtag,**extra_kw):
    """
    Uses runsheet to call simfile.files
    """
    extra_kw['simstr']=runtag2mmlsim(runtag=runtag,noinfodict=True)
    fdict=simfile.files(**extra_kw)
    return fdict
def runtag2flist(runtag,method,**extra_kw):
    """
    Uses runsheet to call simfile.get_fiellist
    """
    extra_kw['simstr']=runtag2mmlsim(runtag=runtag,noinfodict=True)
    flist=simfile.get_filelist(method,**extra_kw)
    return flist

    
    
