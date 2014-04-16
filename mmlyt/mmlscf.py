#!/usr/bin/python
####################################################################################################################################
#
# MEAGAN LANG'S SCF METHODS
#
####################################################################################################################################
import sys,os,shutil,glob,copy,pprint,scipy,math
import numpy as np
import matplotlib as mplib
import pNbody
LIST_METHODS=['makefiles','submit','retrieve','clean']
from mmlutils import *
import main as mmlnbody
import simlist,simcalc,simplot,simfile,mmlbuildgal,mmlgadget,mmlgalfit,icgen

def main():
    """
    Provides command line acces to SCF methods
    """
    out = mmlnbody.walk(mtype='scf')
    return out

####################################################################################################################################
####################################################################################################################################
# METHODS FOR GETTING LISTS
def dict_ftype():
    """
    Returns a dictionary of files listed by type 
    """
    fdict={}
    return fdict
def dict_fgroup():
    """
    Returns a dictionary of file types listed by group
    """
    fdict={'clean' :['longin','longout'],
           'input' :['shortin','longin'],
           'output':['shortout','longout']}
    return fdict

####################################################################################################################################
####################################################################################################################################
# HIGH LEVEL METHODS REQUIRING MMLSIM

####################################################################################################################################
# METHOD FOR RUNNING DIFFERENT SCF METHODS
def run(simstr,method,**method_kw):
    """
    Provides interface for running different SCF operations
    """
    # Set constants
    methLIST=LIST_METHODS
    typeLIST=simlist.LIST_PTYPS
    # Pars input
    method=mmlpars.mml_pars(method,list=methLIST)
    # Initialize default output
    out=None
    # Proceed based on method
    if   method=='makefiles':
        method_kw['ptype']=mmlio.askselect('Files for which particle type should be created?',typeLIST)
        method_kw['typList']=[method_kw['ptype']]
        mkfiles(simstr,**method_kw)
    elif method=='submit'   :
        method_kw['ptype']=mmlio.askselect('Run for which particle type should be submitted?',typeLIST)
        subrun(simstr,**method_kw)
    elif method=='retrieve' :
        method_kw['ptype']=mmlio.askselect('Run for which particle type should be retrieved?',typeLIST)
        endrun(simstr,**method_kw)
    elif method=='clean'    :
        method_kw['ftype']='clean'
        method_kw['memcomp']='accre'
        simfile.rmfiles('scf',simstr=simstr,**method_kw)
    else:
        raise Exception('Option for method {} needs to be added to mmlscf.run.'.format(method))
    # Return output
    return out

####################################################################################################################################
# METHOD PERFORM ACTIONS TO SETUP AND RUN SCF ON ACCRE
def subrun(simstr,ptype=None,**extra_kw):
    """
    Sets up and runs SCF on ACCRE
    """
    # Pars input
    ptype=mmlpars.mml_pars(ptype,default='disk',type=str)
    extra_kw['typList']=[ptype]
    files_accre=simstr.mkfiledict(memcomp='accre',checkaccess=True)['scf']
    # Update input files
    if mmlio.yorn('Create new input files?'): mkfiles(simstr,**extra_kw)
    # Move input files to ACCRE
    if mmlio.yorn('Move new input files to ACCRE?'): putfiles(simstr,**extra_kw)
    # Make directories
    if mmlio.yorn('Make SCF output directories?'):
        for ityp in extra_kw['typList']:
            if not mmlio.yorn('    For {} particles?'.format(ityp)): continue
            mmlfiles.mkdirs(files_accre[ityp]['output']['dir'])
            dkeylist=files_accre[ityp]['output']['subdirs']
            dirlist=[files_accre[ityp]['output'][idirkey] for idirkey in dkeylist]
            mmlfiles.mkdirs(dirlist)
    # Run SCF
    if mmlinfo.hostname2compid()=='bender':
        if mmlio.yorn('Run SCF?'):
            pbsfile=files_accre[ptype]['input']['pbs']
            schfile=files_accre[ptype]['input']['sched']
            simstr.subpbs('scf',ptype=ptype,pbsfile=pbsfile,schfile=schfile)
    return

####################################################################################################################################
# METHOD TO PERFORM ACTIONS TO RETRIEVE AND CLEAN-UP A RUN
def endrun(simstr,ptype=None,**extra_kw):
    """
    Retrieves and cleans up results from an SCF run
    """
    # Pars input
    ptype=mmlpars.mml_pars(ptype,default='disk',type=str)
    # Retrieve results
    if mmlio.yorn('Recover SCF results from ACCRE?'): getfiles(simstr,typList=[ptype])
    # Clean up results
    if mmlio.yorn('Clean up SCF results on ACCRE?' ): simfile.rmfiles('scf',ftype='clean',memcomp='accre',simstr=simstr)
    return

####################################################################################################################################
# METHOD TO MOVE INPUT FILES TO ACCRE
def putfiles(simstr,typList=None,**extra_kw):
    """
    Moves files to ACCRE in preparation to running SCF code
    """
    # Pars input
    typListDEF=simlist.LIST_PTYPS
    typList=mmlpars.mml_pars(typList,default=typListDEF,type=list)
    # Move files
    if mmlio.yorn('Move SCF input files to ACCRE?'     ): simfile.mvfiles_accre(simstr,'to','scf','input' ,typList=typList)
    if mmlio.yorn('Move old SCF output files to ACCRE?'): simfile.mvfiles_accre(simstr,'to','scf','output',typList=typList)
    return

####################################################################################################################################
# METHOD TO RETRIEVE OUTPUT FROM AN SCF ANALYSIS RUN
def getfiles(simstr,typList=None):
    """
    Retrieves SCF files from ACCRE
    """
    # Move files
    if mmlio.yorn('Recover SCF input files from ACCRE?') : simfile.mvfiles_accre(simstr,'from','scf','input' ,typList=typList)
    if mmlio.yorn('Recover SCF output files from ACCRE?'): simfile.mvfiles_accre(simstr,'from','scf','output',typList=typList)
    return

####################################################################################################################################
####################################################################################################################################
# FILE CLASSES AND METHODS

####################################################################################################################################
# METHOD FOR PARSING SCF FILE OPTIONS
def pars_fileopt(fopt):
    """
    Returns a parsed dictionary of SCF file options with keys:
        execpath: Str absolute path to location of compiled SCF executable
        parmpath: Str absolute path to location of default SCF par file
        modspath: Str absolute path to location of default SCF mod file
    """
    # Define keys
    form={
        'execpath': mmlpars.parsdict(default='/home/langmm/bin/scf/mpiscf',type=str),
        'parmpath': mmlpars.parsdict(default='/home/langmm/bin/scf/scfpar',type=str),
        'modspath': mmlpars.parsdict(default='/home/langmm/bin/scf/scfmod',type=str)
        }
    # Pars input
    fopt=mmlpars.mml_formpars(fopt,form)
    # Return parsed input
    return fopt
        
####################################################################################################################################
# METHOD TO CREATE DICTIONARY OF SCF FILES
def files(fdict):
    """
    Returns a dicitonary of SCF files & directories with keys:
      dir:      Directory containing all SCF analysis files (input/output)
      exdeft:   Location of static compiled SCF executable
      pmdeft:   Location of static SCF parameter file
      mddeft:   Location of static SCF modification file
      unitfile: Location of static units file (for pNbody stuff)
      input:    Directory containing all SCF simulation input files
        exec:      Compiled SCF executable used to produce run
        par:       Parameter file used by SCF to coordinate analysis
        mod:       SCF mod file
        wrap:      Bash wrapper for looping over SCF snapshots
        pbs:       Submission script to start analysis on ACCRE
        bodsbase:  Input particle information file
        coefbase:  Input coefficient file base string
      output:   Directory containing all SCF simulation output files
        runout:    Runtime output for submitted job
        logbase:   SCF log file
        chkbase:   ?
        outbase:   SCF output file
        bodsbase:  Output particle information file
        coefbase:  Output coefficient file base string
        elbase:    ?
        olilbase:  ?
    """
    # Set constats
    typList=simlist.LIST_PTYPS
    # Pars input
    fdict=mmlpars.mml_pars(fdict,type=dict)
    pfix=fdict['pfix']
    options=pars_fileopt({})
    # Initialize output dictionary
    files=dict(fdict['scf'],
               exdeft=options['execpath'],pmdeft=options['parmpath'],mddeft=options['modspath'])
    files['unitfile']=os.path.join(files['dir'],pfix+'units')
    # Loop over particle types create a directory tree for each
    for ityp in typList:
        files[ityp]={'dir':os.path.join(files['dir'],ityp)}
        # Add input files
        infiles={}
        infiles['dir']=os.path.join(files[ityp]['dir'],'input')
        fkeys_in=['exec','par','mod','wrap','pbs','sched']
        flist_in=['mpiscf','scfpar','scfmod','scfwrap','scfpbs','scfschednum']
        infiles=mmlfiles.filedict(infiles['dir'],filekeys=fkeys_in,filelist=flist_in,pfix=pfix,fdict=infiles)
        # Long input directories
        dkeys_in_long=['snapdir','coefdir']
        dlist_in_long=['snap','coef']
        infiles['subdirs']=dkeys_in_long
        infiles=mmlfiles.filedict(infiles['dir'],dirkeys=dkeys_in_long,dirlist=dlist_in_long,pfix=pfix,fdict=infiles)
        # Long input files
        fkeys_in_long=['snapbase','coefbase']
        flist_in_long=['scfbi','scficoef']
        for ifin in range(len(dkeys_in_long)):
            infiles=mmlfiles.filedict(infiles[dkeys_in_long[ifin]],
                                      filekeys=[fkeys_in_long[ifin]],filelist=[flist_in_long[ifin]],pfix=pfix,fdict=infiles)
        files[ityp]['input']=infiles
        # Add output files
        outfiles={}
        outfiles['dir']=os.path.join(files[ityp]['dir'],'output')
        fkeys_out=['pbsout']
        flist_out=['output']
        outfiles=mmlfiles.filedict(outfiles['dir'],filekeys=fkeys_out,filelist=flist_out,pfix=pfix,fdict=outfiles)
        # Long output directories
        dkeys_out_long=['logdir','outdir','snapdir','coefdir','eldir','chkdir','olildir']
        dlist_out_long=['log','out','snap','coef','el','chk','olil']
        outfiles['subdirs']=dkeys_out_long
        outfiles=mmlfiles.filedict(outfiles['dir'],dirkeys=dkeys_out_long,dirlist=dlist_out_long,pfix=pfix,fdict=outfiles)
        # Long output files
        fkeys_out_long=['logbase','outbase','snapbase','coefbase','elbase','chkbase','olilbase']
        flist_out_long=['scflog','scfout','snap','scfocoef','scfel','scfchkpt','slil']
        for ifout in range(len(dkeys_out_long)):
            outfiles=mmlfiles.filedict(outfiles[dkeys_out_long[ifout]],
                                       filekeys=[fkeys_out_long[ifout]],filelist=[flist_out_long[ifout]],pfix=pfix,fdict=outfiles)
        files[ityp]['output']=outfiles
    # Return file dictionary
    fdict['scf']=files
    return fdict

####################################################################################################################################
# METHOD TO RETURN LIST OF RELAVENT FILES
def get_filelist(fdict,keylist,ptype=None):
    """
    Returns a list of relavant SCF files 
    """
    # Set constants
    ptypeLIST=simlist.LIST_PTYPS
    # Pars input
    fdict=mmlpars.mml_pars(fdict,type=dict)
    keylist=mmlpars.mml_pars(keylist,type=list)
    ptype=mmlpars.mml_pars(ptype,default='disk',list=ptypeLIST)
    # Add files
    filelist={}
    if 'default'  in keylist:
        filelist['default' ]=[fdict['exdeft'],fdict['pmdeft'],fdict['mddeft']]
    if 'static'   in keylist:
        filelist['static'  ]=[fdict['unitfile']]
    if 'shortin'  in keylist:
        filelist['shortin' ]=[fdict[ptype]['input']['exec'],
                              fdict[ptype]['input']['par'],fdict[ptype]['input']['mod'],
                              fdict[ptype]['input']['wrap'],fdict[ptype]['input']['pbs']]
    if 'longin'   in keylist:
        baselist=['snap','coef']
        filelist['longin'  ]=[fdict[ptype]['input'][ibase+'base']+'*' for ibase in baselist]
    if 'shortout' in keylist:
        filelist['shortout']=[fdict[ptype]['output']['pbsout']]
    if 'longout'  in keylist:
        baselist=['snap','coef','log','out','el','chk','olil']
        filelist['longout' ]=[fdict[ptype]['output'][ibase+'base']+'*' for ibase in baselist]
    # Return output
    return filelist

####################################################################################################################################
####################################################################################################################################
# METHODS TO CREATE SCF FILES

####################################################################################################################################
# METHOD TO PARS DICTIONARY OF SCF PARAMETER FILE OPTIONS
def pars_paropt(popt):
    """
    Returns a parsed dictionary of SCF parameter options with keys:
        headline:  Str identifying the run being analyzed
        nsteps:    Int number of timesteps to take
        noutbod:   Int number of steps between output of system state
        noutlog:   Int number of steps between log output
        dteps:     The timestep
        G:         Float gravitational constant
        tfinal:    Float time to evolve system to
        multistep: Bool specifying if multiple steps should be taken
        fixedn:    Bool specifying if n is fixed
        selfgrav:  Bool specifying if self gravity is turned on
        inptcoef:  Bool specifying if expansion coefficients are supplied as input
        outpcoef:  Bool specifying if expansion coefficients should be returned as output
        zeroodd:   Bool specifying if odd coefficients should be zeroed
        zeroeven:  Bool specifying if even coefficients should be zeroed
        fixacc:    Bool specifying if conservation of linear momentum should be forced
        rcrit:     Float specifying ?
        ecrit:     Float specifying ?
        lilout:    Bool specifying ?
        nlilout:   Int specifying ?
    """
    # Define values
    form={
        'headline' : mmlpars.parsdict(default='galaxy',type=str         ),
        'nsteps'   : mmlpars.parsdict(default=0       ,type=int  ,min=0 ),
        'noutbod'  : mmlpars.parsdict(default=1       ,type=int  ,min=1 ),
        'noutlog'  : mmlpars.parsdict(default=1       ,type=int  ,min=1 ),
        'G'        : mmlpars.parsdict(default=1.0     ,type=float,min=0.),
        'tfinal'   : mmlpars.parsdict(default=100.    ,type=float,min=0.),
        'multistep': mmlpars.parsdict(default=False   ,type=bool        ),
        'fixedn'   : mmlpars.parsdict(default=False   ,type=bool        ),
        'selfgrav' : mmlpars.parsdict(default=True    ,type=bool        ),
        'inptcoef' : mmlpars.parsdict(default=False   ,type=bool        ),
        'outpcoef' : mmlpars.parsdict(default=True    ,type=bool        ),
        'zeroodd'  : mmlpars.parsdict(default=False   ,type=bool        ),
        'zeroeven' : mmlpars.parsdict(default=False   ,type=bool        ),
        'fixacc'   : mmlpars.parsdict(default=False   ,type=bool        ),
        'rcrit'    : mmlpars.parsdict(default=0.      ,type=float,min=0.),
        'ecrit'    : mmlpars.parsdict(default=0.      ,type=float,min=0.),
        'lilout'   : mmlpars.parsdict(default=False   ,type=bool        ),
        'nlilout'  : mmlpars.parsdict(default=0       ,type=int,min=0   )
        }
    # Pars input
    popt=mmlpars.mml_formpars(popt,form)
    # Return parsed input
    return popt

####################################################################################################################################
# METHOD FOR PARSING DICTIONARY OF SCF MODIFICATION FILE OPTIONS
def pars_modopt(mopt):
    """
    Returns a parsed dictionary of SCF modification file options with keys:
        iseed:      Int to initialize random number generator
        bhmass:     Float black hole mass
        epsbh:      Float black hole softening length
        tstartbh:   Float time to start black hole growth
        tgrowbh:    Float time to grow black hole
        tlivebh:    Float time black hole lives for
        tdiebh:     Float time to shrink black hole
        xdrag:      Float drag coeff for vx
        ydrag:      Float drag coeff for vy
        zdrag:      Float drag coeff for vz
        tstartdrag: Float time drag starts
        tgrowdrag:  Float time drag grows
        tlivedrag:  Float time drag lasts
        tdiedrag:   Float time drag dies down
        bhgrav:     Bool specifying if black hole is turned on
        usedrag:    Bool specifying if drag is turned on
        stellev:    Bool specifying if there is stellar evolution
    """
    # Define values
    form={
        'iseed'     : mmlpars.parsdict(default=3587 ,type=int         ),
        'bhmass'    : mmlpars.parsdict(default=0.   ,type=float,min=0.),
        'epsbh'     : mmlpars.parsdict(default=0.   ,type=float,min=0.),
        'tstartbh'  : mmlpars.parsdict(default=0.   ,type=float,min=0.),
        'tgrowbh'   : mmlpars.parsdict(default=0.   ,type=float,min=0.),
        'tlivebh'   : mmlpars.parsdict(default=0.   ,type=float,min=0.),
        'tdiebh'    : mmlpars.parsdict(default=0.   ,type=float,min=0.),
        'xdrag'     : mmlpars.parsdict(default=0.   ,type=float,min=0.),
        'ydrag'     : mmlpars.parsdict(default=0.   ,type=float,min=0.),
        'zdrag'     : mmlpars.parsdict(default=0.   ,type=float,min=0.),
        'tstartdrag': mmlpars.parsdict(default=0.   ,type=float,min=0.),
        'tgrowdrag' : mmlpars.parsdict(default=0.   ,type=float,min=0.),
        'tlivedrag' : mmlpars.parsdict(default=0.   ,type=float,min=0.),
        'tdiedrag'  : mmlpars.parsdict(default=0.   ,type=float,min=0.),
        'bhgrav'    : mmlpars.parsdict(default=False,type=bool        ),
        'usedrag'   : mmlpars.parsdict(default=False,type=bool        ),
        'stellev'   : mmlpars.parsdict(default=False,type=bool        )
        }
    # Pars input
    mopt=mmlpars.mml_formpars(mopt,form)
    # Return parsed input
    return mopt

####################################################################################################################################
# METHOD TO CREATE LOCAL COPIES OF SCF FILES
def mkfiles(simstr,paropt=None,ext=None,typList=None,nproc=None,idgal2=None,**extra_kw):
    """
    Creates SCF files
    """
    mmlio.verbose('Creating files for a SCF run locally...',border=True)
    # Set constants
    typListDEF=simlist.LIST_PTYPS
    typDict={typListDEF[ityp]:ityp for ityp in range(len(typListDEF))}
    # Pars input
    paropt=pars_paropt(paropt)
    ext=mmlpars.mml_pars(ext,default='*',type=str)
    typList=mmlpars.mml_pars(typList,default=typListDEF,type=list)
    nproc=mmlpars.mml_pars(nproc,default=simstr['nprocscf'],type=int)
    # Get file names
    files=simstr.fdict
    # Unit file
    if mmlio.yorn('Create units file?'): unitdict=mkunit(simstr,overwrite=True)
    # Loop over particle types creating directories & asking if files should be created
    dict_mkfile={}
    dict_mkexec={}
    dict_mkparm={}
    dict_mksnap={}
    dict_mkcoef={}
    flag_mksnap=False
    flag_mkcoef=False
    print 'Create files for...'
    for iptyp in typList:
        dict_mkfile[iptyp]=mmlio.yorn('    {} particles?'.format(iptyp))
        if dict_mkfile[iptyp]:
            files_typ=files['scf'][iptyp]
            mmlfiles.mkdirs([files_typ['dir'],files_typ['input']['dir']]+files_typ['input']['subdirs'])
            dict_mkexec[iptyp]=mmlio.yorn('      SCF executable?    ({})'.format(iptyp))
            dict_mkparm[iptyp]=mmlio.yorn('      SCF runtime files? ({})'.format(iptyp))
            dict_mksnap[iptyp]=mmlio.yorn('      SCF snapshots?     ({})'.format(iptyp))
            dict_mkcoef[iptyp]=False #mmlio.yorn('      SCF coefficients?  ({})'.format(iptyp))
            if dict_mksnap[iptyp]: flag_mksnap=True
            if dict_mkcoef[iptyp]: flag_mkcoef=True
        else:
            dict_mkexec[iptyp]=False
            dict_mkparm[iptyp]=False
            dict_mksnap[iptyp]=False
            dict_mkcoef[iptyp]=False
    if flag_mksnap: owsnap=mmlio.yorn('Overwrite existing SCF snapshots?')
    if flag_mkcoef: owcoef=mmlio.yorn('Overwrite existing SCF coefficients?')
    # Single files
    for iptyp in typList:
        if not dict_mkfile[iptyp]: continue
        # Create directories
        mmlfiles.mkdirs(files['scf'][iptyp]['input']['dir'])
        # Executable file
        if dict_mkexec[iptyp]:
            if not os.path.isfile(files['scf']['exdeft']):
                raise Exception('Static SCF executable does not exist. Create it first: {}'.format(files['scf']['exdeft']))
            mmlfiles.cpfiles(files['scf']['exdeft'],files['scf'][iptyp]['input']['exec'],overwrite=True)
        # Parameter files
        if dict_mkparm[iptyp]:
            pardict=mkpar(simstr,iptyp,overwrite=True)
            moddict=mkmod(simstr,iptyp,overwrite=True)
            pbsdict=mkpbs(simstr,iptyp,overwrite=True,nproc=nproc)
            wrapdict=mkwrap(simstr,iptyp,overwrite=True,mpipath=pbsdict['mpipath'],nproc=nproc)
    # Snapshot files & icfile
    if flag_mksnap:
        # Create directories
        mmlfiles.mkdirs(files['scf'][iptyp]['input']['snapdir'])
        # Locate gadget snapshots
        gadgsnap=files['gadget']['output']['snapbase']+ext
        srcsnapList=sorted(glob.glob(gadgsnap)+glob.glob(files['gadget']['input']['ic']))
        if len(srcsnapList)==0:
            raise Exception('No GADGET snapshots were found matching the pattern {}'.format(gadgsnap))
        for isrc in srcsnapList:
            # Get extension from source snapshot
            iext=simstr.get_fext(isrc,ftype='gadget')
            flag_readsrc=False
            for iptyp in typList:
                idst=files['scf'][iptyp]['input']['snapbase']+iext
                if not dict_mksnap[iptyp]: continue
                if not owsnap and os.path.isfile(idst): continue
                if not flag_readsrc:
                    flag_readsrc=True
                    inb=simstr.get_pnbody(fname=isrc,ftype='gadget',ptype=iptyp,idgal2=idgal2)
                    typList_INC=inb.get_typlist()
                    for ityp in typList:
                        if ityp in typList_INC:
                            icenter=inb.get_center(centergal=1,centertyp=ityp,move=True,method='dens')
                    inb.set_ftype('scf')
                if inb.npart[typDict[iptyp]]==0: continue
                print 'Creating SCF snapshot: {}'.format(idst)
                index=inb.get_index(ptyp=iptyp)#,pgal=1)
                inb_typ=pNbody.Nbody(p_name=idst,ftype='scf',status='new',unitsfile=files['scf']['unitfile'],
                                     pos=inb.pos[index,:],vel=inb.vel[index,:],mass=inb.mass[index])
#                centerList=inb_typ.get_center(move=True,method='density')
                inb_typ.atime=inb.atime
                inb_typ.nbody=len(inb_typ.mass)
                inb_typ.write()
                del inb_typ
            if flag_readsrc: del inb
    # Coefficient files
    if flag_mkcoef:
        # Create directories
        mmlfiles.mkdirs(files['scf'][iptyp]['input']['coefdir'])
        raise Exception('Creation of SCF coefficient files not currently supported')
    return

####################################################################################################################################
# METHOD TO CREATE SCF PARAMETER FILE
def mkpar(simstr,ptype,overwrite=False):
    """
    Creates an SCF parameter file and returns contents as a dictionary
    """
    # Set constants
    keywidth=29
    headline='**********Basic input parameters**********'
    tailline='******************************************'
    # Pars input
    overwrite=mmlpars.mml_pars(overwrite,type=bool)
    # Get file names
    files_local=simstr.mkfiledict(checkaccess=True)
    files_accre=simstr.mkfiledict(memcomp='accre',checkaccess=False)
    parfile=files_local['scf'][ptype]['input']['par']
    parfile_stat=files_local['scf']['pmdeft']
    # Load existing file if overwrite False
    if not overwrite and os.path.isfile(parfile):
        mmlio.verbose('Loading existing SCF parameter file:',addval=parfile)
        pardict=mmlio.rwdict('R',parfile,style='fortran')
    # Otherwise create parameter file from default
    else:
        mmlio.verbose('Creating new SCF parameter file:',addval=parfile)
        # If static parameter file exists, load it
        if os.path.isfile(parfile_stat):
            pardict=mmlio.rwdict('R',parfile_stat,style='fortran')
        # If the static parameter file dosn't exist, make it and use default
        else:
            pardict=paropt
            pardict['keylist']=['headline','nsteps','noutbod','noutlog','dteps','G',
                                'tfinal','multistep','fixedn','selfgrav','inptcoef','outpcoef',
                                'zeroodd','zeroeven','fixacc','rcrit','ecrit','lilout','nlilout']
            mmlio.rwdict('W',parfile_stat,pardict,width=keywidth,style='fortran',
                         headline=headline,tailline=tailline,overwrite=overwrite)
        # Modify default for this run
        pardict['headline']=simstr['runtag']
        # Write parameter dictionary to file
        mmlio.rwdict('W',parfile,pardict,width=keywidth,style='fortran',
                     headline=headline,tailline=tailline,overwrite=overwrite)
    # Return output
    return pardict


####################################################################################################################################
# METHOD TO CREATE SCF MODIFICATION FILE
def mkmod(simstr,ptype,modopt=None,overwrite=False):
    """
    Create a SCF modifcation file and return results
    """
    # Set constants
    keywidth=29
    headline='**********Basic input parameters**********'
    tailline='******************************************'
    # Pars input
    modopt=pars_modopt(modopt)
    overwrite=mmlpars.mml_pars(overwrite,type=bool)
    # Get file names
    files_local=simstr.mkfiledict(checkaccess=True)
    files_accre=simstr.mkfiledict(memcomp='accre',checkaccess=False)
    modfile=files_local['scf'][ptype]['input']['mod']
    modfile_stat=files_local['scf']['mddeft']
    # Load existing file if overwrite False
    if not overwrite and os.path.isfile(modfile):
        mmlio.verbose('Loading existing SCF modification file:',addval=modfile)
        moddict=mmlio.rwdict('R',modfile,style='fortran')
    # Otherwise create file from default
    else:
        mmlio.verbose('Creating new SCF modification file:',addval=modfile)
        # If static file exists, load it
        if os.path.isfile(modfile_stat):
            moddict=mmlio.rwdict('R',modfile_stat,style='fortran')
        # If the static file dosn't exist, make it and use default
        else:
            moddict=modopt
            moddict['keylist']=['iseed','bhmass','epsbh','tstartbh','tgrowbh','tlivebh','tdiebh',
                                'xdrag','ydrag','zdrag','tstartdrag','tgrowdrag','tlivedrag','tdiedrag',
                                'bhgrav','usedrag','stellev']
            mmlio.rwdict('W',modfile_stat,moddict,width=keywidth,style='fortran',
                         headline=headline,tailline=tailline,overwrite=overwrite)
        # Modify default for this run
        moddict['headline']=simstr['runtag']
        # Write dictionary to file
        mmlio.rwdict('W',modfile,moddict,width=keywidth,style='fortran',
                     headline=headline,tailline=tailline,overwrite=overwrite)
    # Return output
    return moddict


####################################################################################################################################
# METHOD TO CREATE PBS SUBMISSION SCRIPT
def mkpbs(simstr,ptype,overwrite=False,nproc=None,ppn=None):
    """
    Creates an SCF pbs script and returns contents as a dictionary
    """
    # Pars input
    overwrite=mmlpars.mml_pars(overwrite,type=bool)
    nproc=mmlpars.mml_pars(nproc,default=simstr['nprocscf'],type=int,min=1)
    ppn=mmlpars.mml_pars(ppn,default=2,type=int,range=[1,4])
    # Get file names
    files_local=simstr.mkfiledict(checkaccess=True)
    files_accre=simstr.mkfiledict(memcomp='accre',checkaccess=False)
    pbsfile=files_local['scf'][ptype]['input']['pbs']
    parfile=files_local['scf'][ptype]['input']['par']
    pbsdict=mmlio.pbsdict()
    # Load existing file if overwrite False
    if not overwrite and os.path.isfile(pbsfile):
        mmlio.verbose('Loading existing SCF pbs script:',addval=pbsfile)
        pbsdict.rwmpi('R',pbsfile)
    # Otherwise create file from scratch
    else:
        mmlio.verbose('Creating new SCF pbs script:',addval=pbsfile)
        # Update memory
        mcpu=mmlsim2mcpu(simstr)
        tcpu=mmlsim2tcpu(simstr)
        # Variables from simstr
        pbsdict['ppn'      ]=ppn
        pbsdict['nodes'    ]=nproc/ppn
        pbsdict['pmem'     ]=mcpu
        pbsdict['mem'      ]=mcpu*nproc
        pbsdict['walltime' ]=tcpu
        # Varaibles from path
        pbsdict['simout'   ]=files_accre['scf'][ptype]['output']['pbsout']
#        pbsdict['rundir'   ]=os.path.dirname(files_accre['scf'][ptype]['input']['exec'])
#        pbsdict['execpath' ]=os.path.basename(files_accre['scf'][ptype]['input']['exec'])
#        pbsdict['execflag' ]=0
#        pbsdict['execinput']=' '
#        pbsdict.rwmpi('W',pbsfile,overwrite=overwrite)
        pbsdict['command'  ]=['./{}'.format(os.path.basename(files_accre['scf'][ptype]['input']['wrap']))]
        pbsdict.rwpbs('W',pbsfile,overwrite=overwrite)
    # Return output
    return pbsdict

####################################################################################################################################
# METHOD TO CREATE MPISCF WRAPPER
def mkwrap(simstr,ptype,overwrite=None,mpipath=None,nproc=None):
    """
    Creates an SCF wrapper script and returns contents as a dictionary
    """
    # Set constants
    mpipathDEF=mmlio.pbsdict()['mpipath']
    # Pars input
    overwrite=mmlpars.mml_pars(overwrite,default=False,type=bool)
    mpipath=mmlpars.mml_pars(mpipath,default=mpipathDEF,type=str)
    nproc=mmlpars.mml_pars(nproc,default=simstr['nprocscf'],type=int,min=1)
    # Get file names
    files_local=simstr.mkfiledict(checkaccess=True)
    files_accre=simstr.mkfiledict(memcomp='accre',checkaccess=False)
    wrapfile=files_local['scf'][ptype]['input']['wrap']
    # Load existing file if overwrite False
    if not overwrite and os.path.isfile(wrapfile):
        mmlio.verbose('SCF bash wrapper already exists:',addval=wrapfile)
        cmdlist=[]
    # Otherwise create file from scratch
    else:
        mmlio.verbose('Creating new SCF wrapper script:',addval=wrapfile)
        # Initalize command list
        cmdlist=[
            "# Input files",
            "exec      = '{}'".format(files_accre['scf'][ptype]['input']['exec']).replace(' ',''),
            "par       = '{}'".format(files_accre['scf'][ptype]['input']['par']).replace(' ',''),
            "mod       = '{}'".format(files_accre['scf'][ptype]['input']['mod']).replace(' ',''),
            "isnapbase = '{}'".format(files_accre['scf'][ptype]['input']['snapbase']).replace(' ',''),
            "icoefbase = '{}'".format(files_accre['scf'][ptype]['input']['coefbase']).replace(' ',''),
            "execfile  = '{}'".format(os.path.splitext(files_accre['scf'][ptype]['input']['exec'])[1][1:]).replace(' ',''),
            "parfile   = '{}'".format(os.path.splitext(files_accre['scf'][ptype]['input']['par'])[1][1:]).replace(' ',''),
            "modfile   = '{}'".format(os.path.splitext(files_accre['scf'][ptype]['input']['mod'])[1][1:]).replace(' ',''),
            "isnapfile = '{}'".format(os.path.splitext(files_accre['scf'][ptype]['input']['snapbase'])[1][1:]).replace(' ',''),
            "icoeffile = '{}'".format(os.path.splitext(files_accre['scf'][ptype]['input']['coefbase'])[1][1:]).replace(' ',''),
            "# Output files",
            "osnapbase = '{}'".format(files_accre['scf'][ptype]['output']['snapbase']).replace(' ',''),
            "ocoefbase = '{}'".format(files_accre['scf'][ptype]['output']['coefbase']).replace(' ',''),
            "outbase   = '{}'".format(files_accre['scf'][ptype]['output']['outbase']).replace(' ',''),
            "logbase   = '{}'".format(files_accre['scf'][ptype]['output']['logbase']).replace(' ',''),
            "chkbase   = '{}'".format(files_accre['scf'][ptype]['output']['chkbase']).replace(' ',''),
            "elbase    = '{}'".format(files_accre['scf'][ptype]['output']['elbase']).replace(' ',''),
            "olilbase  = '{}'".format(files_accre['scf'][ptype]['output']['olilbase']).replace(' ',''),
            "osnapfile = '{}'".format(os.path.splitext(files_accre['scf'][ptype]['output']['snapbase'])[1][1:]).replace(' ',''),
            "ocoeffile = '{}'".format(os.path.splitext(files_accre['scf'][ptype]['output']['coefbase'])[1][1:]).replace(' ',''),
            "outfile   = '{}'".format(os.path.splitext(files_accre['scf'][ptype]['output']['outbase'])[1][1:]).replace(' ',''),
            "logfile   = '{}'".format(os.path.splitext(files_accre['scf'][ptype]['output']['logbase'])[1][1:]).replace(' ',''),
            "chkfile   = '{}'".format(os.path.splitext(files_accre['scf'][ptype]['output']['chkbase'])[1][1:]).replace(' ',''),
            "elfile    = '{}'".format(os.path.splitext(files_accre['scf'][ptype]['output']['elbase'])[1][1:]).replace(' ',''),
            "olilfile  = '{}'".format(os.path.splitext(files_accre['scf'][ptype]['output']['olilbase'])[1][1:]).replace(' ',''),
            "lensnapbase=${#isnapbase}","tempdir='temprun'"]
        # Command to initialize files
        cmdlist+=[' ',
                  '# Move files to temporary directory w/ built in SCF names',
                  "if [ ! -e $tempdir ]; then","    mkdir $tempdir","fi","cd $tempdir",
                  "cp $exec $execfile",
                  "cp $par $parfile",
                  "cp $mod $modfile",
                  "snaplist=$(ls {}*)".format(files_accre['scf'][ptype]['input']['snapbase'])]
        # Commands to perform loop
        cmdlist+=[' ',
                  '# Loop over SCF input snapshots',
                  'for isnap in $isnapbase*; do',
                  '    # Get file extension',
                  '    iext=${isnap:${lensnapbase}}',
                  '    echo "Starting snapshot: ${isnap}"',
                  '    echo "ext=${iext}"',
                  ' ',
                  '    # Only continue if output is missing',
                  '    if [ ! -e ${osnapbase}${iext} ]; then',
                  ' ',
                  '        # Move necessary input files',
                  '        if [ -e "${isnapbase}${iext}" ]; then','            mv "${isnapbase}${iext}" $isnapfile','        fi',
                  '        if [ -e "${icoefbase}${iext}" ]; then','            mv "${icoefbase}${iext}" $icoeffile','        fi',
                  ' ',
                  '        # Run mpiscf script',
                  '        {} -np {} ./$execfile'.format(mpipath,nproc),
                  ' ',
                  '        # Recover files',
                  '        if [ -e $isnapfile  ]; then','            mv $isnapfile  "${isnapbase}${iext}"','        fi',
                  '        if [ -e $icoeffile  ]; then','            mv $icoeffile  "${icoefbase}${iext}"','        fi',
                  '        if [ -e $osnapfile* ]; then','            mv $osnapfile* "${osnapbase}${iext}"','        fi',
                  '        if [ -e $ocoeffile  ]; then','            mv $ocoeffile  "${ocoefbase}${iext}"','        fi',
                  '        if [ -e $outfile    ]; then','            mv $outfile    "${outbase}${iext}"','        fi',
                  '        if [ -e $logfile    ]; then','            mv $logfile    "${logbase}${iext}"','        fi',
                  '        if [ -e $chkfile    ]; then','            mv $chkfile    "${chkbase}${iext}"','        fi',
                  '        if [ -e $elfile*    ]; then','            mv $elfile*    "${elbase}${iext}"','        fi',
                  '        if [ -e $olilfile*  ]; then','            mv $olilfile*  "${olilbase}${iext}"','        fi',
                  '    fi',
                  'done']
        # Commands to finalize files
        cmdlist+=[' ',
                  '# Move files back to original paths',
                  'rm $execfile',# $exec',
                  'rm $parfile',# $par',
                  'rm $modfile',# $mod',
                  'rmdir $tempdir',
                  'cd ../']
        # Write commands to file
        mmlio.rwbash('W',wrapfile,cmdlist,overwrite=overwrite)
        # Make file executable
        os.chmod(wrapfile,0755)
    # Return output
    return cmdlist

####################################################################################################################################
# METHOD TO CREATE SCF UNITS FILE
def mkunit(simstr,overwrite=None):
    """
    Creates an SCF units file and returns contents as a dicitonary
    """
    # Get file dictionary
    filedict=simstr.fdict
    # Get gadget parameter file & SCF unit file
    parmfile=filedict['gadget']['input']['param']
    unitfile=filedict['scf']['unitfile']
    # Make necessary directory
    mmlfiles.mkdirs(os.path.dirname(unitfile))
    # Copy gadget parameter file to SCF unit file
    mmlfiles.cpfiles(parmfile,unitfile,overwrite=overwrite)
    # Load and return parameter dictionary
    units=mmlio.rwdict('R',unitfile)
    return units

####################################################################################################################################
####################################################################################################################################
# METHODS TO READ SCF FILES

####################################################################################################################################
# METHOD TO READ IN SCF SNAPSHOTS
def read_snapshot(fileDict,ext,flag_out=False):
    """
    Reads in and SCF snapshot
    """
    # Create constants
    typList=['gas','halo','disk','bulge','stars','bndry']
    typDict={typList[ityp]:ityp for ityp in range(len(typList))}
    iokeys=['mass','pos','vel']
    okeys=['pot','adens','dti','mi','potm2']
    # Preallocate
    outDict=dict(atime=0.,npart=len(typList)*[0L],nbody=0L,tpe=np.array([]),
                 mass=np.array([]),pos=np.array([]),vel=np.array([]),
                 pot=np.array([]),adens=np.array([]),dti=np.array([]),
                 mi=np.array([]),potm2=np.array([]))
    # Loop over type reading in each file
    for ityp in typList:
        # Get file contents
        if flag_out:
            ifile=fileDict[ityp]['output']['snap']+ext
            ioutDict=read_obodfile(ifile)
        else:
            ifile=fileDict[ityp]['input']['snap']+ext
            ioutDict=read_ibodfile(ifile)
        # Add contents to total dictionary
        outDict['atime']=ioutDict['time']
        outDict['npart'][typDict[ityp]]=ioutDict['nbod']
        outDict['tpe']=np.concatenate((outDict['tpe'],typDict[ityp]*np.ones((ioutDict['nbod'],1))),axis=0)
        for iokey in iokeys: outDict[iokey]=np.concatenate((outDict[iokey],ioutDict[iokey]),axis=0)
        if flag_out:
            for okey in okeys: outDict[okey]=np.concatenate((outDict[okey],ioutDict[okey]),axis=0)
    # Sum number of particles and return output
    outDict['nbody']=np.array(outDict['npart']).sum()
    return outDict

####################################################################################################################################
# METHOD FOR READING OUTPUT SNAPSHOTS
def read_iobodfile(fname,flag_out=False,nexpect=None):
    """
    NAME:
        mmlscf.read_iobodfile
    PURPOSE:
        To read in output particle info from an SCF obods file.
    CALLING:
        oboddict=read_iobodfile(fname)
    ARGUMENTS:
        fname:    Str absolute path to SCF obods file.
    KEYWORDS:
        flag_out: Bool specifying if the the file is an output snapshot.
        nexpect:  Long number of particles expected in the file.
    OUTPUT:
        oboddict: Dictionary with output particle info in keys:
            nbod:   Long number of particles
            time:   Float time of snapshot output
            bhmass: Float mass of black hole particle
            mass:   Nx1 float array of particle masses
            pos:    Nx3 float array of particle positions
            vel:    Nx3 float array of particle velocities
            pot:    Nx1 float array of particle potentials
            adens:  Nx1 float array of particle densities?
            dti:    Nx1 float array of ?
            mi:     Nx1 float array of ?
            potm2:  Nx1 float array of particle m=2 potentials
    """
    if isinstance(fname,str): fid=open(fname,'r')
    else:                     fid=fname
    # Read first line
    headline=fid.readline()
    if flag_out:
        nbod,time,bhmass=headline.split()
        bhmass=float(bhmass)
#        print nbod,time,bhmass
    else:
        nbod,time=headline.split()
#        print nbod,time
    nbod=long(float(nbod))
    time=float(time)
    # Initialize dictionary
    oboddict=dict(nbod=nbod,time=time)
    if flag_out: oboddict['bhmass']=bhmass
    # Read in bod info
    if flag_out:
        m,x,y,z,vx,vy,vz,p,ad,dti,mi,pm2 = pNbody.io.read_ascii(fid,range(12),skipheader=True)
    else:
        m,x,y,z,vx,vy,vz = pNbody.io.read_ascii(fid,range(7),skipheader=True)
    # Determine amount of padding needed
    if isinstance(nexpect,long) or isinstance(nexpect,int):
        npad=nexpect-nbod
        if   npad == 0: pass
        elif npad >  0: print '[mmlscf.read_iobodfile] There are missing particles, padding data with zeros.'
        else:           raise Exception('There are {} extra particles!?'.format(abs(npad)))
    else:
        npad=0
    arrpad=np.zeros((npad,))
    arrpad3=np.zeros((npad,3))
    # Populate output dictionary
    oboddict['mass']  = np.concatenate((m.astype(np.float32),arrpad),axis=0)
    oboddict['pos']   = np.concatenate((np.transpose(np.array([x,y,z])).astype(np.float32),arrpad3),axis=0)
    oboddict['vel']   = np.concatenate((np.transpose(np.array([vx,vy,vz])).astype(np.float32),arrpad3),axis=0)
    if flag_out:
        oboddict['pot']   = np.concatenate((p.astype(np.float32),arrpad),axis=0)
        oboddict['adens'] = np.concatenate((ad.astype(np.float32),arrpad),axis=0)
        oboddict['dti']   = np.concatenate((dti.astype(np.float32),arrpad),axis=0)
        oboddict['mi']    = np.concatenate((mi.astype(np.float32),arrpad),axis=0)
        oboddict['potm2'] = np.concatenate((pm2.astype(np.float32),arrpad),axis=0)
    # Close file & return output
    fid.close()
    return oboddict
    

####################################################################################################################################
# METHOD TO READ IN SCF COEFFICIENTS
def read_coeffile(fname):
    '''
    NAME:
        mmlscf.read_coeffile
    PURPOSE:
        To read in output coefficients from an SCF ocoef file.
    CALLING:
        coefdict=read_coeffile(fname)
    ARGUMENTS:
        fname:    Str absolute path to SCF ocoef file.
    OUTPUT:
        coefdict: Dictionary with SCF coefficient keys:
            time:   Simulation time that coefficients are written out
            nmax:   Maximum order of n coefficient
            lmax:   Maximum order of l coefficient
            sinsum: [nmax,lmax,lmax] array of n,l,m coefficients
            cossum: [nmax,lmax,lmax] array of n,l,m coefficients
    '''
    fid=open(fname,'r')
    # Read first line
    headline=fid.readline()
    time,nmax,lmax=headline.split()
    time=float(time)
    nmax=int(float(nmax))
    lmax=int(float(lmax))
    print time,nmax,lmax
    # Initialize dictionary
    coefdict=dict(time=time,nmax=nmax,lmax=lmax)
    coefdict['sinsum']=np.zeros((nmax+1,lmax+1,lmax+1))
    coefdict['cossum']=np.zeros((nmax+1,lmax+1,lmax+1))
    # Loop over coefficients reading them in
    for n in range(coefdict['nmax']+1):
        for l in range(coefdict['lmax']+1):
            for m in range(l+1):
                iline=fid.readline()
                isinsum,icossum=iline.split()
                coefdict['sinsum'][n][l][m]=float(isinsum)
                coefdict['cossum'][n][l][m]=float(icossum)
    # Close file & return output
    fid.close()
    return coefdict

####################################################################################################################################
####################################################################################################################################
# METHOD TO CALCULATE SCF STUFF

####################################################################################################################################
# METHOD TO CALCULATE MEMORY REQUIRED TO RUN SCF ON ACCRE
def mmlsim2mcpu(simstr):
    """
    Estimates the amount of memory an SCF run will need
    """
    # Force reasonable values as a safety
    if   simstr['runtyp'].lower()=='galaxy'  : mcpu=1000./simstr['nprocscf'] # mb
    elif simstr['runtyp'].lower()=='interact': mcpu=2000./simstr['nprocscf'] # mb
    # Return rounded output
    return mmlmath.oom(mcpu,nsig=1,method='CEIL')

####################################################################################################################################
# METHOD TO CALCULATE TIME REQUIRED TO RUN SCF ON ACCRE
def mmlsim2tcpu(simstr):
    """
    Estimates the amount of time required to run SCF
    """
    # Force reasonable values as a safety
    if   simstr['runtyp'].lower()=='galaxy'  : tcpu=1. # hr
    elif simstr['runtyp'].lower()=='interact': tcpu=2. # hr
    # Return rounded output
    return tcpu

####################################################################################################################################
# METHOD TO CALCULATE POTENTIAL FROM COEFFICIENTS
def calc_pot(coefdict,pos,nlist=None,llist=None,mlist=None):
    """
    Determines the potential from SCF coefficients (bad)
    """
    # Pars input
    nlist=mmlpars.mml_pars(nlist,default=range(coefdict['nmax']+1),type=list)
    llist=mmlpars.mml_pars(llist,default=range(coefdict['lmax']+1),type=list)
    mlist=mmlpars.mml_pars(mlist,default=range(coefdict['lmax']+1),type=list)
    # Convert coordinates
    x=pos[:,0] ; y=pos[:,1] ; z=pos[:,2]
    N=len(x)
    rxyz=np.sqrt(x*x+y*y+z*z)
    rxy=np.sqrt(x*x+y*y)
    xi=(rxyz-1.)/(rxyz+1.)
    costh=z/rxyz
    sinth=rxy/rxyz
    cosmphi=np.zeros((coefdict['lmax']+1,N))
    cosmphi[0,:]=1.
    cosmphi[1,:]=y/rxy
    sinmphi=np.zeros((coefdict['lmax']+1,N))
    sinmphi[0,:]=0.
    sinmphi[1,:]=x/rxy
    for m in range(2,coefdict['lmax']+1):
        cosmphi[m,:]=cosmphi[1,:]*cosmphi[m-1,:] - sinmphi[1,:]*sinmphi[m-1,:]
        sinmphi[m,:]=cosmphi[1,:]*sinmphi[m-1,:] + sinmphi[1,:]*cosmphi[m-1,:]
    # Create polynomials
    rl=np.zeros((coefdict['lmax']+1,N))
    plm=np.zeros((coefdict['lmax']+1,coefdict['lmax']+1,N))
    cnl=np.zeros((coefdict['nmax']+1,coefdict['lmax']+1,N))
    for k in range(len(costh)):
        iplm,idplm=scipy.special.lpmn(coefdict['lmax'],coefdict['lmax'],costh[k])
        plm[:,:,k]=np.transpose(iplm)[:,:]
    for l in range(0,coefdict['lmax']+1):
        rl[l,:]=(rxyz**l)/((1.+rxyz)**(2*l+1))
#        for m in range(0,l+1):
            # Associate Legendre Polynomial of 1st Kind
#            for k in range(len(costh)):
#                plm[l,m,k]=scipy.special.lpmn(m,l,costh[k])
        # Ultra-spherical polynomials
        for n in range(0,coefdict['nmax']+1):
            alpha=2.*float(l)+1.5
            cnl[n,l,:]=scipy.special.gegenbauer(n,alpha)(xi)
    # Calculate potential
    pot=np.zeros(N)
    for l in llist:
        for m in mlist:
            Alm=np.zeros(N) ; Blm=np.zeros(N)
            for n in nlist:
                Alm+=cnl[n,l,:]*coefdict['cossum'][n][l][m]
                Blm+=cnl[n,l,:]*coefdict['sinsum'][n][l][m]
            pot+=rl[l,:]*plm[l,m,:]*(Alm*cosmphi[m,:] + Blm*sinmphi[m,:])
    # Return output
    return pot

####################################################################################################################################
####################################################################################################################################
# PROVIDE COMMAND LINE ACCESS
if __name__ == '__main__': main()
