#!/usr/bin/python
####################################################################################################################################
#
# MEAGAN LANG'S GADGET METHODS
#
####################################################################################################################################
import sys,os,shutil,glob,copy,pprint,scipy,math
import numpy as np
import matplotlib as mplib
import pNbody
LIST_METHODS=['makestatic','makefiles','submit','retrieve','clean']
from mmlutils import *
import main as mmlnbody
import simlist,simcalc,simplot,simfile,mmlbuildgal,mmlscf,mmlgalfit,icgen

def main():
    """
    Provides command line access to GADGET methods
    """
    out = mmlnbody.walk(mtype='gadget')
    return out

####################################################################################################################################
####################################################################################################################################
# METHODS FOR GETTING LISTS
def dict_ftype():
    """
    Returns a dictionary of files listed by type
    """
    fdict={'shortout':['shortout'],
           'shortin' :['exec','param'],
           'longout' :['snapshot','restart','rest_bak'],
           'shortlongout':['shortsnapshot']}
    return fdict
def dict_fgroup():
    """
    Returns a dictionary of file types listed by group
    """
    fdict={'clean' :['snapshot','rest_bak'],
           'input' :['shortin' ,'longin' ],
           'output':['shortout','longout']}
    return fdict

####################################################################################################################################
####################################################################################################################################
# HIGH LEVEL METHODS REQUIRING MMLSIM

####################################################################################################################################
# METHOD FOR RUNNING DIFFERENT GADGET OPERATIONS
def run(simstr,method,**method_kw):
    """
    Provides interface for running different GADGET operations
    """
    # Set constants
    methLIST=LIST_METHODS
    # Pars input
    method=mmlpars.mml_pars(method,list=methLIST)
    # Initialize default output
    out=None
    # Proceed based on method
    if   method=='makestatic':
        print simstr['runtyp']
        if   simstr['runtyp']=='galsim': mmlgadget.mkstatic(simstr,**method_kw)
        elif simstr['runtyp']=='intsim':
            method_kw['snapnum1']=mmlio.askquest('Snapshot # to use from galaxy 1 to create IC file?',default=-1,dtype='int')
            method_kw['snapnum2']=mmlio.askquest('Snapshot # to use from galaxy 2 to create IC file?',default=-1,dtype='int')
            mkstatic(simstr,**method_kw)
        else: raise Exception('Invalid run type: {}'.format(simstr['runtyp']))
    elif method=='makefiles' : mkfiles(simstr,**method_kw)
    elif method=='submit'    : subrun(simstr,**method_kw)
    elif method=='retrieve'  :
        method_kw['snapend']=mmlio.askquest('Number retrieved snapshots should end in?',default=5,dtype='str')
        endrun(simstr,**method_kw)
    elif method=='clean'     :
        method_kw['ftype']='clean'
        method_kw['memcomp']=simstr['runcomp']
        simfile.rmfiles('gadget',simstr=simstr,**method_kw)
    else:
        raise Exception('Option for method {} needs to be added to mmlgadget.run.'.format(method))
    # Return output
    return out

####################################################################################################################################
# METHOD PERFORM ACTIONS TO SETUP AND RUN A GADGET SIMULATION ON ACCRE
def subrun(simstr,**extra_kw):
    """
    Sets up and runs GADGET on ACCRE
    """
    # Get constant
    files_accre=simstr.mkfiledict(memcomp=simstr['runcomp'],checkaccess=True)['gadget']
    # Update input files
    if mmlio.yorn('Create new input files?'): mkfiles(simstr)
    # Move input files to ACCRE
    if mmlio.yorn('Move new input files to ACCRE?'): putfiles(simstr)
    # Make directories
    mmlfiles.mkdirs(files_accre['output']['dir'])
    # Run GADGET
    if mmlinfo.hostname2compid()=='bender':
        if mmlio.yorn('Run GADGET?'):
            icmeth=mmlio.askselect('How should the run be started?',['ic','snapshot','restart'])
            if icmeth=='ic':
                pbsfile=files_accre['input']['pbs']
            elif icmeth=='snapshot':
                pbsfile=files_accre['input']['snpbs']
            elif icmeth=='restart':
                pbsfile=files_accre['input']['repbs']
            schfile=files_accre['input']['sched']
            simstr.subpbs('gadget',pbsfile=pbsfile,schfile=schfile)
    return
        
####################################################################################################################################
# METHOD TO PERFORM ACTIONS TO RETRIEVE AND CLEAN-UP A RUN
def endrun(simstr,snapend=None,**extra_kw):
    """
    Retrieves and cleans up results from a GADGET run
    """
    # Retrieve results
    if mmlio.yorn('Recover GADGET results from run machine?'): getfiles(simstr,snapend=snapend)
    # Clean up results
    if mmlio.yorn('Clean GADGET results on run machine?'    ): simfile.rmfiles('gadget',ftype='clean',memcomp=simstr['runcomp'],simstr=simstr)
    return

####################################################################################################################################
# METHOD TO MOVE FILES TO ACCRE
def putfiles(simstr):
    """
    Moves GADGET files to ACCRE
    """
    if mmlio.yorn('Move GADGET input files to ACCRE?'     ): simfile.mvfiles_accre(simstr,'to','gadget','input' )
    if mmlio.yorn('Move GADGET old output files to ACCRE?'): simfile.mvfiles_accre(simstr,'to','gadget','output')
    return

####################################################################################################################################
# METHOD TO GET FILES FROM ACCRE
def getfiles(simstr,snapend=None):
    """
    Retrieves GADGET files from ACCRE
    """
    # Move files
    if mmlio.yorn('Recover GADGET input files from ACCRE?' ): simfile.mvfiles_accre(simstr,'from','gadget','input' ,snapend=snapend)
    if mmlio.yorn('Recover GADGET output files from ACCRE?'): simfile.mvfiles_accre(simstr,'from','gadget','output',snapend=snapend)

####################################################################################################################################
####################################################################################################################################
# FILE CLASSES AND METHODS

####################################################################################################################################
# METHOD FOR PARSING GADGET FILE OPTIONS
def pars_fileopt(fopt):
    """
    Returns parsed file option dictionary with keys:
        execpath: Str absolute path to location of compiled GADGET executable
        mfilpath: Str absolute path to location of compiled GADGET executable
        parmpath: Str absolute path to location of default GADGET parameter file
    """
    # Define keys
    form={
        'execpath': mmlpars.parsdict(default='/home/langmm/bin/Gadget/Gadget-2.0.6/Gadget2/Gadget2' ,type=str),
        'mfilpath': mmlpars.parsdict(default='/home/langmm/bin/Gadget/Gadget-2.0.6/Gadget2/Makefile',type=str),
        'parmpath': mmlpars.parsdict(default='/home/langmm/bin/Gadget/default.param'                ,type=str)
        }
    # Pars input
    fopt=mmlpars.mml_formpars(fopt,form)
    # Return parsed input
    return fopt

####################################################################################################################################
# METHOD TO CREATE DICTIONARY OF GADGET FILES
def files(fdict):
    """
    Returns a dictionary of GADGET files & directories with keys:
      dir:      Directory containing all GADGET simulation files (input/output)
      exdeft:   Location of default compiled GADGET executable
      mfdeft:   Location of Makefile used to compile default GADGET executable
      pmdeft:   Location of parameter file used to load default values
      icstat:   Location of static IC file for the simulation
      icevol:   Location of evolved IC file for the simulation
      pmstat:   Location of static parameter file for the simulation
      pmevol:   Location of parameter file for evolved IC file
      soften:   Location of particle softenings in cm
      input:    Directory containing all GADGET simulation input files
        exec:      Compiled GADGET used to produce run
        makefile:  Makefile used to compile GADGET
        ic:        Initial conditions files
        param:     Parameter file used by PBS to produce run
        repar:     Parameter file used by REPBS to produce run
        snpar:     Parameter file used by SNPBS to produce run
        pbs:       Submission script to start sim from an IC file
        repbs:     Submission script to start sim from RESTART files
        snpbs:     Submission script to start sim from SNAPSHOT file
        resub:     Script to resubmit a simulation after it finishes
        outlist:   List of times to output snapshots at
      output:   Directory containing all GADGET simulation output files
        usedparam: Parameters used during run
        runout:    Runtime output for submitted job
        snapbase:  Snapshot base string
        restbase:  Restart file base string
        energy:    Energy statistics
        cpu:       CPU statistics
        info:      Info file
        timings:   Timing info
    """
    # Pars input
    fdict=mmlpars.mml_pars(fdict,type=dict)
    rundir=fdict['rundir']
    pfix=fdict['pfix']
    pfix_shrt=fdict['pfix_shrt']
    options=pars_fileopt({})
    # Initialize output dictionary
    files=dict(fdict['gadget'],
               exdeft=options['execpath'],mfdeft=options['mfilpath'],pmdeft=options['parmpath'],
               soften=os.path.join(rundir,pfix_shrt+'softenings'))
    # Add input files
    infiles={}
    infiles['dir']=os.path.join(files['dir'],'input')
    infiles['subdirs']=[]
    fkeys_in=['exec','makefile','ic','param','repar','snpar','pbs','repbs','snpbs','resub','outlist','sched']
    flist_in=['Gadget2','Makefile','ic_gadget','param','repar','snpar','pbs','repbs','snpbs','resub','output_times.txt','schednum']
    infiles=mmlfiles.filedict(infiles['dir'],filekeys=fkeys_in,filelist=flist_in,pfix=pfix,fdict=infiles)
    files['input']=infiles
    # Add output files
    outfiles={}
    outfiles['dir']=os.path.join(files['dir'],'output')
    outfiles['subdirs']=[]
    outfiles['usedparam']=os.path.join(outfiles['dir'],'parameters-usedvalues')
    outfiles['usedparam_param']=infiles['param']+'-usedvalues'
    outfiles['usedparam_repar']=infiles['repar']+'-usedvalues'
    outfiles['usedparam_snpar']=infiles['snpar']+'-usedvalues'
    fkeys_out=['runout','snapbase','restbase','energy','cpu','info','timings']
    flist_out=['output','snapshot','restart','energy','cpu','info','timings']
    outfiles=mmlfiles.filedict(outfiles['dir'],filekeys=fkeys_out,filelist=flist_out,pfix=pfix,fdict=outfiles)
    files['output']=outfiles
    # Return file dictionary
    fdict['gadget']=files
    return fdict

####################################################################################################################################
# METHOD TO RETURN LISTS OF GADGET FILES
def get_filelist(fdict,keylist):
    """
    Returns a list of relavent GADGET files
    """
    # Pars input
    fdict=mmlpars.mml_pars(fdict,type=dict)
    keylist=mmlpars.mml_pars(keylist,type=list)
    # Add files
    filelist={}
    if 'default'   in keylist:
        fkeylist=['exdeft','mfdeft','pmdeft']
        filelist['default' ]=[fdict[ifkey] for ifkey in fkeylist]
    if 'static'    in keylist:
        fkeylist=['icevol','pmevol','soften']
        filelist['static'  ]=[fdict[ifkey] for ifkey in fkeylist]
    if 'longin'    in keylist:
        filelist['longin'  ]=[fdict['input']['ic']]
    if 'shortout'  in keylist:
        fkeylist=['runout','energy','cpu','info','timings',
                  'usedparam','usedparam_param','usedparam_repar','usedparam_snpar']
        filelist['shortout']=[fdict['output'][ifkey] for ifkey in fkeylist]
    # GADGET specific file types
    if 'exec'      in keylist:
        fkeylist=['exec','makefile']
        filelist['exec'    ]=[fdict['input'][ifkey] for ifkey in fkeylist]
    if 'param'     in keylist:
        fkeylist=['param','repar','snpar','pbs','repbs','snpbs','outlist','sched']
        filelist['param'   ]=[fdict['input'][ifkey] for ifkey in fkeylist]
    if 'snapshot'  in keylist:
        filelist['snapshot']=[fdict['output']['snapbase']+'*']
    if 'restart'   in keylist:
        filelist['restart' ]=[fdict['output']['restbase']+'*']
    if 'rest_bak'  in keylist:
        filelist['rest_bak']=[fdict['output']['restbase']+'*.bak']
    if 'shortsnapshot' in keylist:
        snapfiles=glob.glob(fdict['output']['snapbase']+'*')
        if len(snapfiles) > 0: filelist['shortsnapshot']=[sorted(snapfiles)[-1]]
        else                 : filelist['shortsnapshot']=[]
    # Return output
    return filelist
    
####################################################################################################################################
####################################################################################################################################
# METHODS TO CREATE GADGET FILES

####################################################################################################################################
# METHOD TO PARS DICTIONARY OF SIMULATION PARAMETER OPTIONS
def pars_paropt(popt):
    """
    Returns a parsed dictionary of simulation parameter options with keys:
        tsim:      Float length of the simulation in GADGET units (Gyr)
        dtsnap:    Float time between snapshots in GADGET units (Gyr)
        resub:     Bool specifying whether or not a run should be resubmitted.
        icmeth:    Integer specifying what type of file to start from (0:ic,1:restart,2:snapshot)
        outlist:   Bool specifying if output times should be read from a file.
        outspace:  Str specifying what spacing the outputlist should have ("LIN" or "LOG")
        periodic:  Bool specifying if boundary conditions are periodic
        zoom:      Bool specifying if simulation is a zoom
        comovint:  Bool specifying if comoving integration is turned on
        treepm:    Bool specifying if simulation is treepm
        dblfloat:  Bool specifying if double floats should be used
        pmgrid:    Float size of particle-mesh grid
        partalloc: Float particle allocation factor
        treealloc: Float tree allocation factor
        combuff:   Float communication buffer
        nsph:      Integer of gas smoothing neighbors
        gastemp:   Float initial temperature of gas particles
        gasflag:   Bool specifying if gas is included in the simulation
    """
    # Define values
    outspLIST=['LIN','LOG']
    form={
        'tsim'     : mmlpars.parsdict(default=1.0  ,type=float,min=0.0     ),
        'dtsnap'   : mmlpars.parsdict(default=0.01 ,type=float,min=0.0     ),
        'resub'    : mmlpars.parsdict(default=False,type=bool              ),
        'icmeth'   : mmlpars.parsdict(default=0    ,type=int  ,list=[0,1,2]),
        'outlist'  : mmlpars.parsdict(default=False,type=bool              ),
        'outspace' : mmlpars.parsdict(default='LIN',type=str,list=outspLIST),
        'periodic' : mmlpars.parsdict(default=False,type=bool              ),
        'zoom'     : mmlpars.parsdict(default=False,type=bool              ),
        'comovint' : mmlpars.parsdict(default=False,type=bool              ),
        'treepm'   : mmlpars.parsdict(default=False,type=bool              ),
        'dblfloat' : mmlpars.parsdict(default=False,type=bool              ),
        'pmgrid'   : mmlpars.parsdict(default=0    ,type=int  ,min=0       ),
        'partalloc': mmlpars.parsdict(default=1.5  ,type=float,min=0.0     ),
        'treealloc': mmlpars.parsdict(default=0.8  ,type=float,min=0.0     ),
        'combuff'  : mmlpars.parsdict(default=25.0 ,type=float,min=0.0     ),
        'nsph'     : mmlpars.parsdict(default=0    ,type=int  ,min=0       ),
        'gastemp'  : mmlpars.parsdict(default=1.0e4,type=float,min=0.0     ),
        'gasflag'  : mmlpars.parsdict(default=False,type=bool              )
        }
    # Pars input
    popt=mmlpars.mml_formpars(popt,form)
    # Return parsed input
    return popt

####################################################################################################################################
# METHOD TO CREATE STATIC GADGET FILES
def mkstatic(simstr,**options):
    """
    Creates static GADGET files
    """
    # Get file names
    files.simstr.fdict
    mmlio.verbose('Creating static GADGET files...',border=True)
    # Create isolated galaxy files
    if simstr['runtyp']=='galsim':
        # Create unevolved static ic
        if mmlio.yorn('Create new static unevolved IC and parameter files?'):
            mkic(simstr,'galsim',**options)
        # Create evolved static ic
        if mmlio.yorn('Create new static evolved IC and parameter files?'):
            newsnap=files['gadget']['icevol']
            newparm=files['gadget']['pmevol']
            endsnap=sorted(glob.glob(files['gadget']['output']['snapbase']+'*'),reverse=True)[0]
            endparm=files['gadget']['output']['usedparam']
            endmake=files['gadget']['input']['makefile']
            # Snapshot
            nb=pNbody.Nbody(p_name=endsnap,ftype='gadget',unitsfile=endparm,makefile=endmake)
            nb.rename(newsnap)
            nb.update_flags(flag_ic=True)
            nb.write()
            print '    {}'.format(newsnap)
            # Parameter file
            mmlfiles.cpfiles(endparm,newparm,overwrite=True)
            print '    {}'.format(newparm)
    # Create interacting galaxy files
    elif simstr['runtyp']=='intsim':
        if mmlio.yorn('Create new static IC and parameter files?'):
            mkic(simstr,'intsim',**options)
    # Handle error and return control
    else: raise Exception('Invalid runtyp: {}'.format(simstr['runtyp']))
    return

####################################################################################################################################
# METHOD TO CREATE LOCAL COPIES OF GADGET FILES
def mkfiles(simstr,paropt=None,verbose=False,**extra_kw):
    """
    Creates local copies of GADGET files
    """
    if verbose:
        mmlio.verbose('Creating files for a GADGET run locally...',border=True)
    # Pars input
    paropt=pars_paropt(paropt)
    # Get file names
    files=simstr.fdict
    mmlfiles.mkdirs(files['gadget']['input']['dir'])
    # IC file
    if mmlio.yorn('Create new local copy IC file?'):
        if not os.path.isfile(files['icgen']['gdic']):
            raise Exception('Static IC file does not exist. Create it first: {}'.format(files['icgen']['gdic']))
        mmlfiles.cpfiles(files['icgen']['gdic'],files['gadget']['input']['ic'],overwrite=True)
    # Executable file
    if mmlio.yorn('Create new local copy of GADGET executable?'):
        if not os.path.isfile(files['gadget']['exdeft']+'_'+simstr['runcomp']):
            raise Exception('Default executable does not exist. Create it first: {}'.format(files['gadget']['exdeft']))
        if not os.path.isfile(files['gadget']['mfdeft']+'_'+simstr['runcomp']):
            raise Exception('Default Makefile does not exist. Create it first: {}'.format(files['gadget']['mfdeft']))
        mmlfiles.cpfiles(files['gadget']['exdeft']+'_'+simstr['runcomp'],files['gadget']['input']['exec'],overwrite=True)
        mmlfiles.cpfiles(files['gadget']['mfdeft']+'_'+simstr['runcomp'],files['gadget']['input']['makefile'],overwrite=True)
    # Parameter files
    if mmlio.yorn('Create GADGET submission files?'):
        owparam=mmlio.yorn('Overwrite existing files?')
        tcpu=mmlsim2tcpu(simstr)
        for icmeth in range(3):
            pardict=mkparam(simstr,icmeth=icmeth,overwrite=owparam,tcpu=tcpu)
            pbsdict=mkpbs(simstr,icmeth=icmeth,overwrite=owparam,tcpu=tcpu)
        nout=int((pardict['TimeMax']-pardict['TimeBegin'])/pardict['TimeBetSnapshot'])
        outlist=mkoutlist(pardict['TimeBegin'],pardict['TimeMax'],nout,
                          files['gadget']['input']['outlist'],
                          method=paropt['outspace'],overwrite=owparam)

####################################################################################################################################
# METHOD TO CREATE GADGET PARAMETER FILE
def mkparam(simstr,icmeth=0,overwrite=False,parstat=None,parfile=None,
            softdict=None,tcpu=None,**param_kw):
    """
    Creates a GADGET parameter file and returns content as a dictionary
    """
    # Set constants
    parkeylist=['param','repar','snpar']
    keywidth=29
    # Get file names
    files_local=simstr.mkfiledict(checkaccess=True)
    files_accre=simstr.mkfiledict(memcomp=simstr['runcomp'],checkaccess=False)
    # Pars input
    icmeth=mmlpars.mml_pars(icmeth,type=int,list=[0,1,2])
    parfile=mmlpars.mml_pars(parfile,default=files_local['gadget']['input'][parkeylist[icmeth]],type=str)
    parstat=mmlpars.mml_pars(parstat,default=files_local['icgen']['gdpm'],type=str)
    overwrite=mmlpars.mml_pars(overwrite,type=bool)
    softdict=mmlpars.mml_pars(softdict,default={'empty':0},type=dict)
    if not isinstance(tcpu,float): tcpu=mmlsim2tcpu(simstr)
    # Load existing file if overwrite False
    if not overwrite and os.path.isfile(parfile):
        mmlio.verbose('Loading existing GADGET parameter file:',addval=parfile)
        params=mmlio.rwdict('R',parfile)
    # Otherwise create parameter file from default
    else:
        mmlio.verbose('Creating new GADGET parameter file:',addval=parfile)
        # If this is a normal parameter file, write it from scratch
        if icmeth==0:
            params=mmlio.rwdict('R',parstat)
            # Set file names
            params['OutputDir'         ]=files_accre['gadget']['output']['dir']
            params['SnapshotFileBase'  ]=os.path.basename(files_accre['gadget']['output']['snapbase'])
            params['InitCondFile'      ]=files_accre['gadget']['input']['ic']
            params['RestartFile'       ]=os.path.basename(files_accre['gadget']['output']['restbase'])
            params['InfoFile'          ]=os.path.basename(files_accre['gadget']['output']['info'])
            params['TimingsFile'       ]=os.path.basename(files_accre['gadget']['output']['timings'])
            params['CpuFile'           ]=os.path.basename(files_accre['gadget']['output']['cpu'])
            params['EnergyFile'        ]=os.path.basename(files_accre['gadget']['output']['energy'])
            params['ResubmitCommand'   ]=files_accre['gadget']['input']['resub']
            params['OutputListFilename']=files_accre['gadget']['input']['outlist']
            # Set communication buffer
            params['BufferSize'        ]=mmlsim2combuff(simstr)
            # Set timings
            params['TimeLimitCPU'      ]=60.0*60.0*tcpu
            # Read softenings from file if it exists
            if 'empty' in softdict and os.path.isfile(files_local['gadget']['soften']):
                softs=mmlio.rwdict('R',files_local['gadget']['soften'])
                softdict['gas'  ]=softs['Gas'  ]/params['UnitLength_in_cm']
                softdict['halo' ]=softs['Halo' ]/params['UnitLength_in_cm']
                softdict['disk' ]=softs['Disk' ]/params['UnitLength_in_cm']
                softdict['bulge']=softs['Bulge']/params['UnitLength_in_cm']
                softdict['stars']=softs['Stars']/params['UnitLength_in_cm']
                softdict['bndry']=softs['Bndry']/params['UnitLength_in_cm']
                del softdict['empty']
            if 'empty' not in softdict:
                params['SofteningGas'         ]=softdict['gas'  ]
                params['SofteningHalo'        ]=softdict['halo' ]
                params['SofteningDisk'        ]=softdict['disk' ]
                params['SofteningBulge'       ]=softdict['bulge']
                params['SofteningStars'       ]=softdict['stars']
                params['SofteningBndry'       ]=softdict['bndry']
                params['SofteningGasMaxPhys'  ]=softdict['gas'  ]
                params['SofteningHaloMaxPhys' ]=softdict['halo' ]
                params['SofteningDiskMaxPhys' ]=softdict['disk' ]
                params['SofteningBulgeMaxPhys']=softdict['bulge']
                params['SofteningStarsMaxPhys']=softdict['stars']
                params['SofteningBndryMaxPhys']=softdict['bndry']
        # If this a parameter file for a restart file start, copy normal parameter file
        elif icmeth==1:
            params=mmlio.rwdict('R',files_local['gadget']['input']['param'])
        # Otherwise, load normal parameter file and adjust it
        elif icmeth==2:
            params=mmlio.rwdict('R',files_local['gadget']['input']['param'])
            # Set IC file to last snapshot written
            snaplist=glob.glob(os.path.join(params['OutputDir'],params['SnapshotFileBase']+'*'))
            if len(snaplist)==0:
                params['InitCondFile']=os.path.join(params['OutputDir'],params['SnapshotFileBase']+'_000')
            else:
                params['InitCondFile']=snaplist[-1]
    # Scale softenings
    softkeyList=['SofteningGas','SofteningHalo','SofteningDisk','SofteningBulge','SofteningStars','SofteningBndry',
                 'SofteningGasMaxPhys','SofteningHaloMaxPhys','SofteningDiskMaxPhys','SofteningBulgeMaxPhys',
                 'SofteningStarsMaxPhys','SofteningBndryMaxPhys']
#    for isoftkey in softkeyList: params[isoftkey]*=softfact
    # Override based on input keys
    for ikey in param_kw.keys():
        if params.has_key(ikey):
            print 'Overriding GADGET parameter file option {} with value {}'.format(ikey,param_kw[ikey])
            params[ikey]=param_kw[ikey]
    # Write parameter dictionary to file
    mmlio.rwdict('W',parfile,params,width=keywidth,overwrite=overwrite)
    # Return output
    return params


####################################################################################################################################
# METHOD TO CREATE PBS SUBMISSION SCRIPT
def mkpbs(simstr,icmeth=0,overwrite=False,tcpu=None,ppn=None):
    """
    Creates a GADGET pbs script and returns content as a dictionary
    """
    # Set constants
    pbskeylist=['pbs','repbs','snpbs']
    parkeylist=['param','repar','snpar']
    # Pars input
    icmeth=mmlpars.mml_pars(icmeth,type=int,list=[0,1,2])
    overwrite=mmlpars.mml_pars(overwrite,type=bool)
    ppn=mmlpars.mml_pars(ppn,default=2,type=int,range=[1,4])
    # Get file names
    files_local=simstr.mkfiledict(checkaccess=True)
    files_accre=simstr.mkfiledict(memcomp=simstr['runcomp'],checkaccess=False)
    pbsfile=files_local['gadget']['input'][pbskeylist[icmeth]]
    parfile=files_local['gadget']['input'][parkeylist[icmeth]]
    pbsdict=mmlio.pbsdict()
    # Load existing file if overwrite False
    if not overwrite and os.path.isfile(pbsfile):
        mmlio.verbose('Loading existing GADGET pbs script:',addval=pbsfile)
        pbsdict.rwmpi('R',pbsfile)
    # Otherwise create file from scratch
    else:
        mmlio.verbose('Creating new GADGET pbs script:',addval=pbsfile)
        # Update memory
#        simstr.get_memory('gadget',update=True)
        if not isinstance(tcpu,float): tcpu=mmlsim2tcpu(simstr)
        mcpu=mmlsim2mcpu(simstr)
        # Variables from simstr
        pbsdict['ppn'      ]=ppn
        pbsdict['nodes'    ]=simstr['nprocgd']/ppn
        pbsdict['pmem'     ]=mcpu
        pbsdict['mem'      ]=mcpu*simstr['nprocgd']
        pbsdict['walltime' ]=tcpu
        # Varaibles from path
        pbsdict['simout'   ]=files_accre['gadget']['output']['runout']
        pbsdict['rundir'   ]=os.path.dirname(files_accre['gadget']['input']['exec'])
        pbsdict['execpath' ]=os.path.basename(files_accre['gadget']['input']['exec'])
        pbsdict['execflag' ]=icmeth
        pbsdict['execinput']=os.path.basename(parfile)
        pbsdict.rwmpi('W',pbsfile,overwrite=overwrite)
    # Return output
    return pbsdict

####################################################################################################################################
# METHOD TO CREATE LIST OF OUTPUT TIMES
def mkoutlist(tbeg,tend,nout,fname,method=None,overwrite=False):
    '''
    NAME:
        mmlgadget.mkoutlist
    PURPOSE:
        To create a file containing a list of output times for a GADGET sim.
    CALLING:
        outlist=mkoutlist(tbeg,tend,nout,fname[,method=None,overwrite=False])
    ARGUMENTS:
        tbeg:      Float time that simulation starts
        tend:      Float time that simulation ends
        nout:      Int number of simulation output times between tbeg & tend
        fname:     Str absolute path of file to write output list to
    KEYWORDS:
        method:    Str specifying if list should be linear ("LIN") or log ("LOG")
        overwrite: Bool specifying if existing file should be overwritten
    OUTPUT:
        outlist:   List of float output times
    '''
    # Pars input
    tbeg=mmlpars.mml_pars(tbeg,type=float,min=0.0)
    tend=mmlpars.mml_pars(tend,type=float,min=tbeg)
    nout=mmlpars.mml_pars(nout,type=int,min=1)
    fname=mmlpars.mml_pars(fname,type=str)
    method=mmlpars.mml_pars(method,type=str,list=['LIN','LOG'],default='LIN')
    overwrite=mmlpars.mml_pars(overwrite,type=bool)
    # If overwrite not set and file exists, load data from it
    if not overwrite and os.path.isfile(fname):
        fid=open(fname,'r')
        outlist=[]
        for line in fid:
            outlist.append(float(line.split()))
        fid.close()
    # Otherwise write data to file
    else:
        if os.path.isfile(fname): os.remove(fname)
        if method=='LIN':
            outlist=list(np.linspace(tbeg,tend,nout))
        elif method=='LOG':
            outlist=list(np.logspace(tbeg,tend,nout))
        fid=open(fname,'w')
        for itime in outlist:
            fid.write(repr(itime)+'\n')
        fid.close()
    # Return output
    return outlist

####################################################################################################################################
####################################################################################################################################
# METHOD TO READ GADGET FILES

####################################################################################################################################
# METHOD TO READ SOFTENINGS FROM GADGET PARAMETER FILE
def read_softenings(fname):
    """
    Reads and returns softening lengths from a GADGET parameter file as a dictionary with type keys
    """
    # Read parameters from parameter file
    params=mmlio.rwdict('R',fname)
    # Add softenings to dictionary
    soft={'gas'  : params['SofteningGas'         ],
          'halo' : params['SofteningHalo'        ],
          'disk' : params['SofteningDisk'        ],
          'bulge': params['SofteningBulge'       ],
          'stars': params['SofteningStars'       ],
          'bndry': params['SofteningBndry'       ]
          }
    # Return dictionary of softening lengths
    return soft

####################################################################################################################################
# METHOD TO READ GADGET MAKE FILE
def read_makefile(fname):
    '''
    NAME:
        mmlgadget.read_makefile
    PURPOSE:
        To read options from a GADGET makefile
    CALLING:
        makeopt=read_makefile(fname)
    ARGUMENTS:
        fname:   Str absolute path to makefile that should be read
    OUTPUT:
        makeopt: Dictionary of makefile options with option/value pairs.
    '''
    # Check if file is valid
    if not os.path.isfile(fname):
        raise Exception('Makefile name is not valid: {}'.format(fname))
    # Open file
    fid=open(fname,'r')
    # Check each line for options
    makeopt={}
    for line in fid:
        if line.startswith('OPT'):
            if line.startswith('OPTIMIZE') or line.startswith('OPTIONS'): continue
            vars=line.split('-D')
            keyval=vars[-1].split('=')
            if len(keyval) == 2:
                key=keyval[0]
                val=mmlstring.str2val(keyval[1].split('\n')[0])
            else:
                key=keyval[0].split('\n')[0]
                val=True
            makeopt[key]=val
#            print '{} = {}'.format(key,val)
        elif line.startswith('#OPT'):
            if line.startswith('#OPTIMIZE') or line.startswith('#OPTIONS'): continue
            vars=line.split('-D')
            keyval=vars[-1].split('=')
            if len(keyval) == 2:
                key=keyval[0]
            else:
                key=keyval[0].split('\n')[0]
            val=False
            makeopt[key]=val
#            print '{} = {}'.format(key,val)
        else:
            continue
    fid.close()

    return makeopt

####################################################################################################################################
# METHOD TO READ/PLOT ENERGY FILES
def read_energy(fname):
    '''
    NAME:
        mmlgadget.read_energy
    PURPOSE:
        To read in data from a GADGET energy file.
    CALLING:
        edict=read_energy(fname)
    ARGUMENTS:
        fname: Str absolute path to the GADGET energy file that should be read.
    OUTPUT:
        Dictionary of information read from FNAME with key/value pairs:
            t:  List of times
            ei: List of internal energies
            ek: List of kinetic energies
            ep: List of potential energies
    '''
    if not os.path.isfile(fname):
        raise Exception('Energy file does not exist: {}'.format(fname))
    # Initialize
    colKeys=['t','ei','ep','ek',
             'eiGas','epGas','ekGas',
             'eiHalo','epHalo','ekHalo',
             'eiDisk','epDisk','ekDisk',
             'eiBulge','epBulge','ekBulge',
             'eiStars','epStars','ekStars',
             'eiBndry','epBndry','ekBndry',
             'massGas','massHalo','massDisk','massBulge','massStars','massBndry']
    edict={ikey:[] for ikey in colKeys}
    # Open file
    fid=open(fname,'r')
    # Read input
    for iline in fid:
        colList=iline.split() ; ncols=len(colList)
        if ncols != 28:
            raise Exception('The number of columns ({}) is not the number expected ({}).'.format(ncols,28))
        for icol in range(ncols):
            edict[colKeys[icol]].append(float(colList[icol]))
    # Return output
    return edict


####################################################################################################################################
####################################################################################################################################
# METHODS FOR PLOTTING GADGET STUFF

####################################################################################################################################
# METHOD FOR PLOTTING GADGET ENERGY STATS
def plot_energy(fname,plotfile=None,overwrite=False,**plot_kw):
    '''
    NAME:
        mmlgadget.plot_energy
    PURPOSE:
        To plot data from a GADGET energy file
    CALLING:
        plot_energy(fname[,plotfile=None,overwrite=False,**plot_kw])
    ARGUMENTS:
        fname:     Str absolute path to GADGET energy file
    KEYWORDS:
        plotfile:  Str absolute path to file where plot should be saved
        overwrite: Bool specifying if existing plot should be overwritten
        Addtional keywords are passed on to the plotting routine (simplot.plot_energy)
    '''
    # Pars input
    if os.path.isfile(plotfile):
        if overwrite:
            os.remove(plotfile)
        else:
            print 'Plot of GADGET energy statistics already exists: {}'.format(plotfile)
            return
    # Read in data
    edict=read_energy(fname)
    # Plot
    t=np.array(edict['t'],float)
    ek=np.array(edict['ek'],float)
    ep=np.array(edict['ep'],float)
    simplot.plot_energy(t,ek,ep,fname=plotfile,overwrite=overwrite,**plot_kw)

####################################################################################################################################
####################################################################################################################################
# METHODS FOR CALCULATING THINGS

####################################################################################################################################
# METHOD TO DETERMINE COMMUNICATION BUFFER
def mmlsim2combuff(simstr):
    """
    Estimates the communication buffer a GADGET simulation will need
    """
    # Get file dictionary
    files=simstr.fdict
    icfile=files['icgen']['gdic']
    # Get buffer size from IC file if it exists
    if os.path.isfile(icfile):
        statinfo=os.stat(icfile)
        combuff=float(statinfo.st_size)/(1.0e6)
    # Otherwise get it by scaling previous results
    else:
        combuff=float(simstr.ntot)*(24.0/560000.0)
    # Return rounded output
    return int(mmlmath.oom(combuff,nsig=2,method='CEIL'))
    
####################################################################################################################################
# METHOD TO DETERMINE MEMORY FROM MMLSIM OBJECT
def mmlsim2mcpu(simstr):
    """
    Estimates the amount of memory a GADGET simulation will need
    """
    # Get file dictionary
    files=simstr.fdict
    parfile=files['icgen']['gdpm']
    # Get parameters from parameter file if it exists
    if os.path.isfile(parfile):
        params=mmlio.rwdict('R',parfile)
        partalloc=float(params['PartAllocFactor'])
        treealloc=float(params['TreeAllocFactor'])
        combuff=int(params['BufferSize'])
        perflag=params['PeriodicBoundariesOn']
        if perflag==1: periodic=True
        else:          periodic=False
    # Otherwise force defaults
    else:
        partalloc=None
        treealloc=None
        combuff=mmlsim2combuff(simstr)
        periodic=None
    # Calculate memory based on GADGET user manual
    mcpu=calcmcpu(simstr.ntot,nproc=simstr['nprocgd'],
                  partalloc=partalloc,treealloc=treealloc,
                  combuff=combuff,periodic=periodic)
    # Force reasonable values as a safety
    if   simstr['runtyp']=='galsim': mcpu=1000. # mb
    elif simstr['runtyp']=='intsim': mcpu=2000. # mb
    else: raise Exception('Invalid run type: {}'.format(simstr['runtyp']))
    # Return rounded output
    return mmlmath.oom(mcpu,nsig=1,method='CEIL')

####################################################################################################################################
# METHOD TO DETERMINE TIME REQUIRED TO RUN GADGET
def mmlsim2tcpu(simstr):
    """
    Estimates the amount of time required to run GADGET on ACCRE
    """
    # Force reasonable values as a safety
    tcpu=mmlio.askquest('What should be the overall time limit?',default=48.,dtype='float')
    # Return output
    return tcpu
    
####################################################################################################################################
# METHOD TO DETERMINE MEMORY REQUIRED TO RUN GADGET
def calcmcpu(ntot,nproc=None,nsph=None,partalloc=None,treealloc=None,combuff=None,
             periodic=None,double=None,treepm=None,pmgrid=None,zoom=None):
    '''
    NAME:
        mmlgadget.calcmcpu
    PURPOSE:
        To determine the amound of memory a GADGET simulation would need.
    CALLING:
        mcpu=calcmcpu(ntot[,nproc=,nsph=,partalloc=,treealloc=,combuff=,
                      periodic=,double=,treepm=,pmgrid=,zoom)
    ARGUMENTS:
        ntot:      Long number of particles
    KEYWORDS:
        nproc:     Int number of processors [DEFAULT=1]
        nsph:      Long number of SPH particles [DEFAULT=0]
        partalloc: Float particle allocation factor [DEFAULT=1.5]
        treealloc: Float tree allocation factor [DEFAULT=0.8]
        combuff:   Int ommunication buffer [DEFAULT=25]
        periodic:  Bool specifying if boundaries are periodic
        double:    Bool specifying if calculations are double floats
        treepm:    Bool specifying if simulation is particle-mesh
        pmgrid:    Int size of particle-mesh grid [DEFAULT=0]
        zoom:      Bool specifying if simulation is a zoom
    OUTPUT:
        mcpu:      Int memory required on each processor (MB)
    '''
    # Pars input
    ntot=mmlpars.mml_pars(ntot,type=long,min=0L)
    nproc=mmlpars.mml_pars(nproc,default=1,type=int,min=0)
    nsph=mmlpars.mml_pars(nsph,default=0L,type=long,min=0L)
    partalloc=mmlpars.mml_pars(partalloc,default=1.5,type=float,min=0)
    treealloc=mmlpars.mml_pars(treealloc,default=0.8,type=float,min=0)
    combuff=mmlpars.mml_pars(combuff,default=25,type=int,min=0)
    periodic=mmlpars.mml_pars(periodic,default=False,type=bool)
    double=mmlpars.mml_pars(double,default=False,type=bool)
    treepm=mmlpars.mml_pars(treepm,default=False,type=bool)
    pmgrid=mmlpars.mml_pars(pmgrid,default=0,type=int)
    zoom=mmlpars.mml_pars(zoom,default=False,type=bool)
    # Additional factors
    if double:
        dbl=2.0
    else:
        dbl=1.0
    # General particle stuff
    M1=dbl*float(ntot)*partalloc*(68.0+64.0*treealloc)
    # SPH particle stuff
    M2=dbl*float(nsph)*partalloc*84.0
    # Communication buffer
    M3=float(nproc)*float(combuff)*(1024.0**2)
    # TreePM stuff
    if treepm:
        if periodic:
            M4=dbl*((float(pmgrid)**3)*12.0-16.0)
            if zoom:
                M5=dbl*((2.0*float(pmgrid))**3)*16.0
            else:
                M5=0.0
        else:
            M4=0.0
            if zoom:
                M5=dbl*((2.0*float(pmgrid))**3)*20.0
            else:
                M5=dbl*((2.0*float(pmgrid))**3)*16.0
    else:
        M4=0.0
        M5=0.0
    # Summ memory & return in MB
    mcpu=(M1+M2+M3+M4+M5)/(1024.0**2)
    # Divide accros processors
    return mcpu/float(nproc)


####################################################################################################################################
####################################################################################################################################
# PROVIDE COMMAND LINE ACCESS
if __name__ == '__main__': main()
