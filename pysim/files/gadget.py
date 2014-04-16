#!/usr/bin/python
####################################################################################################################################
#
# MEAGAN LANG'S GADGET METHODS
#
####################################################################################################################################
import sys,os,shutil,glob,copy,pprint,scipy,math,struct,fnmatch
import numpy as np
import matplotlib as mplib
from mmlutils import *

LIST_PTYPS=['gas','halo','disk','bulge','stars','bndry']
N_PTYP=len(LIST_PTYPS)

LIST_BLOCKS_STD=['HEAD','POS','VEL','ID','MASS']
LIST_BLOCKS_SPH=['U','RHO','HSML']
LIST_BLOCKS_OPT=['POT','ACCE','ENDT','TSTP']
LIST_BLOCKS=LIST_BLOCKS_STD+LIST_BLOCKS_SPH+LIST_BLOCKS_OPT

LIST_METHODS=['mk_ic','mk_files','setup','retrieve','clean']
LIST_METH_MKIC=['buildgal','snapshot','restart','interact','timetest']

FILELIST_SETUP=['exec','input','ic','restart']
FILELIST_OUTPT=['output','snapshot','restart','input']
FILELIST_CLEAN=['snapshot','restbak']

FILELIST_SPLIT={'stat_ic':'.{:d}','ic':'.{:d}','snapshot':'.{:d}'}

METHLIST_SETUP=[]
METHLIST_OUTPT=[]
METHLIST_CLEAN=[]

VALIDINITKEY='ic'
VALIDSNAPKEY='snapshot'
VALIDSNAPEXT=['000','*5']

INCLEXT_SETUP={}
INCLEXT_OUTPT={'inclext_end' :{'snapshot':'5'},
               'inclext_list':{'snapshot':VALIDSNAPEXT}}
INCLEXT_CLEAN={}
EXCLEXT_SETUP={}
EXCLEXT_OUTPT={}#'exclext_end':{'restart':'.bak'}}
EXCLEXT_CLEAN={}

FEXT_MAKEFILE=['Makefile','makefile','gdmk']
FEXT_PARMFILE=['param','gdpar']

def runtyp2runtime(runtyp):
    if   runtyp=='galsim': runtime=3.
    elif runtyp=='intsim': runtime=5.
    else: raise Exception('Unsupported runtyp: {}'.format(runtyp))
    return runtime

import pysim

####################################################################################################################################
####################################################################################################################################
# METHODS FOR REPORTING STUFF FROM FILES

####################################################################################################################################
# METHOD FOR DETERMINING HOW MANY FILES SHOULD BE USED
def get_nsplit(simstr,key,**exkw):
    """Determine how many files a file class is split into"""
    # Handle file that is singular
    if key not in FILELIST_SPLIT: nfile=1
    # Handle files that can be plit
    else:
        if   key in ['stat_ic','ic']: nfile=nsplit_ic(simstr,**exkw)
        elif key=='snapshot'        : nfile=nsplit_snaps(simstr,**exkw)
        else: raise Exception('Key {} is in FILELIST_SPLIT, but is not supported.'.format(key))
    # Return nfile (0 if singular and no ext added)
    if nfile>1: return nfile
    else      : return 0
def nsplit_ic(simstr,npartlim=2.0e7):
    """Determine how many files IC files are split into"""
    # Get stat_ic path
    fic=simstr.fdict['gadget']['icgen']['stat_ic']
    # Count files
    if   os.path.isfile(fic     ): nfile=1
    elif os.path.isfile(fic+'.0'): nfile=len(glob.glob(fic+'.*'))
    else: raise Exception('stat_ic file(s) do not exist. Why do you want to know?')
    # Return count
    return nfile
def nsplit_snaps(simstr,npartlim=2.0e7):
    """Determine how many files snapshots should be split into"""
    # Check if snapshots already exist
    fsnap=simstr.fdict['gadget']['output']['snapshot'].format(0)
    if   os.path.isfile(fsnap     ): nfile=1
    elif os.path.isfile(fsnap+'.0'): nfile=len(glob.glob(fsnap+'.*'))
    # Otherwise determine this based on number of particles per processor
    else:
        # Get variables
        npart=float(simstr.ntot)
        nproc=float(simstr['nprocev'])
        # Determine number per processor
        nfile=int(np.ceil(npart/npartlim))
        if nfile==1: return nfile
        # Ensure number is multiple of processors
        while nproc%nfile!=0: nfile+=1
    # Return number
    return nfile

####################################################################################################################################
# METHOD FOR GETTING EXECUTION TIME
def performance(simstr,stepmin=128,reqmin=False,unit='Gyr',**exkw):
    """
    Return execution time for 1 Gyr
    """
    # Pars input
    if unit=='Gyr': tscl=(3.085678e21/1.0e5)/(3.15569e16)
    else: raise Exception('Invalid unit: {}'.format(unit))
    # Load timing stats from files
    fout=simstr.fdict['gadget']['output']['cpu']
    print '    '+fout
    if not os.path.isfile(fout): raise Exception('Timings file does not exist.')
    cpustats=rw_cpustat('R',fout)
    # Identify time of last output
    steparr=np.array(cpustats['step'])
    tsimarr=np.array(cpustats['tsim'])
    if reqmin and steparr.max()<stepmin:
        raise Exception('Only {} steps, but {} required.'.format(steparr.max(),stepmin))
    idxend=np.arange(len(steparr))[steparr<=stepmin].max()
    if idxend==len(tsimarr): idxend-=1
    # Get cpu time 
    time=cpustats['tot'][idxend]*simstr['nprocev']
    print simstr.ntot,simstr['nprocev'],len(cpustats['tsim']),idxend,cpustats['tsim'][idxend]
    # Convert to per Gyr
    time/=(cpustats['tsim'][idxend]/tscl)
    # Return time
    return time,'s/'+unit

####################################################################################################################################
####################################################################################################################################
# METHODS TO CREATE GADGET FILES
 
####################################################################################################################################
# METHOD TO CREATE GADGET EXECUTABLE FILES
def mk_exec(simstr,**exkw):
    """
    Creates a local copy of executable files
    """
    # Select default files
    exdeft=simstr.finfo['gadget']['exec']['default']+'_'+simstr['runcomp']
    mfdeft=simstr.finfo['gadget']['makefile']['default']+'_'+simstr['runcomp']
    excopy=simstr.finfo['gadget']['exec']['path']
    mfcopy=simstr.finfo['gadget']['makefile']['path']
    # Check existence
    if not os.path.isfile(exdeft): raise Exception('Executable does not exist: {}'.format(exdeft))
    if not os.path.isfile(mfdeft): raise Exception('Makefile does not exist: {}'.format(mfdeft))
    # Copy files to local directory
    mmlfiles.cpfiles(exdeft,excopy,overwrite=True)
    mmlfiles.cpfiles(mfdeft,mfcopy,overwrite=True)
    # Return
    return excopy,mfcopy

####################################################################################################################################
# METHOD TO CREATE GADGET PARAMETER FILE
def mk_param(simstr=None,runtag=None,owmake=False,overwrite=None,askuser=False,mkw={}):
    """
    Creates a parameter file
    """
    # Pars input
    if simstr==None: simstr=pysim.loadsim(runtag)
    # Get file names
    fpar_stat=simstr.finfo['gadget']['stat_pm']['path']
    fpar_copy=simstr.fdict['gadget']['input']['param']
    # Check for overwrite
    if not isinstance(overwrite,bool):
        if os.path.isfile(fpar_stat):
            if askuser: overwrite=mmlio.yorn('Overwrite existing static parameter file?')
            else      : overwrite=mmlpars.mml_pars(overwrite,type=bool,default=False)
        else: overwrite=False
    # If the static parameter file dosn't exist, create it
    if overwrite or not os.path.isfile(fpar_stat):
        # Initialize make log
        mklog=mk_makelog(simstr,overwrite=owmake,**mkw)
        # Get default parameter file
        if mklog['baserun']==simstr['runtag']: 
            fpar_deft=simstr.finfo['gadget']['stat_pm']['default']
        else:
            fpar_deft=pysim.loadsim(mklog['baserun']).finfo['gadget']['stat_pm']['path']
        # Initialize static parameter file
        params=rw_param('R',fpar_deft)
        # Set file names
        files_run=simstr.fdict_run
        params['OutputDir'         ]=files_run['gadget']['output']['dir']
        params['SnapshotFileBase'  ]=os.path.basename(files_run['gadget']['output']['snapshot'].rstrip('_'+pysim.files.get_extfmt(files_run['gadget']['output']['snapshot'])))
        params['InitCondFile'      ]=files_run['gadget']['input']['ic']
        params['RestartFile'       ]=os.path.basename(files_run['gadget']['output']['restart'].rstrip('.'+pysim.files.get_extfmt(files_run['gadget']['output']['restart'])))
        params['InfoFile'          ]=os.path.basename(files_run['gadget']['output']['info'])
        params['TimingsFile'       ]=os.path.basename(files_run['gadget']['output']['timings'])
        params['CpuFile'           ]=os.path.basename(files_run['gadget']['output']['cpu'])
        params['EnergyFile'        ]=os.path.basename(files_run['gadget']['output']['energy'])
        params['ResubmitCommand'   ]=files_run['gadget']['input']['resub']
        params['OutputListFilename']=files_run['gadget']['input']['outlist']
        # Number of files per snapshot
        params['NumFilesPerSnapshot']=nsplit_snaps(simstr)
        # Run time
        if mklog['method']=='timetest':
            params['TimeMax']=1.01*params['TimeBetSnapshot']
            params['TimeLimitCPU']=3600.*1.
        else:
            time2Gyr=(params['UnitLength_in_cm']/params['UnitVelocity_in_cm_per_s'])/((3.15569e7)*(1.0e9))
            params['TimeMax']=runtyp2runtime(simstr['runtyp'])*time2Gyr
        # Set communication buffer
        if mklog['method']!='timetest':
            params['BufferSize']=mmlsim2combuff(simstr)
            params['PartAllocFactor']=mmlsim2partalloc(simstr)
            params['TreeAllocFactor']=mmlsim2treealloc(simstr)
        # Set softenings
        for ityp in LIST_PTYPS:
            params['Softening'+ityp.capitalize()]=mklog['softs'][ityp]
            params['Softening'+ityp.capitalize()+'MaxPhys']=mklog['softs'][ityp]
        # Save parameters
        rw_param('W',fpar_stat,params,overwrite=overwrite)
        # Check parameter file
        print fpar_stat
        mmlio.yorn('Make changes to the above parameter file now!')
    # Make local copy of parameter file
    mmlfiles.cpfiles(fpar_stat,fpar_copy,overwrite=True)
    # Return output
    return fpar_copy

####################################################################################################################################
# METHOD TO CREATE SUBMISSION SCRIPTS
def mk_pbs(simstr,**kwargs): return mk_submit(simstr,'pbs',**kwargs)
def mk_sbatch(simstr,**kwargs): return mk_submit(simstr,'sbatch',**kwargs)
def mk_submit(simstr,ftype,param=None,verbose=False,**kwargs):
    """
    Creates a submission script
    """
    # Pars input
    ftype=mmlpars.mml_pars(ftype,list=mmlio.DICT_SUBOPT.keys())
    # Get file names
    files_run=simstr.fdict_run
    subfile=simstr.fdict['gadget']['input'][ftype]
    # Get parameter file info
    if not isinstance(param,dict):
        param=mmlpars.mml_pars(param,default=simstr.fdict['gadget']['input']['param'],type=str)
        param=rw_param('R',param)
    # Get run command info
    exefile=os.path.basename(files_run['gadget']['input']['exec'])
    parfile=os.path.basename(files_run['gadget']['input']['param'])
    mpicmd=mmlinfo.COMPDICT_MPI[simstr['runcomp']]
    # Variables from simstr
    if 'timetest' in simstr['runtag']: name='gd'+simstr['subtag']+'_'+simstr['runtag']
    else                             : name='gd_'+simstr['runtag']+simstr['subtag']
    subdict={'jobname' :name,
             'nproc'   :int(simstr['nprocev']),
             'ppn'     :int(mmlinfo.COMPDICT_PPN[simstr['runcomp']]),
             'pmem'    :int(mmlsim2mcpu(simstr)),
             'twall'   :float(param['TimeLimitCPU'])/3600.,
             'outfile' :files_run['gadget']['output']['runout'],
             'cmdlist' :['mkdir -p {}'.format(files_run['gadget']['output']['dir']),
                         '{} -np {} ./{} {} 0'.format(mpicmd,simstr['nprocev'],exefile,parfile)]}
    # Steampede stuff
    if simstr['runcomp']=='stampede':
        subdict['account']=files.STAMPEDE_ACCT
    # Write submission script
    mmlio.rw_submit('W',subfile,ftype,subdict,verbflag=verbose,overwrite=True)
    # Create output directory
    mmlfiles.mkdirs(files_run['gadget']['output']['dir'],host=simstr['runcomp'])
    # Return output
    return subfile

####################################################################################################################################
def mk_resub(simstr,**kwargs):
    """
    Creates a resubmission script
    """
    # Get file names
    fname=simstr.fdict['gadget']['input']['resub']
    fname=False
    # Return output
    return fname

####################################################################################################################################
# METHODS TO CREATE LIST OF OUTPUT TIMES
def mk_outlist(simstr,method=None,tbeg=None,tend=None,askuser=False,**kwargs):
    """
    Creates a file containing a list of output times
    """
    methLIST=['lin','log'] 
    methDEF='lin'
    tbegDEF=0.
    tendDEF=10.0
    noutDEF=100
    # Pars input
    fname=simstr.fdict['gadget']['input']['outlist']
    parfile=simstr.fdict['gadget']['input']['param']
    # Get info from parameter file
    if os.path.isfile(parfile):
        param=rw_param('R',parfile)
        # Discontinue if not set
        if param['OutputListOn']==0: return False
        # Get parameters
        tbegDEF=param['TimeBegin']
        tendDEF=param['TimeMax']
    # Ask for info on output list
    if askuser:
        if method not in methLIST: method=mmlio.askselect('How should the timing of output be scaled?',methLIST)
        if not isinstance(tbeg,float): 
            tbeg=mmlio.askquest('What time (in Gyr) should the snapshots start at?',dtype='float',default=tbegDEF)
        if not isinstance(tend,float): 
            tend=mmlio.askquest('What time (in Gyr) should the snapshots end at?',dtype='float',default=tendDEF)
        if not isinstance(nout,int  ): 
            nout=mmlio.askquest('How many snapshot should be output between {} and {} Gyrs?'.format(tbeg,tend),dtype=int,default=noutDEF)
    else:
        method=mmlpars.mml_pars(method,list=methLIST,default=methDEF)
        tbeg=mmlpars.mml_pars(tbeg,type=float,default=tbegDEF)
        tend=mmlpars.mml_pars(tend,type=float,default=tendDEF)
        nout=mmlpars.mml_pars(nout,type=int,default=noutDEF)
    # Create outlist
    if   method=='lin': outlist=list(np.linspace(tbeg,tend,nout))
    elif method=='log': outlist=list(np.logspace(np.log10(tbeg),np.log10(tend),nout))
    else: raise Exception('Invalid method: {}'.format(method))
    # Write it to file
    fid=open(fname,'w')
    for itime in outlist:
        fid.write(repr(itime)+'\n')
    fid.close()
    # Return output
    return fname

####################################################################################################################################
# METHOD TO CREATE IC FILES
def mk_ic(simstr=None,runtag=None,askuser=False,overwrite=None,owmake=False,mkw={},use_pNbody=None,**kwargs):
    """
    Creates IC files
    """
    # Pars input
    if simstr==None: simstr=pysim.loadsim(runtag)
    # Get filenames
    fic_stat=simstr.fdict['gadget']['icgen']['stat_ic']
    fic_copy=simstr.fdict['gadget']['input']['ic']
    funitfile=simstr.fdict['gadget'  ]['icgen']['stat_pm']
    fic_flag=os.path.isfile(fic_stat) or os.path.isfile(fic_stat+'.0')
    # Check for overwrite
    if not isinstance(overwrite,bool):
        if fic_flag:
            if askuser: overwrite=mmlio.yorn('Overwrite existing static IC file?')
            else      : overwrite=mmlpars.mml_pars(overwrite,type=bool,default=False)
        else: overwrite=False
    # Check for pNbody
    if askuser and not isinstance(use_pNbody,bool):
        use_pNbody=mmlio.yorn('Use pNbody?')
    else:
        use_pNbody=mmlpars.mml_pars(use_pNbody,type=bool,default=False)
    if use_pNbody: import pNbody
    # If the static IC file dosn't exist, create it
    print '[gadget.mk_ic]',overwrite,fic_flag
    if overwrite or not fic_flag:
        # Initialize make log
        mklog=mk_makelog(simstr,overwrite=owmake,**mkw)
        method=mklog['method']
        # Create IC file
        if   method=='buildgal': 
            iunitfile=simstr.fdict['buildgal']['input']['units'  ]
            iudict=pysim.files.rw_units('R',iunitfile)
            fudict=pysim.files.rw_units('R',funitfile)
            # Load buildgal IC file
            pnb0=simstr.loadsnap(fname=simstr.fdict['buildgal']['icgen']['stat_ic'],ftype='buildgal',
                                 unitfile=iunitfile,pNbody=use_pNbody)
            # Convert format
            if use_pNbody:
                pnb=pNbody.Nbody(p_name=fic_stat,ftype='gadget',status='new',unitsfile=funitfile,
                                 mass=pnb0.mass*iudict['M'].ratio(fudict['M']),
                                 pos =pnb0.pos *iudict['L'].ratio(fudict['L']),
                                 vel =pnb0.vel *iudict['V'].ratio(fudict['V']),
                                 num =pnb0.num                                ,
                                 tpe =pnb0.tpe                                )
            else: pnb=pnb0
            # Reset time
            time0=pysim.nbody.array.SimArray(0.,'Gyr')
        elif method=='timetest':
            isimstr=pysim.loadsim(mklog['baserun'])
            mmlio.verbose(mklog['baserun'])
            # Get counts
            ntot0=isimstr.ntot
            ntot=mklog['ntot']
            # Load snapshot file
            pnb0=isimstr.loadsnap(fext=mklog['snapext'],ftype='gadget',pNbody=use_pNbody)
            # Reset time
            time0=pysim.nbody.array.SimArray(0.,'Gyr')
            # Downsample
            print 'old'
            print pnb0.header.npart
            print pnb0.header.mass
            print len(pnb0),'=',ntot0
            print ''
            pnb=pnb0.downsample(ntot)
            print 'new'
            print pnb.header.npart
            print pnb.header.mass
            print len(pnb),'=',ntot
            print ''
        elif method=='snapshot':
            isimstr=pysim.loadsim(mklog['baserun'])
            # Load snapshot file
            pnb=isimstr.loadsnap(fext=mklog['snapext'],ftype='gadget')
            mklog['snaptime']=float(pnb.properties['time'].in_units('Gyr'))
            # Reset time
            if mklog['resettime']: time0=pysim.nbody.array.SimArray(mklog['time0'],'Gyr')
            else:
                if use_pNbody: time0=pysim.nbody.array.SimArray(pnb.atime*(pnb.localsystem_of_units.convertionFactorTo(pNbody.units.Unit_Myr)/1e3),'Gyr')
                else         : time0=pnb.properties['time']
            # Heat disk
            pnb.heatdisk(mklog['heatfact'])
        elif method=='restart' : raise Exception('Not currently supported.')
        elif method=='interact':
            from pysim.nbody.analysis import langmm,halo
            # Load snapshots
            isimstr1=pysim.loadsim(simstr.infodict['gal1'])
            isimstr2=pysim.loadsim(simstr.infodict['gal2'])
            pnb1=isimstr1.loadsnap(fext=mklog['snapext1'],ftype='gadget',pNbody=use_pNbody)
            pnb2=isimstr2.loadsnap(fext=mklog['snapext2'],ftype='gadget',pNbody=use_pNbody)
            if use_pNbody:
                mklog['snaptime1']=pnb1.atime*pnb1.localsystem_of_units.convertionFactorTo(pNbody.units.Unit_Myr)/1e3
                mklog['snaptime2']=pnb2.atime*pnb2.localsystem_of_units.convertionFactorTo(pNbody.units.Unit_Myr)/1e3
            else:
                mklog['snaptime1']=float(pnb1.properties['time'].in_units('Gyr'))
                mklog['snaptime2']=float(pnb2.properties['time'].in_units('Gyr'))
            # Center galaxies
            if use_pNbody:
                pnb1.center(move=True)
                pnb2.center(move=True)
            else:
                halo.center(pnb1)
                halo.center(pnb2)
            # Get total mass
            if use_pNbody:
                M1=pnb1.mass_tot*pnb1.localsystem_of_units.convertionFactorTo(pNbody.units.Unit_Msol)
                M2=pnb2.mass_tot*pnb2.localsystem_of_units.convertionFactorTo(pNbody.units.Unit_Msol)
            else:
                M1=float(pnb1['mass'].sum().in_units('Msol'))
                M2=float(pnb2['mass'].sum().in_units('Msol'))
            R1=None ; R2=None
            # Add orbit to second galaxy
            if use_pNbody: pnb2.add_orbit(M1,M2,R1_pc=R1,R2_pc=R2,**mklog['mkinfo'])
            else         : langmm.add_orbit(pnb2,M1,M2,R1_pc=R1,R2_pc=R2,**mklog['mkinfo'])
            # Add galaxies
            pnb=pnb1+pnb2
            if use_pNbody: print pnb.nbody
            else         : print len(pnb)
            # Reset time
            time0=pysim.nbody.array.SimArray(mklog['time0'],'Gyr')
        else: raise Exception('Invalid method: {}'.format(method))
        # Reset time
        mklog['time0']=float(time0.in_units('Gyr'))
        if use_pNbody: pnb.atime=time0.in_units('Gyr')/(pnb.localsystem_of_units.convertionFactorTo(pNbody.units.Unit_Myr)/1e3)
        else         : pnb.properties['time']=time0 
        # Save gadget IC file
        mmlfiles.mkdirs(os.path.dirname(fic_stat))
        if overwrite and fic_flag: os.remove(fic_stat)
        if use_pNbody:
            pnb.rename(fic_stat)
            pnb.update_flags(flag_ic=True)
            pnb.write()
        else:
            pnb.write(filename=fic_stat,ftype='gadget',overwrite=overwrite,unitfile=funitfile,flag_ic=True)
        del pnb
        #mmlio.yorn(fic_stat)
    # Make local copy of IC file
    if os.path.isfile(fic_stat):
        mmlfiles.cpfiles(fic_stat,fic_copy,overwrite=True)
    elif os.path.isfile(fic_stat+'.0'):
        for ific_stat in glob.glob(fic_stat+'.*'):
            iext=ific_stat.split(fic_stat)[-1]
            if iext.endswith('~'): continue
            ific_copy=fic_copy+iext
            mmlfiles.cpfiles(ific_stat,ific_copy,overwrite=True)
    # Return
    return fic_copy

####################################################################################################################################
# METHOD TO CREATE MAKELOG
def mk_makelog(simstr,method=None,overwrite=False,testopt={},**exkw):
    """
    Creates a make log
    """
    fmklog=simstr.fdict['gadget']['icgen']['makelog']
    # Read make log
    if os.path.isfile(fmklog) and not overwrite: mklog=pysim.files.rw_makelog('R',fmklog)
    # Create it
    else:
        # Pars input
        if method not in LIST_METH_MKIC: method=mmlio.askselect('Pick a method for making ICs',LIST_METH_MKIC)
        # Initialize makelog 
        mklog=dict(ftype ='gadget',
                   method=method,
                   runtag=simstr['runtag'],
                   runtyp=simstr['runtyp'],
                   addtag=simstr.infodict['addtag'],
                   ntot  =simstr.get_ntot())
        mklog.update(**exkw)
        if   method=='buildgal': mklog.update(baserun=simstr['runtag'],
                                              nproc  =simstr['nprocmk'])
        elif method=='snapshot': mklog.update(baserun=simstr['runtag'].split(mklog['addtag']))
        elif method=='restart' : mklog.update(nproc  =simstr['nprocev'])
        elif method=='interact': mklog.update(baserun=simstr['runtag'])
        elif method=='timetest': 
            ifodict=copy.deepcopy(simstr.infodict) ; ifotag=ifodict.pop('tagstr',None)
            # if simstr['runtyp']=='galsim':
            #     if   simstr.ntot > long(1e7): ifodict['ntot']=long(1e8)
            #     elif simstr.ntot > long(1e6): ifodict['ntot']=long(1e7)
            #     else                        : ifodict['ntot']=long(1e6)
            mklog.update(baserun=mmlparam.par2tag(pysim.simlist,simstr['runtyp'],inpar=ifodict))
        else: raise Exception('Invalid method: {}'.format(method))
        if method=='timetest': simstr0=pysim.loadsim(mklog['parmrun'])
        else                 : simstr0=pysim.loadsim(mklog['baserun'])
        mklog['mkinfo']=simstr0.infodict
        # Set defaults
        mklogDEF={}
        if   method=='buildgal': 
            import buildgal
            mdict=buildgal.model2param(simstr0.infodict)
            mklogDEF.update(softs={ityp:mdict.get(ityp+'_eps',0.)/1000. for ityp in LIST_PTYPS})
        elif method=='snapshot': mklog.update(softs=rw_softenings(simstr0.fdict['gadget']['icgen']['stat_pm']))
        elif method=='interact':
            mklog['snaptime1']=0.
            mklog['snaptime2']=0.
            mklog['mkinfo']['incl']=simstr0.infodict['incl']*np.pi/180.
            mklogDEF.update(softs={})
            simstr1=pysim.loadsim(simstr.infodict['gal1'])
            simstr2=pysim.loadsim(simstr.infodict['gal2'])
            sdict1=rw_softenings(simstr1.fdict['gadget']['icgen']['stat_pm'])
            sdict2=rw_softenings(simstr2.fdict['gadget']['icgen']['stat_pm'])
            for ityp in LIST_PTYPS:
                print ityp,sdict1[ityp],sdict2[ityp]
                if sdict1[ityp]==sdict2[ityp]: mklogDEF['softs'][ityp]=sdict1[ityp]
                else:
                    mklogDEF['softs'][ityp]=mmlio.askselect('Select softening length:',[sdict1[ityp],sdict2[ityp]])
            #else: raise Exception('The softening lengths are not equal. Decide how to treat this.')
        elif method=='timetest':
            mklog.update(softs=rw_softenings(simstr0.fdict['gadget']['icgen']['stat_pm']))
            if testopt.get('adjsoft',False):
                mfrac=float(simstr0.ntot)/float(mklog['ntot'])
                mmlio.verbose('Adjusting softening by cube root of mass ratio ({})...'.format(mfrac))
                print '    Baserun: {}'.format(simstr0['runtag'])
                pprint.pprint(mklog['softs'])
                for ityp in LIST_PTYPS: mklog['softs'][ityp]*=(mfrac**(1./3.))
                print '    Testrun: {}'.format(simstr['runtag'])
                pprint.pprint(mklog['softs'])
            #print simstr.ntot
        else: raise Exception('Invalid method: {}'.format(method))
        # Pars makelog
        mklog=mmlparam.parspar('pysim.files','mklog_{}_{}'.format('gadget',method),inpar=mklog,defpar=mklogDEF,
                               askuser=True,load=False,save=False)
        # Save to file
        pysim.files.rw_makelog('W',fmklog,mklog,overwrite=overwrite)
    # Return 
    return mklog

####################################################################################################################################
####################################################################################################################################
# METHOD TO TEST FILES
def test_snap(fname):
    """
    Tests validity of a snapshot
    """
    # Try to open file
    try: fd=open(fname)
    except:
        # Try to open first in a set with file extension
        try: fd=open(fname+'.0')
        except: return False
    # Get fileinfo
    out=info_snap(fd)['valid']
    # Close file & return output
    fd.close()
    return out
def info_snap(fd):
    """
    Gets info on the format of the snapshot
    """
    # Initalize output
    out={'valid':False,'version':'None','byteswap':False,'endian':'='}
    # Read first 4 bytes & close file
    (r,) = struct.unpack('=I',fd.read(4))
    # Interpret first integer
    #     Gadget 1: 256 (or byteswapped equivalent)
    #     Gadget 2: 8   (or byteswapped equivalent)
    v={'Gadget1':np.int32(256),'Gadget2':np.int32(8)}
    for iv in v:
        if   r==v[iv]           : out.update(valid=True,version=iv)
        elif r==v[iv].byteswap(): out.update(valid=True,version=iv,byteswap=True)
    # Get endianness of file
    if out['byteswap']:
        if sys.byteorder=='little': out['endian']='>'
        else                      : out['endian']='<'
    # Reset file position & return info
    fd.seek(0,0)
    return out

####################################################################################################################################
# METHOD TO RETURN BLOCK LIST FOR A SNAPSHOT
def get_blocklist(fname,makefile=None):
    """
    Returns a list of blocks for version 1 snapshots
    """
    # Initialize blocklist
    blocklist=copy.deepcopy(LIST_BLOCKS)
    mkopt2block=dict(OUTPUTPOTENTIAL      ='POT' ,
                     OUTPUTACCELERATION   ='ACCE',
                     OUTPUTCHANGEOFENTROPY='ENDT',
                     OUTPUTTIMESTEP       ='TSTP')
    # Check for makefile
    if not isinstance(makefile,str): makefile=mmlfiles.search(fname,FEXT_MAKEFILE,nlvlup=1)
    # Read options from a snapshot
    if isinstance(makefile,str):
        makeopt=rw_makefile('R',makefile)
        for iopt in mkopt2block:
            if not makeopt.has_key(iopt): makeopt[iopt]=False
            if not makeopt[iopt]: blocklist.remove(mkopt2block[iopt])
    # Otherwise remove all optional blocks
    else:
        for iopt in LIST_BLOCKS_OPT: blocklist.remove(iopt)
    return blocklist

####################################################################################################################################
####################################################################################################################################
# METHODS TO READ/WRITE FILES

####################################################################################################################################
# METHOD TO READ/WRITE SNAPSHOTS
def rw_snap(rwid,fname,fdict=None,overwrite=False,makefile=None):
    """
    Reads/writes snapshot files
    """
    if not test_snap(fname): raise Exception('Invalid file type.')
    # Read
    if rwid=='R':
        # Open file
        fd=open(fname)
        # Get info from first 4 bytes
        info=info_snap(fd)
        # Get block list for version 1 snapshots
        if info['version']=='Gadget1': 
            blockidx=0
            blocklist=get_blocklist(fname,makefile=makefile)
        # Continue loading blocks until entire file read
        outdict={}
        while True:
            # Read name & size from block header
            if info['version']=='Gadget2':
                head=fd.read(5*4)
                if len(head)!=5*4: name,blocklen="    ",0
                else:
                    head=struct.unpack(info['endian']+'I4sIII',head)
                    if head[0] != 8 or head[3] != 8 or head[4] != head[2]-8 :
                        raise IOError, "Corrupt header record. Possibly incorrect file format"
                    name,blocklen=head[1],head[2]-8
            # Read size from block header & get name from list
            elif info['version']=='Gadget1':
                head=fd.read(4)
                if len(head)!=4: name,blocklen="    ",0
                else:
                    (blocklen,)=struct.unpack(info['endian']+'I',head)
                    try: name=blocklist[blockidx]
                    except IndexError:
                        warnings.warn("Run out of block names. Using fallbacks: UNK*",RuntimeWarning)
                        name='UNK'+str(blockidx-len(blocklist))
                blockidx+=1
            else: raise IOError,'Unsupprted snapshot version: {}'.format(info['version'])
            tag=name.strip()
            # Break if block is empty
            if blocklen==0: break
            # Header block
            if tag=='HEAD':
                if blocklen != 256: raise IOError, "Mis-sized {} block in {}".format(tag,fname)
                outdict[tag]=rw_snaphead(fd,info['endian'])
                # Update block list based on particles
                if outdict[tag]['npart'][0]==0:
                    for isph in LIST_BLOCKS_SPH: 
                        blocklist.remove(isph)
                if outdict[tag]['nmass'].sum()==0: blocklist.remove('MASS')
            # Body block
            else:
                if   tag=='MASS'           : t_part=outdict['HEAD']['nmass'].sum()
                elif tag in LIST_BLOCKS_SPH: t_part=outdict['HEAD']['npart'][0]
                else                       : t_part=outdict['HEAD']['npart'].sum()
                datatyp=False
                partlen=blocklen/t_part
                if tag=='ID':
                    if   partlen==4: datatyp=np.int32
                    elif partlen==8: datatyp=np.int64
                elif tag in ['POS','VEL','ACCE']:
                    #partlen/=3
                    if   partlen==12: datatyp=np.float32
                    elif partlen==24: datatyp=np.float64
                else:
                    if   partlen==4: datatyp=np.float32
                    elif partlen==8: datatyp=np.float64
                if isinstance(datatyp,bool): raise IOError,"No valid type for {} (partlen={})".format(tag,partlen)
                # Record position & move to end of block
                blockstart=fd.tell()
                extra_len=t_part*partlen
                if extra_len >= 2**32: blocklen=extra_len
                fd.seek(blocklen,1)
                # Add block to output dictionary
                outdict[tag]=dict(start=blockstart,
                                  length=blocklen,
                                  partlen=partlen,
                                  datatyp=datatyp)
            # Read block footer
            foot=fd.read(4)
            if len(foot)!=4: raise IOError, "Could not read block footer"
            (footlen,)=struct.unpack(info['endian']+'I',foot)
            if footlen != blocklen: raise IOError, 'Bad record size for {} in {}'.format(name,fname)
        # Close file
        fd.close()
        # Add mass block
        if 'MASS' not in outdict:
            outdict['MASS']=dict(start=0,length=0,partlen=8,datatyp=np.float64)
        # Return dictionary
        return outdict
    # Write
    elif rwid=='W':
        pass
    # Error
    else: raise Exception('Invalid rwid: {}'.format(rwid))
    return


####################################################################################################################################
# METHOD TO READ/WRITE SNAPSHOT HEADERS
def rw_snaphead(fd,endian='='):
    """
    Reads/writes snapshot header
    """
    # Read data
    if fd.mode=='r':
        data=fd.read(256)
        if len(data) != 256: raise IOError, "Could not read HEAD block in "+fname
        # Construct header structure
        h=dict(npart    = np.zeros(N_PTYP,dtype=np.uint32),
               massarr  = np.zeros(N_PTYP),
               time=0., redshift = 0.,
               flag_sfr=0,flag_feedback=0,
               npartAll = np.zeros(N_PTYP,dtype=np.uint32),
               flag_cooling=0, Nfiles=0,
               BoxSize=0., Omega0=0., OmegaLambda=0., HubbleParam=0.,
               flag_stellarage=0, flag_metals=0,
               NallHW   = np.zeros(N_PTYP,dtype=np.uint32),
               flag_entropy_instead_u=0, flag_doubleprecision=0, 
               flag_ic_info=0, lpt_scalingfactor=0., fill=0)
        # Return if empty
        if data=='': return
        # Construct format string
        fmt=endian+"IIIIIIddddddddiiIIIIIIiiddddiiIIIIIIiiif48s"
        if struct.calcsize(fmt)!=256: raise Exception('Bug: the header format string is not 256 bytes')
        # Unpack data
        (h['npart'][0],h['npart'][1],h['npart'][2],h['npart'][3],h['npart'][4],h['npart'][5],
         h['massarr'][0],h['massarr'][1],h['massarr'][2],h['massarr'][3],h['massarr'][4],h['massarr'][5],
         h['time'],h['redshift'],h['flag_sfr'],h['flag_feedback'],
         h['npartAll'][0],h['npartAll'][1],h['npartAll'][2],h['npartAll'][3],h['npartAll'][4],h['npartAll'][5],
         h['flag_cooling'],h['Nfiles'],h['BoxSize'],h['Omega0'],h['OmegaLambda'],h['HubbleParam'],h['flag_stellarage'],h['flag_metals'],
         h['NallHW'][0],h['NallHW'][1],h['NallHW'][2],h['NallHW'][3],h['NallHW'][4],h['NallHW'][5],
         h['flag_entropy_instead_u'],h['flag_doubleprecision'],h['flag_ic_info'],h['lpt_scalingfactor'],h['fill'])=struct.unpack(fmt,data)
        # Add simply calculated values
        h['nbody']=h['npart'].sum()
        h['nmass']=copy.deepcopy(h['npart'])
        h['nmass'][h['massarr']!=0]=0.
        # Return header dictionary
        return h
    # Write
    elif fd.mode=='w': pass
    # Error
    else: raise Exception('Invalid mode: {}'.format(fd.mode))
    return

####################################################################################################################################
# METHOD TO READ/WRITE CPU FILES
def rw_cpustat(rwid,fname,fdict=None,overwrite=False):
    """
    Reads/writes cpu statistics files
    Line 1:
    Line 2:
        1. Total consumption
        2. Gravitational force computation
        3. Hydrodynamics
        4. Domain decomposition
        5. Potential energy computation
        6. Drifting the particle set
        7. Timestep determination and 'kicking' of particles
        8. Writing of snapshot files
        9. Pure time for tree walks
        10. Tree construction
        11. Communication and summation
        12. Imbalance of gravitational tree
        13. SPH computations
        14. Communication and summation in the SPH routines
        15. Losses due to work-load imbalance in the SPH part
        16. Neighbour adjustment
        17. Time for PM computation
        18. Time needed to establish Peano-Hilbert order
        19. Time spent in cooling and star formation routines
    """
    ln1_keys=['step','tsim','ncpu']
    ln2_keys=['tot','grav','hydro','domain_decomp','pot','drift','kick','snap','tree_walk','tree_construct',
              'comm','tree_imbalance','sph','sph_comm','sph_imbalance','neighbour','pm_comp','order','extraphys']
    tot_keys=ln1_keys+ln2_keys
    # Read
    if    rwid=='R':
        fdict={k:[] for k in tot_keys}
        fid=open(fname,'r')
        for line in fid:
            # Head line
            if line.startswith('Step'):
                for v in line.split(', '):
                    if   v.startswith('Step '): fdict['step'].append(int(float(v.split()[-1])))
                    elif v.startswith('Time:'): fdict['tsim'].append(float(v.split()[-1]))
                    elif v.startswith('CPUs:'): fdict['ncpu'].append(int(float(v.split()[-1])))
                    else: raise Exception('Unknown variable: {}'.format(v))
            # Body line
            else: 
                for k,v in zip(ln2_keys,line.split()): fdict[k].append(float(v))
        fid.close()
        return fdict
    # Write
    elif rwid=='W':
        if not overwrite and os.path.isfile(fname): return
        fid=open(fname,'w')
        for idx in range(len(fdict['step'])):
            ifdict={k:fdict[k][idx] for k in tot_keys}
            # Head line
            line1='Step {step}, Time: {tsim}, CPUs: {ncpu}\n'.format(**ifdict)
            fid.write(line1)
            # Body line
            line2=['{'+k+': >8.2f}' for k in ln2_keys]
            fid.write('  '+line2.join('   ')+'\n')
        fid.close()
        return
    # Error
    else: raise Exception('Invalid rwid: {}'.format(rwid))
    return

####################################################################################################################################
# METHOD TO READ/WRITE PARAMETER FILES
def rw_param(rwid,fname,fdict=None,overwrite=False):
    """
    Reads/writes parameter files
    """
    # Read
    if   rwid=='R': 
        fdict=mmlio.rwdict('R',fname)
        return fdict
    # Write
    elif rwid=='W': 
        mmlio.rwdict('W',fname,fdict,width=29,overwrite=overwrite)
        return
    # Error
    else: raise Exception('Invalid rwid: {}'.format(rwid))
    return

####################################################################################################################################
# METHOD TO READ/WRITE LIST OF OUTPUT TIMES
def rw_outlist(rwid,fname,outlist=None,overwrite=None):
    """
    Reads/writes a list of output times to file
    """
    # Read
    if   rwid=='R': 
        if not os.path.isfile(fname): raise Exception('File name is not valid: {}'.format(fname))
        fid=open(fname,'r')
        outlist=[]
        for line in fid: outlist.append(float(line.strip()))
        fid.close()
    # Write
    elif rwid=='W':
        outlist=mmlpars.mml_pars(outlist,type=list)
        if os.path.isfile(fname) and not overwrite: return outlist
        fid=open(fname,'w')
        for itime in outlist: fid.write(repr(itime)+'\n')
        fid.close()
    # Return output
    return outlist

####################################################################################################################################
# METHOD TO READ SOFTENINGS FROM GADGET PARAMETER FILE
def rw_softenings(fname):
    """
    Reads and returns softening lengths from a GADGET parameter file as a dictionary with type keys
    """
    # Read parameters from parameter file
    params=rw_param('R',fname)
    # Add softenings to dictionary
    soft={ityp:params['Softening'+ityp.capitalize()] for ityp in LIST_PTYPS}
    # Return dictionary of softening lengths
    return soft

####################################################################################################################################
# METHOD TO READ GADGET MAKE FILE
def rw_makefile(rwid,fname):
    """
    Reads options from a GADGET makefile
    """
    # Read
    if rwid=='R':
        # Check if file is valid
        if not os.path.isfile(fname):
            raise Exception('Makefile name is not valid: {}'.format(fname))
        # Open file
        fid=open(fname,'r')
        # Check each line for options
        makeopt={}
        for line in fid:
            # Handle uncommented options
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
            # Handle commented options
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
            else:
                continue
        fid.close()
        return makeopt
    # Write
    elif rwid=='W': raise Exception('Why do you want to make a Makefile?')
    # Error
    else: raise Exception('Invalid rwid: {}'.format(rwid))
    return

####################################################################################################################################
# METHOD TO READ/PLOT ENERGY FILES
def rw_energy(fname):
    """
    Reads data from gadget energy files
    """
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
             'massGas','massHalo','massDisk',
             'massBulge','massStars','massBndry']
    edict={ikey:[] for ikey in colKeys}
    # Open file
    fid=open(fname,'r')
    # Read input
    for iline in fid:
        colList=iline.split() ; ncols=len(colList)
        if ncols != len(colKeys):
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
# METHOD TO DETERMINE PARTALLOCFACTOR
def mmlsim2partalloc(simstr):
    """
    Estimates the PartAllocFactor a GADGET simulation will need using the minimum amount of memory
    """
    # Set fiducial values
    ntot0=2.0e8
    partalloc0=1.5
    # Scale directly (PartAlloc increases and npart/processor decreases)
    partalloc=partalloc0*(ntot0/float(simstr.ntot))
    # Round to two sig figs
    partalloc=mmlmath.oom(partalloc,nsig=2,method='CEIL')
    # Return
    return partalloc

####################################################################################################################################
# METHOD TO DETERMINE TREEALLOCFACTOR
def mmlsim2treealloc(simstr):
    """
    Estimates the TreeAllocFactor a GADGET simulation will need
    """
    treealloc=1.0
    return treealloc

####################################################################################################################################
# METHOD TO DETERMINE COMMUNICATION BUFFER
def mmlsim2combuff(simstr):
    """
    Estimates the communication buffer a GADGET simulation will need
    """
    minbuff=15
    # Get file dictionary
    icfile=simstr.fdict['gadget']['icgen']['stat_ic']
    # Get buffer size from IC file if it exists
    if os.path.isfile(icfile):
        statinfo=os.stat(icfile)
        combuff=float(statinfo.st_size)/(1.0e6)
    # Otherwise get it by scaling previous results
    else:
        combuff=float(simstr.ntot)*(24.0/560000.0)
    # Round and set minimum
    combuff=int(mmlmath.oom(combuff,nsig=2,method='CEIL'))
    combuff=max(combuff,minbuff)
    # Return rounded output
    return combuff
    
####################################################################################################################################
# METHOD TO DETERMINE MEMORY FROM MMLSIM OBJECT
def mmlsim2mcpu(simstr):
    """
    Estimates the amount of memory a GADGET simulation will need
    """
    # Get file dictionary
    parfile=simstr.fdict['gadget']['icgen']['stat_pm']
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
    mcpu=calcmcpu(simstr.ntot,nproc=simstr['nprocev'],
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
# PARAMETER LISTS
def listpar(partyp,**exkw):
    partyp=mmlpars.mml_pars(partyp,type=str)
    if   partyp=='mklog_buildgal': par={'keylist':['nproc','softs','time0'],
                                        'nproc'    :{'form':int  ,'def':0       },
                                        'softs'    :{'form':dict ,'def':{ityp:0. for ityp in LIST_PTYPS}},
                                        'time0'    :{'form':float,'def':0.}}
    elif partyp=='mklog_snapshot': par={'keylist':['snapext','snaptime','resettime','time0','heatfact','softs'],
                                        'snapext'  :{'form':str  ,'def':'_000ic'},
                                        'resettime':{'form':bool ,'def':True    },
                                        'heatfact' :{'form':float,'def':1.0     },
                                        'snaptime' :{'form':float,'def':0.},
                                        'softs'    :{'form':dict ,'def':{ityp:0. for ityp in LIST_PTYPS}},
                                        'time0'    :{'form':float,'def':0.}}
    elif partyp=='mklog_restart' : par={'keylist':['nproc'],
                                        'nproc'    :{'form':int  ,'def':0       }}
    elif partyp=='mklog_interact': par={'keylist':['snapext1','snapext2','snaptime1','snaptime2','time0','softs'],
                                        'snapext1' :{'form':str  ,'def':'_000ic'},
                                        'snapext2' :{'form':str  ,'def':'_000ic'},
                                        'snaptime1':{'form':float,'def':0.},
                                        'snaptime2':{'form':float,'def':0.},
                                        'softs'    :{'form':dict ,'def':{ityp:0. for ityp in LIST_PTYPS}},
                                        'time0'    :{'form':float,'def':0.}}
    elif partyp=='mklog_timetest': par={'keylist':['snapext','parmrun','softs'],
                                        'snapext'  :{'form':str  ,'def':'ic'},
                                        'softs'    :{'form':dict ,'def':{ityp:0. for ityp in LIST_PTYPS}},
                                        'parmrun'  :{'form':str  ,'def':''}}
    elif partyp=='param': par={'OutputDir'                :{'def':'output'          ,'form':str    },
                               'SnapshotFileBase'         :{'def':'snapshot'        ,'form':str    },
                               'SnapFormat'               :{'def':1                 ,'form':[1,2,3]},
                               'NumFilesPerSnapshot'      :{'def':1                 ,'form':int    },
                               'InitCondFile'             :{'def':'galaxy.dat'      ,'form':str    },
                               'ICFormat'                 :{'def':1                 ,'form':[1,2,3]},
                               'RestartFile'              :{'def':'restart'         ,'form':str    },
                               'InfoFile'                 :{'def':'info.txt'        ,'form':str    },
                               'TimingsFile'              :{'def':'timings.txt'     ,'form':str    },
                               'CpuFile'                  :{'def':'cpu.txt'         ,'form':str    },
                               'EnergyFile'               :{'def':'energy.txt'      ,'form':str    },
                               # CPU-time limit and restart options
                               'TimeLimitCPU'             :{'def':172800            ,'form':int    },
                               'ResubmitCommand'          :{'def':resub.sh          ,'form':str    },
                               'ResubmitOn'               :{'def':0                 ,'form':[0,1]  },
                               'CpuTimeBetRestartFile'    :{'def':10800             ,'form':int    },
                               # Simulation specific parameters
                               'TimeBegin'                :{'def':0.0               ,'form':float  },
                               'TimeMax'                  :{'def':10.2269           ,'form':float  },
                               'BoxSize'                  :{'def':1000.0            ,'form':float  },
                               'PeriodicBoundariesOn'     :{'def':0                 ,'form':[0,1]  },
                               'ComovingIntegrationOn'    :{'def':0                 ,'form':[0,1]  },
                               # Cosmological parameters
                               'HubbleParam'              :{'def':0.7               ,'form':float  },
                               'Omega0'                   :{'def':0.3               ,'form':float  },
                               'OmegaLambda'              :{'def':0.7               ,'form':float  },
                               'OmegaBaryon'              :{'def':0.04              ,'form':float  },
                               # Memory allocation
                               'BufferSize'               :{'def':29                ,'form':int    },
                               'PartAllocFactor'          :{'def':3.0               ,'form':float  },
                               'TreeAllocFactor'          :{'def':2.0               ,'form':float  },
                               # Gravitational force accuracy
                               'TypeOfOpeningCriterion'   :{'def':0                 ,'form':[0,1]  },
                               'ErrTolTheta'              :{'def':0.5               ,'form':float  },
                               'ErrTolForceAcc'           :{'def':0.002             ,'form':float  },
                               # Time integration accuracy
                               'MaxSizeTimestep'          :{'def':0.0102269         ,'form':float  },
                               'MinSizeTimestep'          :{'def':0.0               ,'form':float  },
                               'TypeOfTimestepCriterion'  :{'def':0                 ,'form':[0]    },
                               'ErrTolIntAccuracy'        :{'def':0.00125           ,'form':float  },
                               'TreeDomainUpdateFrequency':{'def':0.1               ,'form':float  },
                               'MaxRMSDisplacementFac'    :{'def':0.2               ,'form':float  },
                               # Output of snapshot files
                               'OutputListOn'             :{'def':0                 ,'form':[0,1]  },
                               'OutputListFilename'       :{'def':'output_times.txt','form':str    },
                               'TimeOfFirstSnapshot'      :{'def':0.0               ,'form':float  },
                               'TimeBetSnapshot'          :{'def':0.0102269         ,'form':float  },
                               'TimeBetStatistics'        :{'def':0.00102269        ,'form':float  },
                               'NumFilesWrittenInParallel':{'def':1                 ,'form':int    },
                               # System of units
                               'UnitVelocity_in_cm_per_s' :{'def':1.0e5             ,'form':float  },
                               'UnitLength_in_cm'         :{'def':3.085678e21       ,'form':float  },
                               'UnitMass_in_g'            :{'def':1.989e43          ,'form':float  },
                               'GravityConstantInternal'  :{'def':0.                ,'form':float  },
                               # SPH parameters
                               'DesNumNgb'                :{'def':64                ,'form':int    },
                               'MaxNumNgbDeviation'       :{'def':2                 ,'form':int    },
                               'ArtBulkViscCons'          :{'def':1.0               ,'form':float  },
                               'CourantFac'               :{'def':0.15              ,'form':float  },
                               'InitGasTemp'              :{'def':10000.            ,'form':float  },
                               'MinGasTemp'               :{'def':20.0              ,'form':float  },
                               'MinGasHsmlFractional'     :{'def':0.1               ,'form':float  },
                               # Softenings
                               'SofteningGas'             :{'def':0.0               ,'form':float  },
                               'SofteningHalo'            :{'def':0.1               ,'form':float  },
                               'SofteningDisk'            :{'def':0.05              ,'form':float  },
                               'SofteningBulge'           :{'def':0.05              ,'form':float  },
                               'SofteningStars'           :{'def':0.0               ,'form':float  },
                               'SofteningBndry'           :{'def':0.0               ,'form':float  },
                               'SofteningGasMaxPhys'      :{'def':0.0               ,'form':float  },
                               'SofteningHaloMaxPhys'     :{'def':0.1               ,'form':float  },
                               'SofteningDiskMaxPhys'     :{'def':0.05              ,'form':float  },
                               'SofteningBulgeMaxPhys'    :{'def':0.05              ,'form':float  },
                               'SofteningStarsMaxPhys'    :{'def':0.0               ,'form':float  },
                               'SofteningBndryMaxPhys'    :{'def':0.0               ,'form':float  }}
    else: raise Exception('Invalid parameter type: {}'.format(partyp))
    return par
