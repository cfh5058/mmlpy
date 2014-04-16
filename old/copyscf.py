####################################################################################################################################
#
# MEAGAN LANG'S SCF METHODS
#
####################################################################################################################################
import mmlclass,mmlpars,mmlfiles,os,shutil,mmlio,glob,mmlstring
import mmlsimplot as msplot
import mmlmath as mmath
import numpy as np
import pNbody

####################################################################################################################################
# CLASS CONTAINING GADGET FILE OPTIONS
####################################################################################################################################
class fileopt(mmlclass.mydict):
    '''
    NAME:
        mmlscf.fileopt
    PURPOSE:
        Dictionary for passing around SCF file options.
    FORMAT:
        A dictionary subclass with keys:
            execpath: Str absolute path to location of compiled SCF executable
            parmpath: Str absolute path to location of default SCF par file
            modspath: Str absolute path to location of default SCF mod file
        See fileopt.get_form() dictionary for additional details.
    '''

    def __init__(self,execpath=None,parmpath=None):
        '''
        NAME:
            mmlscf.fileopt.__init__
        PURPOSE:
            To initialize fileopt objects.
        CALLING:
            obj=fileopt(keys)
        KEYWORDS:
            Keywords include all keys listed for the fileopt object.
        OUTPUT:
            New instance of a mmlscf.fileopt object.
        '''
        mmlclass.mydict.__init__(self,execpath=execpath,parmpath=parmpath)
        self.pars()

    def get_form(self):
        '''
        NAME:
            mmlscf.fileopt.get_form
        PURPOSE:
            To return a dictionary defining fileopt key/value pairs.
        CALLING:
            form=fileopt.get_form()
        OUTPUT:
            Dictionary defining fileopt key/value pairs.
        '''
        form={
            'execpath': mmlpars.parsdict(default='/scratch/langmm/scf/mpiscf',type=str),
            'parmpath': mmlpars.parsdict(default='/scratch/langmm/scf/scfpar',type=str),
            'modspath': mmlpars.parsdict(default='/scratch/langmm/scf/scfmod',type=str)
            }
        return form

    def pars(self):
        '''
        NAME:
            mmlscf.fileopt.pars
        PURPOSE:
            To parse fileopt objects.
        CALLING:
            fileopt.pars()
        '''
        self=mmlpars.mml_formpars(self,self.get_form())
        
####################################################################################################################################
# CLASS CONTAINING SCF PARAMETER OPTIONS
####################################################################################################################################
class paropt(mmlclass.mydict):
    '''
    NAME:
        mmlscf.paropt
    PURPOSE:
        Dictionary for passing around SCF parameter options.
    FORMAT:
        A dictionary subclass with keys:
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
        See paropt.get_form() dictionary for additional details.
    '''

    def __init__(self,headline=None,nsteps=None,noutbod=None,noutlog=None,dteps=None,G=None,
                 tfinal=None,multistep=None,fixedn=None,selfgrav=None,inptcoef=None,outpcoef=None,
                 zeroodd=None,zeroeven=None,fixacc=None,rcrit=None,ecrit=None,lilout=None,nlilout=None):
        '''
        NAME:
            mmlscf.paropt.__init__
        PURPOSE:
            To initialize mmlscf.paropt objects.
        CALLING:
            obj=paropt(keys)
        KEYWORDS:
            Keywords include all keys listed for the paropt object class.
        OUTPUT:
            New instance of a mmlscf.paropt object.
        '''
        mmlclass.mydict.__init__(self,headline=headline,nsteps=nsteps,noutbod=noutbod,noutlog=noutlog,
                                 G=G,tfinal=tfinal,multistep=multistep,fixedn=fixedn,selfgrav=selfgrav,
                                 inptcoef=inptcoef,outpcoef=outpcoef,zeroodd=zeroodd,zeroeven=zeroeven,
                                 fixacc=fixacc,rcrit=rcrit,ecrit=ecrit,lilout=lilout,nlilout=nlilout)
        self.pars()

    def get_form(self):
        '''
        NAME:
            mmlscf.paropt.get_form
        PURPOSE:
            To return a dictionary that defines paropt key/value pairs.
        CALLING:
            form=paropt.get_form()
        OUTPUT:
            Dictionary defining paropt key/value pairs. Used by paropt.pars to parse paropt objects.
        '''
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
        return form

    def pars(self):
        '''
        NAME:
            mmlscf.paropt.pars
        PURPOSE:
            To parse paropt objects.
        CALLING:
            paropt.pars()
        '''
        self=mmlpars.mml_formpars(self,self.get_form())

####################################################################################################################################
# CLASS CONTAINING SCF MODIFICATION OPTIONS
####################################################################################################################################
class modopt(mmlclass.mydict):
    '''
    NAME:
        mmlscf.modopt
    PURPOSE:
        Dictionary for passing around SCF modification options.
    FORMAT:
        A dictionary subclass with keys:
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
        See paropt.get_form() dictionary for additional details.
    '''

    def __init__(self,iseed=None,bhmass=None,epsbh=None,tstartbh=None,tgrowbh=None,tlivebh=None,
                 tdiebh=None,xdrag=None,ydrag=None,zdrag=None,tstartdrag=None,tgrowdrag=None,tlivedrag=None,
                 tdiedrag=None,bhgrav=None,usedrag=None,stellev=None):
        '''
        NAME:
            mmlscf.modopt.__init__
        PURPOSE:
            To initialize mmlscf.modopt objects.
        CALLING:
            obj=modopt(keys)
        KEYWORDS:
            Keywords include all keys listed for the modopt object class.
        OUTPUT:
            New instance of a mmlscf.modopt object.
        '''
        mmlclass.mydict.__init__(self,iseed=iseed,bhmass=bhmass,epsbh=epsbh,tstartbh=tstartbh,
                                 tgrowbh=tgrowbh,tlivebh=tlivebh,tdiebh=tdiebh,xdrag=xdrag,ydrag=ydrag,zdrag=zdrag,
                                 tstartdrag=tstartdrag,tgrowdrag=tgrowdrag,tlivedrag=tlivedrag,
                                 tdiedrag=tdiedrag,bhgrav=bhgrav,usedrag=usedrag,stellev=stellev)
        self.pars()

    def get_form(self):
        '''
        NAME:
            mmlscf.modopt.get_form
        PURPOSE:
            To return a dictionary that defines paropt key/value pairs.
        CALLING:
            form=modopt.get_form()
        OUTPUT:
            Dictionary defining modopt key/value pairs. Used by modopt.pars to parse modopt objects.
        '''
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
        return form

    def pars(self):
        '''
        NAME:
            mmlscf.modopt.pars
        PURPOSE:
            To parse modopt objects.
        CALLING:
            modopt.pars()
        '''
        self=mmlpars.mml_formpars(self,self.get_form())

####################################################################################################################################
# METHOD TO CREATE DICTIONARY OF GADGET FILES
####################################################################################################################################
def files(rundir,runid,addid='',options=fileopt(),accre=False):
    '''
    NAME:
        mmlscf.files
    PURPOSE:
        To return files associated with an SCF N-body simulation analysis.
    CALLING:
        fdict=files(rundir,runid,addid="",options=fileopt,accre=False)
    ARGUMENTS:
        rundir:   Directory containing all files for SCF analysis
        runid:    String uniquely identifying a N-body simulation/analysis
    KEYWORDS:
        addid:    String to add to the end of runid to create prefix for file names
        options:  mmlscf.fileopt object containing file options
        accre:    Bool specifying if files are on ACCRE
    OUTPUT:
        fdict:    Dictionary containing the following key/value pairs:
          dir:      Directory containing all SCF analysis files (input/output)
          exstat:   Location of static compiled SCF executable
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
    '''
    # Pars input
    if options == None: options=fileopt()
    else: options.pars()
    # Initialize output dictionary
    pfix=runid+addid+'.'
    files={'dir':os.path.join(rundir,'scf'+addid),'exstat':options['execpath'],'pmstat':options['parmpath'],'mdstat':options['modspath']}
    # Add input files
    infiles={}
    infiles['dir']=os.path.join(files['dir'],'input')
    fkeys_in=['exec','par','mod','wrap','pbs']
    flist_in=['mpiscf','scfpar','scfmod','scfwrap','scfpbs']
    infiles=mmlfiles.filedict(infiles['dir'],filekeys=fkeys_in,filelist=flist_in,pfix=pfix,fdict=infiles)
    # Long input directories
    dkeys_in_long=['snapdir','coefdir']
    dlist_in_long=['snap','coef']
    infiles=mmlfiles.filedict(infiles['dir'],dirkeys=dkeys_in_long,dirlist=dlist_in_long,pfix=pfix,fdict=infiles)
    # Long input files
    fkeys_in_long=['snapbase','coefbase']
    flist_in_long=['scfbi','scficoef']
    for ifin in range(len(dkeys_in_long)):
        infiles=mmlfiles.filedict(infiles[dkeys_in_long[ifin]],
                                  filekeys=[fkeys_in_long[ifin]],filelist=[flist_in_long[ifin]],pfix=pfix,fdict=infiles)
    files['input']=infiles
    # Add output files
    outfiles={}
    outfiles['dir']=os.path.join(files['dir'],'output')
    fkeys_out=['runout']
    flist_out=['output']
    outfiles=mmlfiles.filedict(outfiles['dir'],filekeys=fkeys_out,filelist=flist_out,pfix=pfix,fdict=outfiles)
    # Long output directories
    dkeys_out_long=['logdir','outdir','snapdir','coefdir','eldir','chkdir','olildir']
    dlist_out_long=['log','out','snap','coef','el','chk','olil']
    outfiles=mmlfiles.filedict(outfiles['dir'],dirkeys=dkeys_out_long,dirlist=dlist_out_long,pfix=pfix,fdict=outfiles)
    # Long output files
    fkeys_out_long=['logbase','outbase','snapbase','coefbase','elbase','chkbase','olilbase']
    flist_out_long=['scflog','scfout','snap','scfocoef','scfel','scfchkpt','slil']
    for ifout in range(len(dkeys_out_long)):
        outfiles=mmlfiles.filedict(outfiles[dkeys_out_long[ifout]],
                                   filekeys=[fkeys_out_long[ifout]],filelist=[flist_out_long[ifout]],pfix=pfix,fdict=outfiles)
    files['output']=outfiles
    # Return generic names if on ACCRE (prevents need to modify and recompile tmhscf.h for every run)
#    if accre:
#        for ikey in files['input'].keys():
#            if ikey != 'dir':
#                ifullfile=files['input'][ikey]
#                ishrtfile=os.path.join(os.path.dirname(ifullfile),os.path.splitext(ifullfile)[1][1:])
#                files['input'][ikey]=ishrtfile
#        for ikey in files['output'].keys():
#            if ikey != 'dir':
#                ifullfile=files['output'][ikey]
#                ishrtfile=os.path.join(os.path.dirname(ifullfile),os.path.splitext(ifullfile)[1][1:])
#                files['output'][ikey]=ishrtfile
    # Return file dictionary
    return files

####################################################################################################################################
# METHOD TO PERFORM PREP WORK FOR RUNNING SCF CODE ON ACCRE
####################################################################################################################################
def run(simstr):
    '''
    NAME:
        mmlscf.run
    PURPOSE:
        To preform prep work for running SCF analysis code on ACCRE.
    CALLING:
        run(simstr)
    ARGUMENTS:
        simstr:    icgen.mmlsim object describing a simulation
    '''
    if mmlio.yorn('Move SCF input files to ACCRE?'):      mvfiles(simstr,'TO','INPUT')
    if mmlio.yorn('Move old SCF output files to ACCRE?'): mvfiles(simstr,'TO','OUTPUT')
    filedict=simstr.mkfiledict(accre=True)
    mmlfiles.mkdirs(filedict['scf']['output']['dir'])
    dkeylist=['snapdir','coefdir','outdir','logdir','chkdir','eldir','olildir']
    dirlist=[filedict['scf']['output'][idirkey] for idirkey in dkeylist]
    mmlfiles.mkdirs(dirlist)

####################################################################################################################################
# METHOD TO MOVE SCF FILES TO/FROM ACCRE
####################################################################################################################################
def mvfiles(simstr,method,filetype,overwrite=False,**options):
    '''
    NAME:
        mmlscf.mvfiles
    PURPOSE:
        To move SCF analysis files to/from ACCRE.
    CALLING:
        mvfiles(simstr,method,filetype,overwrite=False)
    ARGUMENTS:
        simstr:    icgen.mmlsim object describing a simulation
        method:    string specifying whether to move files "to"/"from" ACCRE
        filetype:  type of file to move to/from ACCRE ("input" or "output")
    KEYWORDS:
        overwrite: Boolean specifying whether or not to overwite existing files.
    '''
    # Get file dictionary
    files_local=simstr.mkfiledict(accre=False,**options)
    files_accre=simstr.mkfiledict(accre=True ,**options)
    # Select source/destination files
    if method.upper()=='TO':
        files_beg=files_local['scf']
        files_end=files_accre['scf']
    elif method.upper()=='FROM':
        files_beg=files_accre['scf']
        files_end=files_local['scf']
    else:
        raise Exception('Keyword method ({}) must be "TO" or "FROM".'.format(method))
    # Select specific files
    filekeydict={'executable':['exec'],
                 'parameter':['par','mod','pbs','wrap'],
                 'shortout':['runout'],'longout':['outbase','logbase','chkbase','elbase','olilbase'],
                 'snapshot':['snapbase'],'coefficient':['coefbase']}
    if filetype.upper()=='INPUT':
        mmlfiles.mkdirs([files_end['input']['dir'],files_end['input']['snapdir'],files_end['input']['coefdir']])
        subtypelist=['snapshot','executable','parameter']
    elif filetype.upper()=='OUTPUT':
        mmlfiles.mkdirs([files_end['output']['dir'],files_end['output']['snapdir'],files_end['output']['coefdir'],
                         files_end['output']['outdir'],files_end['output']['logdir'],files_end['output']['chkdir'],
                         files_end['output']['eldir'],files_end['output']['olildir']])
        subtypelist=['genoutput','coefficient','longout']
    else:
        raise Exception('Keyword filetype ({}) must be "INPUT" or "OUTPUT",'.format(filetype))
    # Loop over list moving files
    for itype in subtypelist:
        if mmlio.yorn('Move '+itype+' file(s)?'):
            if itype in ['snapshot','coefficient','longout']:
                iowflag=mmlio.yorn('Overwrite existing '+itype+' '+filetype.lower()+' file(s)?')
                for itype_file in filekeydict[itype]:
                    isrclist=glob.glob(files_beg[filetype.lower()][itype_file]+'*')
                    idstlist=[]
                    for isrc in isrclist:
                        iext=isrc.split(files_beg[filetype.lower()][itype_file])[1]
                        idstlist.append(files_end[filetype.lower()][itype_file]+iext)
            else:
                iowflag=True
                isrclist=[] ; idstlist=[]
                for ifilekey in filekeydict[itype]:
                    isrclist.append(files_beg[filetype.lower()][ifilekey])
                    idstlist.append(files_end[filetype.lower()][ifilekey])
            print itype
            for isrc,idst in zip(isrclist,idstlist):
                print '{} ---> {}'.format(isrc,idst)
            mmlfiles.cpfiles(isrclist,idstlist,overwrite=iowflag)
        else: continue


####################################################################################################################################
# METHOD TO CREATE LOCAL COPIES OF SCF FILES
####################################################################################################################################
def mkfiles(simstr,fileopt=None,paropt=paropt(),ext='*'):
    '''
    NAME:
        mmlscf.mkfiles
    PURPOSE:
        To create files required to run a SCF simulation analysis.
    CALLING:
        mkfiles(simstr,fileopt=None,paropt=mmlscf.paropt)
    ARGUMENTS:
        simstr:  icgen.mmlsim object containing info on an Nbody simulation.
    KEYWORDS:
        fileopt: Dictionary of options for simstr.mkfiledict
        paropt:  mmlgadget.paropt object containing options for mmlgadget.mkparam
    '''
    mmlio.verbose('Creating files for a SCF run locally...',border=True)
    # Get file names
    files=simstr.mkfiledict(accre=False,options=fileopt)
    mmlfiles.mkdirs(files['scf']['input']['dir'])
    # Snapshot files
    if mmlio.yorn('Create local copies of SCF snapshots?'):
        owsnap=mmlio.yorn('Overwrite existing SCF snapshots?')
        gadgsnap=files['gadget']['output']['snapbase']+ext
        srcsnapList=glob.glob(gadgsnap)
        if len(srcsnapList)==0:
            raise Exception('No GADGET snapshots were found matching the pattern {}'.format(gadgsnap))
        mmlfiles.mkdirs(os.path.dirname(files['scf']['input']['snapbase']))
        for isrc in srcsnapList:
            iext=isrc.split(files['gadget']['output']['snapbase'])[-1]
            idst=files['scf']['input']['snapbase']+iext
            if not owsnap and os.path.isfile(idst):
                pass
            else:
                inb=pNbody.Nbody(isrc,ftype='gadget',unitsfile=files['gadget']['input']['param'],makefile=files['gadget']['input']['makefile'])
                inb=inb.set_ftype('scf')
                inb.rename(idst)
                inb.write()
                del inb
    # Executable file
    if mmlio.yorn('Create new local copy of SCF executable?'):
        if not os.path.isfile(files['scf']['exstat']):
            raise Exception('Static SCF executable does not exist. Create it first: {}'.format(files['scf']['exstat']))
        mmlfiles.cpfiles(files['scf']['exstat'],files['scf']['input']['exec'],overwrite=True)
    # Parameter files
    if mmlio.yorn('Create new SCF submission files?'):
        pardict=mkpar(simstr,fileopt=fileopt,overwrite=True)
        moddict=mkmod(simstr,fileopt=fileopt,overwrite=True)
        pbsdict=mkpbs(simstr,fileopt=fileopt,overwrite=True)
        wrapdict=mkwrap(simstr,fileopt=fileopt,overwrite=True,mpipath=pbsdict['mpipath'])

####################################################################################################################################
# METHOD TO CREATE SCF PARAMETER FILE
####################################################################################################################################
def mkpar(simstr,fileopt=None,overwrite=False):
    '''
    NAME:
        mmlscf.mkpar
    PURPOSE:
        To generate a SCF parameter file.
    CALLING:
        pardict=mkpar(simstr[,fileopt=,overwrite=])
    ARGUMENTS:
        simstr:    icgen.mmlsim obj containing basic info on a Nbody simulation
    KEYWORDS:
        fileopt:   Dictionary containing file naming options.
        overwrite: Bool determining if existing parameter file is overwritten.
    OUTPUT:
        pardict:   Dictionary containing SCF parameters
    '''
    # Set constants
    keywidth=29
    headline='**********Basic input parameters**********'
    tailline='******************************************'
    # Pars input
    overwrite=mmlpars.mml_pars(overwrite,type=bool)
    # Get file names
    files_local=simstr.mkfiledict(accre=False,options=fileopt)
    files_accre=simstr.mkfiledict(accre=True ,options=fileopt)
    parfile=files_local['scf']['input']['par']
    parfile_stat=files_local['scf']['pmstat']
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
            pardict=paropt()
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
####################################################################################################################################
def mkmod(simstr,fileopt=None,overwrite=False):
    '''
    NAME:
        mmlscf.mkmod
    PURPOSE:
        To generate a SCF modification file.
    CALLING:
        moddict=mkmod(simstr[,fileopt=,overwrite=])
    ARGUMENTS:
        simstr:    icgen.mmlsim obj containing basic info on a Nbody simulation
    KEYWORDS:
        fileopt:   Dictionary containing file naming options.
        overwrite: Bool determining if existing modification file is overwritten.
    OUTPUT:
        moddict:   Dictionary containing SCF modifications
    '''
    # Set constants
    keywidth=29
    headline='**********Basic input parameters**********'
    tailline='******************************************'
    # Pars input
    overwrite=mmlpars.mml_pars(overwrite,type=bool)
    # Get file names
    files_local=simstr.mkfiledict(accre=False,options=fileopt)
    files_accre=simstr.mkfiledict(accre=True ,options=fileopt)
    modfile=files_local['scf']['input']['mod']
    modfile_stat=files_local['scf']['mdstat']
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
            moddict=modopt()
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
####################################################################################################################################
def mkpbs(simstr,fileopt=None,overwrite=False):
    '''
    NAME:
        mmlscf.mkpbs
    PURPOSE:
        To create a .pbs script for submitting a SCF analysis to ACCRE.
    CALLING:
        pbsdict=mkpbs(simstr[,fileopt=None,overwrite=False])
    ARGUMENTS:
        simstr:    icgen.mmlsim obj containing basic info on a Nbody simulation
    KEYWORDS:
        fileopt:   Dictionary containing file naming options
        overwrite: Bool determining if existing parameter file is overwritten
    OUTPUT:
        pbsdict:   Dictionary containing parameters necessary to write a pbs script
    '''
    # Pars input
    overwrite=mmlpars.mml_pars(overwrite,type=bool)
    # Get file names
    files_local=simstr.mkfiledict(accre=False,options=fileopt)
    files_accre=simstr.mkfiledict(accre=True ,options=fileopt)
    pbsfile=files_local['scf']['input']['pbs']
    parfile=files_local['scf']['input']['par']
    pbsdict=mmlio.pbsdict()
    # Load existing file if overwrite False
    if not overwrite and os.path.isfile(pbsfile):
        mmlio.verbose('Loading existing SCF pbs script:',addval=pbsfile)
        pbsdict.rwmpi('R',pbsfile)
    # Otherwise create file from scratch
    else:
        mmlio.verbose('Creating new SCF pbs script:',addval=pbsfile)
        # Update memory
        simstr.get_memory('scf',update=True)
        # Variables from simstr
        pbsdict['ppn'      ]=simstr['ppn']
        pbsdict['nodes'    ]=simstr['nproc']/simstr['ppn']
        pbsdict['pmem'     ]=simstr['mcpu']
        pbsdict['mem'      ]=simstr['mcpu']*simstr['nproc']
        pbsdict['walltime' ]=simstr['tcpu']
        # Varaibles from path
        pbsdict['simout'   ]=files_accre['scf']['output']['runout']
#        pbsdict['rundir'   ]=os.path.dirname(files_accre['scf']['input']['exec'])
#        pbsdict['execpath' ]=os.path.basename(files_accre['scf']['input']['exec'])
#        pbsdict['execflag' ]=0
#        pbsdict['execinput']=' '
#        pbsdict.rwmpi('W',pbsfile,overwrite=overwrite)
        pbsdict['command'  ]=['./{}'.format(os.path.basename(files_accre['scf']['input']['wrap']))]
        pbsdict.rwpbs('W',pbsfile,overwrite=overwrite)
    # Return output
    return pbsdict

####################################################################################################################################
# METHOD TO CREATE MPISCF WRAPPER
####################################################################################################################################
def mkwrap(simstr,fileopt=None,overwrite=None,mpipath=None):
    '''
    NAME:
        mmlscf.mkwrap
    PURPOSE:
        To create a bash script to handle loop over mpiscf cases
    CALLING:
        wrapdict=mkwrap(simstr[,fileopt=None,overwrite=False])
    ARGUMENTS:
        simstr:    icgen.mmlsim obj containing basic info on a Nbody simulation
    KEYWORDS:
        fileopt:   Dictionary containing file naming options
        overwrite: Bool determining if existing parameter file is overwritten
    OUTPUT:
        pbsdict:   Dictionary containing parameters necessary to write a pbs script
    '''
    # Set constants
    mpipathDEF=mmlio.pbsdict()['mpipath']
    # Pars input
    overwrite=mmlpars.mml_pars(overwrite,default=False,type=bool)
    mpipath=mmlpars.mml_pars(mpipath,default=mpipathDEF,type=str)
    # Get file names
    files_local=simstr.mkfiledict(accre=False,options=fileopt)
    files_accre=simstr.mkfiledict(accre=True ,options=fileopt)
    wrapfile=files_local['scf']['input']['wrap']
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
            "exec      = '{}'".format(files_accre['scf']['input']['exec']).replace(' ',''),
            "par       = '{}'".format(files_accre['scf']['input']['par']).replace(' ',''),
            "mod       = '{}'".format(files_accre['scf']['input']['mod']).replace(' ',''),
            "isnapbase = '{}'".format(files_accre['scf']['input']['snapbase']).replace(' ',''),
            "icoefbase = '{}'".format(files_accre['scf']['input']['coefbase']).replace(' ',''),
            "execfile  = '{}'".format(os.path.splitext(files_accre['scf']['input']['exec'])[1][1:]).replace(' ',''),
            "parfile   = '{}'".format(os.path.splitext(files_accre['scf']['input']['par'])[1][1:]).replace(' ',''),
            "modfile   = '{}'".format(os.path.splitext(files_accre['scf']['input']['mod'])[1][1:]).replace(' ',''),
            "isnapfile = '{}'".format(os.path.splitext(files_accre['scf']['input']['snapbase'])[1][1:]).replace(' ',''),
            "icoeffile = '{}'".format(os.path.splitext(files_accre['scf']['input']['coefbase'])[1][1:]).replace(' ',''),
            "# Output files",
            "osnapbase = '{}'".format(files_accre['scf']['output']['snapbase']).replace(' ',''),
            "ocoefbase = '{}'".format(files_accre['scf']['output']['coefbase']).replace(' ',''),
            "outbase   = '{}'".format(files_accre['scf']['output']['outbase']).replace(' ',''),
            "logbase   = '{}'".format(files_accre['scf']['output']['logbase']).replace(' ',''),
            "chkbase   = '{}'".format(files_accre['scf']['output']['chkbase']).replace(' ',''),
            "elbase    = '{}'".format(files_accre['scf']['output']['elbase']).replace(' ',''),
            "olilbase  = '{}'".format(files_accre['scf']['output']['olilbase']).replace(' ',''),
            "osnapfile = '{}'".format(os.path.splitext(files_accre['scf']['output']['snapbase'])[1][1:]).replace(' ',''),
            "ocoeffile = '{}'".format(os.path.splitext(files_accre['scf']['output']['coefbase'])[1][1:]).replace(' ',''),
            "outfile   = '{}'".format(os.path.splitext(files_accre['scf']['output']['outbase'])[1][1:]).replace(' ',''),
            "logfile   = '{}'".format(os.path.splitext(files_accre['scf']['output']['logbase'])[1][1:]).replace(' ',''),
            "chkfile   = '{}'".format(os.path.splitext(files_accre['scf']['output']['chkbase'])[1][1:]).replace(' ',''),
            "elfile    = '{}'".format(os.path.splitext(files_accre['scf']['output']['elbase'])[1][1:]).replace(' ',''),
            "olilfile  = '{}'".format(os.path.splitext(files_accre['scf']['output']['olilbase'])[1][1:]).replace(' ',''),
            "lensnapbase=${#isnapbase}","tempdir='temprun'"]
        # Command to initialize files
        cmdlist+=[' ',
                  '# Move files to temporary directory w/ built in SCF names',
                  "if [ ! -e $tempdir ]; then","    mkdir $tempdir","fi","cd $tempdir",
                  "cp $exec $execfile",
                  "cp $par $parfile",
                  "cp $mod $modfile",
                  "snaplist=$(ls {}*)".format(files_accre['scf']['input']['snapbase'])]
        # Commands to count input files
#        cmdlist+=[' ',
#                  '# Count input files',
#                  "snaplist=$(ls {}*)".format(files_accre['scf']['input']['snapbase']),
#                  "coeflist=$(ls {}*)".format(files_accre['scf']['input']['coefbase']),
#                  ]
        # Commands to perform loop
        cmdlist+=[' ',
                  '# Loop over SCF input snapshots',
                  'for isnap in $isnapbase*; do',
                  '    # Get file extension',
                  '    iext=${isnap:${lensnapbase}}',
                  '    echo "Starting snapshot: ${isnap}"',
                  '    echo "ext=${iext}"',
                  ' ',
                  '    # Move necessary input files',
                  '    if [ -e "${isnapbase}${iext}" ]; then','        mv "${isnapbase}${iext}" $isnapfile','    fi',
                  '    if [ -e "${icoefbase}${iext}" ]; then','        mv "${icoefbase}${iext}" $icoeffile','    fi',
                  ' ',
                  '    # Run mpiscf script',
                  '    {} -np {} ./$execfile'.format(mpipath,simstr['nproc']),
                  ' ',
                  '    # Recover files',
                  '    if [ -e $isnapfile  ]; then','        mv $isnapfile  "${isnapbase}${iext}"','    fi',
                  '    if [ -e $icoeffile  ]; then','        mv $icoeffile  "${icoefbase}${iext}"','    fi',
                  '    if [ -e $osnapfile* ]; then','        mv $osnapfile* "${osnapbase}${iext}"','    fi',
                  '    if [ -e $ocoeffile  ]; then','        mv $ocoeffile  "${ocoefbase}${iext}"','    fi',
                  '    if [ -e $outfile    ]; then','        mv $outfile    "${outbase}${iext}"','    fi',
                  '    if [ -e $logfile    ]; then','        mv $logfile    "${logbase}${iext}"','    fi',
                  '    if [ -e $chkfile    ]; then','        mv $chkfile    "${chkbase}${iext}"','    fi',
                  '    if [ -e $elfile*    ]; then','        mv $elfile*    "${elbase}${iext}"','    fi',
                  '    if [ -e $olilfile*  ]; then','        mv $olilfile*  "${olilbase}${iext}"','    fi',
                  'done']
        # Commands to finalize files
        cmdlist+=[' ',
                  '# Move files back to original paths',
                  'rm $execfile $exec',
                  'rm $parfile $par',
                  'rm $modfile $mod',
                  'rmdir $tempdir',
                  'cd ../']
        # Write commands to file
        mmlio.rwbash('W',wrapfile,cmdlist,overwrite=overwrite)
        # Make file executable
        os.chmod(wrapfile,0755)
    # Return output
    return cmdlist

####################################################################################################################################
# METHOD TO READ IN SCF COEFFICIENTS
####################################################################################################################################
def read_coeffile(fname):
    '''
    NAME:
        mmlscf.read_coeffile
    PURPOSE:
        To read in output coefficients from an SCF ocoef file.
    CALLING:
        coefdict=read_coeffile(fname)
    ARGUMENTS:
        fname:
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
    print time,nmax,lmax
    # Initialize dictionary
    coefdict=dict(time=float(time),nmax=int(float(nmax)),lmax=int(float(lmax)))
    print coefdict
    coefdict['sinsum']=np.zeros((nmax,lmax,lmax))
    coefdict['cossum']=np.zeros((nmax,lmax,lmax))
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
