#!/usr/bin/python
####################################################################################################################################
#
# MEAGAN LANG'S SCF METHODS
#
####################################################################################################################################
import sys,os,shutil,glob,copy,pprint,scipy,math,re,subprocess
import numpy as np
from mmlutils import *
from pysim import files

FLAG_RUNSERIAL=True

FILELIST_SETUP=['exec','input','snapin']#,'coefin']                                                                                                                                  
FILELIST_OUTPT=['output']
FILELIST_CLEAN=['output']

FILELIST_SPLIT={}

METHLIST_SETUP=[]
METHLIST_OUTPT=[]
METHLIST_CLEAN=[]

VALIDSNAPKEY='osnap'

INCLEXT_SETUP={}#'inclext_fchk':{'isnap':'osnap'}}
INCLEXT_OUTPT={}
INCLEXT_CLEAN={}
EXCLEXT_SETUP={'exclext_fchk':{'isnap':'osnap'}}
EXCLEXT_OUTPT={}
EXCLEXT_CLEAN={}

FEXT_UNITFILE=['units','scfunits']


####################################################################################################################################
####################################################################################################################################
# TESTING METHODS
def test_performance(compid='stampede'):
    """
    Plot performance info for the SCF code
    """
    # Data
    if compid=='stampede':
        npart=[long(1e4),long(1e5),long(1e6),long(1e7)]
        nproc=[1,10,50]
        data=[[ 5., 6.,28.,67.],
              [ 6., 7.,22.,53.],
              [20.,13.,27.,55.]]
    else: raise Exception('No performance data for {}'.format(compid))
    # Plot
    plotfile='/home/langmm/code/python/mine/pysim/performance/perf_scf.png'
    mmlparallel.plot_performance(nproc,npart,data,fname=plotfile,timelabel='Time/particle/snap (s)')
    print '    '+plotfile
    # Return
    return

####################################################################################################################################
####################################################################################################################################
# METHODS TO CREATE SCF FILES

####################################################################################################################################
# METHOD TO CREATE SCF SNAPSHOTS
def mk_snap(simstr,fext=None,overwrite=None,snaptype0='gadget',**exkw):
    """
    Creates SCF snapshots
    """
    nbodsmax=2048000
    from pysim import files,nbody
    typlist=['halo','disk','bulge','gas']
    #typlist=['disk']
    gallist=[1,2]
    # Get files
    srckey,srcpath,srclist=files.listsnap(simstr,snaptype0,fext=fext)
    dstpath=simstr.fdict['scf']['input']['isnap']
    extfmt=files.get_extfmt(srcpath)
    if len(srclist)==0:
        raise Exception('No {} snapshots were found matching the pattern {}'.format(snaptype0.upper(),srcpath))
    # Create directories
    mmlfiles.mkdirs(os.path.dirname(dstpath))
    # Ask about overwrite
    if overwrite is None: overwrite=mmlio.yorn('Overwrite existing snapshots?')
    # Loop over source snapshots
    dstlist=[]
    for src in srclist:
        snapObj=None
        ext=int(float(files.get_ext(src,srcpath)))
        # Loop over particle types
        for typ in typlist:
            f=nbody.family.get_family(typ)
            # Loop over galaxies
            for gal in gallist:
                g=nbody.galaxy.get_galaxy(str(gal))
                # Check for destination file
                dst=dstpath.format(ext,typ+str(gal))
                #mmlio.verbose(str((dst,os.path.isfile(dst),overwrite)))
                if not overwrite and os.path.isfile(dst): continue
                # Read in source snapshot/select particle type
                if not snapObj: snapObj=simstr.loadsnap(fext=ext,ftype=snaptype0)
                if g not in snapObj.galaxy_names(): continue
                try:
                    sObj0=snapObj[f][g]
                except KeyError:
                    # print f,snapObj.families()
                    # print g,snapObj.galaxies()
                    continue
                if len(sObj0)==0: continue
                # Down sample if too many particles
                if len(sObj0)>nbodsmax:
                    mfact=float(len(sObj0))/float(nbodsmax)
                    sObj=sObj0[np.sort(np.random.random_integers(0,len(sObj0)-1,nbodsmax))]
                else: 
                    mfact=1.0
                    sObj=sObj0
                sObj.properties['time']=snapObj.properties['time']
                # Center
                nbody.analysis.halo.center(sObj,mode='com')
                # Write
                sObj.write(filename=dst,ftype='scf',overwrite=overwrite,massmult=mfact)
                dstlist.append(dst)
                print '    '+dst
    # Return control
    return dstlist

####################################################################################################################################
# METHOD TO CREATE SCF EXEC FILES
def mk_exec(simstr,**exkw):
    """
    Creates a local copy of executable files
    """
    # Select default files
    exdeft=simstr.finfo['scf']['exec']['default']+'_'+simstr['scfcomp']
    excopy=simstr.finfo['scf']['exec']['path']
    # Check existence
    if not os.path.isfile(exdeft): raise Exception('Executable does not exist: {}'.format(exdeft))
    # Copy files to local directory
    mmlfiles.cpfiles(exdeft,excopy,overwrite=True)
    # Return
    return excopy

####################################################################################################################################
# METHOD TO CREATE SCF PARAMETER FILE
def mk_param(simstr=None,runtag=None,overwrite=None,askuser=False,**exkw):
    """
    Creates a parameter file
    """
    import pysim
    # Pars input
    if simstr==None: simstr=pysim.loadsim(runtag)
    if not isinstance(overwrite,bool):
        if askuser: overwrite=mmlio.yorn('Overwrite existing static parameter file?')
        else      : overwrite=mmlpars.mml_pars(overwrite,type=bool,default=False)
    # Get file names
    fpar_stat=simstr.fdict['scf']['icgen']['stat_pm']
    fpar_copy=simstr.fdict['scf']['input']['param'  ]
    # If static parameter file dosn't exist, create it
    if overwrite or not os.path.isfile(fpar_stat):
        # Initialize static parameter file from default
        pardict=rw_param('R',simstr.finfo['scf']['stat_pm']['default'])
        # Add simstr info
        pardict['headline']=simstr['runtag']
        # Save parameters
        rw_param('W',fpar_stat,pardict,overwrite=overwrite)
        # Check parameter file
        print fpar_stat
        mmlio.yorn('Make changes to the above parameter file now!')
    # Make local copy of the parameter file
    mmlfiles.cpfiles(fpar_stat,fpar_copy,overwrite=True)
    # Return output
    return fpar_copy

####################################################################################################################################
# METHOD TO CREATE SCF MODIFICATION FILE
def mk_mods(simstr=None,runtag=None,overwrite=None,askuser=False,**exkw):
    """
    Creates modification files
    """
    import pysim
    # Pars input
    if simstr==None: simstr=pysim.loadsim(runtag)
    if not isinstance(overwrite,bool):
        if askuser: overwrite=mmlio.yorn('Overwrite existing static parameter file?')
        else      : overwrite=mmlpars.mml_pars(overwrite,type=bool,default=False)
    # Get file names
    fmod_stat=simstr.fdict['scf']['icgen']['stat_md']
    fmod_copy=simstr.fdict['scf']['input']['mods'   ]
    # If static modification file dosn't exist, create it
    if overwrite or not os.path.isfile(fmod_stat):
        # Initialize static modification file from default
        moddict=rw_mods('R',simstr.finfo['scf']['stat_md']['default'])
        # Save modifications
        rw_mods('W',fmod_stat,moddict,overwrite=overwrite)
        # Check modification file
        print fmod_stat
        mmlio.yorn('Make changes to the above modification file now!')
    # Make local copy of the parameter file
    mmlfiles.cpfiles(fmod_stat,fmod_copy,overwrite=True)
    # Return output
    return fmod_copy

####################################################################################################################################
# METHODS TO CREATE SUBMISSION SCRIPTS
def mk_pbs(simstr,**kwargs): return mk_submit(simstr,'pbs',**kwargs)
def mk_sbatch(simstr,**kwargs): return mk_submit(simstr,'sbatch',**kwargs)
def mk_submit(simstr,ftype,verbose=False,runserial=FLAG_RUNSERIAL,**kwargs):
    """
    Creates submission scripts
    """
    # Pars input
    ftype=mmlpars.mml_pars(ftype,list=mmlio.DICT_SUBOPT.keys())
    # Get file names
    files_run=files.get_fdict(compid=simstr['scfcomp'],**simstr)
    subfile=simstr.fdict['scf']['input'][ftype]
    # Get run command info
    wrpfile=os.path.basename(files_run['scf']['input']['wrap'])
    mpicmd=mmlinfo.COMPDICT_MPI[simstr['scfcomp']]
    # Variables from simstr
    if 'timetest' in simstr['runtag']: name='scf'+simstr['subtag']+'_'+simstr['runtag']
    else                             : name='scf_'+simstr['runtag']+simstr['subtag']
    subdict={'jobname' :name,
             'nproc'   :int(simstr['nprocscf']),
             'ppn'     :int(mmlinfo.COMPDICT_PPN[simstr['scfcomp']]),
             'pmem'    :int(mmlsim2mcpu(simstr)),
             'twall'   :float(mmlsim2tcpu(simstr)),
             'outfile' :files_run['scf']['output']['runout'],
             'cmdlist' :['mkdir -p {}'.format(files_run['scf']['output']['dir'])]}
    # Stampede stuff
    if simstr['scfcomp']=='stampede': 
        subdict['account']=files.STAMPEDE_ACCT
        runserial=False
    # Run command
    if runserial: subdict['cmdlist'].append('{} -np {} -comm none ./{} 0'.format(mpicmd,subdict['nproc'],wrpfile))
    else        : subdict['cmdlist'].append('./{}'.format(wrpfile))
    # Exit command
    # subdict['cmdlist']+=['','wait','echo finished mpiexec','exit 0']
    # Write submission script
    mmlio.rw_submit('W',subfile,ftype,subdict,verbflag=verbose,overwrite=True)
    # Return output
    return subfile

####################################################################################################################################
# METHOD TO CREATE MPISCF WRAPPER
def mk_wrapper(simstr,runserial=FLAG_RUNSERIAL,**exkw):
    """
    Creates a wrapper script
    """
    extvar='${iext}'
    cmdlist=['start_time=`date +%s`',
             'nproc={}'.format(simstr['nprocscf']),
             'rank=${MPIEXEC_RANK}',
             'if (($nproc == 1)); then',
             '    rank=0',
             'fi']
    # Get file names
    wrpfile=simstr.fdict['scf']['input']['wrap']
#    files_run,finfo_run,ftab_run=files.get_fdict(compid=simstr['scfcomp'],**simstr)
    flist_run=files.get_flist('scf',compid=simstr['scfcomp'],**simstr)
    # Get size sorted lists of files
    flist={'ismall':[],'ilarge':[],'osmall':[],'olarge':[]}
    for igrp in flist_run['keylist']:
        for finfo in copy.deepcopy(flist_run[igrp]):
            if finfo['alias']=='None': continue
            # Extra fields
            finfo['base']=finfo['path'].split('{')[0]
            finfo['ext' ]=extvar
            # Input files
            if   finfo['class']=='input' :
                if   finfo['size']=='small': flist['ismall'].append(finfo)
                elif finfo['size']=='large': flist['ilarge'].append(finfo)
            # Output files
            elif finfo['class']=='output':
                if   finfo['size']=='small': flist['osmall'].append(finfo)
                elif finfo['size']=='large': flist['olarge'].append(finfo)
            # Special files directly referenced
            if finfo['key']=='isnap'  : loopbase=finfo['base' ]
            if finfo['key']=='osnap'  : chckbase=finfo['base' ]
            if finfo['key']=='exec'   : execbase=finfo['alias']
            if finfo['key']=='procout': poutbase=finfo['alias']
    # Command to initialize directory
    if runserial: tempdir='"temprun${rank}"'
    else        : tempdir='"temprun"'
    cmdlist+=[' ',
              '# Move files to temporary directory w/ built in SCF names',
              "tempdir={}".format(tempdir),
              "echo $rank",
              "echo $tempdir",
              "if [ ! -e $tempdir ]; then","    mkdir $tempdir","fi","cd $tempdir",
              "if [ -e 'scflog' ]; then","    rm '*'","fi"]
    # Copy input files for the entire run
    for finfo in flist['ismall']: cmdlist+=['cp {} {}'.format(finfo['path'],finfo['alias'])]
    # Commands to search for files
    cmdlist+=[' ',
              'loopbase="{}"'.format(loopbase),
              #'lenloopbase=${#loopbase}',
              'filelist=($loopbase*)']
    # Run serial code on each snapshot in parallel
    if runserial: 
        # Commands to identify subset of variables
        cmdlist+=[' ',
                  'neach=$((${#filelist[@]}/$nproc))',
                  'if (($((${#filelist[@]} % $nproc)) > 0)); then',
                  '    neach=$((${#filelist[@]}/($nproc - 1)))',
                  'fi',
                  'if (($neach == 0)); then',
                  '    neach=1',
                  'fi',
                  'idx_beg=$(($rank * $neach))',
                  'idx_end=$(($idx_beg + $neach - 1))']
    # Run parallel code in series
    else:
        # Commands to identify subset of variables
        cmdlist+=[' ',
                  'idx_beg=0',
                  'idx_end=$((${#filelist[@]} - 1))']
    # Commands to perform loop
    cmdlist+=[' ',
              '# Loop over SCF input snapshots',
              'let i=$idx_beg',
              'while ((i<=$idx_end)) && ((i<${#filelist[@]})); do',
              '    # Select proper file',
              '    iloop=${filelist[i]}']
    # Commands to get extension and print info
    cmdlist+=['    # Get file extension',
              '    iext="${iloop##$loopbase}"']
    # Commands to check that output dosn't exist
    cmdlist+=[' ',
              '    # Only continue if output is missing',
              '    if [ -e "{}{}" ]; then'.format(chckbase,extvar),
              '        rm "{}{}"'.format(loopbase,extvar),
              '        let i++',
              #'        echo "Skipping snapshot ${iext}"',
              '        continue',
              '    fi']
    if runserial: cmdlist.append('    echo "Starting snapshot ${iext} on ${rank} ${iloop}"')
    else        : cmdlist.append('    echo "Starting snapshot ${iext}: ${iloop}"')
    # Move per snapshot input files
    cmdlist+=[' ','    # Move necessary input files']
    for finfo in flist['ilarge']:
        cmdlist+=['    if [ -e "{base}{ext}" ]; then'.format(**finfo),
                  '        mv "{base}{ext}" "{alias}"'.format(**finfo),
                  '    fi']
    # Run mpiscf
    cmdlist+=[' ','    # Run mpiscf script']
    if runserial: cmdlist+=['    ./{}'.format(execbase)]
    else        : cmdlist+=['    {} -np $nproc ./{}'.format(mmlinfo.COMPDICT_MPI[simstr['scfcomp']],execbase)]
    # Move persnapshot output files
    cmdlist+=[' ','    # Recover files']
    for finfo in flist['ilarge']+flist['olarge']:
        cmdlist+=['    if [ -e "{alias}" ]; then'.format(**finfo),
                  '        mv "{alias}" "{base}{ext}"'.format(**finfo),
                  '    fi']
        #mmlfiles.mkdirs(os.path.dirname(finfo['path']),host=simstr['scfcomp'])
    # Remove snapshot input files if output exists
    cmdlist+=[' ',
              '    # Remove input if output exists',
              '    if [ -e "{}{}" ]; then'.format(chckbase,extvar),
              '        rm "{}{}"'.format(loopbase,extvar),
              '    fi']
    # Close for loop
    cmdlist+=['','    # Advance file count','    let i++',
              '    end_time=`date +%s`',
              '    echo execution time = `expr $end_time - $start_time` s',
              'done']
    # Commands to finalize files
    cmdlist+=[' ','# Move files back to original paths']
    for finfo in flist['ismall']: cmdlist.append('rm {alias}'.format(**finfo))
    cmdlist+=['cd ../','wait','rmdir $tempdir','exit 0']
    # Write commands to file
    rw_wrapper('W',wrpfile,linelist=cmdlist,overwrite=True)
    # Make file executable
    subprocess.call(['chmod','755',wrpfile])
    # Return output
    return wrpfile

####################################################################################################################################
# METHOD TO CREATE SCF UNITS FILE
def mk_units(simstr,snaptype0='gadget',**exkw):
    """
    Creates a units file
    """
    # Get original units file & SCF unit file
    parmfile=files.unitsfile(simstr.fdict,snaptype0)
    if 'timetest' in simstr['runtag'] and not os.path.isfile(parmfile):
        parmfile=parmfile.replace('on{}'.format(simstr['nprocscf']),'on{}'.format(256))
    unitfile=simstr.fdict['scf']['input']['units']
    # Make necessary directory
    mmlfiles.mkdirs(os.path.dirname(unitfile))
    # Copy original units file to SCF unit file
    mmlfiles.cpfiles(parmfile,unitfile,overwrite=True)
    # Return output
    return unitfile

####################################################################################################################################
####################################################################################################################################
# METHODS TO READ SCF FILES

####################################################################################################################################
# METHOD TO READ/WRITE PARAMETER FILES
def rw_param(rwid,fname,fdict=None,overwrite=False):
    """
    Reads/writes parameter files
    """
    # Set constants
    keywidth=29
    headline='**********Basic input parameters**********'
    tailline='******************************************'
    # Read
    if   rwid=='R': 
        fdict=mmlio.rwdict('R',fname,style='fortran')
        fdict=mmlparam.parspar('pysim.files.scf','param',inpar=fdict)
        return fdict
    # Write
    elif rwid=='W': 
        fdict=mmlparam.parspar('pysim.files.scf','param',inpar=fdict)
        fdict['keylist']=mmlparam.listpar('pysim.files.scf','param')['keylist']
        mmlio.rwdict('W',fname,fdict,width=keywidth,style='fortran',
                     headline=headline,tailline=tailline,overwrite=overwrite)
        return
    # Error
    else: raise Exception('Invalid rwid: {}'.format(rwid))
    return

####################################################################################################################################
# METHOD TO READ/WRITE MODIFICATION FILE
def rw_mods(rwid,fname,fdict=None,overwrite=False):
    """
    Reads/writes modification files
    """
    # Set constants
    keywidth=29
    headline='**********Basic input parameters**********'
    tailline='******************************************'
    # Read
    if   rwid=='R':
        fdict=mmlio.rwdict('R',fname,style='fortran')
        fdict=mmlparam.parspar('pysim.files.scf','mods',inpar=fdict)
        return fdict
    # Write
    elif rwid=='W':
        fdict=mmlparam.parspar('pysim.files.scf','mods',inpar=fdict)
        fdict['keylist']=mmlparam.listpar('pysim.files.scf','mods')['keylist']
        mmlio.rwdict('W',fname,fdict,width=keywidth,style='fortran',
                     headline=headline,tailline=tailline,overwrite=overwrite)
        return
    # Error
    else: raise Exception('Invalid rwid: {}'.format(rwid))
    return

####################################################################################################################################
# METHOD TO READ/WRITE WRAPPER SCRIPT
def rw_wrapper(*args,**kwargs):
    """
    Reads/writes wrapper scripts
    """
    return mmlio.rw_bash(*args,**kwargs)

####################################################################################################################################
# METHOD TO READ/WRITE UNITS FILES
def rw_units(*args,**kwargs):
    """
    Reads/writes units files
    """
    return mmlio.rwdict(*args,**kwargs)

####################################################################################################################################
# METHOD FOR READING OUTPUT SNAPSHOTS
def rw_isnap(*args,**kws): 
    kws['flag_out']=False
    return rw_snap(*args,**kws)
def rw_osnap(*args,**kws):
    kws['flag_out']=True
    return rw_snap(*args,**kws)
def rw_snap(rwid,fname,fdict=None,flag_out=False,overwrite=False,massmult=1.0):
    """
    Reads/writes snapshot files
    """
    namedict={'out':('mass','x','y','z','vx','vy','vz','pot','ad','dti','mi','pm2'),
              'in' :('mass','x','y','z','vx','vy','vz')}
    frmtdict={'out':('f','f','f','f','f','f','f','f','f','f','f','f'),
              'in' :('f','f','f','f','f','f','f')}
    # Read
    if   rwid=='R':
        # Header
        fid=open(fname,'r')
        headline=fid.readline()
        header={}
        headout=headline.split()
        flag_out=(len(headout)==3)
        if flag_out: nbody,time,bhmass=headout ; header['bhmass']=float(bhmass)
        else       : nbody,time=headout
        header['nbody']=long(float(nbody))
        header['time' ]=float(time)
        # Body
        if flag_out: 
            names=namedict['out']
            frmts=frmtdict['out']
        else       : 
            names=namedict['in' ]
            frmts=frmtdict['in' ]
        fdata=np.loadtxt(fid,dtype={'names':names,'formats':frmts})
        fid.close()
        fdict={nm:fdata[nm] for nm in names}
        fdict['pos']=np.vstack(tuple(fdict[nm] for nm in [ 'x', 'y', 'z'])).T
        fdict['vel']=np.vstack(tuple(fdict[nm] for nm in ['vx','vy','vz'])).T
        for nm in ['x','y','z','vx','vy','vz']: del fdict[nm]
        # Add header and return
        fdict['header']=header
        return fdict
    # Write
    elif rwid=='W':
        if not overwrite and os.path.isfile(fname): return
        # Header
        flag_out=('bhmass' in fdict['header'])
        if flag_out: headline='{nbody}    {time}    {bhmass}\n'.format(**fdict['header'])
        else       : headline='{nbody}    {time}\n'.format(**fdict['header'])
        fid=open(fname,'w')
        fid.write(headline)
        # Body
        if flag_out: 
            names=namedict['out']
            frmts=frmtdict['out']
        else       : 
            names=namedict['in' ]
            frmts=frmtdict['in' ]
        fdict[ 'x']=fdict['pos'][:,0] ; fdict[ 'y']=fdict['pos'][:,1] ; fdict[ 'z']=fdict['pos'][:,2]
        fdict['vx']=fdict['vel'][:,0] ; fdict['vy']=fdict['vel'][:,1] ; fdict['vz']=fdict['vel'][:,2]
        fdict['mass']=copy.deepcopy(fdict['mass'])*massmult
        frmts=tuple('%'+ifmt for ifmt in frmts)
        farry=np.vstack(tuple(fdict[inam] for inam in names)).T
        np.savetxt(fid,farry,fmt=frmts)#,header=headline)
        # Return control
        fid.close()
        return
    # Error
    else: raise Exception('Invalid rwid: {}'.format(rwid))
    return

####################################################################################################################################
# METHOD TO READ IN SCF COEFFICIENTS
def rw_icoef(*args,**kws):
    kws['flag_out']=False
    return rw_coef(*args,**kws)
def rw_ocoef(*args,**kws):
    kws['flag_out']=True
    return rw_coef(*args,**kws)
def rw_coef(rwid,fname,fdict=None,flag_out=False,overwrite=False):
    """
    coefdict: Dictionary with SCF coefficient keys:
        time:   Simulation time that coefficients are written out
        nmax:   Maximum order of n coefficient
        lmax:   Maximum order of l coefficient
        sinsum: [nmax,lmax,lmax] array of n,l,m coefficients
        cossum: [nmax,lmax,lmax] array of n,l,m coefficients
    """
    # Read
    if rwid=='R':
        fid=open(fname,'r')
        # Read first line
        headline=fid.readline()
        time,nmax,lmax=headline.split()
        time=float(time)
        nmax=int(float(nmax))
        lmax=int(float(lmax))
        print time,nmax,lmax
        # Initialize dictionary
        fdict=dict(time=time,nmax=nmax,lmax=lmax)
        fdict['sinsum']=np.zeros((nmax+1,lmax+1,lmax+1))
        fdict['cossum']=np.zeros((nmax+1,lmax+1,lmax+1))
        # Loop over coefficients reading them in
        for n in range(fdict['nmax']+1):
            for l in range(fdict['lmax']+1):
                for m in range(l+1):
                    iline=fid.readline()
                    isinsum,icossum=iline.split()
                    fdict['sinsum'][n][l][m]=float(isinsum)
                    fdict['cossum'][n][l][m]=float(icossum)
        # Close file & return output
        fid.close()
        return fdict
    # Write
    elif rwid=='W':
        if not overwrite and os.path.isfile(fname): return
        fid=open(fname,'w')
        # Header
        headline='{time}    {nmax}    {lmax}'.format(**fdict)
        fid.write(headline)
        # Body
        for n in range(fdict['nmax']+1):
            for l in range(fdict['lmax']+1):
                for m in range(l+1):
                    isinsum=str(fdict['sinsum'][n][l][m])
                    icossum=str(fdict['cossum'][n][l][m])
                    iline='    '.join(isinsum,icossum)
                    fid.write(iline)
        # Close and return control
        fid.close()
        return
    # Error
    else: raise Exception('Invalid rwid: {}'.format(rwid))
    return

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
    mcpu=float(simstr.ntot)*(100000./float(2e7))
    # Return rounded output
    return mmlmath.oom(mcpu/simstr['nprocscf'],nsig=1,method='CEIL')

####################################################################################################################################
# METHOD TO CALCULATE TIME REQUIRED TO RUN SCF ON ACCRE
def mmlsim2tcpu(simstr):
    """
    Estimates the amount of time required to run SCF
    """
    # Force reasonable values as a safety
    if   simstr['runtyp'].lower()=='galsim': tcpu=1. # hr
    elif simstr['runtyp'].lower()=='intsim': tcpu=2. # hr
    if 'timetest' in simstr['runtag']: tcpu=24.
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
def listpar(partyp,**exkw):
    partyp=mmlpars.mml_pars(partyp,type=str)
    if   partyp=='param': par={'keylist':['headline','nsteps','noutbod','noutlog','dteps','G',
                                          'tfinal','multistep','fixedn','selfgrav','inptcoef','outpcoef',
                                          'zeroodd','zeroeven','fixacc','rcrit','ecrit','lilout','nlilout'],
                               'headline'  :{'def':'galaxy','form':str  },
                               'nsteps'    :{'def':0       ,'form':int  },
                               'noutbod'   :{'def':1       ,'form':int  },
                               'noutlog'   :{'def':1       ,'form':int  },
                               'dteps'     :{'def':0.      ,'form':float},
                               'G'         :{'def':1.0     ,'form':float},
                               'tfinal'    :{'def':0.      ,'form':float}, # 100.
                               'multistep' :{'def':False   ,'form':bool },
                               'fixedn'    :{'def':False   ,'form':bool },
                               'selfgrav'  :{'def':True    ,'form':bool },
                               'inptcoef'  :{'def':False   ,'form':bool },
                               'outpcoef'  :{'def':True    ,'form':bool },
                               'zeroodd'   :{'def':False   ,'form':bool },
                               'zeroeven'  :{'def':False   ,'form':bool },
                               'fixacc'    :{'def':False   ,'form':bool },
                               'rcrit'     :{'def':0.      ,'form':float},
                               'ecrit'     :{'def':0.      ,'form':float},
                               'lilout'    :{'def':False   ,'form':bool },
                               'nlilout'   :{'def':0       ,'form':int  }}
    elif partyp=='mods' : par={'keylist':['iseed','bhmass','epsbh','tstartbh','tgrowbh','tlivebh','tdiebh',
                                          'xdrag','ydrag','zdrag','tstartdrag','tgrowdrag','tlivedrag','tdiedrag',
                                          'bhgrav','usedrag','stellev'],
                               'iseed'     :{'def':3587    ,'form':int  },
                               'bhmass'    :{'def':0.      ,'form':float},
                               'epsbh'     :{'def':0.      ,'form':float},
                               'tstartbh'  :{'def':0.      ,'form':float},
                               'tgrowbh'   :{'def':0.      ,'form':float},
                               'tlivebh'   :{'def':0.      ,'form':float},
                               'tdiebh'    :{'def':0.      ,'form':float},
                               'xdrag'     :{'def':0.      ,'form':float},
                               'ydrag'     :{'def':0.      ,'form':float},
                               'zdrag'     :{'def':0.      ,'form':float},
                               'tstartdrag':{'def':0.      ,'form':float},
                               'tgrowdrag' :{'def':0.      ,'form':float},
                               'tlivedrag' :{'def':0.      ,'form':float},
                               'tdiedrag'  :{'def':0.      ,'form':float},
                               'bhgrav'    :{'def':False   ,'form':bool },
                               'usedrag'   :{'def':False   ,'form':bool },
                               'stellev'   :{'def':False   ,'form':bool }}
    else: raise Exception('Invalid parameter type: {}'.format(partyp))
    return par

