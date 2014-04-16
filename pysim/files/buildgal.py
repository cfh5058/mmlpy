#!/usr/bin/python
####################################################################################################################################
#
# MEAGAN LANG'S BUILDGAL METHODS
#
####################################################################################################################################
import sys,os,shutil,glob,copy,pprint,scipy,math,subprocess,random,re
import numpy as np

LIST_PTYPS=['halo','disk','bulge','gas']
N_PTYP=len(LIST_PTYPS)

LIST_PARMETH=['basepar','diskpar','gaspar','gaspar_flatdisk','bulgepar','bulgepar_selfgrav','bulgepar_nonspher','bulgepar_rotate',
              'halopar','halopar_selfgrav','halopar_NF','halopar_IS','halopar_LH','satpar','satpar_selfgrav','twomodpar']
LIST_HALOPROFS=['LH','NF','IS']
from mmlutils import *
from mmlastro import mmlprofiles,mmlcosmo,mmlconst

FILELIST_SETUP=['exec','input','longin']
FILELIST_OUTPT=['ic','output','longout']
FILELIST_CLEAN=['ic','output','longout']

FILELIST_SPLIT={}

METHLIST_SETUP=[]
#METHLIST_OUTPT=['mk_ic']
METHLIST_OUTPT=['mk_ic2']
METHLIST_CLEAN=[]

VALIDINITKEY='stat_ic'
VALIDSNAPKEY='ic'

INCLEXT_SETUP={}
INCLEXT_OUTPT={}
INCLEXT_CLEAN={}
EXCLEXT_SETUP={}
EXCLEXT_OUTPT={}
EXCLEXT_CLEAN={}

FEXT_UNITFILE=['units','bgunits']

####################################################################################################################################
# METHOD TO MAKE BUILDGAL IC WITHOUT REPEATS
def mk_ic(simstr,pNbody=None,**exkw):
    """
    Creates merged IC file with repeats removed
    """
    from pysim import files,nbody
    # Ask about using pNbody
    if pNbody is None:
        pNbody=mmlio.yorn('Use pNbody to remove repeats?')
    if pNbody:
        from pysim.files import gadget
    # Get filenames
    funit=simstr.fdict['buildgal']['input']['units']
    fin=simstr.fdict['buildgal']['icgen']['stat_mg']
    fout=simstr.fdict['buildgal']['icgen']['stat_ic']
    # Load parameters
    par0=mmlio.rwdict('R',simstr.fdict['buildgal']['icgen']['stat_pm'],style='fortran')
    par0=model2param(simstr.infodict,scale=True,par0=par0)
    # Determine minimum separation
    sftmax=0.
    for ityp in LIST_PTYPS:
        if ityp+'_n' not in par0: continue
        sftmax=max(sftmax,par0[ityp+'_eps'])
    sftmax=mmlio.askquest('What should the minimum separation be?',dtype='float',default=sftmax)
    # Create merged IC if it does not exist
    if not os.path.isfile(fin): pnb=merge_ic(simstr,pNbody=pNbody,**exkw)
    else                      : pnb=simstr.loadsnap(fname=fin,ftype='buildgal',pNbody=pNbody)
    mmlio.verbose('N (initial)={}'.format(len(pnb)))
    # Remove repeats
    if pNbody: 
        idx_rep=pnb.index_norepeats(cut=sftmax)
        pnb_rep=pnb.selecti(idx_rep)
    else:
        idx_rep=nbody.analysis.langmm.index_norepeats(pnb,cut=sftmax)
        pnb_rep=pnb[idx_rep]
    mmlio.verbose('N (no reps)={}'.format(len(pnb_rep)))
    # Remove extra particles
    idx_ext=np.array([],dtype=int)
    for ityp in LIST_PTYPS:
        if ityp+'_n' not in par0: continue
        if par0[ityp+'_n']==0: continue
        ifam=nbody.family.get_family(ityp)
        # Select particles by type
        if pNbody: itypslc=pnb_rep.get_index(ptyp=ityp)
        else     : itypslc=pnb_rep._get_family_slice(ityp)
        itypidx=np.arange(len(pnb_rep),dtype=long)[itypslc]
        iN=len(itypidx) ; fN=par0[ityp+'_n']
        if fN==0: continue
        mmlio.verbose('    {}: Ni={}, Nf={}'.format(ityp,iN,fN))
        # Add random sample of particles
        if fN>iN: raise Exception('There are too few {} particles.'.format(ityp))
        np.random.shuffle(itypidx)
        idx_ext=np.append(idx_ext,itypidx[:fN])
        # Fix masses
        if pNbody:
            pnb_rep.massarr[gadget.LIST_PTYPS.index(ityp)]=par0[ityp+'_mp']
            pnb_rep.mass[itypidx]=par0[ityp+'_mp']
        else: pnb_rep['mass'][itypidx]=par0[ityp+'_mp']
    # Remove extra particles
    if pNbody: pnb_ext=pnb_rep.selecti(np.sort(idx_ext))
    else     : pnb_ext=pnb_rep[np.sort(idx_ext)]
    mmlio.verbose('N (no extr)={}'.format(len(pnb_ext)))
    # Fix masses
    for ityp in LIST_PTYPS:
        if pNbody: itypslc=pnb_ext.get_index(ptyp=ityp)
        else     : itypslc=pnb_ext._get_family_slice(ityp)
        itypidx=np.arange(len(pnb_rep),dtype=long)[itypslc]
        iN=len(itypidx)
        if iN==0: continue
        if pNbody: Mtot=pnb_ext.mass[itypidx].sum()
        else     : Mtot=pnb_ext['mass'][itypidx].sum()
        mmlio.verbose('    {} ({}): Mi={}, Mf={}'.format(ityp,iN,Mtot,par0[ityp+'_mtot']))
    # Save
    if pNbody:
        pnb_ext.rename(fout)
        pnb_ext.write()
    else: pnb_ext.write(filename=fout,ftype='buildgal')
    # Return control
    return pnb


def mk_ic2(simstr,pNbody=None,proclist=None,Nsave=100,**exkw):
    """
    Creates merged IC file with repeats removed
    """
    from pysim import files,nbody
    # Ask about using pNbody
    if pNbody is None:
        pNbody=mmlio.yorn('Use pNbody to remove repeats?')
    if pNbody:
        from pysim.files import gadget
    # Get filenames
    funit=simstr.fdict['buildgal']['input']['units']
    fin_extfmt=files.get_extfmt(simstr.fdict['buildgal']['output']['ic'])
    fin_base=re.sub(fin_extfmt,'*',simstr.fdict['buildgal']['output']['ic'])
    fout=simstr.fdict['buildgal']['icgen']['stat_ic']
    # Count files
    fin_list=glob.glob(fin_base)
    if proclist is None: proclist=range(len(fin_list))
    nproc=len(proclist)
    # Check number of processors
    if nproc==0: raise Exception('There are not any IC files to merge.')
    if nproc!=simstr['nprocmk']: 
        if not mmlio.yorn('{} files were expected but there are {}. OK?'.format(simstr['nprocmk'],nproc)): return
    # Load parameters
    par0=mmlio.rwdict('R',simstr.fdict['buildgal']['icgen']['stat_pm'],style='fortran')
    par0=model2param(simstr.infodict,scale=True,par0=par0)
    # Determine minimum separation
    sftmax=0. ; typList=[] ; typNum=[]
    for ityp in LIST_PTYPS:
        if ityp+'_n' not in par0: continue
        if par0[ityp+'_n']==0: continue
        typNum.append(par0[ityp+'_n'])
        typList.append(ityp)
        sftmax=max(sftmax,par0[ityp+'_eps'])
    sftmax=mmlio.askquest('What should the minimum separation be?',dtype='float',default=sftmax)
    # Check if file should be loaded
    pnb=None ; Ntot0=0 ; Nreq0=par0['total_n'] ; sproc=0
    if os.path.isfile(fout):
        if mmlio.yorn('Load existing stat_ic file as a starting point?'):
            sproc=mmlio.askquest('Enter processor number to start at:',dtype='int',default=min(proclist))
            pnb=simstr.loadsnap(fname=fout,ftype='buildgal',pNbody=pNbody)
            Ntot0=len(pnb)
    proclist=[ip for ip in proclist if ip>=sproc]
    Nrem0=Nreq0-Ntot0
    NPmin=Nreq0/min(typNum)
    # Load each file and merge them
    if Nrem0==0: mmlio.verbose('All particles present.')
    else:
        # Nperproc=Nrem0/nproc
        for iproc in proclist:
            mmlio.verbose('[Nproc={:04d}/{:04d}] Remaining: {: 8d} ({:3.1f} %)'.format(iproc,nproc,Nrem0,100.*float(Nrem0)/float(Nreq0)))
            inb = simstr.loadsnap(fext=iproc,ftype='buildgal',pNbody=pNbody)
            if pnb==None: pnb=inb
            else        : pnb=nbody.analysis.langmm.merge(pnb,inb,N1min=NPmin*iproc,N2min=NPmin,typList=typList,
                                                          pNbody=pNbody,repcut=sftmax,keepTypeRatio=True,warn=True)
            Ntot0=len(pnb)
            Nrem0=Nreq0-Ntot0
            mmlio.verbose('    [Added {: 6d} (min {: 6d})] Ntot (current) = {: 8d}/{: 8d}'.format(len(inb),NPmin,Ntot0,Nreq0))
            del inb
            # Save
            if (iproc!=proclist[0] and iproc%Nsave==0) or iproc==proclist[-1]:
                mmlio.verbose('Saving after processor: {:04d}'.format(iproc))
                if pNbody:
                    pnb.rename(fout)
                    pnb.write()
                else: pnb.write(filename=fout,ftype='buildgal')
                print '    '+fout
            
    # Check numbers
    mmlio.verbose('Overall: {} found, {} required, {} remain'.format(Ntot0,Nreq0,Nrem0))
    idxtot=np.zeros([],dtype=int)
    rng=np.arange(Ntot0)
    for ityp in LIST_PTYPS:
        if ityp+'_n' not in par0: continue
        Nreq=par0[ityp+'_n']
        if Nreq==0: continue
        # Select particles by type
        if pNbody: itypslc=pnb.get_index(ptyp=ityp)
        else     : itypslc=pnb._get_family_slice(ityp)
        itypidx=rng[itypslc]
        # Fix masses
        if pNbody:
            Ntot=len(pnb.mass[itypslc])
            pnb.massarr[gadget.LIST_PTYPS.index(ityp)]=par0[ityp+'_mp']
            pnb.mass[itypslc]=par0[ityp+'_mp']
            Mtot=pnb.mass[itypslc].sum()
        else: 
            Ntot=len(pnb['mass'][itypslc])
            pnb['mass'][itypslc]=par0[ityp+'_mp']
            Mtot=pnb['mass'][itypslc].sum()
        mmlio.verbose('    {}: {} found, {} required, {} remain'.format(ityp,Ntot,Nreq,Nreq-Ntot))
        #mmlio.verbose('    {} ({}): Mi={}, Mf={}'.format(ityp,Ntot,Mtot,par0[ityp+'_mtot']))
        # Remove extra particles
        if Ntot>Nreq: np.random.shuffle(itypidx)
        idxtot=np.append(idxtot,itypidx[:min(Ntot,Nreq)])
    # Downsample to correct numbers
    if len(idxtot)<Ntot0: 
        if pNbody: pnb=pnb.selecti(idxtot)
        else     : pnb=pnb[idxtot]
    # Save
    if not mmlio.yorn('Write to {}?'.format(fout)): return pnb
    if pNbody:
        pnb.rename(fout)
        pnb.write()
    else: pnb.write(filename=fout,ftype='buildgal')
    print '    '+fout
    # Return control
    return

####################################################################################################################################
# METHOD TO MERGE BUILDGAL IC OUTPUT
def merge_ic(simstr,pNbody=None,**exkw):
    """
    Merges files into single IC file
    """
    from pysim import files
    # Ask about using pNbody
    if pNbody is None:
        pNbody=mmlio.yorn('Use pNbody to merge snapshots?')
    # Get file names
    fin_extfmt=files.get_extfmt(simstr.fdict['buildgal']['output']['ic'])
    fin_base=re.sub(fin_extfmt,'*',simstr.fdict['buildgal']['output']['ic'])
    fout=simstr.fdict['buildgal']['icgen']['stat_mg']
    # Count files
    fin_list=glob.glob(fin_base)
    nproc=len(fin_list)
    # Check number of processors
    if nproc==0: raise Exception('There are not any IC files to merge.')
    if nproc!=simstr['nprocmk']: raise Exception('{} files were expected but there are {}'.format(simstr['nprocmk'],nproc))
    # Load each files and merge them
    pnb=None
    for iproc in range(nproc):
        inb = simstr.loadsnap(fext=iproc,ftype='buildgal',pNbody=pNbody)
        if pNbody: inb.mass/=float(nproc)
        else     : inb['mass']/=float(nproc)
        if pnb==None: pnb=inb
        else        : pnb+=inb
        mmlio.verbose('[Nproc={:04d}] N (current)={}'.format(iproc,len(pnb)))
    mmlio.verbose('N (initial)={}'.format(len(pnb)))
    # Save
    if pNbody:
        pnb.rename(fout)
        pnb.write()
    else: pnb.write(filename=fout,ftype='buildgal')
    print '    '+fout
    # Return control
    return pnb

####################################################################################################################################
####################################################################################################################################
# METHODS TO CREATE BUILDGAL FILES

####################################################################################################################################
# METHOD TO CREATE BUILDGAL EXECUTABLE FILES
def mk_exec(simstr,**exkw):
    """
    Creates a local copy of executable files
    """
    # Select default files
    exdeft=simstr.finfo['buildgal']['exec']['default']+'_'+simstr['mkcomp']
    mfdeft=simstr.finfo['buildgal']['makefile']['default']
    excopy=simstr.finfo['buildgal']['exec']['path']
    mfcopy=simstr.finfo['buildgal']['makefile']['path']
    # Check existence
    if not os.path.isfile(exdeft): raise Exception('Executable does not exist: {}'.format(exdeft))
    if not os.path.isfile(mfdeft): raise Exception('Makefile does not exist: {}'.format(mfdeft))
    # Copy files to local directory
    mmlfiles.cpfiles(exdeft,excopy,overwrite=True)
    mmlfiles.cpfiles(mfdeft,mfcopy,overwrite=True)
    # Return
    return excopy,mfcopy

####################################################################################################################################
# METHOD TO CREATE SINGLE PARAMETER FILE
def mk_param(simstr=None,overwrite=None,parstat=None,askuser=False,multfact=2.,**param_kw):
    """
    Creates a buildgal parameter file
    """
    import pysim
    # Pars input
    if simstr==None: simstr=pysim.loadsim(runtag)
    if not isinstance(overwrite,bool):
        if askuser: overwrite=mmlio.yorn('Overwrite existing static parameter file?')
        else      : overwrite=mmlpars.mml_pars(overwrite,type=bool,default=False)
    # Get file names
    fpar_stat=simstr.finfo['buildgal']['stat_pm']['path']
    fpar_copy=simstr.fdict['buildgal']['input']['param']
    # If the static parameter file dosn't exist, create it
    if overwrite or not os.path.isfile(fpar_stat):
        # Get parameters
        if isinstance(parstat,dict): params=getparam('ask',inpar=parstat)
        else                       : params=model2param(simstr.infodict)
        # Save to file
        params['keylist']=getparam('list',inpar=params)
        rw_param('W',fpar_stat,params,overwrite=overwrite)
        # Check the parameter file
        print fpar_stat
        mmlio.yorn('Make changes to the above parameter file now!')
    # Check that particles evenly divisible across nproc
    if (simstr.ntot%simstr['nprocmk'])!=0:
        raise Exception('# of particles ({}) not evenly divided across # of processors ({}).'.format(simstr.ntot,simstr['nprocmk']))
    # Load parameters & scale them
    par0=rw_param('R',fpar_stat)
    par0=model2param(simstr.infodict,scale=True,par0=par0)
    par0['keylist']=getparam('list',inpar=par0)
    # Adjust particle number
    for ityp in LIST_PTYPS+['sat']: 
        if ityp+'_n' in par0: par0[ityp+'_n']=long(par0[ityp+'_n']*multfact/float(simstr['nprocmk']))
    # Loop over processers creating a file for each
    random.seed(par0['randseed'])
    for iproc in range(simstr['nprocmk']):
        ipar=copy.deepcopy(par0)
        ipar['randseed']=random.randint(0,20000000)
#        ipar['randseed']+=iproc
        rw_param('W',fpar_copy.format(iproc),ipar,overwrite=True)
    # Return output
    return fpar_copy
    
####################################################################################################################################
# METHODS TO CREATE SUBMISSION SCRIPTS
def mk_pbs(simstr,**kwargs): return mk_submit(simstr,'pbs',**kwargs)
def mk_sbatch(simstr,**kwargs): return mk_submit(simstr,'sbatch',**kwargs)
def mk_submit(simstr,ftype,param=None,verbose=False,memlim=500000,**kwargs):
    """
    Creates a submission script
    """
    # Pars input
    ftype=mmlpars.mml_pars(ftype,list=mmlio.DICT_SUBOPT.keys())
    # Get file names
    files_run=simstr.fdict_run
    subfile=simstr.fdict['buildgal']['input'][ftype]
    # Get run command info
    exefile=os.path.basename(files_run['buildgal']['input']['wrap'])
    mpicmd=mmlinfo.COMPDICT_MPI[simstr['mkcomp']]
    # Variables from simstr
    if 'timetest' in simstr['runtag']: name='bgal'+simstr['subtag']+'_'+simstr['runtag']
    else                             : name='bgal_'+simstr['runtag']+simstr['subtag']
    subdict={'jobname' :name,
             'nproc'   :int(simstr['nprocmk']),
             'ppn'     :int(mmlinfo.COMPDICT_PPN[simstr['mkcomp']]),
             'pmem'    :int(mmlsim2mcpu(simstr)),
             'twall'   :float(mmlsim2tcpu(simstr)),
             'outfile' :files_run['buildgal']['output']['runout'],
             'cmdlist' :['mkdir -p {}'.format(files_run['buildgal']['output']['dir'])]}
    # Check for too much memory
    memtot=subdict['pmem']*subdict['nproc']
    if memtot>memlim:
        nsub=int(memtot/memlim)
        subdict['nproc']/=nsub
        print nsub,subdict['nproc']
        for isub in range(nsub):
            print isub*subdict['nproc']
            subdict['cmdlist'].append('{} -comm none -np {} ./{} {}'.format(mpicmd,subdict['nproc'],exefile,isub*subdict['nproc']))
    else:
        subdict['cmdlist'].append('{} -comm none -np {} ./{} 0'.format(mpicmd,subdict['nproc'],exefile))
    # Exit command
    # subdict['cmdlist']+=['','wait',"echo 'finished mpiexec'"]
    # Write submission script
    mmlio.rw_submit('W',subfile,ftype,subdict,verbflag=verbose,overwrite=True)
    # Return output
    return subfile

####################################################################################################################################
def mk_wrapper(simstr,**exkw):
    """
    Creates a BUILDGAL wrapper file
    """
    import pysim
    rankfmt='${rank}'
    # Get file names
    wrpfile=simstr.fdict['buildgal']['input']['wrap']
    mmlfiles.mkdirs(os.path.dirname(wrpfile))
    flist_run=pysim.files.get_flist('buildgal',compid=simstr['mkcomp'],**simstr)
    #bgdir=os.path.join(simstr.fdict_run['buildgal']['input']['dir'],'bg'+rankfmt)
    bgdir='bg'+rankfmt
    # Create lines
    linelist=['rank=`expr $1 + $MPIEXEC_RANK`']
    rmlist=[]
    # Setup lines
    linelist+=['echo '+bgdir,'mkdir -p '+bgdir,'cd '+bgdir]
    for ftype in FILELIST_SETUP:
        for finfo in flist_run[ftype]:
            if finfo['alias']=='None' : continue
            else                      : ffile=finfo['alias'] 
            if finfo['size' ]=='large': ifile=re.sub(r'({.*})',rankfmt,finfo['path'])
            else                      : ifile=finfo['path']
            linelist.append('cp {} {}'.format(ifile,ffile))
            rmlist.append(ffile)
    # Execution line
    linelist.append('./{} > {}'.format(simstr.finfo_run['buildgal']['exec'   ]['alias'],
                                       simstr.finfo_run['buildgal']['procout']['alias']))
    # Recovery lines
    for ftype in FILELIST_OUTPT:
        for finfo in simstr.flist_run['buildgal'][ftype]:
            if finfo['alias']=='None' : continue
            else                      : ifile=finfo['alias']
            if finfo['size' ]=='large': ffile=re.sub(r'({.*})',rankfmt,finfo['path'])
            else                      : ffile=finfo['path']
            linelist.append('mv {} {}'.format(ifile,ffile))
            mmlfiles.mkdirs(os.path.dirname(ffile),host=simstr['mkcomp'])
    # Clean up lines
    linelist+=['rm {}'.format(rmfile) for rmfile in rmlist]
    linelist+=['cd '+simstr.fdict_run['buildgal']['input']['dir'],
               'wait',
               'rmdir '+bgdir,
               'exit 0']
    # Write file
    rw_wrapper('W',wrpfile,linelist=linelist,overwrite=True)
    # Convert file to executable
    subprocess.call(['chmod','755',wrpfile])
    # Return output
    return wrpfile

####################################################################################################################################
# METHOD TO CREATE UNITS FILE
def mk_units(simstr,**exkw):
    """
    Creates a units file
    """
    # Get file name
    unifile=simstr.fdict['buildgal']['input']['units']
    # Get parameters
    param=model2param(simstr.infodict)
    units=mmlconst.const_units('buildgal')
    physG=mmlconst.const_phys('cgs')['G']
    # Create dictionary
    unidict={'UnitLength_in_cm':units['L']*param['disk_rscl'],
             'UnitMass_in_g'   :units['M']*param['disk_mtot']}
    unidict['UnitVelocity_in_cm_per_s']=np.sqrt(physG*unidict['UnitMass_in_g']/unidict['UnitLength_in_cm'])
    #unidict['UnitMass_in_g']/=float(simstr['nprocmk'])
    # Save dictionary
    rw_units('W',unifile,unidict,overwrite=True)
    # Return dictionary
    return unifile

####################################################################################################################################
####################################################################################################################################
# METHODS TO READ/WRITE BUILDGAL FILES

####################################################################################################################################
# METHOD TO READ/WRITE SNAPSHOT FILES
def rw_snap(rwid,fname,fdict=None,overwrite=False):
    """
    Reads/writes snapshot files
    """
    names=('mass','x','y','z','vx','vy','vz','pot','tag')
    frmts=('f','f','f','f','f','f','f','f','i')
    # Read
    if   rwid=='R':
        fdata=np.loadtxt(fname,dtype={'names':names,'formats':frmts})
        fdict={nm:fdata[nm] for nm in names}
        fdict['pos']=np.vstack(tuple(fdict[nm] for nm in [ 'x', 'y', 'z'])).T
        fdict['vel']=np.vstack(tuple(fdict[nm] for nm in ['vx','vy','vz'])).T
        for nm in ['x','y','z','vx','vy','vz']: del fdict[nm]
        fdict['header']={'nbody':len(fdict['mass'])}
        return fdict
    # Write
    elif rwid=='W':
        if not overwrite and os.path.isfile(fname): return
        fdict[ 'x']=fdict['pos'][:,0] ; fdict[ 'y']=fdict['pos'][:,1] ; fdict[ 'z']=fdict['pos'][:,2]
        fdict['vx']=fdict['vel'][:,0] ; fdict['vy']=fdict['vel'][:,1] ; fdict['vz']=fdict['vel'][:,2]
        frmts=tuple('%'+ifmt for ifmt in frmts)
        farry=np.vstack(tuple(fdict[inam] for inam in names)).T
        np.savetxt(fname,farry,fmt=frmts)
        return
    # Error
    else: raise Exception('Invalid rwid: {}'.format(rwid))
    return

####################################################################################################################################
# METHOD TO READ/WRITE PARAMETER FILES
def rw_param(rwid,fname,fdict=None,overwrite=None):
    """
    Reads/writes parameter files
    """
    # Read
    if   rwid=='R':
        fdict=mmlio.rwdict('R',fname,style='fortran')
        return fdict
    # Write
    elif rwid=='W':
        mmlio.rwdict('W',fname,fdict,style='fortran',width=29,overwrite=overwrite)
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
# METHOD TO CONVERT MODEL PARAMETERS TO BUILDGAL PARAMETERS
def model2param(sim,delvir=None,randseed=None,minfracn=None,scale=None,par0={},rsoft=None,
                verbose=False):
    """
    Converts model parameters to buildgal parameters
    """
    # Pars input
    delvir=mmlpars.mml_pars(delvir,default=178.,type=float)
    randseed=mmlpars.mml_pars(randseed,default=123,type=int)
    minfracn=mmlpars.mml_pars(minfracn,default=0.01,type=float,min=0)
    scale=mmlpars.mml_pars(scale,default=False,type=bool)
    sim=mmlparam.parspar('pysim.simlist','galsim',inpar=sim)
    # Get galaxy model info
    mod=mmlparam.parspar('pysim.simlist','galmod',askuser=True,tagstr=sim['model'])
    rsoft=mmlpars.mml_pars(rsoft,default=mod['disk_rscl'],type=float,min=0.)
    # Populate buildgal parameter structure
    par=par0
    par['randseed']=randseed
    # Profiles
    par['halo_prof']=sim['haloprof']
    par['disk_prof']='HIYASHI'
    par['bulge_prof']='LH'
    par['gas_prof']=par['disk_prof']
    # Get scale masses
    for ityp in LIST_PTYPS: 
        if ityp!='halo': mod[ityp+'_prof']=par[ityp+'_prof']
        mod[ityp+'_mscl']=mmlprofiles.menc2ms(mod[ityp+'_prof'],mod[ityp+'_mtot'],
                                              r=mod[ityp+'_rmax'],rs=mod[ityp+'_rscl'],
                                              z=mod[ityp+'_zmax'],zs=mod[ityp+'_zscl'])
    if mod['mvir']==0:
        if mod['concen']==0:
            mod['mvir']=mmlprofiles.ms2mvir(mod['halo_prof'],mod['halo_mscl'],c=sim['concen'])
            mod['rvir']=sim['concen']*mmlprofiles.convert_rs(mod['halo_prof'],'NFW',mod['halo_rscl'])
        else:
            mod['mvir']=mmlprofiles.ms2mvir(mod['halo_prof'],mod['halo_mscl'],c=mod['concen'])
            mod['rvir']=mod['concen']*mmlprofiles.convert_rs(mod['halo_prof'],'NFW',mod['halo_rscl'])
    else:
        mod['rvir']=mmlprofiles.mvir2rvir(mod['mvir'],units='buildgal',deltavir=delvir)
    sim['rvir']=mmlprofiles.mvir2rvir(sim['mvir'],units='buildgal',deltavir=delvir)
#    print 'Mvir:'
#    print mod['mvir']
#    print sim['mvir']
#    print 'Rvir:'
#    print mod['rvir']
#    print sim['rvir']
    if scale:
        Mscl=1.0/mod['disk_mtot']
        Rscl=1.0/mod['disk_rscl']
    else:
        Mscl=sim['mvir']/mod['mvir']
        Rscl=sim['rvir']/mod['rvir']
    # Masses and radii
    minM=sim['mvir'] ; minMtyp='all'
    for ityp in LIST_PTYPS:
        if ityp=='halo':
            # Scale mass
            par['halo_mscl']=Mscl*mmlprofiles.convert_ms(mod['halo_prof'],par['halo_prof'],mod['halo_mscl'])
            par['halo_mtot']=Mscl*mmlprofiles.convert_ms(mod['halo_prof'],par['halo_prof'],mod['halo_mtot'])
            # Distances
            par['halo_rscl']=Rscl*mmlprofiles.convert_rs(mod['halo_prof'],par['halo_prof'],mod['halo_rscl'])
            par['halo_rmax']=Rscl*mod['halo_rmax'] 
            par['halo_zscl']=par['halo_rscl']
            par['halo_zmax']=par['halo_rmax']
        else:
            # Scale mass
            if ityp=='gas': par[ityp+'_mscl']=sim['fgas']*Mscl*mod['disk_mscl']
            else          : par[ityp+'_mscl']=Mscl*mod[ityp+'_mscl']
            if ityp=='gas': par[ityp+'_mtot']=sim['fgas']*Mscl*mod['disk_mtot']
            else          : par[ityp+'_mtot']=Mscl*mod[ityp+'_mtot']
            # Distances
            par[ityp+'_rscl']=Rscl*mod[ityp+'_rscl']
            par[ityp+'_zscl']=Rscl*mod[ityp+'_zscl']
            par[ityp+'_rmax']=Rscl*mod[ityp+'_rmax']
            par[ityp+'_zmax']=Rscl*mod[ityp+'_zmax']
        # Set include flag
        if par[ityp+'_mscl']==0: par[ityp]='N'
        else                   : par[ityp]='Y'
        par[ityp+'_selfgrav']='Y'
        # Minimum mass component
        if par[ityp]=='Y' and par[ityp+'_mtot']<minM:
            minM=par[ityp+'_mtot']
            minMtyp=ityp
    # Particle number & mass
    minN=minfracn*float(sim['ntot'])
    Mtot_equ=par['disk_mtot']+par['bulge_mtot']+par['gas_mtot'] # These particles are forced to have the same mass
    if minMtyp=='halo': m0=Mtot_equ/float(sim['ntot']-minN)
    else              : m0=par[minMtyp+'_mtot']/float(minN)
    Ntot_equ=Mtot_equ/m0
    for ityp in LIST_PTYPS:
        if ityp=='halo': par[ityp+'_n']=float(sim['ntot'])-Ntot_equ
        else           : par[ityp+'_n']=par[ityp+'_mtot']/m0
        if par[ityp+'_n']!=0: par[ityp+'_mp']=par[ityp+'_mtot']/par[ityp+'_n']
        else                : par[ityp+'_mp']=0.
    # Fractions
    Mpar=['mscl','mtot','n']
    Rpar=['rscl','rmax','zscl','zmax']
    for ityp in LIST_PTYPS:
        if ityp=='gas': continue
        if sim['f'+ityp]!=1:
            fMscl=sim['f'+ityp]
            fRscl=fMscl**(1./3.)
            for iMpar in Mpar: par['{}_{}'.format(ityp,iMpar)]*=fMscl
            for iRpar in Rpar: par['{}_{}'.format(ityp,iRpar)]*=fRscl
            if par[ityp+'_mtot']==0: par[ityp]='N'
    # Softenings (not needed? - yes they are)
    par['disk_eps']=0.0246689
    par['bulge_eps']=0.115122
    par['halo_eps']=0.0139791
    for ityp in LIST_PTYPS: 
        if ityp+'_eps' not in par: par[ityp+'_eps']=0.01
#        if par[ityp]=='N': par[ityp+'_eps']=0.0
#        else             : par[ityp+'_eps']=(par[ityp+'_mp']/mmlprofiles.main(par[ityp+'_prof'],method='rho',ms=par[ityp+'_mscl'],
#                                                                              rs=par[ityp+'_rscl'],r=rsoft,
#                                                                              zs=par[ityp+'_zscl'],z=0.))**(1./3.)
    par['outsoft']='N'
    # Disk parameters
    par['disk_solarR']=2.428571429*par['disk_rscl']
    if sim['toomreQ']==0: par['disk_solarQ']=1.1
    else                : par['disk_solarQ']=sim['toomreQ']
    # Gas parameters
    par['gas_rmin']=0.0
    # Bulge parameters (nonspher:nsimp)
    if par['bulge_rscl']==par['bulge_zscl']: par['bulge_nonspher']='N'
    else                                   : par['bulge_nonspher']='Y'
    if par['bulge_nonspher']=='Y':
        par['bulge_fracrev']=sim['bulgespin']
        if par['bulge_fracrev']!=0: par['bulge_rotate']='Y'
    # Halo parameters
    if   par['halo_prof']=='NF': 
        par['halo_NF_concen']=sim['concen']
        par['halo_NF_delvir']=delvir
    elif par['halo_prof']=='LH':
        par['halo_LH_rscl']=par['halo_rscl']
    elif par['halo_prof']=='IS':
        par['halo_IS_rcore']=par['halo_rscl']
        par['halo_IS_rtide']=par['halo_rmax']
    # Turn off satellite and twomod
    par['sat']='N'
    par['twomod']='N'
    # Secure number
    par['total_n']=0L
    for ityp in LIST_PTYPS: 
        par[ityp+'_n']=long(mmlmath.oom(par[ityp+'_n'],nsig=2))
        par['total_n']+=par[ityp+'_n']
    # Fill in missing values with user/default input
    par=getparam('ask',inpar=par)
    if verbose:
        print 'mscl={disk_mtot}'.format(**par)
        print 'rscl={disk_rscl}'.format(**par)
    # Return parameters
    return par

####################################################################################################################################
# METHOD TO ASK FOR BUILDGAL PARAMETERS
def getparam(method,inpar=None,**exkw):
    """
    Asks user for buildgal parameters
    """
    # Pars input
    method=mmlpars.mml_pars(method,list=['ask','pars','list'])
    inpar=mmlpars.mml_pars(inpar,default={},type=dict)
    if   method=='ask' : parkw=dict(inclextra=True,askuser=True ,load=False,save=False,**exkw)
    elif method=='pars': parkw=dict(inclextra=True,askuser=False,load=False,save=False,**exkw)
    # Base parameters
    if method=='list': flist=['randseed','outsoft']
    else             : inpar=mmlparam.parspar('pysim.files.buildgal','basepar',inpar=inpar,**parkw)
    # Disk parameters
    if inpar['disk']=='Y':
        if method=='list': flist=['disk_mtot']+flist+['disk_n','disk_zscl','disk_solarR','disk_solarQ','disk_eps','disk_zmax','disk_rmax','gas']
        else             : inpar=mmlparam.parspar('pysim.files.buildgal','diskpar',inpar=inpar,**parkw)
        # Gas parameters
        if inpar['gas']=='Y':
            if method=='list': flist+=['gas_n','gas_mtot','gas_T','gas_zscl','gas_zmax','gas_rmax','gas_rmin','gas_flatdisk']
            else             : inpar=mmlparam.parspar('pysim.files.buildgal','gaspar',inpar=inpar,**parkw)
            if inpar['gas_flatdisk']=='Y':
                if method=='list': flist+=['gas_flt_r','gas_flt_mint','gas_flt_mext','gas_flt_nint','gas_flt_next','gas_flt_mstar']
                else             : inpar=mmlparam.parspar('pysim.files.buildgal','gaspar_flatdisk',inpar=inpar,**parkw)
            if method=='list': flist+=['gas_selfgrav']
    # Bulge parameters
    if method=='list': flist+=['bulge']
    if inpar['bulge']=='Y':
        if method=='list': flist+=['bulge_mtot','bulge_rscl','bulge_selfgrav']
        else             : inpar=mmlparam.parspar('pysim.files.buildgal','bulgepar',inpar=inpar,**parkw)
        if inpar['bulge_selfgrav']=='Y':
            if method=='list': flist+=['bulge_n','bulge_rmax','bulge_eps','bulge_nonspher']
            else             : inpar=mmlparam.parspar('pysim.files.buildgal','bulgepar_selfgrav',inpar=inpar,**parkw)
            if inpar['bulge_nonspher']=='Y':
                if method=='list': flist+=['bulge_zscl','bulge_zmax','bulge_nsimp','bulge_rotate']
                else             : inpar=mmlparam.parspar('pysim.files.buildgal','bulgepar_nonspher',inpar=inpar,**parkw)
                if inpar['bulge_rotate']=='Y':
                    if method=='list': flist+=['bulge_fracrev']
                    else             : inpar=mmlparam.parspar('pysim.files.buildgal','bulgepar_rotate',inpar=inpar,**parkw)
    # Halo parameters
    if method=='list': flist+=['halo']
    if inpar['halo']=='Y':
        if method=='list': flist+=['halo_selfgrav','halo_rmax']
        else             : inpar=mmlparam.parspar('pysim.files.buildgal','halopar',inpar=inpar,**parkw)
        if inpar['halo_selfgrav']=='Y':
            if method=='list': flist+=['halo_n','halo_eps']
            else             : inpar=mmlparam.parspar('pysim.files.buildgal','halopar_selfgrav',inpar=inpar,**parkw)
        if method=='list':
            flist+=['halo_prof','halo_mtot']
            if   inpar['halo_prof']=='IS': flist+=['halo_IS_rcore','halo_IS_rtide']
            elif inpar['halo_prof']=='LH': flist+=['halo_LH_rscl']
            elif inpar['halo_prof']=='NF': flist+=['halo_NF_concen','halo_NF_delvir']
        else: inpar=mmlparam.parspar('pysim.files.buildgal','halopar_{}'.format(inpar['halo_prof']),inpar=inpar,**parkw)
    # Satellite parameters
    if method=='list': flist+=['sat']
    if inpar['sat']=='Y':
        if method=='list': flist+=['sat_mtot','sat_rscl','sat_x','sat_y','sat_z','sat_rmax','sat_selfgrav']
        else             : inpar=mmlparam.parspar('pysim.files.buildgal','satpar',inpar=inpar,**parkw)
        if inpar['sat_selfgrav']=='Y':
            if method=='list': flist+=['sat_n','sat_eps']
            else             : inpar=mmlparam.parspar('pysim.files.buildgal','satpar_selfgrav',inpar=inpar,**parkw)
        if method=='list': flist+=['sat_vx','sat_vy','sat_vz']
    # Two model parameters
    if method=='list': flist+=['twomod']
    if inpar['twomod']=='Y':
        if method=='list': flist+=['twomod_rperi','twomod_rinit','twomod_theta1','twomod_phi1','twomod_theta2','twomod_phi2']
        else             : inpar=mmlparam.parspar('pysim.files.buildgal','twomodpar',inpar=inpar,**parkw)
    # Return output
    if method=='list': return flist
    else             : return inpar

####################################################################################################################################
####################################################################################################################################
# METHODS TO CALCULATE BUILDGAL THINGS

####################################################################################################################################
# METHOD TO DETERMINE MEMORY REQUIRED TO RUN BUILDGAL
def mmlsim2mcpu(simstr):
    """
    Estimates the amount of memory required to run BUILDGAL
    """
    # Force reasonable values as a safety
    mcpu=500000.*(float(simstr.ntot)/float(1e7))
    # Return rounded value
    return mmlmath.oom(mcpu/simstr['nprocmk'],nsig=1,method='CEIL')

####################################################################################################################################
# METHOD TO DETERMINE TIME REQUIRED TO RUN BUILDGAL
def mmlsim2tcpu(simstr):
    """
    Estimates the amount of time required to run BUILDGAL
    """
    # Force reasonable values as a safety
    tcpu=48.
    return tcpu


####################################################################################################################################
####################################################################################################################################
# PARAMETER METHODS
def listpar(method,**exkw):
    """
    Returns a list of supported parameter lists
    """
    list_yorn=['Y','N']
    method=mmlpars.mml_pars(method,list=LIST_PARMETH)
    if   method=='basepar'          : fpar={'randseed'      :{'def':123 ,'form':int              },
                                            'outsoft'       :{'def':'N' ,'form':list_yorn        },
                                            'disk'          :{'def':'Y' ,'form':list_yorn        },
                                            'bulge'         :{'def':'Y' ,'form':list_yorn        },
                                            'halo'          :{'def':'Y' ,'form':list_yorn        },
                                            'sat'           :{'def':'N' ,'form':list_yorn        },
                                            'twomod'        :{'def':'N' ,'form':list_yorn        }}
    elif method=='diskpar'          : fpar={'disk_mtot'     :{'def':1.0 ,'form':float            },
                                            'disk_n'        :{'def':0L  ,'form':long             },
                                            'disk_rscl'     :{'def':1.0 ,'form':float            },
                                            'disk_zscl'     :{'def':0.1 ,'form':float            },
                                            'disk_eps'      :{'def':0.01,'form':float            },
                                            'disk_rmax'     :{'def':10.0,'form':float            },
                                            'disk_zmax'     :{'def':1.0 ,'form':float            },
                                            'disk_solarR'   :{'def':2.428571429,'form':float            },
                                            'disk_solarQ'   :{'def':1.1 ,'form':float            },
                                            'gas'           :{'def':'N' ,'form':list_yorn        }}
    elif method=='gaspar'           : fpar={'gas_n'         :{'def':0L  ,'form':long             },
                                            'gas_mtot'      :{'def':0.0 ,'form':float            },
                                            'gas_T'         :{'def':1.0e4,'form':float            },
                                            'gas_zscl'      :{'def':0.1 ,'form':float            },
                                            'gas_rmax'      :{'def':10.0,'form':float            },
                                            'gas_rmin'      :{'def':0.0 ,'form':float            },
                                            'gas_zmax'      :{'def':1.0 ,'form':float            },
                                            'gas_flatdisk'  :{'def':'N' ,'form':list_yorn        },
                                            'gas_selfgrav'  :{'def':'N' ,'form':list_yorn        }}
    elif method=='gaspar_flatdisk'  : fpar={'gas_flt_r'     :{'def':0.0 ,'form':float            },
                                            'gas_flt_mint'  :{'def':0.0 ,'form':float            },
                                            'gas_flt_mext'  :{'def':0.0 ,'form':float            },
                                            'gas_flt_nint'  :{'def':0L  ,'form':long             },
                                            'gas_flt_next'  :{'def':0L  ,'form':long             },
                                            'gas_flt_mstar' :{'def':0.0 ,'form':float            }}
    elif method=='bulgepar'         : fpar={'bulge_mtot'    :{'def':0.0 ,'form':float            },
                                            'bulge_rscl'    :{'def':0.0 ,'form':float            },
                                            'bulge_selfgrav':{'def':'N' ,'form':list_yorn        },
                                            'bulge_nonspher':{'def':'N' ,'form':list_yorn        }}
    elif method=='bulgepar_selfgrav': fpar={'bulge_n'       :{'def':0L  ,'form':long             },
                                            'bulge_rmax'    :{'def':0.0 ,'form':float            },
                                            'bulge_eps'     :{'def':0.0 ,'form':float            }}
    elif method=='bulgepar_nonspher': fpar={'bulge_zscl'    :{'def':0.0 ,'form':float            },
                                            'bulge_zmax'    :{'def':0.0 ,'form':float            },
                                            'bulge_nsimp'   :{'def':0L  ,'form':long             },
                                            'bulge_rotate'  :{'def':'N' ,'form':list_yorn        }}
    elif method=='bulgepar_rotate'  : fpar={'bulge_fracrev' :{'def':0.0 ,'form':float            }}
    elif method=='halopar'          : fpar={'halo_mtot'     :{'def':0.0 ,'form':float            },
                                            'halo_rmax'     :{'def':0.0 ,'form':float            },
                                            'halo_selfgrav' :{'def':'N' ,'form':list_yorn        },
                                            'halo_prof'     :{'def':'NF','form':LIST_HALOPROFS   }}
    elif method=='halopar_selfgrav' : fpar={'halo_n'        :{'def':0L  ,'form':long             },
                                            'halo_eps'      :{'def':0.0 ,'form':float            }}
    elif method=='halopar_NF'       : fpar={'halo_NF_concen':{'def':0.0 ,'form':float            },
                                            'halo_NF_delvir':{'def':0.0 ,'form':float            }}
    elif method=='halopar_IS'       : fpar={'halo_IS_rcore' :{'def':0.0 ,'form':float            },
                                            'halo_IS_rtide' :{'def':0.0 ,'form':float            }}
    elif method=='halopar_LH'       : fpar={'halo_LH_rscl'  :{'def':0.0 ,'form':float            }}
    elif method=='satpar'           : fpar={'sat_mtot'      :{'def':0.0 ,'form':float            },
                                            'sat_rscl'      :{'def':0.0 ,'form':float            },
                                            'sat_x'         :{'def':0.0 ,'form':float            },
                                            'sat_y'         :{'def':0.0 ,'form':float            },
                                            'sat_z'         :{'def':0.0 ,'form':float            },
                                            'sat_vx'        :{'def':0.0 ,'form':float            },
                                            'sat_vy'        :{'def':0.0 ,'form':float            },
                                            'sat_vz'        :{'def':0.0 ,'form':float            },
                                            'sat_rmax'      :{'def':0.0 ,'form':float            },
                                            'sat_selfgrav'  :{'def':'N' ,'form':list_yorn        }}
    elif method=='satpar_selfgrav'  : fpar={'sat_n'         :{'def':0L  ,'form':long             },
                                            'sat_eps'       :{'def':0.0 ,'form':float            }}
    elif method=='twomodpar'        : fpar={'twomod_rinit'  :{'def':0.0 ,'form':float            },
                                            'twomod_rperi'  :{'def':0.0 ,'form':float            },
                                            'twomod_phi1'   :{'def':0.0 ,'form':float            },
                                            'twomod_theta1' :{'def':0.0 ,'form':float            },
                                            'twomod_phi2'   :{'def':0.0 ,'form':float            },
                                            'twomod_theta2' :{'def':0.0 ,'form':float            }}
    else: raise Exception('Invalid parameter method: {}'.format(method))
    return fpar
def par2tag(method,inpar=None):
    """
    Returns a tag given the set of parameters
    """
    method=mmlpars.mml_pars(method,list=LIST_PARMETH)
    if method in LIST_PARMETH:
        import datetime
        ftag=datetime.date.isoformat(datetime.date.today())
    else:
        outpar=mmlparam.parspar('pysim.files.buildgal',method,inpar=inpar)
        raise Exception('Add method for other tags')
    return ftag

####################################################################################################################################
####################################################################################################################################
# PROVIDE COMMAND LINE ACCESS
if __name__ == '__main__': main()
