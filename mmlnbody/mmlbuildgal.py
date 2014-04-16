#!/usr/bin/python
####################################################################################################################################
#
# MEAGAN LANG'S BUILDGAL METHODS
#
####################################################################################################################################
import sys,os,shutil,glob,copy,pprint,scipy,math,subprocess,random
import numpy as np
import matplotlib as mplib
import pNbody
LIST_METHODS=['makefiles','submit','retrieve','merge','clean']
LIST_PARMETH=['basepar','diskpar','gaspar','gaspar_flatdisk','bulgepar','bulgepar_selfgrav','bulgepar_nonspher','bulgepar_rotate',
              'halopar','halopar_selfgrav','halopar_NF','halopar_IS','halopar_LH','satpar','satpar_selfgrav','twomodpar']
#LIST_GALMODELS=['MW','M81']
LIST_HALOPROFS=['LH','NF','IS']
LIST_PTYPS=['halo','disk','bulge','gas']
from mmlutils import *
from mmlastro import mmlprofiles,mmlcosmo,mmlconst
import main as mmlnbody
import simlist,simcalc,simplot,simfile,mmlgadget,mmlscf,mmlgalfit,icgen

def main():
    """
    Provides command line access to BUILDGAL methods
    """
    out = mmlnbody.walk(mtype='buildgal')
    return out

####################################################################################################################################
####################################################################################################################################
# METHODS FOR GETTING LISTS
def list_fpar(method):
    """
    Returns a list of supported parameter lists
    """
    method=mmlpars.mml_pars(method,list=LIST_PARMETH)
    if   method=='basepar'          : fpar={'randseed'      :int              ,
                                            'outsoft'       :simlist.LIST_YORN,
                                            'disk'          :simlist.LIST_YORN,
                                            'bulge'         :simlist.LIST_YORN,
                                            'halo'          :simlist.LIST_YORN,
                                            'sat'           :simlist.LIST_YORN,
                                            'twomod'        :simlist.LIST_YORN}
    elif method=='diskpar'          : fpar={'disk_mtot'     :float            ,
                                            'disk_n'        :long             ,
                                            'disk_rscl'     :float            ,
                                            'disk_zscl'     :float            ,
                                            'disk_eps'      :float            ,
                                            'disk_rmax'     :float            ,
                                            'disk_zmax'     :float            ,
                                            'disk_solarR'   :float            ,
                                            'disk_solarQ'   :float            ,
                                            'gas'           :simlist.LIST_YORN}
    elif method=='gaspar'           : fpar={'gas_n'         :long             ,
                                            'gas_mtot'      :float            ,
                                            'gas_T'         :float            ,
                                            'gas_zscl'      :float            ,
                                            'gas_rmax'      :float            ,
                                            'gas_rmin'      :float            ,
                                            'gas_zmax'      :float            ,
                                            'gas_flatdisk'  :simlist.LIST_YORN}
    elif method=='gaspar_flatdisk'  : fpar={'gas_flt_r'     :float            ,
                                            'gas_flt_mint'  :float            ,
                                            'gas_flt_mext'  :float            ,
                                            'gas_flt_nint'  :long             ,
                                            'gas_flt_next'  :long             ,
                                            'gas_flt_mstar' :float            }
    elif method=='bulgepar'         : fpar={'bulge_mtot'    :float            ,
                                            'bulge_rscl'    :float            ,
                                            'bulge_selfgrav':simlist.LIST_YORN,
                                            'bulge_nonspher':simlist.LIST_YORN}
    elif method=='bulgepar_selfgrav': fpar={'bulge_n'       :long             ,
                                            'bulge_rmax'    :float            ,
                                            'bulge_eps'     :float            }
    elif method=='bulgepar_nonspher': fpar={'bulge_zscl'    :float            ,
                                            'bulge_zmax'    :float            ,
                                            'bulge_nsimp'   :long             ,
                                            'bulge_rotate'  :simlist.LIST_YORN}
    elif method=='bulgepar_rotate'  : fpar={'bulge_fracrev' :float            }
    elif method=='halopar'          : fpar={'halo_mtot'     :float            ,
                                            'halo_rmax'     :float            ,
                                            'halo_selfgrav' :simlist.LIST_YORN,
                                            'halo_prof'     :LIST_HALOPROFS   }
    elif method=='halopar_selfgrav' : fpar={'halo_n'        :long             ,
                                            'halo_eps'      :float            }
    elif method=='halopar_NF'       : fpar={'halo_NF_concen':float            ,
                                            'halo_NF_delvir':float            }
    elif method=='halopar_IS'       : fpar={'halo_IS_rcore' :float            ,
                                            'halo_IS_rtide' :float            }
    elif method=='halopar_LH'       : fpar={'halo_LH_rscl'  :float            }
    elif method=='satpar'           : fpar={'sat_mtot'      :float            ,
                                            'sat_rscl'      :float            ,
                                            'sat_x'         :float            ,
                                            'sat_y'         :float            ,
                                            'sat_z'         :float            ,
                                            'sat_vx'        :float            ,
                                            'sat_vy'        :float            ,
                                            'sat_vz'        :float            ,
                                            'sat_rmax'      :float            ,
                                            'sat_selfgrav'  :simlist.LIST_YORN}
    elif method=='satpar_selfgrav'  : fpar={'sat_n'         :long             ,
                                            'sat_eps'       :float            }
    elif method=='twomodpar'        : fpar={'twomod_rinit'  :float            ,
                                            'twomod_rperi'  :float            ,
                                            'twomod_phi1'   :float            ,
                                            'twomod_theta1' :float            ,
                                            'twomod_phi2'   :float            ,
                                            'twomod_theta2' :float            }
    else: raise Exception('Invalid parameter method: {}'.format(method))
    return fpar
def list_fparDEF(method):
    """
    Returns a dictionary with default parameter values
    """
    method=mmlpars.mml_pars(method,list=LIST_PARMETH)
    if   method=='basepar'          : fpar={'randseed'      :123 ,
                                            'outsoft'       :'N' ,
                                            'disk'          :'Y' ,
                                            'bulge'         :'Y' ,
                                            'halo'          :'Y' ,
                                            'sat'           :'N' ,
                                            'twomod'        :'N' }
    elif method=='diskpar'          : fpar={'disk_mtot'     :1.0 ,
                                            'disk_n'        :0L  ,
                                            'disk_rscl'     :1.0 ,
                                            'disk_zscl'     :0.1 ,
                                            'disk_eps'      :0.01,
                                            'disk_rmax'     :10.0,
                                            'disk_zmax'     :1.0 ,
                                            'disk_solarR'   :2.428571429,
                                            'disk_solarQ'   :1.1 ,
                                            'gas'           :'N' }
    elif method=='gaspar'           : fpar={'gas_n'         :0L  ,
                                            'gas_mtot'      :0.0 ,
                                            'gas_T'         :1.0e4,
                                            'gas_zscl'      :0.1 ,
                                            'gas_rmax'      :10.0,
                                            'gas_rmin'      :0.0 ,
                                            'gas_zmax'      :1.0 ,
                                            'gas_flatdisk'  :'N' ,
                                            'gas_selfgrav'  :'N' }
    elif method=='gaspar_flatdisk'  : fpar={'gas_flt_r'     :0.0 ,
                                            'gas_flt_mint'  :0.0 ,
                                            'gas_flt_mext'  :0.0 ,
                                            'gas_flt_nint'  :0L  ,
                                            'gas_flt_next'  :0L  ,
                                            'gas_flt_mstar' :0.0 }
    elif method=='bulgepar'         : fpar={'bulge_mtot'    :0.0 ,
                                            'bulge_rscl'    :0.0 ,
                                            'bulge_selfgrav':'N' ,
                                            'bulge_nonspher':'N' }
    elif method=='bulgepar_selfgrav': fpar={'bulge_n'       :0L  ,
                                            'bulge_rmax'    :0.0 ,
                                            'bulge_eps'     :0.0 }
    elif method=='bulgepar_nonspher': fpar={'bulge_zscl'    :0.0 ,
                                            'bulge_zmax'    :0.0 ,
                                            'bulge_nsimp'   :0L  ,
                                            'bulge_rotate'  :'N' }
    elif method=='bulgepar_rotate'  : fpar={'bulge_fracrev' :0.0 }
    elif method=='halopar'          : fpar={'halo_mtot'     :0.0 ,
                                            'halo_rmax'     :0.0 ,
                                            'halo_selfgrav' :'N' ,
                                            'halo_prof'     :'NF'}
    elif method=='halopar_selfgrav' : fpar={'halo_n'        :0L  ,
                                            'halo_eps'      :0.0 }
    elif method=='halopar_NF'       : fpar={'halo_NF_concen':0.0 ,
                                            'halo_NF_delvir':0.0 }
    elif method=='halopar_IS'       : fpar={'halo_IS_rcore' :0.0 ,
                                            'halo_IS_rtide' :0.0 }
    elif method=='halopar_LH'       : fpar={'halo_LH_rscl'  :0.0 }
    elif method=='satpar'           : fpar={'sat_mtot'      :0.0 ,
                                            'sat_rscl'      :0.0 ,
                                            'sat_x'         :0.0 ,
                                            'sat_y'         :0.0 ,
                                            'sat_z'         :0.0 ,
                                            'sat_vx'        :0.0 ,
                                            'sat_vy'        :0.0 ,
                                            'sat_vz'        :0.0 ,
                                            'sat_rmax'      :0.0 ,
                                            'sat_selfgrav'  :'N' }
    elif method=='satpar_selfgrav'  : fpar={'sat_n'         :0L  ,
                                            'sat_eps'       :0.0 }
    elif method=='twomodpar'        : fpar={'twomod_rinit'  :0.0 ,
                                            'twomod_rperi'  :0.0 ,
                                            'twomod_phi1'   :0.0 ,
                                            'twomod_theta1' :0.0 ,
                                            'twomod_phi2'   :0.0 ,
                                            'twomod_theta2' :0.0 }
    else: raise Exception('Invalid parameter method: {}'.format(method))
    return fpar
def fpar2tag(method,inpar=None):
    """
    Returns a tag given the set of parameters
    """
    method=mmlpars.mml_pars(method,list=LIST_PARMETH)
    if method in LIST_PARMETH:
        import datetime
        ftag=datetime.date.isoformat(datetime.date.today())
    else:
        outpar=simfile.parsfpar('buildgal',method,fpar=inpar)
        raise Exception('Add method for other tags')
    return ftag
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
    fdict={'clean' :['longout'],
           'input' :['shortin' ,'longin' ],
           'output':['shortout','longout']}
    return fdict

####################################################################################################################################
####################################################################################################################################
# HIGH LEVEL METHODS REQUIRING MMLSIM

####################################################################################################################################
# METHOD FOR RUNNING DIFFERENT BUILDGAL METHODS
def run(simstr,method,**method_kw):
    """
    Provides interface for running different BUILDGAL operations
    """
    # Set constants
    methLIST=LIST_METHODS
    # Pars input
    method=mmlpars.mml_pars(method,list=methLIST)
    # Initialize default output
    out=None
    # Proceed based on method
    if   method=='makefiles': mkfiles(simstr,**method_kw)
    elif method=='submit'   : subrun(simstr,**method_kw)
    elif method=='retrieve' : endrun(simstr,**method_kw)
    elif method=='merge'    : mergefiles(simstr,**method_kw)
    elif method=='clean'    :
        method_kw['ftype']='clean'
        method_kw['memcomp']='accre'
        simfile.rmfiles('buildgal',simstr=simstr,**method_kw)
    else:
        raise Exception('Option for method {} needs to be added to mmlbuildgal.run.'.format(method))
    # Return output
    return out

####################################################################################################################################
# METHOD PERFORM ACTIONS TO SETUP AND RUN BUILDGAL ON ACCRE
def subrun(simstr,**extra_kw):
    """
    Sets up and runs BUILDGAL on ACCRE
    """
    # Get file names
    files_accre=simstr.mkfiledict(memcomp='accre',checkaccess=True)['buildgal']
    # Update input files
    if mmlio.yorn('Create new input files?'): mkfiles(simstr,**extra_kw)
    # Move input files to ACCRE
    if mmlio.yorn('Move new input files to ACCRE?'): putfiles(simstr,**extra_kw)
    # Make directories
    mmlfiles.mkdirs(files_accre['output']['dir'])
    # Run BUILDGAL
    if os.environ['HOSTNAME'].lower()=='bender':
        if mmlio.yorn('Run BUILDGAL?'):
            pbsfile=files_accre['input']['pbs']
            schfile=files_accre['input']['sched']
            simstr.subpbs('buildgal',pbsfile=pbsfile,schfile=schfile)
    return

####################################################################################################################################
# METHOD TO PERFORM ACTIONS TO RETRIEVE AND CLEAN-UP A RUN
def endrun(simstr,**exkw):
    """
    Retrieves and cleans up results from a BUILDGAL run
    """
    # Retrieve results
    if mmlio.yorn('Recover BUILDGAL results from ACCRE?'): getfiles(simstr)
    # Clean up results
    if mmlio.yorn('Clean up BUILDGAL results on ACCRE?' ): simfile.rmfiles('buildgal',ftype='clean',memcomp='accre',simstr=simstr)
    # Merge files
    if mmlio.yorn('Merge output files into single IC?'  ): mergefiles(simstr)
    return

####################################################################################################################################
# METHOD TO MOVE INPUT FILES TO ACCRE
def putfiles(simstr,**exkw):
    """
    Moves files to ACCRE in preparation to running BUILDGAL
    """
    # Move files
    if mmlio.yorn('Move BUILDGAL input files to ACCRE?'     ): simfile.mvfiles_accre(simstr,'to','buildgal','input' )
#    if mmlio.yorn('Move BUILDGAL old output files to ACCRE?'): simfile.mvfiles_accre(simstr,'to','buildgal','output')
    return

####################################################################################################################################
# METHOD TO GET FILES FROM ACCRE
def getfiles(simstr,**exkw):
    """
    Retrieves BUILDGAL files from ACCRE
    """
    # Move files
#    if mmlio.yorn('Recover BUILDGAL input files from ACCRE?' ): simfile.mvfiles_accre(simstr,'from','buildgal','input' )
    if mmlio.yorn('Recover BUILDGAL output files from ACCRE?'): simfile.mvfiles_accre(simstr,'from','buildgal','output')
    return

####################################################################################################################################
def mergefiles(simstr,**exkw):
    """
    Merges files into single IC file
    """
    # Get file names
    fin_base=simstr.fdict['buildgal']['output']['icbase']
    if '*' not in fin_base: fin_base+='*'
    fout=simstr.fdict['icgen']['bgic']
    # Count files
    fin_list=glob.glob(fin_base)
    nproc=len(fin_list)
    if nproc==0: raise Exception('There are not any IC files to merge.')
    if nproc!=simstr['nprocbg']: raise Exception('{} files were expected but there are {}'.format(simstr['nprocbg'],nproc))
    # Load parameters
    par0=mmlio.rwdict('R',simstr.fdict['icgen']['bgpm'],style='fortran')
    par0=model2param(simstr.infodict,scale=True,par0=par0)
    # Load each files and merge them
    pnb=None
    for iproc in range(nproc):
        inb=simstr.get_pnbody(fname=fin_base.replace('*',str(iproc)),ftype='buildgal',noidgal2=True)
        inb.mass/=float(nproc)
        if pnb==None: pnb=inb
        else        : pnb+=inb
    mmlio.verbose('N (initial)={}'.format(pnb.nbody))
    # Remove repeats
    sftmax=0.
    for ityp in LIST_PTYPS:
        if ityp+'_n' not in par0: continue
        sftmax=max(sftmax,par0[ityp+'_eps'])
    sftmax=mmlio.askquest('What should the minimum separation be?',dtype='float',default=sftmax)
    idx_rep=pnb.index_norepeats(cut=sftmax)
    pnb_rep=pnb.selecti(idx_rep)
    mmlio.verbose('N (no reps)={}'.format(pnb_rep.nbody))
    # Remove extra particles
    idx_ext=np.array([],dtype=int)
    for ityp in LIST_PTYPS:
        if ityp+'_n' not in par0: continue
        itypidx=np.arange(pnb_rep.nbody,dtype=int)[pnb_rep.get_index(ptyp=ityp)]
        iN=len(itypidx) ; fN=par0[ityp+'_n']
        if fN==0: continue
        mmlio.verbose('    {}: Ni={}, Nf={}'.format(ityp,iN,fN))
        # Add random sample of particles
        if fN>iN: raise Exception('There are too few {} particles.'.format(ityp))
        np.random.shuffle(itypidx)
        idx_ext=np.append(idx_ext,itypidx[:fN])
        # Fix masses
        pnb_rep.massarr[simlist.LIST_PTYPBASE.index(ityp)]=par0[ityp+'_mp']
        pnb_rep.mass[itypidx]=par0[ityp+'_mp']
    pnb_ext=pnb_rep.selecti(np.sort(idx_ext))
    mmlio.verbose('N (no extr)={}'.format(pnb_ext.nbody))
    # Fix masses
    for ityp in LIST_PTYPS:
        itypidx=pnb_ext.get_index(ptyp=ityp)
        mmlio.verbose('    {}: Mi={}, Mf={}'.format(ityp,pnb_ext.mass[itypidx].sum(),par0[ityp+'_mtot']))
    # Save
    pnb_ext.rename(fout)
    pnb_ext.write()
    # Return control
    return

####################################################################################################################################
####################################################################################################################################
# FILE CLASSES AND METHODS

####################################################################################################################################
# METHOD FOR PARSING BUILDGAL FILE OPTIONS
def pars_fileopt(fopt):
    """
    Returns a parsed dictionary of BUILDGAL file options
        execpath: Str absolute path to location of compiled BUILDGAL executable
        parmpath: Str absolute path to location of default BUILDGAL parameter file
        mfilpath: Str absolute path to location of default BUILDGAL compile file
    """
    # Define values
    form={
        'execpath': mmlpars.parsdict(default='/home/langmm/bin/buildgalquiet/buildgal.stall'   ,type=str),
        'parmpath': mmlpars.parsdict(default='/home/langmm/bin/buildgalquiet/buildgal.in'      ,type=str),
        'mfilpath': mmlpars.parsdict(default='/home/langmm/bin/buildgalquiet/buildgal.compfile',type=str)
        }
    # Pars input
    fopt=mmlpars.mml_formpars(fopt,form)
    # Return parsed input
    return fopt

####################################################################################################################################
# METHOD TO CREATE DICTIONARY OF BUILDGAL FILES
def files(fdict):
    """
    Returns a dictionary of BUILDGAL files & directories
    """
    # Pars input
    fdict=mmlpars.mml_pars(fdict,type=dict)
    rundir=fdict['rundir']
    pfix=fdict['pfix']
    pfix_shrt=fdict['pfix_shrt']
    options=pars_fileopt({})
    # Initialize output dictionary
    files=dict(fdict['buildgal'],exdeft=options['execpath'],mfdeft=options['mfilpath'],pmdeft=options['parmpath'])
    # Add input files
    infiles={}
    infiles['dir']=os.path.join(files['dir'],'input')
    fkeys_in=['exec','makefile','sub','wrap','pbs','sched','units']
    flist_in=['buildgal','bgMakfile','bgsub','bgwrap','bgpbs','bgschednum','bgunits']
    infiles=mmlfiles.filedict(infiles['dir'],filekeys=fkeys_in,filelist=flist_in,pfix=pfix,fdict=infiles)
    # Long input keys
    longin_keys=['par']
    longin_list=['bgpar']
    longin_accr=['buildgal.in']
    # Long input directories
    if fdict['memcomp']=='accre':
        long_dir=files['dir']
        long_pfx=''
        dkeys_in_long=[ikey+'dir' for ikey in longin_keys]
        dlist_in_long=['bg*' for ikey in longin_keys]
    else:
        long_dir=infiles['dir']
        long_pfx=pfix
        dkeys_in_long=[ikey+'dir' for ikey in longin_keys]
        dlist_in_long=longin_list
    infiles['subdirs']=dkeys_in_long
    infiles=mmlfiles.filedict(long_dir,dirkeys=dkeys_in_long,dirlist=dlist_in_long,pfix=long_pfx,fdict=infiles)
    # Long input files
    if fdict['memcomp']=='accre':
        fkeys_in_long=[ikey+'base' for ikey in longin_keys]
        flist_in_long=longin_accr
    else:
        fkeys_in_long=[ikey+'base' for ikey in longin_keys]
        flist_in_long=longin_list
    for ifin in range(len(dkeys_in_long)):
        infiles=mmlfiles.filedict(infiles[dkeys_in_long[ifin]],filekeys=[fkeys_in_long[ifin]],filelist=[flist_in_long[ifin]],pfix=long_pfx,fdict=infiles)
    files['input']=infiles
    # Add output files
    outfiles={}
    outfiles['dir']=os.path.join(files['dir'],'output')
    fkeys_out=['runout']
    flist_out=['bgrunout']
    outfiles=mmlfiles.filedict(outfiles['dir'],filekeys=fkeys_out,filelist=flist_out,pfix=pfix,fdict=outfiles)
    # Long output keywords
    longout_keys=['out','ic','DDAT','STAT','TREEBI','rotc','sigr','sigt','sigz','toom','scflog','ticklog']
    longout_list=['bgout','ic_buildgal','DDAT','STAT','TREEBI','rotc','sigr','sigt','sigz','toom','scflog','ticklog']
    longout_accr=['bgout','ic_buildgal.out','DDATA','STATS','TREEBI','rotcurve.dat','sigratio.dat','sigt.dat','sigzero.dat','toomreqave.dat','scflog','tickerlog']
    # Long output directories
    if fdict['memcomp']=='accre':
        long_dir=files['dir']
        dkeys_out_long=[ikey+'dir' for ikey in longout_keys]
        dlist_out_long=['bg*' for ikey in longout_keys]
    else:
        long_dir=outfiles['dir']
        dkeys_out_long=[ikey+'dir' for ikey in longout_keys]
        dlist_out_long=longout_list
    outfiles['subdirs']=dkeys_out_long
    outfiles=mmlfiles.filedict(long_dir,dirkeys=dkeys_out_long,dirlist=dlist_out_long,pfix=long_pfx,fdict=outfiles)
    # Long output files
    if fdict['memcomp']=='accre':
        fkeys_out_long=[ikey+'base' for ikey in longout_keys]
        flist_out_long=longout_accr
    else:
        fkeys_out_long=[ikey+'base' for ikey in longout_keys]
        flist_out_long=longout_list
    for ifout in range(len(dkeys_out_long)):
        outfiles=mmlfiles.filedict(outfiles[dkeys_out_long[ifout]],
                                   filekeys=[fkeys_out_long[ifout]],filelist=[flist_out_long[ifout]],pfix=long_pfx,fdict=outfiles)
    files['output']=outfiles
    # Return dictionary
    fdict['buildgal']=files
    return fdict

####################################################################################################################################
# METHOD TO RETURN LIST OF RELAVENT FILES
def get_filelist(fdict,keylist):
    """
    Returns a list of relevant BUILDGAL files
    """
    # Pars input
    fdict=mmlpars.mml_pars(fdict,type=dict)
    keylist=mmlpars.mml_pars(keylist,type=list)
    # Add files
    filelist={}
    if 'default'   in keylist:
        fkeylist=['exdeft','mfdeft','pmdeft']
        filelist['default' ]=[fdict[ifkey] for ifkey in fkeylist]
    if 'shortin'   in keylist:
        shrtList=['exec','makefile','sub','wrap','pbs','sched','units']
        filelist['shortin' ]=[fdict['input'][ishrt] for ishrt in shrtList]
    if 'longin'    in keylist:
        longList=['par']
        filelist['longin'  ]=[]
        for ilong in longList:
            ifile=fdict['input'][ilong+'base']
            if '*' not in ifile: ifile+='*'
            filelist['longin'  ].append(ifile)
    if 'shortout'  in keylist:
        shrtList=['runout']
        filelist['shortout']=[fdict['output'][ishrt] for ishrt in shrtList]
    if 'longout'   in keylist:
        longList=['out','ic','DDAT','STAT','TREEBI','rotc','sigr','sigt','sigz','toom','scflog','ticklog']
        filelist['longout' ]=[]
        for ilong in longList:
            ifile=fdict['output'][ilong+'base']
            if '*' not in ifile: ifile+='*'
            filelist['longout' ].append(ifile)
    # Return filelist
    return filelist

####################################################################################################################################
####################################################################################################################################
# METHODS TO CREATE BUILDGAL FILES

####################################################################################################################################
# METHOD TO CREATE BUILDGAL FILES
def mkfiles(simstr,nproc=None,tcpu=None,ppn=None,verbose=False,**exkw):
    """
    Creates BUILDGAL files
    """
    if verbose: mmlio.verbose('Creating files for a BUILDGAL run locally...',border=True)
    # Executable file
    if mmlio.yorn('Create new local copy of BUILDGAL executable?'):
        if not os.path.isfile(simstr.fdict['buildgal']['exdeft']):
            raise Exception('Default executable does not exist. Create it first: {}'.format(simstr.fdict['buildgal']['exdeft']))
        if not os.path.isfile(simstr.fdict['buildgal']['mfdeft']):
            raise Exception('Default Makefile does not exist. Create it first: {}'.format(simstr.fdict['buildgal']['mfdeft']))
        mmlfiles.mkdirs(simstr.fdict['buildgal']['input']['dir'])
        mmlfiles.cpfiles(simstr.fdict['buildgal']['exdeft'],simstr.fdict['buildgal']['input']['exec'],overwrite=True)
        mmlfiles.cpfiles(simstr.fdict['buildgal']['mfdeft'],simstr.fdict['buildgal']['input']['makefile'],overwrite=True)
    # Submission files
    if mmlio.yorn('Create BUILDGAL submission files?'):
        owparam=mmlio.yorn('Overwrite existing files?')
        pbsdict=mkpbs(simstr,tcpu=tcpu,ppn=ppn,overwrite=owparam,verbose=verbose)
        wrpdict=mkwrap(simstr,overwrite=owparam,verbose=verbose)
    # Parameter files
    if mmlio.yorn('Create BUILDGAL parameter files for each processor?'):
        owparam=mmlio.yorn('Overwrite existing files?')
        pardict=mkparam_proc(simstr,nproc=nproc,overwrite=owparam)
    # Units file
    if mmlio.yorn('Create BUILDGAL units file?'):
        owunits=mmlio.yorn('Overwrite existing file?')
        unidict=mkunits(simstr,overwrite=owunits,verbose=verbose)
    return

####################################################################################################################################
# METHOD TO CREATE UNITS FILE
def mkunits(simstr,overwrite=None,verbose=None):
    """
    Creates buildgal units file
    """
    # Pars input
    overwrite=mmlpars.mml_pars(overwrite,type=bool,default=False)
    verbose=mmlpars.mml_pars(verbose,type=bool,default=False)
    # Get file name
    unifile=simstr.fdict['buildgal']['input']['units']
    # Load file if it exists and overwrite not set
    if os.path.isfile(unifile) and not overwrite:
        if verbose: mmlio.verbose('Loading existing BUILDGAL units file:',addval=unifile)
        unidict=mmlio.rwdict('R',unifile)
    # Otherwise create it from scratch
    else:
        if verbose: mmlio.verbose('Creating BUILDGAL units file:',addval=unifile)
        # Get parameters
        param=model2param(simstr.infodict)
        units=mmlconst.const_units('buildgal')
        physG=mmlconst.const_phys('cgs')['G']
        # Create dictionary (assuming G=1)
        unidict={'UnitLength_in_cm':units['L']*param['disk_rscl'],
                 'UnitMass_in_g'   :units['M']*param['disk_mtot']}
        unidict['UnitVelocity_in_cm_per_s']=np.sqrt(physG*unidict['UnitMass_in_g']/unidict['UnitLength_in_cm'])
#        unidict['UnitMass_in_g']/=float(simstr['nprocbg'])
        # Save dictionary
        mmlio.rwdict('W',unifile,unidict,overwrite=overwrite)
    # Return dictionary
    return unidict

####################################################################################################################################
# METHOD TO CREATE PARAMETER FILES FOR EACH PROCESSOR
def mkparam_proc(simstr,nproc=None,overwrite=None,verbose=None,multfact=None):
    """
    Creates buildgal parameter files for each processor
    """
    # Pars input
    nproc=mmlpars.mml_pars(nproc,default=simstr['nprocbg'],type=int,min=0)
    multfact=mmlpars.mml_pars(multfact,type=float,default=2.,min=0.)
    overwrite=mmlpars.mml_pars(overwrite,type=bool,default=False)
    verbose=mmlpars.mml_pars(verbose,type=bool,default=False)
    # Check that particles evenly divisible across nproc
    if (simstr.ntot%nproc)!=0:
        raise Exception('The number of particles ({}) must be evenly divided across the number of processors ({}).'.format(simstr.ntot,nproc))
    # Load parameters
    par0=mmlio.rwdict('R',simstr.fdict['icgen']['bgpm'],style='fortran')
    par0=model2param(simstr.infodict,scale=True,par0=par0)
    # Adjust particle number
    for ityp in LIST_PTYPS+['sat']: 
        if ityp+'_n' in par0: par0[ityp+'_n']=long(par0[ityp+'_n']*multfact/float(nproc))
    # Loop over processers creating a file for each
    parbase=simstr.fdict['buildgal']['input']['parbase']
    random.seed(par0['randseed'])
    for iproc in range(nproc):
        ipar=copy.deepcopy(par0)
        ipar['randseed']=random.randint(0,20000000)
#        ipar['randseed']+=iproc
        ifile=parbase+str(iproc)
        ipar=mkparam(simstr,overwrite=overwrite,parfile=ifile,parstat=ipar)
    return par0
    
####################################################################################################################################
# METHOD TO CREATE SINGLE PARAMETER FILE
def mkparam(simstr,overwrite=None,parfile=None,parstat=None,verbose=None,**param_kw):
    """
    Creates a buildgal parameter file
    """
    # Set constants
    keywidth=29
    # Get file names
    files_local=simstr.mkfiledict(checkaccess=True)
    # Pars input
    parfile=mmlpars.mml_pars(parfile,default=files_local['icgen']['bgpm'],type=str)
    overwrite=mmlpars.mml_pars(overwrite,type=bool,default=False)
    verbose=mmlpars.mml_pars(verbose,type=bool,default=False)
    # Load existing file if overwrite False
    if not overwrite and os.path.isfile(parfile):
        if verbose: mmlio.verbose('Loading existing BUILDGAL parameter file:',addval=parfile)
        params=mmlio.rwdict('R',parfile,style='fortran')
    else:
        if isinstance(parstat,dict): params=getparam('ask',fpar=parstat)
        else                       : params=model2param(simstr.infodict)
    # Override based on input keys
    for ikey in param_kw.keys():
        if params.has_key(ikey):
            if verbose: mmlio.verbose('Overriding BUILDGAL parameter file option {} with value {}'.format(ikey,param_kw[ikey]))
            params[ikey]=param_kw[ikey]
    # Write parameter dictionary to file
    params['keylist']=getparam('list',fpar=params)
    mmlio.rwdict('W',parfile,params,style='fortran',width=keywidth,overwrite=overwrite)
    # Return output
    return params
    
####################################################################################################################################
# METHOD TO CREATE PBS SUBMISSION SCRIPT 
def mkpbs(simstr,tcpu=None,ppn=None,overwrite=None,verbose=None):
    """
    Create a BUILDGAL pbs script and returns content as a dictionary
    """
    # Pars input
    ppn=mmlpars.mml_pars(ppn,default=2,type=int,range=[1,4])
    overwrite=mmlpars.mml_pars(overwrite,type=bool,default=False)
    verbose=mmlpars.mml_pars(verbose,type=bool,default=False)
    # Get file names
    files_local=simstr.mkfiledict(checkaccess=True)
    files_accre=simstr.mkfiledict(memcomp='accre',checkaccess=False)
    pbsfile=files_local['buildgal']['input']['pbs']
    pbsdict=mmlio.pbsdict()
    # Load existing file if overwrite false
    if not overwrite and os.path.isfile(pbsfile):
        if verbose: mmlio.verbose('Loading existing BUILDGAL pbs script:',addval=pbsfile)
        pbsdict.rwmpi('R',pbsfile)
    # Otherwise create file from scratch
    else:
        if verbose: mmlio.verbose('Creating new BUILDGAL pbs script:',addval=pbsfile)
        # Update memory
        if not isinstance(tcpu,float): tcpu=mmlsim2tcpu(simstr)
        mcpu=mmlsim2mcpu(simstr)
        # Variables from simstr
        pbsdict['ppn'      ]=ppn
        pbsdict['nodes'    ]=simstr['nprocbg']/ppn
        pbsdict['pmem'     ]=mcpu
        pbsdict['mem'      ]=mcpu*simstr['nprocbg']
        pbsdict['walltime' ]=tcpu
        # Variables from path
        pbsdict['simout'   ]=files_accre['buildgal']['output']['runout']
        pbsdict['rundir'   ]=files_accre['buildgal']['input']['dir']
        pbsdict['execpath' ]=os.path.basename(files_accre['buildgal']['input']['wrap'])
        pbsdict['mpiopt'   ]='-comm none'
        pbsdict['execflag' ]=0
        pbsdict['execinput']=''
        pbsdict.rwmpi('W',pbsfile,overwrite=overwrite)
    # Return output
    return pbsdict

####################################################################################################################################
def mkwrap(simstr,overwrite=None,verbose=None):
    """
    Creates a BUILDGAL wrapper file
    """
    # Pars input
    overwrite=mmlpars.mml_pars(overwrite,default=False,type=bool)
    verbose=mmlpars.mml_pars(verbose,default=False,type=bool)
    # Get file names
    files_local=simstr.mkfiledict(checkaccess=True)
    files_accre=simstr.mkfiledict(memcomp='accre',checkaccess=False)
    wrpfile=files_local['buildgal']['input']['wrap']
    # Load file if it exists
    if os.path.isfile(wrpfile) and not overwrite:
        if verbose: mmlio.verbose('Loading existing BUILDGAL wrapper script:',addval=wrpfile)
        linelist=mmlio.rwbash('R',wrpfile)
    # Otherwise create it from scratch
    else:
        mmlio.verbose('Creating new BUILDGAL wrapper script:',addval=wrpfile)
        # Create lines
        linelist=['rank=`expr $1 + $MPIEXEC_RANK`',
                  'cd '+files_accre['buildgal']['input']['pardir'].replace('*','${rank}'),
                  'cp {} {}'.format(files_accre['buildgal']['input']['exec'],'buildgal.stall'),
                  './buildgal.stall > bgout']
        # Write file
        mmlio.rwbash('W',wrpfile,linelist=linelist,overwrite=overwrite)
        # Convert file to executable
        subprocess.call(['chmod','755',wrpfile])
#        os.chmod(wrpfile,stat.S_IXUSR)
#        os.chmod(wrpfile,stat.S_IXGRP)
#        os.chmod(wrpfile,stat.S_IXOTH)
    # Return output
    return linelist

####################################################################################################################################
# METHOD TO CONVERT MODEL PARAMETERS TO BUILDGAL PARAMETERS
def model2param(sim,delvir=None,randseed=None,minfracn=None,scale=None,par0={},rsoft=None):
    """
    Converts model parameters to buildgal parameters
    """
    # Pars input
    delvir=mmlpars.mml_pars(delvir,default=178.,type=float)
    randseed=mmlpars.mml_pars(randseed,default=123,type=int)
    minfracn=mmlpars.mml_pars(minfracn,default=0.01,type=float,min=0)
    scale=mmlpars.mml_pars(scale,default=False,type=bool)
    sim=simfile.parsfpar('icgen','galsim',fpar=sim)
    # Get galaxy model info
    mod=simfile.parsfpar('icgen','galmod',askuser=True,tagstr=sim['model'])
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
    for ityp in LIST_PTYPS: par[ityp+'_n']=long(mmlmath.oom(par[ityp+'_n'],nsig=2))
    # Fill in missing values with user/default input
    par=getparam('ask',fpar=par)
    print 'mscl={disk_mtot}'.format(**par)
    print 'rscl={disk_rscl}'.format(**par)
    # Return parameters
    return par

####################################################################################################################################
# METHOD TO ASK FOR BUILDGAL PARAMETERS
def getparam(method,fpar=None,**exkw):
    """
    Asks user for buildgal parameters
    """
    # Pars input
    method=mmlpars.mml_pars(method,list=['ask','pars','list'])
    fpar=mmlpars.mml_pars(fpar,default={},type=dict)
    if   method=='ask' : parkw=dict(inclextra=True,askuser=True ,load=False,save=False,**exkw)
    elif method=='pars': parkw=dict(inclextra=True,askuser=False,load=False,save=False,**exkw)
    # Base parameters
    if method=='list': flist=['randseed','outsoft']
    else             : fpar=simfile.parsfpar('buildgal','basepar',fpar=fpar,**parkw)
    # Disk parameters
    if fpar['disk']=='Y':
        if method=='list': flist=['disk_mtot']+flist+['disk_n','disk_zscl','disk_solarR','disk_solarQ','disk_eps','disk_zmax','disk_rmax','gas']
        else             : fpar=simfile.parsfpar('buildgal','diskpar',fpar=fpar,**parkw)
        # Gas parameters
        if fpar['gas']=='Y':
            if method=='list': flist+=['gas_n','gas_mtot','gas_T','gas_zscl','gas_zmax','gas_rmax','gas_rmin','gas_flatdisk']
            else             : fpar=simfile.parsfpar('buildgal','gaspar',fpar=fpar,**parkw)
            if fpar['gas_flatdisk']=='Y':
                if method=='list': flist+=['gas_flt_r','gas_flt_mint','gas_flt_mext','gas_flt_nint','gas_flt_next','gas_flt_mstar']
                else             : fpar=simfile.parsfpar('buildgal','gaspar_flatdisk',fpar=fpar,**parkw)
            if method=='list': flist+=['gas_selfgrav']
    # Bulge parameters
    if method=='list': flist+=['bulge']
    if fpar['bulge']=='Y':
        if method=='list': flist+=['bulge_mtot','bulge_rscl','bulge_selfgrav']
        else             : fpar=simfile.parsfpar('buildgal','bulgepar',fpar=fpar,**parkw)
        if fpar['bulge_selfgrav']=='Y':
            if method=='list': flist+=['bulge_n','bulge_rmax','bulge_eps','bulge_nonspher']
            else             : fpar=simfile.parsfpar('buildgal','bulgepar_selfgrav',fpar=fpar,**parkw)
            if fpar['bulge_nonspher']=='Y':
                if method=='list': flist+=['bulge_zscl','bulge_zmax','bulge_nsimp','bulge_rotate']
                else             : fpar=simfile.parsfpar('buildgal','bulgepar_nonspher',fpar=fpar,**parkw)
                if fpar['bulge_rotate']=='Y':
                    if method=='list': flist+=['bulge_fracrev']
                    else             : fpar=simfile.parsfpar('buildgal','bulgepar_rotate',fpar=fpar,**parkw)
    # Halo parameters
    if method=='list': flist+=['halo']
    if fpar['halo']=='Y':
        if method=='list': flist+=['halo_selfgrav','halo_rmax']
        else             : fpar=simfile.parsfpar('buildgal','halopar',fpar=fpar,**parkw)
        if fpar['halo_selfgrav']=='Y':
            if method=='list': flist+=['halo_n','halo_eps']
            else             : fpar=simfile.parsfpar('buildgal','halopar_selfgrav',fpar=fpar,**parkw)
        if method=='list':
            flist+=['halo_prof','halo_mtot']
            if   fpar['halo_prof']=='IS': flist+=['halo_IS_rcore','halo_IS_rtide']
            elif fpar['halo_prof']=='LH': flist+=['halo_LH_rscl']
            elif fpar['halo_prof']=='NF': flist+=['halo_NF_concen','halo_NF_delvir']
        else: fpar=simfile.parsfpar('buildgal','halopar_{}'.format(fpar['halo_prof']),fpar=fpar,**parkw)
    # Satellite parameters
    if method=='list': flist+=['sat']
    if fpar['sat']=='Y':
        if method=='list': flist+=['sat_mtot','sat_rscl','sat_x','sat_y','sat_z','sat_rmax','sat_selfgrav']
        else             : fpar=simfile.parsfpar('buildgal','satpar',fpar=fpar,**parkw)
        if fpar['sat_selfgrav']=='Y':
            if method=='list': flist+=['sat_n','sat_eps']
            else             : fpar=simfile.parsfpar('buildgal','satpar_selfgrav',fpar=fpar,**parkw)
        if method=='list': flist+=['sat_vx','sat_vy','sat_vz']
    # Two model parameters
    if method=='list': flist+=['twomod']
    if fpar['twomod']=='Y':
        if method=='list': flist+=['twomod_rperi','twomod_rinit','twomod_theta1','twomod_phi1','twomod_theta2','twomod_phi2']
        else             : fpar=simfile.parsfpar('buildgal','twomodpar',fpar=fpar,**parkw)
    # Return output
    if method=='list': return flist
    else             : return fpar

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
    mcpu=500000./simstr['nprocbg']
    # Return rounded value
    return mmlmath.oom(mcpu,nsig=1,method='CEIL')

####################################################################################################################################
# METHOD TO DETERMINE TIME REQUIRED TO RUN BUILDGAL
def mmlsim2tcpu(simstr):
    """
    Estimates the amount of time required to run BUILDGAL
    """
    # Force reasonable values as a safety
    tcpu=150.
    return tcpu

####################################################################################################################################
####################################################################################################################################
# PROVIDE COMMAND LINE ACCESS
if __name__ == '__main__': main()
