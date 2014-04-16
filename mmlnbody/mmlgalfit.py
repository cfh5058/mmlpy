#!/usr/bin/python
####################################################################################################################################
#
# MEAGAN LANG'S GALFIT METHODS
#
####################################################################################################################################
import sys,os,shutil,glob,copy,pprint,scipy,math
import numpy as np
import matplotlib as mplib
import pNbody
LIST_METHODS=['makefiles','submit','retrieve','clean']
from mmlutils import *
import main as mmlnbody
import simlist,simcalc,simplot,simfile,mmlbuildgal,mmlgadget,mmlscf,icgen

def main():
    """
    Provides command line interfact to GALFIT methods
    """
    out = mmlnbody.walk(mtype='galfit')
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
           'input' :['shortin' ,'longin' ],
           'output':['shortout','longout']}
    return fdict

####################################################################################################################################
####################################################################################################################################
# HIGH LEVEL METHODS REQUIRING MMLSIM

####################################################################################################################################
# METHOD FOR RUNNING DIFFERENT GALFIT METHODS
def run(simstr,method,**method_kw):
    """
    Provides interface for running different GALFIT operations
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
        simfile.rmfiles('galfit',simstr=simstr,**method_kw)
    else:
        raise Exception('Option for method {} needs to be added to mmlgalfit.run.'.format(method))
    # Return output
    return out

####################################################################################################################################
# METHOD PERFORM ACTIONS TO SETUP AND RUN GALFIT ON ACCRE
def subrun(simstr,ptype=None,**extra_kw):
    """
        Sets up and runs GALFIT on ACCRE
            """
    # Pars input
    ptype=mmlpars.mml_pars(ptype,default='disk',type=str)
    extra_kw['typList']=[ptype]
    files_accre=simstr.mkfiledict(memcomp='accre',checkaccess=True)['galfit']
    # Update input files
    if mmlio.yorn('Create new input files?'): mkfiles(simstr,**extra_kw)
    # Move input files to ACCRE
    if mmlio.yorn('Move new input files to ACCRE?'): putfiles(simstr,**extra_kw)
    # Make directories
    if mmlio.yorn('Make GALFIT output directories?'):
        for ityp in typList:
            if not mmlio.yorn('    For {} particles?'.format(ityp)): continue
            mmlfiles.mkdirs(files_accre[ityp]['output']['dir'])
            dkeylist=files_accre[ityp]['output']['subdirs']
            dirlist=[files_accre[ityp]['output'][idirkey] for idirkey in dkeylist]
            mmlfiles.mkdirs(dirlist)
    # Run GALFIT
    if os.environ['HOSTNAME'].lower()=='bender':
        if mmlio.yorn('Run GALFIT?'):
            pbsfile=files_accre[ptype]['input']['pbs']
            schfile=files_accre[ptype]['input']['sched']
            simstr.subpbs('galfit',ptype=ptype,pbsfile=pbsfile,schfile=schfile)
    return

####################################################################################################################################
# METHOD TO PERFORM ACTIONS TO RETRIEVE AND CLEAN-UP A RUN
def endrun(simstr,ptype=None):
    """
    Retrieves and cleans up results from an GALFIT run
    """
    # Pars input
    ptype=mmlpars.mml_pars(ptype,default='disk',type=str)
    # Retrieve results
    if mmlio.yorn('Recover GALFIT results from ACCRE?'): getfiles(simstr,typList=[ptype])
    # Clean up results
    if mmlio.yorn('Clean up GALFIT results on ACCRE?' ): simfile.rmfiles('galfit',ftype='clean',memcomp='accre',simstr=simstr)
    return
                                                                            
####################################################################################################################################
# METHOD TO MOVE INPUT FILES TO ACCRE
def putfiles(simstr,typList=None):
    """
    Moves files to ACCRE in preparation to running GALFIT code
    """
    # Pars input
    typListDEF=simlist.LIST_PTYPS
    typList=mmlpars.mml_pars(typList,default=typListDEF,type=list)
    # Move files
    if mmlio.yorn('Move GALFIT input files to ACCRE?'     ): simfile.mvfiles_accre(simstr,'to','galfit','input' ,typList=typList)
    if mmlio.yorn('Move old GALFIT output files to ACCRE?'): simfile.mvfiles_accre(simstr,'to','galfit','output',typList=typList)
    return

####################################################################################################################################
# METHOD TO RETRIEVE OUTPUT FROM AN GALFIT ANALYSIS RUN
def getfiles(simstr,typList=None):
    """
    Retrieves GALFIT files from ACCRE
    """
    # Move files
    if mmlio.yorn('Recover GALFIT input files from ACCRE?') : simfile.mvfiles_accre(simstr,'from','galfit','input' ,typList=typList)
    if mmlio.yorn('Recover GALFIT output files from ACCRE?'): simfile.mvfiles_accre(simstr,'from','galfit','output',typList=typList)
    return

####################################################################################################################################
####################################################################################################################################
# FILE CLASSES AND METHODS

####################################################################################################################################
# METHOD FOR PARSING GALFIT FILE OPTIONS
def pars_fileopt(fopt):
    """
    Returns parsed GALFIT file option dictionary with keys:
        execpath: Str absolute path to location of compiled GALFIT executable
        parmpath: Str absolute path to location of default GALFIT executable
        temppath: Str absolute path to temporary directory for running GALFIT locally
    """
    # Define values
    form={
        'execpath': mmlpars.parsdict(default='/home/langmm/bin/galfit/bin/galfit'                  ,type=str),
        'parmpath': mmlpars.parsdict(default='/home/langmm/bin/galfit/galfit-example/EXAMPLE.INPUT',type=str),
        'temppath': mmlpars.parsdict(default='/home/langmm/gfittemp'                               ,type=str)
        }
    # Pars input
    fopt=mmlpars.mml_formpars(fopt,form)
    # Return parsed input
    return fopt

####################################################################################################################################
# METHOD TO RETURN DICTIONARY DEFINING GALFIT FILES & DIRECTORIES
def files(fdict):
    """
    Returns a dictionary defining Galfit files & directories with keys:
    """
    # Set constants
    typList=simlist.LIST_PTYPS
    # Pars input
    fdict=mmlpars.mml_pars(fdict,type=dict)
    pfix=fdict['pfix']
    options=pars_fileopt({})
    # Initialize directory and default/static files
    files=dict(fdict['galfit'],
               exdeft=options['execpath'],pmdeft=options['parmpath'],tempdir=options['temppath'])
    # Loop over particle types creating a directory tree for each
    for ityp in typList:
        files[ityp]={'dir':os.path.join(files['dir'],ityp)}
        # Input directory
        infiles={}
        infiles['dir']=os.path.join(files[ityp]['dir'],'input')
        # Short input files
        fkeys_in=['exec','par','wrap','pbs']
        flist_in=['galfit','gfitpar','gfitwrap','gfitpbs']
        infiles=mmlfiles.filedict(infiles['dir'],filekeys=fkeys_in,filelist=flist_in,pfix=pfix,fdict=infiles)
        # Long input directories
        dkeys_in_long=['fitsdir','sigdir','psfdir','maskdir','cstrdir']
        dlist_in_long=['fits','sigma','psf','mask','constrain']
        infiles['subdirs']=dkeys_in_long
        infiles=mmlfiles.filedict(infiles['dir'],dirkeys=dkeys_in_long,dirlist=dlist_in_long,pfix=pfix,fdict=infiles)
        # Long input files
        fkeys_in_long=['fitsbase','sigbase','psfbase','maskbase','cstrbase']
        flist_in_long=['gfitfits','gfitsig','gfitpsf','gfitmask','gfitconstr']
        for ifin in range(len(dkeys_in_long)):
            infiles=mmlfiles.filedict(infiles[dkeys_in_long[ifin]],filekeys=[fkeys_in_long[ifin]],
                                      filelist=[flist_in_long[ifin]],pfix=pfix,fdict=infiles)
        files[ityp]['input']=infiles
        # Output directory
        outfiles={}
        outfiles['dir']=os.path.join(files[ityp]['dir'],'output')
        # Short output files
        fkeys_out=['pbsout']
        flist_out=['gfitout']
        outfiles=mmlfiles.filedict(outfiles['dir'],filekeys=fkeys_out,filelist=flist_out,pfix=pfix,fdict=outfiles)
        # Long output directories
        dkeys_out_long=['logdir','pardir','imgdir']
        dlist_out_long=['log','par','img']
        outfiles['subdirs']=dkeys_out_long
        outfiles=mmlfiles.filedict(outfiles['dir'],dirkeys=dkeys_out_long,dirlist=dlist_out_long,pfix=pfix,fdict=outfiles)
        # Long output files
        fkeys_out_long=['logbase','parbase','imgbase']
        flist_out_long=['gfitlog','gfitparout','gfitimgout']
        for ifout in range(len(dkeys_out_long)):
            outfiles=mmlfiles.filedict(outfiles[dkeys_out_long[ifout]],filekeys=[fkeys_out_long[ifout]],
                                       filelist=[flist_out_long[ifout]],pfix=pfix,fdict=outfiles)
        files[ityp]['output']=outfiles
    # Return file dictionary
    fdict['galfit']=files
    return fdict

####################################################################################################################################
# METHOD TO RETURN LIST OF RELAVENT FILES
def get_filelist(fdict,keylist,ptype=None):
    """
    Returns a list of relavent GALFIT files
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
        filelist['default' ]=[fdict['exdeft'],fdict['pmdeft']]
    if 'static'   in keylist:
        filelist['static'  ]=[]
    if 'shortin'  in keylist:
        filelist['shortin' ]=[fdict[ptype]['input']['exec'],
                              fdict[ptype]['input']['par'],fdict[ptype]['input']['wrap'],
                              fdict[ptype]['input']['pbs']]
    if 'longin'   in keylist:
        baselist=['fits','sig','psf','mask','cstr']
        filelist['longin'  ]=[fdict[ptype]['input'][ibase+'base']+'*' for ibase in baselist]
    if 'shortout' in keylist:
        filelist['shortout']=[fdict[ptype]['output']['pbsout']]
    if 'longout'  in keylist:
        baselist=['log','par','img']
        filelist['longout' ]=[fdict[ptype]['output'][ibase+'base']+'*' for ibase in baselist]
    # Return filelist
    return filelist

####################################################################################################################################
####################################################################################################################################
# METHODS FOR CREATING GALFIT FILES

####################################################################################################################################
# METHOD TO PARS PARAMETER FILE OPTION DICTIONARY
def pars_paropt(popt):
    """
    Returns a parsed dictionary containing GALFIT parameter file options
    """
    # Define values
    form={
        }
    # Pars input
    popt=mmlpars.mml_formpars(popt,form)
    # Return parsed input
    return popt

####################################################################################################################################
# METHOD TO CREATE GALFIT FILES LOCALLY
def mkfiles(simstr,paropt=None,ext=None,typList=None,nproc=None,idgal2=None):
    """
    Creates GALFIT files
    """
    mmlio.verbose('Creating files for GALFIT run locally...',border=True)
    # Set constants
    typListDEF=simlist.LIST_PTYPS
    typDict={typListDEF[ityp]:ityp for ityp in range(len(typListDEF))}
    # Pars input
    paropt=pars_paropt(paropt)
    ext=mmlpars.mml_pars(ext,default='*',type=str)
    typList=mmlpars.mml_pars(typList,default=typListDEF,type=list)
    nproc=mmlpars.mml_pars(nproc,default=simstr['nprocgfit'],type=int)
    # Get file names
    files=simstr.fdict
    # Loop over particle types creating directories & asking if files should be created
    # NOTE: create sigma,psf,mask,cstr files?
    dict_mkfile={}
    dict_mkexec={}
    dict_mkparm={}
    dict_mksnap={}
    flag_mksnap=False
    for iptyp in typList:
        dict_mkfile[iptyp]=mmlio.yorn('    {} particles?'.format(iptyp))
        if dict_mkfile[iptyp]:
            files_typ=files['galfit'][iptyp]
            mmlfiles.mkdirs([files_typ['dir'],files_typ['input']['dir']]+files_typ['input']['subdirs'])
            dict_mkexec[iptyp]=mmlio.yorn('      GALFIT executable?    ({})'.format(iptyp))
            dict_mkparm[iptyp]=mmlio.yorn('      GALFIT runtime files? ({})'.format(iptyp))
            dict_mksnap[iptyp]=mmlio.yorn('      GALFIT snapshots?     ({})'.format(iptyp))
            if dict_mksnap[iptyp]: flag_mksnap=True
        else:
            dict_mkexec[iptyp]=False
            dict_mkparm[iptyp]=False
            dict_mksnap[iptyp]=False
        if flag_mksnap: owsnap=mmlio.yorn('Overwrite existing GALFIT snapshots?')
    # Single files
    for iptyp in typList:
        if not dict_mkfile[iptyp]: continue
        # Executable file
        if dict_mkexec[iptyp]:
            if not os.path.isfile(files['galfit']['exdeft']):
                raise Exception('Static GALFIT executable does not exist. Create it first: {}'.format(files['galfit']['exdeft']))
            mmlfiles.cpfiles(files['galfit']['exdeft'],files['galfit'][iptyp]['input']['exec'],overwrite=True)
        # Parameter files
        if dict_mkparm[iptyp]:
            pardict=mkpar(simstr,iptyp,overwrite=True)
            pbsdict=mkpbs(simstr,iptyp,overwrite=True,nproc=nproc)
            wrapdict=mkwrap(simstr,iptyp,overwrite=True,mpipath=pbsdict['mpipath'],nproc=nproc)
    # Snapshots & IC file

    raise Exception('mmlgalfit.mkfiles is a work in progress')
    return

####################################################################################################################################
# METHOD TO CREATE GALFIT PARAMETER FILE
def mkpar(simstr,ptype,overwrite=None):
    """
    Creates a GALFIT parameter file and returns contents as a dictionary
    """

####################################################################################################################################
# METHOD TO CREATE GALFIT PBS SCRIPT
def mkpbs(simstr,ptype,overwrite=None,nproc=None):
    """
    Creates a GALFIT pbs script and returns contents as a dictionary
    """

####################################################################################################################################
# METHOD TO CREATE GALFIT WRAPPER SCRIPT
def mkwrap(simstr,ptype,overwrite=None,mpipath=None,nproc=None):
    """
    Creates a GALFIT wrapper script and returns contents as a dictionary
    """

####################################################################################################################################
####################################################################################################################################
# PROVIDE COMMAND LINE ACCESS
if __name__ == '__main__': main()
