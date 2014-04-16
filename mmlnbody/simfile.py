#!/usr/bin/python
####################################################################################################################################
#
# MEAGAN LANG'S SIMFILE METHODS
#
####################################################################################################################################
import sys,os,shutil,glob,copy,pprint,scipy,math,itertools,re
import numpy as np
import matplotlib as mplib
import pNbody
LIST_METHODS=['getfdict','getflist','mvrun','mvtree','ziprun','unziprun','procnum','cleantree']
from mmlutils import *
import main as mmlnbody
import simlist,simcalc,simhist,simprof,simstat,simplot,mmlbuildgal,mmlgadget,mmlscf,mmlgalfit,icgen,compare
LIST_FILEMETHODS=copy.deepcopy(simlist.LIST_METHTYP)
for imtyp in simlist.LIST_METHTYP_NOF: LIST_FILEMETHODS.remove(imtyp)
DIR_FPAR='/home/langmm/analysis/'

def main():
    """
    Provides command line access to FILE methods
    """
    out = mmlnbody.walk(mtype='file')
    return out

####################################################################################################################################
####################################################################################################################################
# METHODS FOR GETTING LISTS
def dict_ftype(fmeth=None):
    """
    Returns a dictionary of files listed by type
    """
    if   fmeth=='all'     : fdict={}
    elif fmeth=='icgen'   : fdict=icgen.dict_ftype()
    elif fmeth=='buildgal': fdict=mmlbuildgal.dict_ftype()
    elif fmeth=='gadget'  : fdict=mmlgadget.dict_ftype()
    elif fmeth=='scf'     : fdict=mmlscf.dict_ftype()
    elif fmeth=='galfit'  : fdict=mmlgalfit.dict_ftype()
    elif fmeth=='calc'    : fdict=simcalc.dict_ftype()
    elif fmeth=='stat'    : fdict=simstat.dict_ftype()
    elif fmeth=='hist'    : fdict=simhist.dict_ftype()
    elif fmeth=='plot'    : fdict=simplot.dict_ftype()
    elif fmeth=='compare' : fdict=compare.dict_ftype()
    else: raise Exception('Invalid file method: {}'.format(fmeth))
#    ftyps=list_filetypes(fmeth=fmeth,baseonly=True)
#    for iftyp in ftyps:
#        if iftyp not in fdict: fdict[iftyp]=[]
    return fdict
def dict_fgroup(fmeth=None):
    """
    Returns a dictionary of files types list by group
    """
    if   fmeth=='all'     : fdict={}
    elif fmeth=='icgen'   : fdict=icgen.dict_fgroup()
    elif fmeth=='buildgal': fdict=mmlbuildgal.dict_fgroup()
    elif fmeth=='gadget'  : fdict=mmlgadget.dict_fgroup()
    elif fmeth=='scf'     : fdict=mmlscf.dict_fgroup()
    elif fmeth=='galfit'  : fdict=mmlgalfit.dict_fgroup()
    elif fmeth=='calc'    : fdict=simcalc.dict_fgroup()
    elif fmeth=='stat'    : fdict=simstat.dict_fgroup()
    elif fmeth=='hist'    : fdict=simhist.dict_fgroup()
    elif fmeth=='plot'    : fdict=simplot.dict_fgroup()
    elif fmeth=='compare' : fdict=compare.dict_fgroup()
    else: raise Exception('Invalid file method: {}'.format(fmeth))
#    fgrps=list_filegroups(fmeth=fmeth,baseonly=True)
#    for ifgrp in fgrps:
#        if ifgrp not in fdict: fdict[ifgrp]=[]
    return fdict
def list_filetypes(fmeth=None,flength=None,baseonly=None):
    """
    Returns a list of supported file types
    """
    # Add default files based on file length
    flength=mmlpars.mml_pars(flength,default='all',list=list_filelengths())
    baseonly=mmlpars.mml_pars(baseonly,default=False,type=bool)
    if   flength=='all'  : typeLIST=['default','static',
                                     'shortin','longin',
                                     'shortout','longout']
    elif flength=='long' : typeLIST=['longin','longout']
    elif flength=='short': typeLIST=[]
    else: raise Exception('Invalid file length {}'.format(flength))
    # Add listed long and short files
    if not baseonly:
        ftypes=dict_ftype(fmeth=fmeth)
        for itype in typeLIST:
            if itype in ftypes: typeLIST+=ftypes[itype]
    return typeLIST
def list_filegroups(fmeth=None,baseonly=None):
    """
    Returns a list of supported file groups
    """
    fmeth=mmlpars.mml_pars(fmeth,default='all',list=['all']+LIST_FILEMETHODS)
    baseonly=mmlpars.mml_pars(baseonly,default=False,type=bool)
    groupLIST=['all','clean','input','output']
    if not baseonly:
        fgrps=dict_fgroup(fmeth=fmeth)
        for igrp in groupLIST:
            if igrp in fgrps: groupLIST+=fgrps[igrp]
    return groupLIST
def list_filelengths():
    """
    Returns a list of strings identifying length of file list
    """
    flenLIST=['all','long','short']
    return flenLIST
def list_fpar(fmeth,method,plot=None):
    """
    Returns a dictionary of required parameter lists
    """
    fmeth=mmlpars.mml_pars(fmeth,list=LIST_FILEMETHODS)
    if   fmeth=='prof'    : fpar=simprof.list_fpar(method)
    elif fmeth=='stat'    : fpar=simstat.list_fpar(method,plot=plot)
    elif fmeth=='calc'    : fpar=simcalc.list_fpar(method)
    elif fmeth=='hist'    : fpar=simhist.list_fpar(method,plot=plot)
    elif fmeth=='plots'   : fpar=simplot.list_fpar(method)
    elif fmeth=='icgen'   : fpar=icgen.list_fpar(method)
    elif fmeth=='buildgal': fpar=mmlbuildgal.list_fpar(method)
    else: raise Exception('Invalid file method: {}'.format(fmeth))
    if 'keylist' in fpar:
        keylist=copy.deepcopy(fpar['keylist'])
        for ikey in fpar.keys():
            if ikey=='keylist': continue
            if ikey not in keylist:
                mmlio.verbose('Key {} not in {} {} parameter keylist.'.format(ikey,fmeth,method))
                fpar['keylist'].append(ikey)
        for ikey in keylist:
            if ikey not in fpar:
                mmlio.verbose('Key {} in {} {} parameter keylist, but not fpar.'.format(ikey,fmeth,method))
                fpar['keylist'].remove(ikey)
    else: fpar['keylist']=sorted(fpar.keys())
    return fpar
def list_fparDEF(fmeth,method,plot=None):
    """
    Returns a dictionary of default parameters
    """
    fmeth=mmlpars.mml_pars(fmeth,list=LIST_FILEMETHODS)
    if   fmeth=='prof'    : fpar=simprof.list_fparDEF(method)
    elif fmeth=='stat'    : fpar=simstat.list_fparDEF(method,plot=plot)
    elif fmeth=='calc'    : fpar=simcalc.list_fparDEF(method)
    elif fmeth=='hist'    : fpar=simhist.list_fparDEF(method,plot=plot)
    elif fmeth=='plots'   : fpar=simplot.list_fparDEF(method)
    elif fmeth=='icgen'   : fpar=icgen.list_fparDEF(method)
    elif fmeth=='buildgal': fpar=mmlbuildgal.list_fparDEF(method)
    else: raise Exception('Invalid file method: {}'.format(fmeth))
    return fpar
def iter_fpar(fmeth,method,fpar=None):
    """
    Returns a list of dictionaries iterating over parameters
    """
    if fpar==None: fpar=list_fpar(fmeth,method)
    pnam=list_fpar(fmeth,method)['keylist']
    if len(pnam)==0: return []
    pval=[fpar[inam] for inam in pnam]
    piter=itertools.product(*pval)
    iterpar=[]
    for iplist in piter:
        ifpar={ipnam:ipval for ipnam,ipval in zip(pnam,iplist)}
        iterpar.append(ifpar)
    return iterpar

####################################################################################################################################
####################################################################################################################################
# HIGH LEVEL METHODS REQUIRING MMLSIM

####################################################################################################################################
# METHOD FOR RUNNING DIFFERENT FILE METHODS
def run(simstr,method,*method_arg,**method_kw):
    """
    Provides interface for running different FILE methods
    """
    # Pars input
    method=mmlpars.mml_pars(method,list=LIST_METHODS)
    # Initialize default output
    out=None
    # Proceed based on method
    if   method == 'getfdict' : out=files(simstr=simstr,**method_kw)
    elif method == 'getflist' : out=get_filelist(simstr=simstr,*method_arg,**method_kw)
    elif method == 'mvrun'    : mvrun(simstr,**method_kw)
    elif method == 'mvtree'   : mvtree(simstr,**method_kw)
    elif method == 'ziprun'   : ziprun(simstr,**method_kw)
    elif method == 'unziprun' : unziprun(simstr,**method_kw)
    elif method == 'procnum'  : set_procnum(simstr,**method_kw)
    elif method == 'cleantree': cleantree(simstr=simstr,**method_kw)
    else:
        raise Exception('Option for method {} needs to be added to simfile.run.'.format(method))
    # Return output
    return out

####################################################################################################################################
# METHOD TO TURN PARAMETERS INTO A TAG
def fpar2tag(fmeth,method,inpar,plot=None,**exkw):
    """
    Returns a tag for a given set of parameters and method
    """
    fmeth=mmlpars.mml_pars(fmeth,list=LIST_FILEMETHODS)
    inpar=mmlpars.mml_pars(inpar,type=dict,default={})
    if   fmeth=='prof'    : tagstr=simprof.fpar2tag(method,inpar)
    elif fmeth=='stat'    : tagstr=simstat.fpar2tag(method,inpar,plot=plot)
    elif fmeth=='calc'    : tagstr=simcalc.fpar2tag(method,inpar)
    elif fmeth=='hist'    : tagstr=simhist.fpar2tag(method,inpar,plot=plot)
    elif fmeth=='plots'   : tagstr=simplot.fpar2tag(method,inpar)
    elif fmeth=='icgen'   : tagstr=icgen.fpar2tag(method,inpar)
    elif fmeth=='buildgal': tagstr=mmlbuildgal.fpar2tag(method,inpar)
    else: raise Exception('Invalid file method: {}'.format(fmeth))
    return tagstr

####################################################################################################################################
# METHOD TO RETURN LIST OF FILES
def fpar_taglist(fmeth,method):
    """
    Returns a list of existing fpar file tags
    """
    fmeth=mmlpars.mml_pars(fmeth,list=LIST_FILEMETHODS)
    method=mmlpars.mml_pars(method,type=str)
    fdir=os.path.join(DIR_FPAR,fmeth,method)
    flist=sorted(glob.glob(fdir+'/*'))
    tlist=[]
    for ifile in flist: 
        if not ifile.endswith('~'): tlist.append(fname2ftag(fmeth,method,ifile))
    return tlist
    
####################################################################################################################################
# METHOD TO RETURN FILE NAME FROM TAGSTR
def ftag2fname(fmeth,method,tagstr):
    """
    Converts tagstr to a file name
    """
    fmeth=mmlpars.mml_pars(fmeth,list=LIST_FILEMETHODS)
    method=mmlpars.mml_pars(method,type=str)
    tagstr=mmlpars.mml_pars(tagstr,type=str)
    tagstr=tagstr.split('{}_'.format(method))[-1]
    fname=os.path.join(DIR_FPAR,fmeth,method,'{}_{}'.format(method,tagstr))
    return fname
def fname2ftag(fmeth,method,fname):
    """
    Converts file name to tagstr
    """
    fmeth=mmlpars.mml_pars(fmeth,list=LIST_FILEMETHODS)
    method=mmlpars.mml_pars(method,type=str)
    fname=mmlpars.mml_pars(fname,type=str)
    tagstr=os.path.basename(fname).split('{}_'.format(method))[-1]
    return tagstr
def fpar2fname(fmeth,method,inpar=None,tagstr=None,**exkw):
    """
    Convers parameters to a file name
    """
    fmeth=mmlpars.mml_pars(fmeth,list=LIST_FILEMETHODS)
    if not isinstance(tagstr,str): tagstr=fpar2tag(fmeth,method,inpar,**exkw)
    if tagstr==os.path.basename(tagstr): fname=ftag2fname(fmeth,method,tagstr)
    else                               : fname=tagstr
    return fname

####################################################################################################################################
# METHOD TO SAVE FILE PARAMETERS
def savefpar(fmeth,method,inpar,tagstr=None,overwrite=None,**exkw):
    """
    Saves parameters to a file 
    """
    fmeth=mmlpars.mml_pars(fmeth,list=LIST_FILEMETHODS)
    outpar=parsfpar(fmeth,method,fpar=inpar,**exkw)
    fname=fpar2fname(fmeth,method,inpar=outpar,tagstr=tagstr,**exkw)
    overwrite=mmlpars.mml_pars(overwrite,default=False,type=bool)
    outpar['keylist']=list_fpar(fmeth,method)['keylist']
    if not os.path.isfile(fname) or overwrite:
        mmlio.rwdict('W',fname,outpar,overwrite=overwrite)
        mmlio.verbose('Saved file parameters to file:')
        print '    '+fname
    return outpar

####################################################################################################################################
# METHOD TO LOAD PARAMETERS FROM FILE
def loadfpar(fmeth,method,fname,**exkw):
    """
    Loads parameters from a file
    """
    fmeth=mmlpars.mml_pars(fmeth,list=LIST_FILEMETHODS)
    fname=fpar2fname(fmeth,method,tagstr=fname,**exkw)
    if not os.path.isfile(fname): raise Exception('Invalid file name: {}'.format(fname))
    inpar=mmlio.rwdict('R',fname)
    outpar=parsfpar(fmeth,method,fpar=inpar,**exkw)
    return outpar

####################################################################################################################################
# METHOD TO ASK USER FOR PARAMETERS
def parsfpar(fmeth,method,fpar=None,tagstr=None,fpar_deft=None,plot=None,overwrite=None,
             askuser=None,save=None,load=None,init=None,inclextra=None,**extrakw):
    """
    Asks user for info required to create a histogram
    """
    # Set constants
    fpar_list=list_fpar(fmeth,method,plot=plot)
    fpar_deft0=mmlpars.mml_pars(fpar_deft,type=dict,default={})
    fpar_deft=list_fparDEF(fmeth,method,plot=plot)
    fpar_deft.update(**fpar_deft0)
    fpar_keys=fpar_list['keylist']
    # Pars input 
    fpar=copy.deepcopy(mmlpars.mml_pars(fpar,default={},type=dict))
    askuser=mmlpars.mml_pars(askuser,default=False,type=bool)
    if askuser: loadDEF=True  ; saveDEF=True  ; initDEF=True
    else      : loadDEF=False ; saveDEF=False ; initDEF=False
    load=mmlpars.mml_pars(load,default=loadDEF,type=bool)
    save=mmlpars.mml_pars(save,default=saveDEF,type=bool)
    init=mmlpars.mml_pars(init,default=initDEF,type=bool)
    inclextra=mmlpars.mml_pars(inclextra,default=False,type=bool)
    # Get fpar dictionary from keys 
    for ikey in fpar_keys:
        if ikey in extrakw : fpar[ikey]=extrakw[ikey]
        if ikey not in fpar: 
            if init: fpar[ikey]=None
            else   : raise Exception('Parameter {} must be provided for {} method {}.'.format(ikey,fmeth,method))
    # Get fpar dictionary from file 
    if load:
        if not isinstance(tagstr,str) and askuser:
            tagstr=mmlio.askquest('Enter file containing {} {} parameters'.format(fmeth,method),default='None',dtype='str')
            if tagstr=='None': tagstr=None
        if isinstance(tagstr,str): 
            if os.path.isfile(fpar2fname(fmeth,method,tagstr=tagstr)):
                ofpar=loadfpar(fmeth,method,tagstr)
                for ikey in ofpar.keys(): fpar[ikey]=ofpar[ikey]
            else:
                mmlio.verbose('File does not exist. Creating it...')
    # Ask for lists
    for ikey in fpar_keys:
        if isinstance(fpar_list[ikey],list):
            if fpar[ikey] not in fpar_list[ikey]:
                if askuser: fpar[ikey]=mmlio.askselect('Select valid list member for {}.'.format(ikey),fpar_list[ikey])
                else      : fpar[ikey]=mmlpars.mml_pars(fpar[ikey],list=fpar_list[ikey],default=fpar_deft[ikey])
        elif isinstance(fpar_list[ikey],type):
            if not isinstance(fpar[ikey],fpar_list[ikey]):
                if askuser: fpar[ikey]=mmlio.askquest('Enter valid value for {}.'.format(ikey),dtype=str(fpar_list[ikey]),default=fpar_deft[ikey])
                else      : fpar[ikey]=mmlpars.mml_pars(fpar[ikey],type=fpar_list[ikey],default=fpar_deft[ikey])
        else: raise Exception('Invalid {} fpar list type: {}'.format(ikey,type(fpar_list)))
    # Remove extra tags
    if not inclextra:
        for ikey in fpar.keys():
            if ikey not in fpar_keys: del fpar[ikey]
    # Keylist
#    fpar['keylist']=fpar_keys
    # Save file
    if save: ofpar=savefpar(fmeth,method,fpar,inclextra=inclextra,tagstr=tagstr,overwrite=overwrite)
    # Return dictionaries
    return fpar

####################################################################################################################################
# METHOD TO RETURN FILES FOR A GIVEN METHOD
def fpar2file(simstr,fmeth,method,inpar,fext=None,test=None,plot=None,sing=None,
              ext=None,teststr0=None,plotstr0=None,animstr0=None):
    """
    Returns a dictionary of file info
    """
    # Set constants
    plotopt=mmlplot.parsopt()
    plotext0=plotopt['plotext']
    animext0=plotopt['animext']
    # Pars input
    if fmeth not in LIST_FILEMETHODS:
        print LIST_FILEMETHODS
        mmlio.yorn('{}'.format(fmeth))
    fmeth=mmlpars.mml_pars(fmeth,list=LIST_FILEMETHODS)
    ftag=fpar2tag(fmeth,method,inpar,plot=plot)
    fext=mmlpars.mml_pars(fext,type=str,default='_*')
    test=mmlpars.mml_pars(test,type=bool,default=False)
    plot=mmlpars.mml_pars(plot,type=bool,default=False)
    sing=mmlpars.mml_pars(sing,type=bool,default=False)
    ext=mmlpars.mml_pars(ext,type=str,default='')
    teststr0=mmlpars.mml_pars(teststr0,default='test_',type=str)
    plotstr0=mmlpars.mml_pars(plotstr0,default='plot',type=str)
    animstr0=mmlpars.mml_pars(animstr0,default='anim',type=str)
    # Handle flags
    if test: teststr=teststr0
    else   : teststr=''
    if plot: plotstr=plotstr0
    else   : plotstr=''
    if plot: ext=plotext0
    if sing: fext=''
    # Get file info
    try:
        fdict=simstr.fdict
        methdir=fdict[fmeth]['dir']
        pfix=fdict['pfix']
    except:
        fdict=simstr
        if fmeth in fdict and 'pfix' in fdict:
            methdir=fdict[fmeth]['dir']
            pfix=fdict['pfix']
        elif 'dir' in fdict and 'pfix' in fdict:
            methdir=fdict['dir']
            pfix=fdict['pfix']
        else: raise
    # Initialize dictionary
    snfdict={fmeth+'dir':methdir}
    # Tag strings
    begtag='{}{}'.format(plotstr,method)
    if len(ftag)==0: endtag=''
    else           : endtag='_{}'.format(ftag)
    # Directory
    if sing:
        snfdict['dir']=os.path.join(methdir,begtag)
    else:
        dir='{}{}'.format(begtag,endtag)
        if len(ftag)==0: snfdict['dir']=os.path.join(methdir,dir)
        else           : snfdict['dir']=os.path.join(methdir,begtag,dir)
    # File base
    base='{}{}{}{}'.format(teststr,pfix,begtag,endtag)
    snfdict['base']=os.path.join(snfdict['dir'],base)
    # File name
    snfdict['file']='{}{}'.format(snfdict['base'],fext)
    if len(ext)>0: snfdict['file']+='.{}'.format(ext)
    # Animation
    if plot:
        # Extensions
        snfdict['animext']=animext0
        snfdict['plotext']=plotext0
        # Directory
        snfdict['animdir']=os.path.join(methdir,animstr0+method)
        # File
        anim='{}{}{}{}{}.{}'.format(teststr,pfix,animstr0,method,endtag,animext0)
        snfdict['anim']=os.path.join(snfdict['animdir'],anim)
    # Return dictionary
    return snfdict

####################################################################################################################################
# METHOD TO RETURN SIMULATION FILES
def files(fmeth=None,**input_kw):
    """
    Generates a dictionary of simulation file names
    """
    # Pars input
    fdict=pars_fileinput(**input_kw)
    fdict['profpar']=os.path.join('rundir',fdict['pfix_shrt']+'profpar')
    # Create directories
    for ifmeth in LIST_FILEMETHODS:
        fdict[ifmeth]={'dir' :os.path.join(fdict['rundir'],ifmeth+fdict['subtag']),
                       'pfix':fdict['pfix'],'pfix_shrt':fdict['pfix_shrt']}
    # Create ICGEN file names
    fdict=icgen.files(fdict)
    # Create BUILDGAL file names
    fdict=mmlbuildgal.files(fdict)
    # Create GADGET file names
    fdict=mmlgadget.files(fdict)
    # Create SCF file names
    fdict=mmlscf.files(fdict) # Forcing SCF to use default options
    # Create GALFIT file names
    fdict=mmlgalfit.files(fdict)
    # Create CALC file names
    fdict=simcalc.files(fdict)
    # Create HIST file names
#    fdict=simhist.files(fdict)
    # Create STAT file names
    fdict=simstat.files(fdict)
    # Create PLOT file names
    fdict=simplot.files(fdict)
    # Create COMPARE file names
    fdict=compare.files(fdict,**input_kw)
    # Select subset of files as necessary
    if fmeth in LIST_FILEMETHODS:
        return fdict[fmeth]
    else:
        return fdict

####################################################################################################################################
# METHOD TO RETURN LIST OF RELEVANT FILES
def get_filelist(fmeth,ftype=None,ptype=None,**input_kw):
    """
    Returns a list of relavent files
    """
    # Pars input
    fmeth=mmlpars.mml_pars(fmeth,list=LIST_FILEMETHODS)
    # Get groups and list
    ftypDICT=dict_ftype(fmeth=fmeth)
    fgrpDICT=dict_fgroup(fmeth=fmeth)
    ftypLIST=ftypDICT.keys()
    fgrpLIST=fgrpDICT.keys()
    # Pars file type & get files
    ftype=mmlpars.mml_pars(ftype,default='all',list=['all']+fgrpLIST+ftypLIST)
    genkeys=dict(ftype=ftype,**input_kw)
    fdict=files(**genkeys)
    # Determine what files to add
    if   ftype=='all'     : keylist=ftypLIST.remove('default')
    else:
        if   ftype in fgrpLIST: keylist=fgrpDICT[ftype]
        elif ftype in ftypLIST: keylist=[ftype]
        else: raise Exception('Invalid file type: {}'.format(ftype))
    # Add files grouped by type
    keylist0=copy.deepcopy(keylist)
    for ikey in keylist0:
        if ikey in ftypDICT and len(ftypDICT[ikey])>0:
            keylist+=ftypDICT[ikey]
            keylist.remove(ikey)
    # Call correct funciton for fmeth
    inlist=[fdict[fmeth],keylist]
    if   fmeth=='icgen'   : filelist=icgen.get_filelist(*inlist)
    elif fmeth=='buildgal': filelist=mmlbuildgal.get_filelist(*inlist)
    elif fmeth=='gadget'  : filelist=mmlgadget.get_filelist(*inlist)
    elif fmeth=='scf'     : filelist=mmlscf.get_filelist(*inlist,ptype=ptype)
    elif fmeth=='galfit'  : filelist=mmlgalfit.get_filelist(*inlist,ptype=ptype)
    elif fmeth=='calc'    : filelist=simcalc.get_filelist(*inlist)
    elif fmeth=='stat'    : filelist=simstat.get_filelist(*inlist)
    elif fmeth=='hist'    : filelist=simhist.get_filelist(*inlist)
    elif fmeth=='plots'   : filelist=simplot.get_filelist(*inlist)
    elif fmeth=='compare' : filelist=simplot.get_filelist(*inlist)
    else: raise Exception('Invalid fmeth: {}'.format(fmeth))
    # Fill in missing values
    for ikey in keylist:
        if ikey not in filelist: filelist[ikey]=[]
    filelist['keylist']=keylist
    # Return output
    return filelist

####################################################################################################################################
####################################################################################################################################
# METHODS FOR MOVING/RENAMING RUN FILES/DIRECTORIES

####################################################################################################################################
# METHOD FOR MOVING/RENAMING A RUN
def mvrun(isimstr,overwrite=None,**options):
    """
    Move/rename files assosiated with a run
    """
    # Get sim objects
    print 'Getting info on new file names...'
    fsimstr=mmlnbody.asksimobj(infodict=isimstr.infodict,runtyp=isimstr['runtyp'])
    fname=fpar2fname('icgen','mmlsim',fpar=dict(fsimstr))
    # Rename old run files
    fname_new=fname+'_new'
    fname_old=fname+'_old'
    if os.path.isfile(fname): os.rename(fname,fname_old)
    if os.path.isfile(fname_new):
        print 'Using info from new file: {}'.format(fname_new)
        shutil.copy2(fname_new,fname)
    # Get sim objects
    shutil.copy2(fname,fname_new)
    # Move files
    try:
        mvfiles(isimstr,fsimstr,fmeth='all',overwrite=overwrite,**options)
        cleantree(simstr=isimstr,**options)
        if os.path.isfile(fname_old): os.remove(fname_old)
        os.remove(fname_new)
    except:
        print "Error moving files:", sys.exc_info()[0]
        print "Moving run files back..."
        os.rename(fname,fname_new)
        if os.path.isfile(fname_old): os.rename(fname_old,fname)
        raise
    # Return control
    return

####################################################################################################################################
# METHOD FOR MOVING AN ENTIRE SIMULATION DIRECTORY TREE
def mvtree(simstr,comp_end=None):
    """
    Moves an entire run tree without changing file names.
    """
    # Set constants
    usrcomp=mmlinfo.hostname2compid()
    complist=mmlinfo.complist()
    usrindx=complist['compid'].index(usrcomp)
    # Get location of sim
    comp_beg=simstr['memcomp']
    # Get computer memory info
    if comp_beg not in complist['memlist'][usrindx]:
        raise Exception('Cannot access computer {} from the current one {}. Switch to bender.'.format(comp_beg,usrcomp))
    print 'Current tree location: {}'.format(comp_beg)
    if comp_end not in complist['memlist'][usrindx]:
        comp_end=mmlio.askselect('What computer should the tree be moved to?',complist['memlist'][usrindx])
    # Get filenames
    files_beg=simstr.mkfiledict(memcomp=comp_beg,checkaccess=True)
    files_end=simstr.mkfiledict(memcomp=comp_end,checkaccess=True)
    # Copy the tree
    print 'Source      run tree: '+files_beg['rundir']
    print 'Destination run tree: '+files_end['rundir']
    if mmlio.yorn('Continue moving tree?'):
        mmlfiles.mkdirs(files_end['simdir'])
        shutil.copytree(files_beg['rundir'],files_end['rundir'])
    else:
        return
    # Remove the old tree
    if mmlio.yorn('Remove old run tree?'):
        shutil.rmtree(files_beg['rundir'])
    # Update the run files
    if mmlio.yorn('Update the run files?'):
        simstr['memcomp']=comp_end
        savefpar('icgen','galsim',dict(simstr),overwrite=True)
    return

####################################################################################################################################
# METHOD FOR MOVING/RENAMING RUN FILES
def mvfiles(isimstr,fsimstr,fmeth=None,memcomp1=None,memcomp2=None,
            overwrite=None,typList=None,snapend=None,**options):
    """
    Moves/renames simulation files
    """
    # Pars input
    fmeth=mmlpars.mml_pars(fmeth,default='all',list=['all']+LIST_FILEMETHODS)
    fmeth=fmeth.lower()
    # Create dictionary containing general keys
    rmsrc=mmlio.yorn('[mvfiles] Remove source files?')
    genkeys=dict(overwrite=overwrite,rmsrc=rmsrc,memcomp1=memcomp1,memcomp2=memcomp2,
                 snapend=snapend,typList=typList,**options)
    # Loop over list
    if fmeth=='all': fmethlist=LIST_FILEMETHODS
    else           : fmethlist=[fmeth]
    for ifmeth in fmethlist:
        if mmlio.yorn('Move {} files?'.format(ifmeth)):
            mvfiles_fmeth(isimstr,fsimstr,ifmeth,'all',**genkeys)
    return

####################################################################################################################################
# METHOD FOR MOVING/RENAMING RUN FILES BY TYPE
def mvfiles_fmeth(simstr1,simstr2,fmeth,ftype,overwrite=None,rmsrc=None,snapend=None,typList=None,
                  memcomp1=None,memcomp2=None,newstyletag1=None,newstyletag2=None,**extra_kw):
    """
    Method to move/rename run files
    """
    # Set constants
    typListDEF=simlist.LIST_PTYPS
    # Pars input
    fmeth=mmlpars.mml_pars(fmeth,type=str)
    ftype=mmlpars.mml_pars(ftype,type=str)
    memcomp1=mmlpars.mml_pars(memcomp1,default=simstr1['memcomp'],type=str)
    memcomp2=mmlpars.mml_pars(memcomp2,default=simstr2['memcomp'],type=str)
    overwrite=mmlpars.mml_pars(overwrite,default=True,type=bool)
    fmeth=fmeth.lower()
    ftype=ftype.lower()
    memcomp1=memcomp1.lower()
    memcomp2=memcomp2.lower()
    # Create file keys
    fkeys_beg=dict(ftype=ftype,ptype=None,memcomp=memcomp1,checkaccess=True,newstyletag=newstyletag1)
    fkeys_end=dict(ftype=ftype,ptype=None,memcomp=memcomp2,checkaccess=True,newstyletag=newstyletag2)
    # Get file dictionary
    files_beg=simstr1.mkfiledict(**fkeys_beg)[fmeth]
    files_end=simstr2.mkfiledict(**fkeys_end)[fmeth]
    filelist_beg=simstr1.get_filelist(fmeth,**fkeys_beg)
    filelist_end=simstr2.get_filelist(fmeth,**fkeys_end)
    # Get source/destination hosts
    srchost = mmlinfo.computers('CFHMAC')['host'] if memcomp1=='cfhmac' else ''
    dsthost = mmlinfo.computers('CFHMAC')['host'] if memcomp2=='cfhmac' else ''
    # Set fmeth dependent things
    longlist=list_filetypes(fmeth=fmeth,flength='long' )
    shrtlist=list_filetypes(fmeth=fmeth,flength='short')
    if   fmeth=='buildgal': pass
    if   fmeth=='gadget'  : snapend=mmlpars.mml_pars(snapend,default='',type=str)
    if   fmeth=='scf'     : typList=mmlpars.mml_pars(typList,default=typListDEF,type=list)
    elif fmeth=='galfit'  : typList=mmlpars.mml_pars(typList,default=typListDEF,type=list)
    else                  : typList=['all']
    # Loop over particle type
    flag_static=True
    for iptyp in typList:
        fkeys_beg['ptype']=iptyp
        fkeys_end['ptype']=iptyp
        # Ask if files should be moved for particle type (scf only)
        if fmeth=='scf':
            if not mmlio.yorn('Move file(s) for {} particles?'.format(iptyp)): continue
        # Select files
        if fmeth=='scf':
            pfilelist_beg=simstr1.get_filelist(fmeth,**fkeys_beg)
            pfilelist_end=simstr2.get_filelist(fmeth,**fkeys_end)
            pfiles_beg=files_beg[iptyp]
            pfiles_end=files_end[iptyp]
        else:
            pfilelist_beg=filelist_beg
            pfilelist_end=filelist_end
            pfiles_beg=files_beg
            pfiles_end=files_end
        # Loop over list of file subtypes
        for iftyp in pfilelist_beg['keylist']:
            if not memcomp1=='cfhmac':
                iexistlist=[]
                for ipattern in pfilelist_beg[iftyp]: iexistlist+=glob.glob(ipattern)
                if len(iexistlist)==0:
                    print '{} {} file(s) for {} particles do not exist. Skipping.'.format(iftyp,fmeth,iptyp)
                    continue
            if not mmlio.yorn('Move {} {} file(s) for {} particles?'.format(iftyp,fmeth,iptyp)): continue
            # Create lists of source/destination files
            isrclist=[] ; idstlist=[]
            # Handle multiple files w/ same base name
            if iftyp in longlist:
                if memcomp1=='cfhmac': raise Exception('Not currently supported.')
                if memcomp2=='cfhmac': raise Exception('Not recommended.')
                iowflag=mmlio.yorn('Overwrite existing {} {} file(s) for {} particles?'.format(iftyp,fmeth,iptyp))
                for ifile_beg,ifile_end in zip(pfilelist_beg[iftyp],pfilelist_end[iftyp]):
                    isrclist_new=sorted(glob.glob(ifile_beg))
                    iextlist=[]
                    # Loop over source files getting extensions
                    for isrc in isrclist_new:
                        iext=re.findall(ifile_beg.replace('*','(.*)'),isrc)[0]
                        iextlist.append(iext)
                        isrclist.append(isrc)
                        idstlist.append(ifile_end.replace('*',iext))
                    # Only move SCF input files that don't have output
                    if fmeth=='scf' and not iowflag and ftype=='longin':
                        for iext in iextlist:
                            if os.path.isfile(pfiles_beg['output']['snapbase']+iext):
                                isrclist.remove(ifile_beg.replace('*',iext))
                                idstlist.remove(ifile_end.replace('*'.iext))
                    # Only move GADGET snapshots that end in snapend
                    if fmeth=='gadget' and iftyp=='snapshot' and len(snapend)>0:
                        for iext in iextlist:
                            if not iext.endswith(snapend):
                                isrclist.remove(ifile_beg.replace('*',iext))
                                idstlist.remove(ifile_end.replace('*',iext))
            elif iftyp in shrtlist:
                iowflag=overwrite
                if memcomp1=='cfhmac':
                    isrcext=mmlio.askquest('Enter extension for file to transfer',default='_000ic',type='str')
                for ifile_beg,ifile_end in zip(pfilelist_beg[iftyp],pfilelist_end[iftyp]):
                    iextlist=[]
                    if memcomp1=='cfhmac': isrc=ifile_beg.replace('*',isrcext)
                    else                 : isrc=sorted(glob.glob(ifile_beg))[-1]
                    iext=re.findall(ifile_beg.replace('*','(.*)'),isrc)[0]
                    iextlist.append(iext)
                    isrclist.append(isrc)
                    idstlist.append(ifile_end.replace('*',iext))
            # Handle single files
            else:
                iowflag=overwrite
                # Handle static file
                if iftyp=='static':
                    if not flag_static: continue
                    flag_static=False
                # Add files to list
                for ifile_beg,ifile_end in zip(pfilelist_beg[iftyp],pfilelist_end[iftyp]):
                    if not memcomp1=='cfhmac':
                        if not os.path.isfile(ifile_beg): continue
                    isrclist.append(ifile_beg)
                    idstlist.append(ifile_end)
            # Move files
            print '({},{},{})'.format(fmeth,iptyp,iftyp)
            for isrc,idst in zip(isrclist,idstlist):
                if isrc==idst: continue
                if not memcomp2=='cfhmac': mmlfiles.mkdirs(os.path.dirname(idst))
                print '{} ---> {}'.format(isrc,idst)
            mmlfiles.cpfiles(isrclist,idstlist,srchost=srchost,dsthost=dsthost,
                             overwrite=iowflag,move=rmsrc)
    return

####################################################################################################################################
# METHOD TO MOVE FILES TO/FROM ACCRE
def mvfiles_accre(simstr,method,fmeth,ftype,overwrite=None,snapend=None,typList=None):
    """
    Moves files to/from ACCRE
    """
    # Pars input
    method=mmlpars.mml_pars(method,list=['to','from'])
    # Select source/destination
    simstr1=copy.deepcopy(simstr)
    simstr2=copy.deepcopy(simstr)
    if method=='to':
        memcomp1=simstr1['memcomp']
        memcomp2='accre'
    elif method=='from':
        memcomp1='accre'
        memcomp2=simstr2['memcomp']
    else: raise Exception('Keyword method ({}) must be "TO" or "FROM".'.format(method))
    # Call function for generic file move
    mvfiles_fmeth(simstr1,simstr2,fmeth,ftype,overwrite=overwrite,
                  snapend=snapend,typList=typList,
                  memcomp1=memcomp1,memcomp2=memcomp2)
    # Return control
    return

####################################################################################################################################
# METHOD TO MOVE FILES TO/FROM CFHMAC
def mvfiles_cfhmac(simstr,method,fmeth,ftype,overwrite=None,snapend=None,typList=None):
    """
    Moves files to/from Meagan's Mac
    """
    # Pars input
    method=mmlpars.mml_pars(method,list=['to','from'])
    # Select source/destination
    simstr1=copy.deepcopy(simstr)
    simstr2=copy.deepcopy(simstr)
    if method=='to':
        memcomp1=simstr['memcomp']
        memcomp2='cfhmac'
    elif method=='from':
        memcomp1='cfhmac'
        memcomp2=simstr['memcomp']
    else: raise Exception('Keyword method ({}) must be "TO" or "FROM".'.format(method))
    # Call function for generic file move
    mvfiles_fmeth(simstr1,simstr2,fmeth,ftype,overwrite=overwrite,
                  snapend=snapend,typList=typList,
                  memcomp1=memcomp1,memcomp2=memcomp2)
    # Return control
    return

####################################################################################################################################
# METHOD TO UPDATE FILES FOR A NEW NUMBER OF PROCESSORS
def set_procnum(simstr):
    """
    Renames files and changes number of processors
    """
    # Ask user for info
    procid=mmlio.askselect('Which processory number should be changed?',['buildgal','gadget','scf'])
    if   procid=='buildgal': nproc_old=simstr['nprocbg' ]
    elif procid=='gadget'  : nproc_old=simstr['nprocgd' ]
    elif procid=='scf'     : nproc_old=simstr['nprocscf']
    nproc_new=mmlio.askquest('The # of {} processors is currently {}. What should it be changed to?'.format(procid,nproc_old),
                             default=nproc_old,dtype='int')
    incltag=mmlio.yorn('Include # of processors in file names?')
    # Stop if the number of processors is the same
    if nproc_new==nproc_old:
        print 'The number of processors is the same. Exiting.'
        return
    # Set new and old mmlsim's
    oldmmlsim=simstr
    newmmlsim=copy.deepcopy(simstr)
    # Proceed based on procid
    if   procid=='buildgal':
        newmmlsim['nprocbg'   ]=nproc_new
        newmmlsim['inclbgtag' ]=incltag
    elif procid=='gadget'  :
        newmmlsim['nprocgd'   ]=nproc_new
        newmmlsim['inclgdtag' ]=incltag
    elif procid=='scf'     :
        newmmlsim['nprocscf'  ]=nproc_new
        newmmlsim['inclscftag']=incltag
    # Rename original files
    print 'Renaming files...'
    mvfiles(oldmmlsim,newmmlsim,ftype='all')
    # Remake files that are nproc dependent
    print 'Remaking files that depend on # of processors...'
    if   procid=='buildgal': mmlbuildgal.mkfiles(newmmlsim)
    elif procid=='gadget'  : mmlgadget.mkfiles(newmmlsim)
    elif procid=='scf'     : mmlscf.mkfiles(newmmlsim)
    return

####################################################################################################################################
####################################################################################################################################
# METHODS FOR CLEANING RUN FILES/DIRECTORIES

####################################################################################################################################
# METHOD TO DELETE SIMULATION FILES
def rmfiles(fmeth,ftype=None,typList=None,**input_kw):
    """
    Deletes files
    """
    # Set constants
    typListDEF=simlist.LIST_PTYPS
    # Pars input
    fmeth=mmlpars.mml_pars(fmeth,list=LIST_FILEMETHODS)
    ftype=mmlpars.mml_pars(ftype,type=str,default='clean')
    if   fmeth in ['scf']   : typList=mmlpars.mml_pars(typList,default=typListDEF,type=list)
    elif fmeth in ['galfit']: typList=mmlpars.mml_pars(typList,default=typListDEF,type=list)
    else                    : typList=['all']
    file_kw=pars_fileinput(**input_kw)
    memcomp=file_kw['memcomp']
    rundir=file_kw['rundir']
    runtag=file_kw['runtag']
    subtag=file_kw['subtag']
    # Double check location of removed files
    if memcomp not in ['accre','cfhmac']:
        if not mmlio.yorn('Are you sure you want to remove {} files from {}? ({})'.format(fmeth,memcomp,runtag)): return
    # Determine which files to remove
    rmlist=[]
    for iptyp in typList:
        # Get list of files
        pfilelist=get_filelist(fmeth,ftype=ftype,ptype=iptyp,**file_kw)
        for iftyp in pfilelist['keylist']:
            # Only move static files once
            if iftyp=='static' and iptyp==typList[0]:
                if mmlio.yorn('Remove {} {} {} file(s)? ({})'.format(memcomp,fmeth,iftyp,runtag)):
                    rmlist+=pfilelist[iftyp]
            # Handle other files
            else:
                if len(glob.glob(pfilelist[iftyp][0]))>0:
                    if mmlio.yorn('Remove {} {} {} file(s)? ({})'.format(memcomp,fmeth,iftyp,runtag)):
                        rmlist+=pfilelist[iftyp]
    # Loop over files removing them
    if not mmlio.yorn('Remove?'): return
    for itemp in rmlist:
        for ifile in glob.iglob(itemp):
            os.remove(ifile)
    # Return control
    return

####################################################################################################################################
# METHOD FOR CLEANING RUN TREE
def cleantree(simstr=None,fmeth=None,typList=None,**options):
    """
    Removes empty directories from the simulation tree
    """
    # Pars input
    fmeth=mmlpars.mml_pars(fmeth,default='all',type=str)
    fmeth=fmeth.lower()
    # Create dictionary containing general keys
    genkeys=dict(simstr=simstr,**options)
    # Add method dependent keys
    if fmeth in ['all','scf'     ]: genkeys['typList']=typList
    # Handle invalid file type
    if fmeth not in ['all']+LIST_FILEMETHODS:
        raise Exception('Invalid file type: {}'.format(fmeth))
    # Loop over file methods
    if fmeth=='all': fmethlist=LIST_FILEMETHODS
    else           : fmethlist=[fmeth]
    for ifmeth in fmethlist:
        cleantree_fmeth(ifmeth,'all',**genkeys)
    return

####################################################################################################################################
# METHOD TO REMOVE EMPTY TREE DIRECTORIES
def cleantree_fmeth(fmeth,ftype,typList=None,**input_kw):
    """
    Method to remove empty tree directories
    """
    # Set constants
    typListDEF=simlist.LIST_PTYPS
    # Pars input
    fmeth=mmlpars.mml_pars(fmeth,list=LIST_FILEMETHODS)
    ftype=mmlpars.mml_pars(ftype,list=['all']+list_filetypes(fmeth))
    fkeys=pars_fileinput(**input_kw)
    # Get file dictionary
    filedict=files(**fkeys)[fmeth]
    filelist=get_filelist(fmeth,ftype=ftype,**fkeys)
    # Set method dependent things
    if   fmeth=='scf'   : typList=mmlpars.mml_pars(typList,default=typListDEF,type=list)
    elif fmeth=='galfit': typList=mmlpars.mml_pars(typList,default=typListDEF,type=list)
    else                : typList=['all']
    # Loop over particle type
    flag_static=True
    for iptyp in typList:
        # Select files
        if fmeth=='scf':
            pfilelist=get_filelist(fmeth,ftype=ftype,ptype=iptyp,**fkeys)
            pfiledict=filedict[iptyp]
        else:
            pfilelist=filelist
            pfiledict=filedict
        # Loop over list of file subtypes
        for iftyp in pfilelist['keylist']:
            for ifile in pfilelist[iftyp]:
                idir=os.path.dirname(ifile)
                if os.path.isdir(idir) and len(glob.glob(os.path.join(idir,'*')))==0:
                    print '{} directory empty. Removing: {}'.format(iftyp,idir)
                    os.removedirs(idir)
                else:
                    pass
    return

####################################################################################################################################
# METHOD TO COMPRESS A RUN
def ziprun(**options):
    """
    Compresses run directories.
    """
    options['checkaccess']=True
    fdict=files(**options)
    zipcmd="tar -zcvf {} {} --remove-files".format(fdict['archiv'],fdict['rundir'])
    print zipcmd
    os.system(zipcmd)
    shutil.rmtree(fdict['rundir'])
    return

####################################################################################################################################
# METHOD TO UNCOMPRESS A RUN
def unziprun(**options):
    """
    Uncompress a run directory.
    """  
    options['checkaccess']=True
    fdict=files(**options)
    zipcmd="tar -zxvf {} -C {}".format(fdict['archiv'],fdict['simdir'])
    print zipcmd
    os.system(zipcmd)
    os.remove(fdict['archiv'])
    return

####################################################################################################################################
####################################################################################################################################
# SUPPORTING METHOD FOR GENERATING FILE NAMES FROM MMLSIM OBJECTS

####################################################################################################################################
# METHOD TO RETURN FILE PREFIX
def get_prefix(runtag,subtag):
    """
    Returns string prefixes for files
    """
    if len(runtag)==0:
        pfix_shrt=''
        pfix=subtag
    else:
        pfix_shrt=runtag+'.'
        pfix=runtag+subtag+'.'
    return pfix,pfix_shrt

####################################################################################################################################
# METHOD TO RETURN INPUT FOR FILE METHODS
def pars_fileinput(simstr=None,rundir=None,runtag=None,subtag=None,
                   topdir=None,memcomp=None,checkaccess=None,
                   simdir=None,runtyp=None,
                   newstyletag=None,
                   inclbgtag =None,nprocbg =None,
                   inclgdtag =None,nprocgd =None,
                   inclscftag=None,nprocscf=None,
                   **extra_kw):
    """
    Returns a dictionary of default file input variables
    """
    # Get memcomp from rundir
    if isinstance(rundir,str):
        memcompDEF=mmlinfo.dir2compid(rundir)
    else:
        memcompDEF=mmlinfo.hostname2compid()
    # Get runtyp from runtag
    if isinstance(runtag,str):
        try:
            runtypDEF=runtag2simtype(runtag)
        except:
            runtypDEF=''
    else:
        runtypDEF=''
    # Pars simstr
    if simstr==None:
        simstr={'memcomp'    :memcompDEF,
                'runtag':'','runtyp':runtypDEF,'subtag':'',
                'newstyletag':True ,
                'inclbgtag'  :False,'nprocbg' :0,
                'inclgdtag'  :False,'nprocgd' :0,
                'inclscftag' :False,'nprocscf':0}
    memcomp=mmlpars.mml_pars(memcomp,default=simstr['memcomp'],type=str)
    runtag=mmlpars.mml_pars(runtag,default=simstr['runtag'],type=str)
    runtyp=mmlpars.mml_pars(runtyp,default=simstr['runtyp'],type=str)
    subtag=mmlpars.mml_pars(subtag,default=simstr['subtag'],type=str)
    newstyletag=mmlpars.mml_pars(newstyletag,default=simstr['newstyletag'],type=bool)
    inclbgtag =mmlpars.mml_pars(inclbgtag ,default=simstr['inclbgtag' ],type=bool)
    inclgdtag =mmlpars.mml_pars(inclgdtag ,default=simstr['inclgdtag' ],type=bool)
    inclscftag=mmlpars.mml_pars(inclscftag,default=simstr['inclscftag'],type=bool)
    nprocbg =mmlpars.mml_pars(nprocbg ,default=simstr['nprocbg' ],type=int,min=0)
    nprocgd =mmlpars.mml_pars(nprocgd ,default=simstr['nprocgd' ],type=int,min=0)
    nprocscf=mmlpars.mml_pars(nprocscf,default=simstr['nprocscf'],type=int,min=0)
    # Get directories
    topdirDEF=get_topdir(memcomp=memcomp,checkaccess=checkaccess)
    topdir=mmlpars.mml_pars(topdir,default=topdirDEF,type=str)
    simdirDEF=get_simdir(runtyp=runtyp,topdir=topdir)
    simdir=mmlpars.mml_pars(simdir,default=simdirDEF,type=str)
    rundirDEF=get_rundir(runtag=runtag,simdir=simdir)
    rundir=mmlpars.mml_pars(rundir,default=rundirDEF,type=str)
    # Get added tag
    subtagproc=get_subtag(subtag=subtag,newstyletag=newstyletag,
                          inclbgtag =inclbgtag ,nprocbg =nprocbg ,
                          inclgdtag =inclgdtag ,nprocgd =nprocgd ,
                          inclscftag=inclscftag,nprocscf=nprocscf)
    # Get prefixes
    pfix,pfix_shrt=get_prefix(runtag,subtagproc)
    # Get other file stuff
    mkinfo=os.path.join(rundir,runtag+'.mkinfo')
    archiv=os.path.join(simdir,runtag+'.tar.gz')
    rapsht=fpar2fname('icgen','mmlsim',tagstr=runtag)
    # Create and return dictionary
    input_kw=dict(rundir=rundir,runtag=runtag,subtag=subtagproc,memcomp=memcomp,
                  topdir=topdir,simdir=simdir,newstyletag=newstyletag,
                  inclbgtag =inclbgtag ,nprocbg =nprocbg ,
                  inclgdtag =inclgdtag ,nprocgd =nprocgd ,
                  inclscftag=inclscftag,nprocscf=nprocscf,
                  pfix=pfix,pfix_shrt=pfix_shrt,
                  mkinfo=mkinfo,archiv=archiv,rapsht=rapsht)
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
def get_simdir(runtyp=None,topdir=None,**topdir_kw):
    """
    Returns the properly assigned simulation type directory
    """
    # Pars input
    topdirDEF=get_topdir(**topdir_kw)
    topdir=mmlpars.mml_pars(topdir,default=topdirDEF,type=str)
    runtyp=mmlpars.mml_pars(runtyp,default='',type=str)
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

####################################################################################################################################
# METHOD TO RETURN DEFAULT ADDED RUN TAG
def get_subtag(**simstr):
    """
    Returns the properly assigned ID tag to append to file names
    """
    # Pars input
    # Create subrunid including # of processors and subtag
    if simstr['newstyletag']:
        subtag=""
        if simstr['inclbgtag' ]: subtag+="_BG%(nprocbg)d"   % simstr
        if simstr['inclgdtag' ]: subtag+="_GD%(nprocgd)d"   % simstr
        if simstr['inclscftag']: subtag+="_SCF%(nprocscf)d" % simstr
    else:
#        subtag="_ON%(nprocgd)d%(subtag)s" % simstr
        subtag="_ON%(nprocgd)d" % simstr
    # Return ID string
    return subtag+simstr['subtag']

####################################################################################################################################
# METHOD TO GENERATE DEFAULT RUNTAG
def get_runtag(simtype,**tagkw):
    """
    Returns the properly assigned runtag for a simulation
    """
    # Pars input
    simtype=mmlpars.mml_pars(simtype,list=simlist.LIST_RUNTYPS)
    tagkw=parsfpar('icgen',simtype,fpar=tagkw,askuser=True)
    # Proceed based on simtype
    runtag=fpar2tag('icgen',simtype,tagkw)
    # Return runtag
    return runtag
def oldget_runtag_gal(newstyletag=None,addruntag=None,**tagkw):
    """
    Returns the propertly assigned runtag for a galaxy simulation
    """
    # Pars input
    newstyletag=mmlpars.mml_pars(newstyletag,default=False,type=bool)
    addruntag=mmlpars.mml_pars(addruntag,default='',type=str)
    tagkw=runsheet.askuser_mkinfo('galaxy',**tagkw)
    # Create short strings
    logm=mmlstring.decimal(np.log10(tagkw['mvir']),formstr='.4g')
    gmod=tagkw['model'].upper()
    logn=mmlstring.decimal(np.log10(tagkw['ntot']))
    gpro=tagkw['haloprof'].upper()
    conc=mmlstring.decimal(tagkw['concen'])
    fgas=mmlstring.decimal(100*tagkw['fgas'])
    # Combine and return runtag
    if newstyletag:
        runtag='{}_{}M{}{}{}FG{}'.format(gmod,logn,logm,gpro,conc,fgas)
    else:
        runtag='{}_{}M{}{}{}FG{}'.format(gmod,logn,logm,gpro,conc,fgas)
    # Return runtag
    return runtag+addruntag
def oldget_runtag_int(newstyletag=None,fromruntag=None,addruntag=None,**tagkw):
    """
    Returns the properly assigned runtag for an interaction simulation
    """
    # Set constants
    gproLIST=mmlbuildgal.LIST_HALOPROFS
    # Pars input
    newstyletag=mmlpars.mml_pars(newstyletag,default=False,type=bool)
    fromruntag=mmlpars.mml_pars(fromruntag,default=False,type=bool)
    addruntag=mmlpars.mml_pars(addruntag,default='',type=str)
    tagkw=runsheet.askuser_mkinfo('interact',**tagkw)
    # Create strings from runtag
    if fromruntag:
        # Create short run ID
        undsep1=tagkw['runtag1'].split('_') ; undsep2=tagkw['runtag2'].split('_')
        gmod1=undsep1[0]                    ; gmod2=undsep2[0]
        imstr1=undsep1[1]                   ; imstr2=undsep2[1]
        # Determine mass ratio
        for igpro in gproLIST:
            ivars1=imstr1.split(igpro)      ; ivars2=imstr2.split(igpro)
            imstr1=ivars1[0]                ; imstr2=ivars2[0]
        ivars1=imstr1.split('M')            ; ivars2=imstr2.split('M')
        instr1=ivars1[0]                    ; instr2=ivars2[0]
        imstr1=ivars1[1]                    ; imstr2=ivars2[1]
        logm1str=imstr1.split('p')          ; logm2str=imstr2.split('p')
        logm1=float(logm1str[0])            ; logm2=float(logm2str[0])
        if len(logm1str) == 2 : logm1+=float(logm1str[1])/100
        if len(logm2str) == 2 : logm2+=float(logm2str[1])/100
        ntot1 = instr1
        ntot2 = instr2
        qmas = mmlstring.decimal(10.0**(logm1-logm2))
    # Create strings from simdict
    else:
        simdict1=runsheet.get_runsheet(runsheet=tagkw['run1'])
        simdict2=runsheet.get_runsheet(runsheet=tagkw['run2'])
        # Create short strings
        gmod1 = simdict1['mkinfo']['model'].upper()
        gmod2 = simdict2['mkinfo']['model'].upper()
        ntot1 = mmlstring.decimal(np.log10(simdict1['mkinfo']['ntot']))
        ntot2 = mmlstring.decimal(np.log10(simdict2['mkinfo']['ntot']))
        qmas  = mmlstring.decimal(simdict1['mkinfo']['mvir']/simdict2['mkinfo']['mvir'])
    # Create other strings
    rper  = mmlstring.decimal(tagkw['rperi'])
    ecc   = mmlstring.decimal(tagkw['ecc'])
    incl  = mmlstring.decimal(tagkw['incl']*180./np.pi) # Angle in degrees
    # Combine to create runtag
    if newstyletag:
        runtag='{}_{}_Q{}RP{}E{}I{}'.format(gmod1+ntot1,gmod2+ntot2,qmas,rper,ecc,incl)
    else:
        runtag='{}_{}_Q{}RP{}E{}I{}'.format(gmod1+ntot1,gmod2+ntot2,qmas,rper,ecc,incl)
    # Return runtag
    return runtag+addruntag

