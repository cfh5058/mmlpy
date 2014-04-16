####################################################################################################################################
#
# MEAGAN LANG'S FILE METHODS
#
####################################################################################################################################
import mmlpars,mmlclass,mmlio
import os,shutil,subprocess,types,copy,itertools,glob,importlib
DIR_PARAM=os.path.join(os.path.expanduser('~'),'pypar')#'/home/langmm/pypar/'

####################################################################################################################################
# METHOD FOR LISTING PARAMETERS
def listpar(module,partyp,default=None,**exkw):
    """
    Returns a dictionary of required parameter lists
    """
    if isinstance(module,str): module=importlib.import_module(module)
    partyp=mmlpars.mml_pars(partyp,type=str)
    default=mmlpars.mml_pars(default,type=bool,default=False)
    fmeth=par_fmeth(module)
    par_list=module.listpar(partyp,default=default,**exkw)
    if not default:
        if 'keylist' in par_list:
            keylist=copy.deepcopy(par_list['keylist'])
            for ikey in par_list.keys():
                if ikey=='keylist': continue
                if ikey not in keylist:
                    mmlio.verbose('Key {} not in {} {} parameter keylist.'.format(ikey,fmeth,partyp))
                    par_list['keylist'].append(ikey)
            for ikey in keylist:
                if ikey not in par_list:
                    mmlio.verbose('Key {} in {} {} parameter keylist, but not dictionary.'.format(ikey,fmeth,partyp))
                    par_list['keylist'].remove(ikey)
        else: par_list['keylist']=sorted(par_list.keys())
        if 'tagstr' in par_list['keylist']:
            raise Exception('{} {}: tagstr is a reserved parameter key.'.format(fmeth,partyp))
    keylist=par_list.get('keylist',par_list.keys())
    if len(keylist)>0:
        key0=keylist[0]
        if isinstance(par_list[key0],dict) and par_list[key0].has_key('def'):
            if default:
                for ikey in keylist: par_list[ikey]=par_list[ikey]['def']
            else:
                for ikey in keylist: par_list[ikey]=par_list[ikey]['form']
    return par_list

####################################################################################################################################
# METHOD FOR ITERATING OVER PARAMETERS
def iterpar(module,partyp,inpar=None,**exkw):
    """
    Returns a list of dictionaries iterating over parameters
    """
    if inpar==None: inpar=listpar(module,partyp,**exkw)
    pnam=listpar(module,partyp,**exkw)['keylist']
    if len(pnam)==0: return []
    pval=[inpar[inam] for inam in pnam]
    piter=itertools.product(*pval)
    iterpar=[]
    for iplist in piter:
        iinpar={ipnam:ipval for ipnam,ipval in zip(pnam,iplist)}
        iterpar.append(iinpar)
    return iterpar

####################################################################################################################################
# METHOD TO TURN PARAMETERS INTO A TAG
def par2tag(module,partyp,inpar,**exkw):
    """
    Returns a tag for a given set of parameters and method
    """
    if isinstance(module,str): module=importlib.import_module(module)
    partyp=mmlpars.mml_pars(partyp,type=str)
    inpar=mmlpars.mml_pars(inpar,type=dict,default={})
    if 'tagstr' in inpar: tagstr=inpar['tagstr']
    else                : tagstr=module.par2tag(partyp,inpar,**exkw)
    return tagstr

####################################################################################################################################
# METHOD TO RETURN LIST OF FILES
def par_taglist(module,partyp):
    """
    Returns a list of existing parameter file tags
    """
    fdir=par_dirname(module,partyp)
    flist=sorted(glob.glob(fdir+'/*'))
    tlist=[]
    for ifile in flist: 
        if not ifile.endswith('~') and not ifile.endswith('.py') and not os.path.isdir(ifile): 
            tlist.append(fname2tag(module,partyp,ifile))
    return tlist

####################################################################################################################################
# METHOD TO GET DIRECTORY FOR PARAMETERS
def par_fpack(module):
    if isinstance(module,types.ModuleType): 
        pack=os.path.basename(os.path.dirname(module.__file__))
    else:
        pack=repr(module.__class__).split()[1].split('.')[0].strip("'>")
#        pack=repr(type(module)).split("'")[1].split('.')[0]
    return pack
def par_fmeth(module):
    if isinstance(module,types.ModuleType): 
        fname=os.path.splitext(os.path.basename(module.__file__))[0]
    else:
        fname=repr(module.__class__).split()[1].split('.')[-1].strip("'>")
#        print fname
#        fname=repr(type(module)).split("'")[1].split('.')[-1]
    return fname
def par_dirname(module,partyp):
    """
    Returns the directory for a given parameter type
    """
    if isinstance(module,str): module=importlib.import_module(module)
    partyp=mmlpars.mml_pars(partyp,type=str)
    fpack=par_fpack(module)
    fmeth=par_fmeth(module)
    fdir=os.path.join(DIR_PARAM,fpack,fmeth,partyp)
    return fdir

####################################################################################################################################
# METHOD TO RETURN FILE NAME FROM TAGSTR
def tag2fname(module,partyp,tagstr):
    """
    Converts tagstr to a file name
    """
    tagstr=mmlpars.mml_pars(tagstr,type=str)
    fdir=par_dirname(module,partyp)
    fname=os.path.join(fdir,tagstr)
    return fname

####################################################################################################################################
# METHOD TO RETURN TAGSTR FROM FILE NAME
def fname2tag(module,partyp,fname):
    """
    Converts file name to tagstr
    """
    tagstr=os.path.basename(fname)
    return tagstr

####################################################################################################################################
# METHOD TO RETURN FILE NAME FROM PARAMETERS
def par2fname(module,partyp,inpar=None,tagstr=None,**exkw):
    """
    Convers parameters to a file name
    """
    if not isinstance(tagstr,str): tagstr=par2tag(module,partyp,inpar,**exkw)
    if tagstr==os.path.basename(tagstr): fname=tag2fname(module,partyp,tagstr)
    else                               : fname=tagstr
    return fname

####################################################################################################################################
# METHOD TO SAVE PARAMETERS TO FILE
def savepar(module,partyp,inpar,tagstr=None,overwrite=None,verbose=True,**exkw):
    """
    Saves parameters to a file
    """
    outpar=parspar(module,partyp,inpar=inpar,**exkw)
    fname=par2fname(module,partyp,inpar=outpar,tagstr=tagstr,**exkw)
    overwrite=mmlpars.mml_pars(overwrite,default=False,type=bool)
    outpar['keylist']=listpar(module,partyp,**exkw)['keylist']
    if not os.path.isfile(fname) or overwrite:
        mmlio.rwdict('W',fname,outpar,overwrite=overwrite)
        if verbose:
            mmlio.verbose('Saved file parameters to file:')
            print '    '+fname
    return outpar

####################################################################################################################################
# METHOD TO LOAD PARAMETERS FROM A FILE
def loadpar(module,partyp,fname,**exkw):
    """
    Loads parameters from a file
    """
    fname=par2fname(module,partyp,tagstr=fname,**exkw)
    if not os.path.isfile(fname): raise Exception('Invalid file name: {}'.format(fname))
    inpar=mmlio.rwdict('R',fname)
    outpar=parspar(module,partyp,inpar=inpar,**exkw)
    outpar['tagstr']=fname2tag(module,partyp,fname)
    return outpar

####################################################################################################################################
# METHOD TO ASK USER FOR PARAMETERS
def parspar(module,partyp,inpar=None,tagstr=None,defpar=None,plot=None,overwrite=None,verbose=True,
            askuser=None,save=None,load=None,init=None,miss=None,inclextra=None,**extrakw):
    """
    Asks user for info required to create a histogram
    """
    # Set constants
    if isinstance(module,str): module=importlib.import_module(module)
    par_list=listpar(module,partyp,plot=plot,**extrakw)
    defpar=mmlpars.mml_pars(defpar,type=dict,default={})
    par_deft=listpar(module,partyp,plot=plot,default=True,**extrakw)
    par_deft.update(defpar)
    par_keys=par_list['keylist']
    fmeth=par_fmeth(module)
    # Pars input
    outpar=copy.deepcopy(mmlpars.mml_pars(inpar,default={},type=dict))
    askuser=mmlpars.mml_pars(askuser,default=False,type=bool)
    if askuser: loadDEF=True  ; saveDEF=True  ; initDEF=True
    else      : loadDEF=False ; saveDEF=False ; initDEF=False
    load=mmlpars.mml_pars(load,default=loadDEF,type=bool)
    save=mmlpars.mml_pars(save,default=saveDEF,type=bool)
    init=mmlpars.mml_pars(init,default=initDEF,type=bool)
    miss=mmlpars.mml_pars(miss,default=False  ,type=bool)
    if init or save: miss=False
    inclextra=mmlpars.mml_pars(inclextra,default=False,type=bool)
    # Get parameter dictionary from keys
    misskeys=[]
    for ikey in par_keys:
        if ikey in extrakw : outpar[ikey]=extrakw[ikey]
        if ikey not in outpar:
            misskeys.append(ikey)
            if   init: outpar[ikey]=None
            elif miss: pass
            else     : raise Exception('Parameter {} must be provided for {} param type {}.'.format(ikey,fmeth,partyp))
    # Get parameter dictionary from file
    if load and len(misskeys)>0:
        if not isinstance(tagstr,str) and askuser:
            tagstr=mmlio.askquest('Enter file containing {} {} parameters'.format(fmeth,partyp),default='None',dtype='str')
            if tagstr=='None': tagstr=None
        if isinstance(tagstr,str):
            if os.path.isfile(par2fname(module,partyp,tagstr=tagstr)):
                intpar=loadpar(module,partyp,tagstr,miss=miss)
                outpar.update(intpar)
            else: 
                if askuser: mmlio.verbose('File {} does not exist. Creating it...'.format(tagstr))
    # Ask for lists
    for ikey in par_keys:
        if miss and ikey in misskeys: continue
        if isinstance(par_list[ikey],list):
            if outpar[ikey] not in par_list[ikey]:
                keynone=(outpar[ikey] is None)
                if   askuser: outpar[ikey]=mmlio.askselect('Select valid list member for {}.'.format(ikey),par_list[ikey],default=par_deft[ikey])
                elif keynone: outpar[ikey]=par_deft[ikey]
                else: raise Exception('Value {} for parameter {} not in list {}'.format(outpar[ikey],ikey,par_list[ikey]))
                #else      : outpar[ikey]=mmlpars.mml_pars(outpar[ikey],list=par_list[ikey],default=par_deft[ikey])
        elif isinstance(par_list[ikey],type):
            if not isinstance(outpar[ikey],par_list[ikey]):
                keynone=(outpar[ikey] is None)
                if   askuser: outpar[ikey]=mmlio.askquest('Enter valid value for {}.'.format(ikey),dtype=str(par_list[ikey]),default=par_deft[ikey])
                elif keynone: outpar[ikey]=par_deft[ikey]
                else: raise Exception('Value {} for parameter {} is not specified type {}'.format(outpar[ikey],ikey,par_list[ikey]))
                #else      : outpar[ikey]=mmlpars.mml_pars(outpar[ikey],type=par_list[ikey],default=par_deft[ikey])
        else: raise Exception('Invalid {} parameter list type: {}'.format(ikey,type(par_list)))
    # Remove extra tags
    if not inclextra:
        for ikey in outpar.keys():
            if ikey not in par_keys+['tagstr']: del outpar[ikey]
    # Add tagstr
#    if isinstance(tagstr,str): outpar['tagstr']=os.path.basename(tagstr)
#    else                     : outpar['tagstr']=par2tag(module,partyp,outpar,plot=plot)
    # Save file
    if save: intpar=savepar(module,partyp,outpar,inclextra=inclextra,tagstr=tagstr,overwrite=overwrite,verbose=verbose)
    # Return output
    return outpar

####################################################################################################################################
# METHOD TO RETURN FILES FOR A GIVEN METHOD
def par2file(simstr,module,partyp,inpar=None,fext=None,test=None,plot=None,sing=None,
             ext=None,teststr0=None,plotstr0=None,animstr0=None,**exkw):
    """
    Returns a dictionary of file info
    """
    # Set constants
    import mmlplot
    plotopt=mmlplot.parsopt()
    plotext0=plotopt['plotext']
    animext0=plotopt['animext']
    # Pars input
    if isinstance(module,str): module=importlib.import_module(module)
    ftag=par2tag(module,partyp,inpar,plot=plot,**exkw)
    fextDEF='' if sing else '*'
    fext=mmlpars.mml_pars(fext,default=fextDEF)
    test=mmlpars.mml_pars(test,type=bool,default=False)
    plot=mmlpars.mml_pars(plot,type=bool,default=False)
    sing=mmlpars.mml_pars(sing,type=bool,default=False)
    ext=mmlpars.mml_pars(ext,type=str,default='')
    teststr0=mmlpars.mml_pars(teststr0,default='test_',type=str)
    plotstr0=mmlpars.mml_pars(plotstr0,default='plot',type=str)
    animstr0=mmlpars.mml_pars(animstr0,default='anim',type=str)
    fmeth=par_fmeth(module)
    # Handle flags
    if test: teststr=teststr0
    else   : teststr=''
    if plot: plotstr=plotstr0
    else   : plotstr=''
    if plot: ext=plotext0
    # Get file info
    try:
        fdict=simstr.fdict
        pfix=fdict['pfix']
        if fmeth in fdict:
            methdir=fdict[fmeth]['dir']
        else:
            methdir=os.path.join(fdict['rundir'],fmeth)
    except:
        fdict=simstr
        if 'pfix' in fdict:
            pfix=fdict['pfix']
            if fmeth in fdict:
                methdir=fdict[fmeth]['dir']
            elif 'dir' in fdict:
                methdir=fdict['dir']
            elif 'rundir' in fdict:
                methdir=os.path.join(fdict['rundir'],fmeth)
            else: raise
    # Initialize dictionary
    snfdict={fmeth+'dir':methdir}
    # Tag strings
    begtag='{}{}'.format(plotstr,partyp)
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
    if '{' in snfdict['base']: snfdict['file']=snfdict['base'].format(fext)
    else                     : snfdict['file']='{}{}'.format(snfdict['base'],fext)
    if len(ext)>0: snfdict['file']+='.{}'.format(ext)
    # Animation
    if plot:
        # Extensions
        snfdict['animext']=animext0
        snfdict['plotext']=plotext0
        # Directory
        snfdict['animdir']=os.path.join(methdir,animstr0+partyp)
        # File
        anim='{}{}{}{}{}.{}'.format(teststr,pfix,animstr0,partyp,endtag,animext0)
        snfdict['anim']=os.path.join(snfdict['animdir'],anim)
    # Return dictionary
    return snfdict
