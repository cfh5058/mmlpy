#!/usr/bin/env python
import inspect,importlib,os,glob,copy,re,subprocess,pprint
from mmlutils import *
from time import gmtime, strftime

STAMPEDE_ACCT='TG-PHY130001'

# Supporting methods and lists
DIR_SNAPFORM=os.path.join(os.path.dirname(inspect.stack()[0][1]),'formats')
def getmodule(name): return importlib.import_module(name)
def ftype2module(ftype): return getmodule('pysim.files.'+ftype)
def ftype2format(ftype): return os.path.join(DIR_SNAPFORM,ftype+'.pm')
def ftype2snaptyp(ftype,pNbody=False): 
    if pNbody: return 'Nbody_'+ftype
    else     : return ftype.capitalize()+'Snap'
def ftype2snapobj(ftype): return getattr(importlib.import_module('pysim.nbody.'+ftype),ftype2snaptyp(ftype))
def ftype2snapkey(ftype): return getattr(ftype2module(ftype),'VALIDSNAPKEY',False)
def ftype2initkey(ftype): return getattr(ftype2module(ftype),'VALIDINITKEY',False)
def ftype2snapext(ftype): return getattr(ftype2module(ftype),'VALIDSNAPEXT',['*'])
def ftype2nproc(ftype,sim):
    if   ftype in ['buildgal']: nproc=sim['nprocmk']
    elif ftype in ['gadget'  ]: nproc=sim['nprocev']
    elif ftype in ['scf'     ]: nproc=sim['nprocscf']
    else: raise Exception('Unsupported ftype: {}'.format(ftype))
    return nproc
def ftype2compid(ftype,sim):
    if   ftype in ['buildgal']: compid=sim['mkcomp']
    elif ftype in ['gadget'  ]: compid=sim['runcomp']
    elif ftype in ['scf'     ]: compid=sim['scfcomp']
    else: raise Exception('Unsupported ftype: {}'.format(ftype))
    return compid
LIST_SNAPFORM=[os.path.basename(ftype).split('.pm')[0] for ftype in sorted(glob.glob(ftype2format('*')))]
LIST_SNAPORDR=['gadget','scf','buildgal']
LIST_FILETYPES=LIST_SNAPFORM
# DICT_SNAPKEY={'order'   :['gadget','scf','buildgal'],
#               'buildgal':['ic'],
#               'gadget'  :['snapshot'],
#               'scf'     :['isnap']}
# DICT_INITKEY={'order'   :['gadget','scf','buildgal'],
#               'buildgal':['stat_ic'],
#               'gadget'  :['ic','stat_ic'],
#               'scf'     :[]}
TESTRUNS=['test5','test6','test7']

####################################################################################################################################
# METHOD FOR GETTING PERFORMANCE
# def performance(method,simstr=None,ftype='gadget',overwrite=False,**exkw):
#     """
#     Return performance
#     """
#     # Record new statistics
#     if method=='record':
#         # Get module
#         fmod=ftype2module(ftype)
#         # Handle module specific performance measure
#         if hasattr(fmod,'performance'): time,unit=fmod.performance(simstr,**exkw)
#         # Handle generic performance measure
#         else:
#             # Get output file name
#             fout=simstr.fdict[ftype]['output']['runout']
#             print '    '+fout 
#             # Get default units
#             if   ftype=='scf'     : unitDEF='s/snapshot'
#             elif ftype=='buildgal': unitDEF='s/ic'
#             else: raise Exception('Not default unit for file type: {}'.format(ftype))
#             unitDEF=exkw.get('unit',unitDEF)
#             # Ask for user input
#             time=mmlio.askquest('Enter {} execution time from above file.'.format(ftype),default=0.,dtype='float')
#             unit=mmlio.askquest('Enter {} units for time just entered.'.format(ftype),default=unitDEF,dtype='string')
#         # Return time and units
#         return time,unit

####################################################################################################################################
# METHOD FOR COUNTING SPLIT FILES
def get_nsplit(simstr,ftype,key,**exkw):
    """Determine if file is split. 0 means it is not."""
    fmodule=ftype2module(ftype)
    if len(fmodule.FILELIST_SPLIT)==0: nfile=0
    else: nfile=fmodule.get_nsplit(simstr,key,**exkw)
    return nfile
def get_splitext(simstr,ftype,key,**exkw):
    """Determine the extension for split files"""
    fmodule=ftype2module(ftype)
    nfile=get_nsplit(simstr,ftype,key,**exkw)
    if nfile: ext=fmodule.FILELIST_SPLIT[key]
    else    : ext=''
    return ext

####################################################################################################################################
# METHOD FOR INFERING FILE OPTIONS
def get_snapopt(simstr,ftype=None,fname=None,fext=None,ptype='disk',pgal=1):
    """
    Infer file options
    """
    # Raise error if fname & fext not provided
    if fname is None and fext is None:
        raise Exception('Either file name or extension must be provided.')
    # File types
    if ftype in LIST_SNAPFORM: ftypelist=[ftype]
    else:
        if ftype!=None: raise Exception('Unrecognized file type: {}'.format(ftype))
        ftypelist=LIST_SNAPORDR
    # Fill in missing file options
    for iftype in ftypelist:
        if ftype and fname and fext: break
        # Get list of possible keys
        if fext=='ic': fkeylist=[ftype2initkey(iftype)]
        else         : fkeylist=[ftype2snapkey(iftype),ftype2initkey(iftype)]
        # Loop over possible keys
        for ifkey in fkeylist:
            if not ifkey: continue
            if ftype and fname and fext: break
            iftemp=simstr.finfo[iftype][ifkey]['path']
            #mmlio.verbose(str((iftemp,fext,type(fext))))
            try:
                # Create a test file name based on fext
                if fext=='ic'            : ifname=iftemp.format(fext,ptype+str(pgal))
                elif isinstance(fext,str): ifname=add_ext(iftemp,fext.lstrip('_'))
                else                     : ifname=iftemp.format(fext,ptype+str(pgal))
                # Test the generated file name
                if fname:
                    iftest=(ifname!=fname)
                    if not iftest and iftype=='gadget': iftest=(ifname!=fname.split('.')[:-1].join('.'))
                else:
                    iftest=os.path.isfile(ifname)
                    if not iftest and iftype=='gadget': iftest=os.path.isfile(ifname+'.0')
                if not iftest: continue
                # Assign values
                fname=ifname
                ftype=iftype
                ftemp=iftemp
                fkey=ifkey
                # Fill in fext if it is missing
                if not fext: fext=get_ext(fname,ftemp)
                if not fext: fext=fkey
            except: continue
    # Create dictionary
    fopt=dict(ftype=ftype,fname=fname,fext=fext,ptype=ptype,pgal=pgal)
    # Raise error if file missing
    if not fname: 
        pprint.pprint(fopt)
        raise Exception('Could not find matching file')
    # Return dictionary of parameters
    return fopt

####################################################################################################################################
# METHOD FOR LISTING SNAPSHOTS
def listsnap(simstr,ftype,fext=None,host=None):
    """
    Returns a list of snapshots on memcomp
    """
    # Get correct module
    fmodule=ftype2module(ftype)
    # Get extension list
    if fext==None:
        fext=getattr(fmodule,'VALIDSNAPEXT',['*'])
    elif isinstance(fext,list): pass
    else: fext=[fext]
    lflag=('ic' in fext and len(fext)>1)
    # Loop over extensions
    templist=[] ; fkey=[] ; ftemp=[]
    for iext in fext:
        # Get correct key
        if iext=='ic': ikey=ftype2initkey(ftype)
        else         : ikey=ftype2snapkey(ftype)
        if not ikey: continue
        # Get template & extension format
        if not isinstance(host,str):
            itemp=simstr.finfo[ftype][ikey]['path']
        else:
            itemp=get_fdict(ret='info',compid=host,**simstr)[ftype][ikey]['path']
        extfmt=get_extfmt(itemp)
        # Add templates to list
        if ikey not in fkey:
            fkey.append(ikey)
            ftemp.append(itemp)
        if isinstance(iext,str): ifile=itemp.replace(extfmt,iext)
        else                   : ifile=itemp.format(iext)
        if ftype=='gadget' and len(glob.glob(ifile))==0: ifile+='.*'
        templist.append(ifile)
    if len(templist)==0: return []
    # List all files
    slist=sorted(mmlfiles.ls(templist,host=host))
    # Only return sets
    if len(fkey)==1:
        fkey =fkey[0]
        ftemp=ftemp[0]
    else: mmlio.verbose('Returning multiple values for fkey and ftemp. Watch behavior with list.')
    # Return
    return fkey,ftemp,slist

####################################################################################################################################
# METHOD FOR DETERMINING UNITS FILE
def unitsfile(fdict,ftype):
    """
    Returns the path of the unitsfile for a given file type
    """
    if ftype=='gadget': ufile=fdict[ftype]['input']['param']
    else              : ufile=fdict[ftype]['input']['units']
    return ufile

####################################################################################################################################
# METHOD FOR LISTING FILE PARAMETERS
def listpar(partyp,**exkw):
    """
    Returns a dictionary for parsing file parameters
    """
    partyp=mmlpars.mml_pars(partyp,type=str)
    if   partyp.startswith('mklog'):
        par={'keylist':['ftype','method','runtyp','runtag','addtag','baserun','mkinfo','ntot'],
             'ftype'    :{'def':'','form':str },
             'method'   :{'def':'','form':str },
             'runtyp'   :{'def':'','form':str },
             'runtag'   :{'def':'','form':str },
             'addtag'   :{'def':'','form':str },
             'baserun'  :{'def':'','form':str },
             'mkinfo'   :{'def':{},'form':dict},
             'ntot'     :{'def':0L,'form':long}}
        mkstr,ftype,method=partyp.split('_')
        module=ftype2module(ftype)
        parmod=module.listpar(mkstr+'_'+method,**exkw)
        par['keylist']+=parmod['keylist'] ; del parmod['keylist']
        par.update(**parmod)
    else: raise Exception('Invalid parameter type: {}'.format(partyp))
    return par

####################################################################################################################################
def get_fdata(ret=None,retall=None,**kwargs):
    """
    Returns file info in specified format
    """
    if isinstance(retall,bool):
        if retall: retDEF='all'
        else     : retDEF='dict'
    else: retDEF='dict'
    ret=mmlpars.mml_pars(ret,type=str,default=retDEF)
    # Get dict/info/tab
    if ret!='list': out=get_fdict(ret=ret,**kwargs)
    # Get list
    if ret in ['all','list']:
        flist={ftype:get_flist(ftype,**kwargs) for ftype in LIST_SNAPFORM}
        if ret=='all': out.append(flist)
        else         : out=flist
    # Return output
    return out

####################################################################################################################################
# METHOD TO RETURN SIMULATION FILES
def get_fdict(ftypes='all',ret=None,retall=None,**input_kw):
    """
    Generates a dictionary of simulation file names
    """
    # Pars input
    fdict=pars_fileinput(**input_kw)
    fdict['profpar']=os.path.join(fdict['rundir'],fdict['pfix_shrt']+'profpar')
    fdict['limits']=os.path.join(fdict['rundir'],fdict['pfix_shrt']+'limits')
    if isinstance(retall,bool):
        if retall: retDEF='all'
        else     : retDEF='dict'
    else: retDEF='dict'
    ret=mmlpars.mml_pars(ret,list=['all','dict','info','tab','dir'],default=retDEF)
    # Pars file types
    if isinstance(ftypes,str): ftypes=[ftypes]
    ftypes=mmlpars.mml_pars(ftypes,type=list,default=copy.deepcopy(LIST_SNAPFORM))
    if 'all' in ftypes: ftypes=copy.deepcopy(LIST_SNAPFORM)
    # Loop over file types
    ftab={} ; finfo={} ; dirlist={}
    for iprog in ftypes:
        # Get file format file
        ifpm=ftype2format(iprog)
        if not os.path.isfile(ifpm): 
            mmlio.verbose('Missing format file for {}'.format(iprog))
            continue
        # Initialize dictionaries
        ifdict={'dir':os.path.join(fdict['rundir'],iprog+fdict['subtag']),
                'pfix':fdict['pfix'],'pfix_shrt':fdict['pfix_shrt']}
        ifinfo={}
        idirlist=[]
        # Read in table 
        iftab=mmlio.rwtable('R',ifpm) ; iftab['path']=len(iftab['key'])*['None'] ; iftab['keylist'].append('path')
        tabkeys=copy.deepcopy(iftab['keylist']) ; tabkeys.remove('key')
        for iclass in set(iftab['class']): ifdict[iclass]={}
        # Add directories
        for i in range(len(iftab['key'])):
            # Parent directory
            if iftab['dir'][i]=='main': ipardir=ifdict['dir']
            else                      : ipardir=iftab['path'][iftab['key'].index(iftab['dir'][i])]
            # Add path to table
            if   iftab['flag'][i]=='dir'    : iftab['path'][i]=os.path.join(ipardir,iftab['name'][i])
            elif iftab['flag'][i]=='selfdir': iftab['path'][i]=os.path.join(ipardir,iftab['name'][i].split('{')[0],fdict['pfix']+iftab['name'][i])
            elif iftab['flag'][i]=='nopfix' : iftab['path'][i]=os.path.join(ipardir,iftab['name'][i])
            elif iftab['flag'][i]=='shpfix' : iftab['path'][i]=os.path.join(ipardir,fdict['pfix_shrt']+iftab['name'][i])
            else                            : iftab['path'][i]=os.path.join(ipardir,fdict['pfix']+iftab['name'][i])
            # Populate dictionary
            if   iftab['flag'][i]=='dir':
                if iftab['class'][i]==iftab['key'][i]: ifdict[iftab['class'][i]]['dir']=iftab['path'][i]
                else                                 : ifdict[iftab['class'][i]][iftab['key'][i]+'dir']=iftab['path'][i]
            else                        : ifdict[iftab['class'][i]][iftab['key'][i]]=iftab['path'][i]
            # Populate info
            ifinfo[iftab['key'][i]]={ikey:iftab[ikey][i] for ikey in tabkeys}
            # Populate direcotry info
            if   iftab['flag'][i]=='dir'    : idirlist.append(iftab['path'][i])
            elif iftab['flag'][i]=='selfdir': idirlist.append(os.path.dirname(iftab['path'][i]))
        # Add ith dictionary
        fdict[iprog]=ifdict ; finfo[iprog]=ifinfo ; ftab[iprog]=iftab ; dirlist[iprog]=idirlist
    if   ret=='all' : return fdict,finfo,ftab
    elif ret=='dict': return fdict
    elif ret=='info': return finfo
    elif ret=='tab' : return ftab
    elif ret=='dir' : return dirlist
    else: raise Exception('Unrecognized return option {}'.format(ret))

####################################################################################################################################
# METHOD TO RETURN LIST OF RELEVANT FILES
def get_flist(ftype,groups=None,**kwargs):
    """
    Returns a list of relavent files
    """
    # Pars input
    ftype=mmlpars.mml_pars(ftype,list=LIST_SNAPFORM)
    # Get files
    kwargs['retall']=True
    fdict,finfo,ftab=get_fdict(ftypes=[ftype],**kwargs)
    fdict,finfo,ftab=fdict[ftype],finfo[ftype],ftab[ftype]
    # Pars
    if isinstance(groups,str): groups=[groups]
    groups=mmlpars.mml_pars(groups,default=['all'],type=list)
    if 'all' in groups: groups=set(ftab['group'])
    # Initialize file list
    tabform={ikey:[] for ikey in ftab['keylist']}
    flist={igrp:[] for igrp in groups}
    flist['keylist']=copy.deepcopy(groups)
    # Loop over files
    for i in range(len(ftab['key'])):
        if ftab['group'][i] not in groups: continue
        flist[ftab['group'][i]].append({ikey:ftab[ikey][i] for ikey in ftab['keylist']})
    # Print info on missing groups
    for grp in groups:
        if len(flist[grp])==0: mmlio.verbose('No files found in group {}'.format(grp))
    # Return output
    return flist

####################################################################################################################################
# METHOD TO SETUP RUN FILES
def subrun(simstr,ftype=None,flag_mk=None,flag_mv=None,flag_run=None,flag_sub=False,**exkw):
    """
    Sets up files for run
    """
    # Get file type & module
    if ftype not in LIST_SNAPFORM:
        ftype=mmlio.askselect('Select a type of run to setup & submit.',LIST_SNAPFORM)
    fmodule=ftype2module(ftype)
    # Create files
    if len(fmodule.FILELIST_SETUP)>0:
        if not isinstance(flag_mk,bool): flag_mk=mmlio.yorn('Create new {} files?'.format(ftype))
        if flag_mk:
            mkfiles(simstr,ftype=ftype,groups=fmodule.FILELIST_SETUP,**exkw)
    # Move files
    if len(fmodule.FILELIST_SETUP)>0:
        if   ftype=='buildgal': dstcomp=simstr['mkcomp']
        elif ftype=='gadget'  : dstcomp=simstr['runcomp']
        elif ftype=='scf'     : dstcomp=simstr['scfcomp']
        if not isinstance(flag_mv,bool): flag_mv=mmlio.yorn('Move new {} files to {}?'.format(ftype,dstcomp))
        if flag_mv:
            extargs=dict(fmodule.INCLEXT_SETUP,**fmodule.EXCLEXT_SETUP)
            mvfiles(simstr,ftype=ftype,groups=fmodule.FILELIST_SETUP,rmsrc=False,
                    compid2=dstcomp,**dict(extargs,**exkw))
    # Run functions
    if len(fmodule.METHLIST_SETUP)>0:
        if not isinstance(flag_run,bool): flag_run=mmlio.yorn('Run {} pre-process methods?'.format(ftype))
        if flag_run:
            for func in fmodule.METHLIST_SETUP: getattr(fmodule,func)(simstr,**exkw)
    # Submit run
    if flag_sub:
        pass
    # Return control
    return

####################################################################################################################################
# METHOD TO RECOVER RUN FILES
def endrun(simstr,ftype=None,flag_rm=None,flag_mv=None,flag_run=None,empty=False,**exkw):
    """
    Retrieves files from run & performs basic post process
    """
    # Get file type & module
    if ftype not in LIST_SNAPFORM:
        ftype=mmlio.askselect('Select a type of file to create.',LIST_SNAPFORM)
    fmodule=ftype2module(ftype)
    # Get remote location of files
    if   ftype=='buildgal': srccomp=simstr['mkcomp']
    elif ftype=='gadget'  : srccomp=simstr['runcomp']
    elif ftype=='scf'     : srccomp=simstr['scfcomp']
    # Retrieve results
    if empty: flag_mv=True
    if len(fmodule.FILELIST_OUTPT):
        if not isinstance(flag_mv,bool): flag_mv=mmlio.yorn('Recover {} files from {}?'.format(ftype,srccomp))
        if flag_mv:
            extargs=dict(fmodule.INCLEXT_OUTPT,**fmodule.EXCLEXT_OUTPT)
            mvfiles(simstr,ftype=ftype,groups=fmodule.FILELIST_OUTPT,rmsrc=False,mvall=empty,
                    compid1=srccomp,**extargs)
    # Clean up files
    cleanrun(simstr,ftype=ftype,compid=srccomp,empty=empty)
    # if len(fmodule.FILELIST_CLEAN):
    #     if not isinstance(flag_rm,bool): flag_rm=mmlio.yorn('Clean up {} files on {}?'.format(ftype,srccomp))
    #     if flag_rm:
    #         rmfiles(ftype,groups=fmodule.FILELIST_CLEAN,compid=srccomp,**simstr)
    #     cleantree(ftype,compid=srccomp,**simstr)
    # Run functions
    if len(fmodule.METHLIST_OUTPT)>0:
        if not isinstance(flag_run,bool): flag_run=mmlio.yorn('Run {} post-process methods?'.format(ftype))
        if flag_run:
            for func in fmodule.METHLIST_OUTPT: getattr(fmodule,func)(simstr)
    # Return control
    return
        
####################################################################################################################################
# METHOD FOR CLEANING UP RUN FILES
def cleanrun(simstr,ftype=None,compid=None,empty=False):
    """
    Method to clean-up run files
    """
    compid=mmlpars.mml_pars(compid,default=simstr['memcomp'],type=str)
    # Get file type & module
    if ftype not in LIST_SNAPFORM:
        ftype=mmlio.askselect('Select a type of file to clean.',LIST_SNAPFORM)
    fmodule=ftype2module(ftype)
    # Clean up files & directories
    if empty: cleangrp=['all']
    else    : cleangrp=fmodule.FILELIST_CLEAN
    if len(cleangrp):
        # Files
        if mmlio.yorn('Clean up {} files on {}?'.format(ftype,compid)):
            rmfiles(ftype,groups=cleangrp,compid=compid,**simstr)
        # Directories
        if mmlio.yorn('Clean up {} directories on {}?'.format(ftype,compid)):
            cleantree(ftype,compid=compid,**simstr)
    # Run functions
    if len(fmodule.METHLIST_CLEAN)>0:
        if mmlio.yorn('Run {} clean up methods?'.format(ftype)):
            for func in fmodule.METHLIST_CLEAN: getattr(fmodule,func)(simstr,compid=compid)
    # Return control
    return

####################################################################################################################################
# METHOD FOR CREATING FILES
def mkfiles(simstr,ftype=None,groups=None,verbose=False,**kwargs):
    """
    Method to create files
    """
    # Get file type & module
    if ftype not in LIST_SNAPFORM:
        ftype=mmlio.askselect('Select a type of file to create.',LIST_SNAPFORM)
    fmodule=ftype2module(ftype)
    # Get info on file type
    flist=get_flist(ftype,groups=groups,**simstr)
    groups=flist['keylist']
    # Print info
    if verbose: mmlio.verbose('Creating files for {} on {}'.format(ftype.upper(),simstr['memcomp']))
    # Loop over groups
    runfunc=[]
    for group in groups:
        # Skip empty groups & askuser for confirmation
        if len(flist[group])==0: continue
        mkfuncs=list(set([iftab['mkfunc'] for iftab in flist[group]]))
        if len(mkfuncs)==1 and mkfuncs[0]=='None': continue
        if not mmlio.yorn('Create new local copies of {} files?'.format(group)): continue
        # Loop over files
        for iftab in flist[group]:
            # Break for directories, files w/o a mkfunc, & mkfuncs previously run
            if iftab['flag'  ]=='dir' : continue
            if iftab['mkfunc']=='None': continue
            if iftab['mkfunc'] in runfunc: continue
            # Run function
            mkfunc=getattr(fmodule,iftab['mkfunc'])
            iout=mkfunc(simstr,**kwargs)
            # Print info
            if iout:
                if   isinstance(iout,str ): print '    '+iout
                elif isinstance(iout,(list,tuple)): 
                    for iiout in iout: print '    '+iiout
                else: raise Exception('Invalid output type: {}'.format(type(iout)))
            # Add function to list of runs
            runfunc.append(iftab['mkfunc'])
            if 'timetest' not in simstr['runtag']: kwargs['owmake']=False
    # Return control
    return

def get_extfmt(path):
    """
    Returns the extension format for a path
    """
    extfmt=re.findall(r'({.*})',path)
    if   len(extfmt)==0: extfmt='{}'
    elif len(extfmt)==1: extfmt=extfmt[0]
    else: raise Exception('>1 format fields encountered: {}'.format(path))
    return extfmt

def get_ext(path,form):
    """
    Returns the extension given a form
    """
    extfmt=get_extfmt(form)
    ext=re.findall(form.replace(extfmt,r'(.*)'),path)
    if   len(ext)==0: ext=None
    elif len(ext)==1: ext=ext[0]
    else: raise Exception('>1 extensions fround in path: {}'.format(path))
    return ext

def add_ext(form,ext):
    """
    Replaces extension format with an extension
    """
    extfmt=get_extfmt(form)
    path=form.replace(extfmt,ext)
    return path

####################################################################################################################################
# METHOD FOR MOVING/RENAMING RUN FILES BY TYPE
def mvfiles(*args,**kwargs):
    """
    Method to move/rename run files
    """
    newmeth=kwargs.pop('newmeth',True)
    if newmeth: return newmvfiles(*args,**kwargs)
    else      : return oldmvfiles(*args,**kwargs)
    
def newmvfiles(simstr1,simstr2=None,ftype='all',groups='all',overwrite=True,rmsrc=None,
               inclext_end={},inclext_list={},inclext_fchk={},mvall=False,
               exclext_end={},exclext_list={},exclext_fchk={},
               compid1=None,compid2=None,**extra_kw):
    """
    Method to move/rename run files (uses rsync)
    """
    # Pars input
    inclext={'end':inclext_end,'list':inclext_list,'fchk':inclext_fchk}
    exclext={'end':exclext_end,'list':exclext_list,'fchk':exclext_fchk}
    fkeys_src=copy.deepcopy(mmlpars.mml_pars(simstr1))
    fkeys_dst=copy.deepcopy(mmlpars.mml_pars(simstr2,default=simstr1))
    ftype=mmlpars.mml_pars(ftype,list=copy.deepcopy(LIST_SNAPFORM))
    compid1=mmlpars.mml_pars(compid1,default=fkeys_src['memcomp'],type=str) ; fkeys_src['memcomp']=compid1
    compid2=mmlpars.mml_pars(compid2,default=fkeys_dst['memcomp'],type=str) ; fkeys_dst['memcomp']=compid2
    if not isinstance(overwrite,bool): overwrite=not mmlio.yorn('[mvfiles] Keep existing files?')
    if not isinstance(rmsrc,bool): rmsrc=mmlio.yorn('[mvfiles] Remove source files?')
    # Get source/destination hosts
    memlist = mmlinfo.computers()['memlist']
    srchost = '' if compid1 in memlist else compid1
    dsthost = '' if compid2 in memlist else compid2
    # Get file dictionary
    outsrc=get_fdict(ftypes=[ftype],retall=True,**fkeys_src)
    outdst=get_fdict(ftypes=[ftype],retall=True,**fkeys_dst)
    fdict_src,finfo_src,ftab_src=(iout[ftype] for iout in outsrc)
    fdict_dst,finfo_dst,ftab_dst=(iout[ftype] for iout in outdst)
    # srcdir=outsrc[0]['topdir']
    # dstdir=outdst[0]['topdir']
    srcdir=fdict_src['dir']
    dstdir=fdict_dst['dir']
    print srcdir
    print dstdir
    # Get file list
    flist_all=get_flist(ftype,**fkeys_src)
    flist_src=get_flist(ftype,groups=groups,**fkeys_src)
    flist_dst=get_flist(ftype,groups=groups,**fkeys_dst)
    # Get directory list
    dirlist=get_fdict(ftypes=[ftype],ret='dir',**fkeys_src)[ftype]
    if not srchost: mmlfiles.mkdirs(dirlist,host=srchost)
    # Loop over particle type
    exclfiles=[]
    inclfiles=[idir.split(srcdir)[-1] for idir in dirlist]
    for grp in flist_all['keylist']:
        # Exclude files not in listed groups or at user request
        inclgrp=True
        if grp not in flist_src: inclgrp=False
        elif not mvall: inclgrp=mmlio.yorn('Move {} {} file(s)?'.format(ftype,grp))
        if not inclgrp:
            for srcinfo in flist_all[grp]:
                # if srcinfo['flag']!='dir':
                srcbase=srcinfo['path'].split(srcdir)[-1]
                srcbase+=add_ext(get_splitext(simstr1,ftype,srcinfo['key']),'*')
                exclfiles.append(add_ext(srcbase,'*'))
            continue
        # Loop over files
        for srcinfo,dstinfo in zip(flist_src[grp],flist_dst[grp]):
            srclist=[] ; dstlist=[]
            # Get relative paths
            splitext=add_ext(get_splitext(simstr1,ftype,srcinfo['key']),'*')
            srcbase=srcinfo['path'].split(srcdir)[-1]+splitext
            dstbase=dstinfo['path'].split(dstdir)[-1]+splitext
            if srcbase!=dstbase:
                print srcbase
                print dstbase
                raise Exception('Source base does not match destination base for {} {} {}'.format(ftype,grp,srcinfo['key']))
            # Initialize restrictions on extensions
            iinclext=dict(end=False,list=False,fchk=False)
            iexclext=dict(end=False,list=False,fchk=False)
            for key in ['end','list','fchk']:
                iinclext[key]=inclext[key].get(srcinfo['key'],False)
                iexclext[key]=exclext[key].get(dstinfo['key'],False)
            if iinclext['fchk']: iinclext['fchk']=finfo_src[iinclext['fchk']]['path']
            if iexclext['fchk']: iexclext['fchk']=finfo_src[iexclext['fchk']]['path']
            # Exclude all files if include rules provided
            if iinclext['fchk'] or iinclext['list'] or iinclext['end']:
                exclfiles.append(add_ext(srcbase,'*'))
            # Restrict based on file check
            if iinclext['fchk']:
                chklist=mmlfiles.ls(add_ext(iinclext['fchk'],'*'),host=srchost)
                inclfiles+=[add_ext(srcbase,get_ext(ichk,iinclext['fchk'])) for ichk in chklist]
            if iexclext['fchk']:
                chklist=mmlfiles.ls(add_ext(iexclext['fchk'],'*'),host=srchost)
                exclfiles+=[add_ext(srcbase,get_ext(ichk,iexclext['fchk'])) for ichk in chklist]
            # Restrict based on list
            if iinclext['end']: inclfiles.append(add_ext(srcbase,'*'+iinclext['end']))
            if iexclext['end']: exclfiles.append(add_ext(srcbase,'*'+iexclext['end']))
            # Restrict based on end of extension
            if iinclext['list']: inclfiles+=[add_ext(srcbase,iext) for iext in iinclext['list']]
            if iexclext['list']: exclfiles+=[add_ext(srcbase,iext) for iext in iexclext['list']]
    # File names
    logdir=os.path.join(os.path.expanduser('~'),'transferlogs')
    mmlfiles.mkdirs(logdir)
    rsync_incl=os.path.join(logdir,'{}_{}_incl'.format(simstr1['runtag'],ftype))
    rsync_excl=os.path.join(logdir,'{}_{}_excl'.format(simstr1['runtag'],ftype))
    rsync_tlog=os.path.join(logdir,'{}_{}_tlog'.format(simstr1['runtag'],ftype))
    rsync_opt=['-avzu','--log-file={}'.format(rsync_tlog),'--itemize-changes']
    if not overwrite: rsync_opt.append('--ignore-existing')
#    rsync_opt.append('--dry-run')
    # Make list of files to include
    if os.path.isfile(rsync_incl): mmlfiles.rm(rsync_incl)
    if len(inclfiles)>0:
        rsync_opt.append('--include-from={}'.format(rsync_incl))
        f_incl=open(rsync_incl,'w')
        for f in inclfiles: f_incl.write("%s\n" % f)
        f_incl.close()
        mmlio.verbose(rsync_incl)
        #if not mmlio.yorn(rsync_incl): return
    # Make list of files to exclude
    if os.path.isfile(rsync_excl): mmlfiles.rm(rsync_excl)
    if len(exclfiles)>0:
        rsync_opt.append('--exclude-from={}'.format(rsync_excl))
        f_excl=open(rsync_excl,'w')
        for f in exclfiles: f_excl.write("%s\n" % f)
        f_excl.close()
        mmlio.verbose(rsync_excl)
        #if not mmlio.yorn(rsync_excl): return
    # Get directories including host names
    if srchost: rsync_srcdir=srchost+':'+srcdir+'/'
    else      : rsync_srcdir=srcdir+'/'
    if dsthost: rsync_dstdir=dsthost+':'+dstdir+'/'
    else      : rsync_dstdir=dstdir+'/'
    # Create command
    cpcmd=['rsync']+rsync_opt+[rsync_srcdir,rsync_dstdir]
    print ' '.join(cpcmd)
    print subprocess.call(cpcmd)
    mmlio.verbose(rsync_tlog)
    return

def oldmvfiles(simstr1,simstr2=None,ftype='all',groups='all',overwrite=None,rmsrc=None,
               inclext_end={},inclext_list={},inclext_fchk={},
               exclext_end={},exclext_list={},exclext_fchk={},
               compid1=None,compid2=None,**extra_kw):
    """
    Method to move/rename run files
    """
    # Pars input
    fkeys_src=copy.deepcopy(mmlpars.mml_pars(simstr1))
    fkeys_dst=copy.deepcopy(mmlpars.mml_pars(simstr2,default=simstr1))
    ftype=mmlpars.mml_pars(ftype,list=copy.deepcopy(LIST_SNAPFORM))
    compid1=mmlpars.mml_pars(compid1,default=fkeys_src['memcomp'],type=str) ; fkeys_src['memcomp']=compid1
    compid2=mmlpars.mml_pars(compid2,default=fkeys_dst['memcomp'],type=str) ; fkeys_dst['memcomp']=compid2
    if not isinstance(overwrite,bool): overwrite=mmlio.yorn('[mvfiles] Overwrite existing files?')
    if not isinstance(rmsrc,bool): rmsrc=mmlio.yorn('[mvfiles] Remove source files?')
    # Get source/destination hosts
    memlist = mmlinfo.computers()['memlist']
    srchost = '' if compid1 in memlist else compid1
    dsthost = '' if compid2 in memlist else compid2
    # Get file dictionary
    outsrc=get_fdict(ftypes=[ftype],retall=True,**fkeys_src)
    outdst=get_fdict(ftypes=[ftype],retall=True,**fkeys_dst)
    fdict_src,finfo_src,ftab_src=(iout[ftype] for iout in outsrc)
    fdict_dst,finfo_dst,ftab_dst=(iout[ftype] for iout in outdst)
    # Get file list
    flist_src=get_flist(ftype,groups=groups,**fkeys_src)
    flist_dst=get_flist(ftype,groups=groups,**fkeys_dst)
    # Loop over particle type
    srclist_tot=[] ; dstlist_tot=[]
    for grp in flist_src['keylist']:
        infolist=flist_src[grp]
        if not mmlio.yorn('Move {} {} file(s)?'.format(ftype,grp)): continue
        # Check for large files
        longgrp=('large' in [info['size'] for info in infolist])
        # Get overwrite options
        if longgrp: owflag=mmlio.yorn('Overwrite existing {} {} file(s)?'.format(ftype,grp))
        else      : owflag=overwrite
        # Loop over files
        for srcinfo,dstinfo in zip(flist_src[grp],flist_dst[grp]):
            srclist=[] ; dstlist=[]
            # Get file template
            extfmt=get_extfmt(srcinfo['path'])
            srctemp=srcinfo['path'].replace(extfmt,'*') ; srcfmt='('+srctemp+')\n'
            dsttemp=dstinfo['path'].replace(extfmt,'*') ; dstfmt='('+dsttemp+')\n'
            # Skip if no files match the template
            isrclist=mmlfiles.ls(srctemp,host=srchost)
            if len(isrclist)==0: continue
            idstlist=mmlfiles.ls(dsttemp,host=dsthost)
            # Create directory if path is a directory
            if srcinfo['flag']=='dir': mmlfiles.mkdir(dstinfo['path'],host=dsthost)
            # Groups with large lists of files
            elif srcinfo['size']=='large':
                # Initialize restrictions on extensions
                inclext=dict(end=False,list=False,fchk=False)
                if srcinfo['key'] in inclext_end : inclext['end' ]=inclext_end[srcinfo['key']]
                if srcinfo['key'] in inclext_list: inclext['list']=inclext_list[srcinfo['key']]
                if srcinfo['key'] in inclext_fchk: inclext['fchk']=finfo_src[inclext_fchk[srcinfo['key']]]['path']
                # Generate lists of source/destination files and extensions
                extlist=[]
                for src in isrclist:
                    # Get the extensions for this ith source
                    ext=re.findall(srcinfo['path'].replace(extfmt,r'(.*)'),src)
                    if   len(ext)==0: continue
                    elif len(ext)==1: ext=ext[0]
                    else: raise Exception('Found to many extensions: {}'.format(src))
                    # Pars restrictions on extensions
                    if inclext['fchk'] and not owflag and os.path.isfile(inclext['fchk'].replace(extfmt,ext)): continue
                    if inclext['end'] or inclext['list']:
                        extvalid=False
                        if inclext['end']  and ext.endswith(inclext['end']): extvalid=True
                        if inclext['list'] and ext in inclext['list']: extvalid=True
                        if not extvalid: continue
                    # Continue if not overwrite
                    dst=dstinfo['path'].replace(extfmt,ext)
                    if not owflag and dst in idstlist: continue
                    # Append lists
                    extlist.append(ext)
                    srclist.append(src)
                    dstlist.append(dst)
                # Check for existing output for SCF
                # inclext_fchk['snapin']='snapout'
                # Only move GADGET snapshots that end in snapend or in the list
                # inclext_end['snapshot']='5'
                # inclext_list['snapshot']='_000'
            # Groups with short lists of files
            else:
                if not owflag and srcinfo['path'] in idstlist: continue
                srclist.append(srcinfo['path'])
                dstlist.append(dstinfo['path'])
            # Add files to total list
            srclist_tot+=srclist
            dstlist_tot+=dstlist
    # Move files
    mmlfiles.cpfiles(srclist_tot,dstlist_tot,srchost=srchost,dsthost=dsthost,
                     overwrite=True,move=rmsrc,verbose=True)
    return

####################################################################################################################################
# METHOD TO DELETE SIMULATION FILES
def rmfiles(ftype,groups=[],**input_kw):
    """
    Deletes files
    """
    # Pars input
    ftype=mmlpars.mml_pars(ftype,list=LIST_SNAPFORM)
    file_kw=pars_fileinput(**input_kw)
    host=file_kw['compid']
    splitkeys=ftype2module(ftype).FILELIST_SPLIT
    # Double check location of removed files
    if host not in ['accre']:
        if not mmlio.yorn('Are you sure you want to remove {} files from {}?'.format(fmeth,host)): return
    # Get list of files
    flist=get_flist(ftype,groups=groups,**file_kw)
    # Loop over file groups
    rmlist=[] ; rmdirlist=[]
    for grp in flist['keylist']:
        infolist=flist[grp]
        rmgrp=None
        for info in infolist:
            temp=add_ext(info['path'],'*')
            filelist=mmlfiles.ls(temp,host=host)
            empty=(len(filelist)==0)
            if empty and info['key'] in splitkeys:
                temp+=add_ext(splitkeys[info['key']],'*')
                filelist=mmlfiles.ls(temp,host=host)
                ftest=(len(filelist)==0)
            if empty: continue
            if rmgrp==None: rmgrp=mmlio.yorn('Remove {} {} {} file(s)?'.format(host,ftype,grp))
            if rmgrp: 
                if os.path.isdir(temp): rmdirlist.append(temp)
                else                  : rmlist.append(temp)
    # Remove files
    pprint.pprint(rmlist)
    if not mmlio.yorn('Remove?'): return
    mmlfiles.rm(rmlist,host=host)
    # Return control
    return

####################################################################################################################################
# METHOD TO CLEAN UP EMPTY DIRECTORIES
def cleantree(ftype,**input_kw):
    """
    Method to remove empty tree directories
    """
    # Get list of files
    flist=get_flist(ftype,**input_kw)
    # Loop over groups
    for grp in flist['keylist']:
        # Loop over files
        for info in flist[grp]:
            idir=os.path.dirname(info['path'])
            # Remove dir if empty
            if os.path.isdir(idir) and len(glob.glob(os.path.join(idir,'*')))==0:
                mmlio.verbose('{} directory empty. Removing: {}'.format(grp,idir))
                os.removedirs(idir)
    return

####################################################################################################################################
# METHOD TO PARS FILE KEYWORDS
def pars_fileinput(compid=None,rundir=None,pfix=None,pfix_shrt=None,subtag=None,**kwargs):
    """
    Method to pars file input
    """
    # Identify computer
    memcomp=kwargs.get('memcomp',mmlinfo.dir2compid(rundir))
    compid=mmlpars.mml_pars(compid,default=memcomp,list=mmlinfo.LIST_COMPUTERS)
    runtyp=kwargs.get('runtyp','')
    runtag=kwargs.get('runtag','')
    # Run directory
    if isinstance(rundir,str): 
        simdir=os.path.dirname(rundir)
        topdir=os.path.dirname(simdir)
    else:
        topdir=mmlinfo.computers(memcomp=compid)['simdir']
        simdir=os.path.join(topdir,runtyp)
        rundir=os.path.join(simdir,runtag)
    # Sub tag
    if subtag is None:
        subtag=''
        if kwargs.get('inclmktag' ,False): subtag+= '_MK{}'.format(kwargs.get('nprocmk' ,0))
        if kwargs.get('inclevtag' ,False): subtag+= '_EV{}'.format(kwargs.get('nprocev' ,0))
        if kwargs.get('inclscftag',False): subtag+='_SCF{}'.format(kwargs.get('nprocscf',0))
        subtag+=kwargs.get('subtag','')
    # Prefixes
    pfix_long=pfix
    if pfix_long is None:
        if pfix_shrt is None: pfix_shrt=runtag
        pfix_long=pfix_shrt+subtag
    else:
        if pfix_shrt is None: pfix_shrt=pfix_long.split(subtag)[0]
    if len(pfix_shrt)>0 and not pfix_shrt.endswith('.'): pfix_shrt+='.'
    if len(pfix_long)>0 and not pfix_long.endswith('.'): pfix_long+='.'
    # Create dictionary
    fdict=dict(compid=compid,topdir=topdir,rundir=rundir,subtag=subtag,pfix=pfix_long,pfix_shrt=pfix_shrt,
               mkinfo=os.path.join(rundir,pfix_shrt+'mkinfo'),
               archiv=os.path.join(simdir,pfix_shrt+'tar.gz'),
               rapsht=mmlparam.par2fname('pysim.simlist','mmlsim',tagstr=runtag))
    # Return file parameters
    return fdict

####################################################################################################################################
# METHOD FOR READ/WRITE OF MAKE LOG
def rw_makelog(rwid,fname,indict=None,overwrite=False):
    """
    Reads/writes makelog files
    """
    # Pars input
    rwid=mmlpars.mml_pars(rwid,list=['R','W'])
    fname=mmlpars.mml_pars(fname,type=str)
    # Read
    if   rwid=='R':
        outdict=mmlio.rwdict('R',fname) ; del outdict['keylist']
        method='{ftype}_{method}'.format(**outdict)
        outdict=mmlparam.parspar('pysim.files','mklog_'+method,inpar=outdict)
        return outdict
    # Write
    elif rwid=='W':
        method='{ftype}_{method}'.format(**indict)
        outdict=mmlparam.parspar('pysim.files','mklog_'+method,inpar=indict)
        outdict['keylist']=['timestamp','makecomp']+mmlparam.listpar('pysim.files','mklog_'+method)['keylist']
        outdict.update(timestamp=strftime("%Y%b%d-%H:%M:%S",gmtime()),
                       makecomp =mmlinfo.hostname2compid())
        mmlio.rwdict('W',fname,outdict,overwrite=overwrite)
        return
    # Error
    else: raise Exception('Invalid rwid: {}'.format(rwid))
    return

def rw_units(rwid,fname,indict=None,overwrite=False):
    """
    Reads/writes unit files
    """
    from pysim.nbody import units
    # Pars input
    rwid=mmlpars.mml_pars(rwid,list=['R','W'])
    fname=mmlpars.mml_pars(fname,type=str)
    # Read
    if   rwid=='R':
        udict=mmlio.rwdict('R',fname)
        outdict={'M':units.Unit(str(udict['UnitMass_in_g'           ])+' g'       ),
                 'L':units.Unit(str(udict['UnitLength_in_cm'        ])+' cm'      ),
                 'V':units.Unit(str(udict['UnitVelocity_in_cm_per_s'])+' cm s**-1')}
        outdict['T']=outdict['L']/outdict['V']
        return outdict
    # Write
    elif rwid=='W':
        if 'T' in indict and indict['T']!=indict['L']/indict['V']:
            mmlio.verbose('T included ({}), but it does not match L/V ({})'.format(indict['T'],indict['L']/indict['V']))
        udict={'UnitMass_in_g'           :indict['M'].ratio('g'       ),
               'UnitLength_in_cm'        :indict['L'].ratio('cm'      ),
               'UnitVelocity_in_cm_per_s':indict['V'].ratio('cm s**-1')}
        mmlio.rwdict('W',fname,udict,overwrite=overwrite)
        return
    # Error
    else: raise Exception('Invalid rwid: {}'.format(rwid))
    return

####################################################################################################################################
####################################################################################################################################
# LIST OF FUNCTIONS
__all__=['listpar','get_fdict','get_flist','subrun','endrun','cleanrun',
         'mkfiles','mvfiles','rmfiles','cleantree','rw_makelog','rw_units']
