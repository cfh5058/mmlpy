#!/usr/bin/env python
####################################################################################################################################
#
# MEAGAN LANG'S SIMFILE METHODS
#
####################################################################################################################################
import sys,os,shutil,glob,copy,pprint,scipy,math,itertools,re,importlib,fnmatch
import numpy as np
LIST_METHODS=['getfdict','getflist','mvrun','mvtree','ziprun','unziprun','procnum','cleantree']
from mmlutils import *
import main as mmlyt
import simlist

DIR_SNAPFORM=os.path.join(os.path.dirname(mmlyt.__file__),'formats')
def ftype2module(ftype): return importlib.import_module('mml'+ftype)
def ftype2format(ftype): return os.path.join(DIR_SNAPFORM,ftype+'.pm')
def ftype2snapobj(ftype): return getattr(ftype2module(ftype),'NbodySnap_'+ftype)

LIST_SNAPFORM=[os.path.basename(ftype).split('.pm')[0] for ftype in sorted(glob.glob(ftype2format('*')))]

LIST_FILETYPES=LIST_SNAPFORM

LIST_FILEMETHODS=copy.deepcopy(simlist.LIST_METHTYP)
for imtyp in simlist.LIST_METHTYP_NOF: LIST_FILEMETHODS.remove(imtyp)

def main(): return mmlnbody.walk(mtype='file')


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
#             'timestamp':{'def':'','form':str },
#             'makecomp' :{'def':'','form':str },
             'baserun'   :{'def':'','form':str },
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
    else: raise Exception('Invalid method: {}'.format(method))
    # Return output
    return out

####################################################################################################################################
# METHOD TO RETURN SIMULATION FILES
def files(ftypes='all',retall=False,**input_kw):
    """
    Generates a dictionary of simulation file names
    """
    # Pars input
    fdict=pars_fileinput(**input_kw)
    fdict['profpar']=os.path.join('rundir',fdict['pfix_shrt']+'profpar')
    # Pars file types
    if isinstance(ftypes,str): ftypes=[ftypes]
    ftypes=mmlpars.mml_pars(ftypes,type=list,default=copy.deepcopy(LIST_FILETYPES))
    if 'all' in ftypes: ftypes=copy.deepcopy(LIST_FILETYPES)
    # Create directories
    ftab={} ; finfo={}
    for iprog in ftypes:
        ifpm=ftype2format(iprog)
        if not os.path.isfile(ifpm): continue
        # Initialize dictionaries
        ifdict={'dir':os.path.join(fdict['rundir'],iprog+fdict['subtag']),
                'pfix':fdict['pfix'],'pfix_shrt':fdict['pfix_shrt']}
        ifinfo={}
        # Read in table
        iftab=mmlio.rwtable('R',ifpm) ; iftab['path']=len(iftab['key'])*['None']
        tabkeys=iftab['keylist'] ; tabkeys.remove('key')
        for iclass in set(iftab['class']): ifdict[iclass]={}
        # Add directories
        for i in range(len(iftab['key'])):
            # Populate info
            ifinfo[iftab['key'][i]]={ikey:iftab[ikey][i] for ikey in tabkeys}
            # Parent directory
            if iftab['dir'][i]=='main': ipardir=ifdict['dir']
            else                      : ipardir=iftab['path'][iftab['key'].index(iftab['dir'][i])]
            # Add path to table
            if   iftab['flag'][i]=='dir'    : iftab['path'][i]=os.path.join(ipardir,iftab['name'][i])
            elif iftab['flag'][i]=='selfdir': iftab['path'][i]=os.path.join(ipardir,iftab['name'][i],fdict['pfix']+iftab['name'][i])
            elif iftab['flag'][i]=='nopfix' : iftab['path'][i]=os.path.join(ipardir,iftab['name'][i])
            elif iftab['flag'][i]=='shpfix' : iftab['path'][i]=os.path.join(ipardir,fdict['pfix_shrt']+iftab['name'][i])
            else                            : iftab['path'][i]=os.path.join(ipardir,fdict['pfix']+iftab['name'][i])
            # Populate dictionary
            if   iftab['flag'][i]=='dir': 
                if iftab['class'][i]==iftab['key'][i]: ifdict[iftab['class'][i]]['dir']=iftab['path'][i]
                else                                 : ifdict[iftab['class'][i]][iftab['key'][i]+'dir']=iftab['path'][i]
            else                        : ifdict[iftab['class'][i]][iftab['key'][i]]=iftab['path'][i]
        # Add ith dictionary
        fdict[iprog]=ifdict ; finfo[iprog]=ifinfo ; ftab[iprog]=iftab
    if retall: return fdict,finfo,ftab
    else     : return fdict

####################################################################################################################################
# METHOD TO RETURN LIST OF RELEVANT FILES
def get_filelist(ftype,groups=None,**file_kw):
    """
    Returns a list of relavent files
    """
    # Pars input
    ftype=mmlpars.mml_pars(ftype,list=LIST_FILETYPES)
    # Get files
    file_kw['retall']=True
    fdict,finfo,ftab=files(ftypes=[ftype],**file_kw)
    fdict,finfo,ftab=fdict[ftype],finfo[ftype],ftab[ftype]
    # Pars
    if isinstance(groups,str): groups=[groups]
    groups=mmlpars.mml_pars(groups,default=['all'],type=list)
    if 'all' in groups: groups=set(ftab['group'])
    # Initialize file list
    tabform={ikey:[] for ikey in ftab['keylist']}
    filelist={igrp:[] for igrp in groups}
    filelist['keylist']=copy.deepcopy(groups)
    # Loop over files
    for i in range(len(ftab['key'])):
        if ftab['group'][i] not in groups: continue
        filelist[ftab['group'][i]].append({ikey:ftab[ikey][i] for ikey in ftab['keylist']})
    # Print info on missing groups
    for grp in groups:
        if len(filelist[grp])==0: mmlio.verbose('No files found in group {}'.format(grp))
    # Return output
    return filelist

####################################################################################################################################
####################################################################################################################################
# FILE UTILITIES
def search_file(topdir,extlist):
    """
    Searchs directory tree for a file
    """
    if os.path.isfile(topdir): topdir=os.path.dirname(os.path.dirname(topdir))
    fname=False
    for root,dirs,files in os.walk(topdir):
        for ext in extlist:
            mkcand=fnmatch.filter(files,'*'+ext+'*')
            if len(mkcand)>0: fname=os.path.join(root,mkcand[0])
            if isinstance(fname,str): break
        if isinstance(fname,str): break
    return fname

####################################################################################################################################
####################################################################################################################################
# METHODS FOR READING/WRITING FILES

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
        outdict=mmlparam.parspar('mmlyt.simfile','mklog_'+method,inpar=outdict)
        return outdict
    # Write
    elif rwid=='W': 
        method='{ftype}_{method}'.format(**indict)
        outdict=mmlparam.parspar('mmlyt.simfile','mklog_'+method,inpar=indict)
        outdict['keylist']=['timestamp','makecomp']+mmlparam.listpar('mmlyt.simfile','mklog_'+method)
        outdict.update(timestamp=strftime("%Y%b%d-%H:%M:%S",gmtime()),
                       makecomp =mmlinfo.hostname2compid())
        mmlio.rwdict('W',fname,outdict,overwrite=overwrite)
        return
    # Error
    else: raise Exception('Invalid rwid: {}'.format(rwid))

####################################################################################################################################
####################################################################################################################################
# METHODS FOR MOVING/RENAMING RUN FILES/DIRECTORIES

####################################################################################################################################
# METHOD FOR MOVING/RENAMING A RUN
def mvrun(isimstr,**options):
    """
    Move/rename files assosiated with a run
    """
    # Get sim objects
    print 'Getting info on new file names...'
    fsimstr=mmlyt.asksimboj(infodict=isimstr.infodict,runtyp=isimstr['runtyp'])
    fname=mmlparam.par2fname('mmlyt.simlist','mmlsim',inpar=dict(fsimstr))
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
        for ftype in LIST_FILETYPES:
            mvfiles(isimstr,fsimstr,ftype,groups='all',**options)
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
        mmlparam.savepar('mmlyt.simlist','galsim',dict(simstr),overwrite=True)
    return

####################################################################################################################################
# METHOD FOR MOVING/RENAMING RUN FILES BY TYPE
def mvfiles(simstr1,simstr2,ftypes='all',groups='all',overwrite=False,rmsrc=None,
            validext_end={},validext_list={},validext_fchk={},
            memcomp1=None,memcomp2=None,**extra_kw):
    """
    Method to move/rename run files
    """
    # Pars input
    ftype=mmlpars.mml_pars(ftype,list=copy.deepcopy(LIST_FILETYPES))
    memcomp1=mmlpars.mml_pars(memcomp1,default=simstr1['memcomp'],type=str)
    memcomp2=mmlpars.mml_pars(memcomp2,default=simstr2['memcomp'],type=str)
    if not isinstance(rmsrc,bool): rmsrc=mmlio.yorn('[mvfiles] Remove source files?')
    # Get source/destination hosts
    srchost = mmlinfo.computers('CFHMAC')['host'] if memcomp1=='cfhmac' else ''
    dsthost = mmlinfo.computers('CFHMAC')['host'] if memcomp2=='cfhmac' else ''
    # Create file keys
    fkeys_src=dict(memcomp=memcomp1,checkaccess=True)
    fkeys_dst=dict(memcomp=memcomp2,checkaccess=True)
    # Get file dictionary
    outsrc=simstr1.mkfiledict(ftypes=[ftype],retall=True,**fkeys_src)
    outdst=simstr2.mkfiledict(ftypes=[ftype],retall=True,**fkeys_dst)
    fdict_src,finfo_src,ftab_src=(iout[ftype] for iout in outsrc)
    fdict_dst,finfo_dst,ftab_dst=(iout[ftype] for iout in outdst)
    flist_src=simstr1.get_filelist(ftype,groups=groups,**fkeys_src)
    flist_dst=simstr2.get_filelist(ftype,groups=groups,**fkeys_dst)
    # Loop over particle type
    for grp,infolist in flist_src.iteritems():
        if not mmlio.yorn('Move {} {} file(s)?'.format(ftype,grp)): continue
        # Check for large files
        longgrp=('large' in [info['size'] for info in infolist])
        # Get overwrite options
        if longgrp: owflag=mmlio.yorn('Overwrite existing {} {} file(s)?'.format(ftype,grp))
        else      : owflag=overwrite
        # Loop over files
        for srcinfo,dstinfo in zip(flist_src[grp],flist_dst[grp]):
            srclist=[] ; dstlist=[]
            temp=srcinfo['path']
            if srcinfo['size']=='large': temp+='*'
            if len(glob.glob(temp))==0 and memcomp1!='cfhmac': continue
            # Groups with large lists of files
            if longgrp: 
                if memcomp1=='cfhmac': raise Exception('Not currently supported.')
                if memcomp2=='cfhmac': raise Exception('Not recommended.')
                # Initialize restrictions on extensions
                validext=dict(end=False,list=False,fchk=False)
                if srcinfo['key'] in validext_end : validext['end' ]=validext_end[srcinfo['key']]
                if srcinfo['key'] in validext_list: validext['list']=validext_list[srcinfo['key']]
                if srcinfo['key'] in validext_fchk: validext['fchk']=finfo_src[validext_fchk[srcinfo['key']]]['path']
                # Generate lists of source/destination files and extensions
                extlist=[]
                for src in sorted(glob.glob(temp)):
                    # Get the extensions for this ith source
                    ext=re.findall(srcinfo['path'].replace('*','(.*)'),src)[0]
                    # Pars restrictions on extensions
                    if validext['end']  and not ext.endswith(validext['end']): continue
                    if validext['list'] and ext not in validext['list']: continue
                    if validext['fchk'] and not owflag and os.path.isfile(validext['fchk']+ext): continue
                    # Append lists
                    extlist.append(ext)
                    srclist.append(src)
                    dstlist.append(dstinfo['path']+ext)
                # Check for existing output for SCF
                # validext_fchk['snapin']='snapout'
                # Only move GADGET snapshots that end in snapend or in the list
                # validext_end['snapshot']='5'
                # validext_list['snapshot']='_000'
            # Groups with short lists of files
            else:
                srclist.append(srcinfo['path'])
                dstlist.append(dstinfo['path'])
            # Move files
            print '({},{},{})'.format(ftype,grp,srcinfo['key'])
            for src,dst in zip(srclist,dstlist):
                if src==dst: continue
                if not memcomp2=='cfmac': mmlfiles.mkdirs(os.path.dirname(dst))
                print '{} ---> {}'.format(src,dst)
            mmlfiles.cpfiles(srclist,dstlist,srchost=srchost,dsthost=dsthost,
                             overwrite=owflag,move=rmsrc)
    return

####################################################################################################################################
# METHOD TO MOVE FILES TO/FROM ACCRE
def mvfiles_accre(simstr,method,fmeth):
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
    file_kw['memcomp1']=memcomp1
    file_kw['memcomp2']=memcomp2
    # Call function for generic file move
    mvfiles(simstr1,simstr2,fmeth,**file_kw)
    # Return control
    return

####################################################################################################################################
# METHOD TO MOVE FILES TO/FROM CFHMAC
def mvfiles_cfhmac(simstr,method,fmeth,**file_kw):
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
    file_kw['memcomp1']=memcomp1
    file_kw['memcomp2']=memcomp2
    # Call function for generic file move
    mvfiles(simstr1,simstr2,fmeth,**file_kw)
    # Return control
    return

####################################################################################################################################
####################################################################################################################################
# METHODS FOR CLEANING RUN FILES/DIRECTORIES

####################################################################################################################################
# METHOD TO DELETE SIMULATION FILES
def rmfiles(ftype,groups=[],**input_kw):
    """
    Deletes files
    """
    # Pars input
    ftype=mmlpars.mml_pars(ftype,list=LIST_FILETYPES)
    file_kw=pars_fileinput(**input_kw)
    memcomp=file_kw['memcomp']
    rundir=file_kw['rundir']
    runtag=file_kw['runtag']
    subtag=file_kw['subtag']
    # Double check location of removed files
    if memcomp not in ['accre','cfhmac']:
        if not mmlio.yorn('Are you sure you want to remove {} files from {}? ({})'.format(fmeth,memcomp,runtag)): return
    # Get list of files
    flist=get_filelist(ftype,groups=groups,**file_kw)
    # Loop over file groups
    rmlist=[]
    for grp,infolist in flist.iteritems():
        rmgrp=None
        for info in infolist:
            temp=info['path']
            if info['size']=='large': temp+='*'
            if len(glob.glob(temp))==0: continue
            if rmgrp==None: rmgrp=mmlio.yorn('Remove {} {} {} file(s)? ({})'.format(memcomp,ftype,grp,runtag))
            if rmgrp: rmlist.append(temp)
    # Loop over files removing them
    if not mmlio.yorn('Remove?'): return
    for itemp in rmlist:
        for ifile in glob.iglob(itemp):
            os.remove(ifile)
    # Return control
    return

####################################################################################################################################
# METHOD TO REMOVE EMPTY TREE DIRECTORIES
def cleantree(ftype,**input_kw):
    """
    Method to remove empty tree directories
    """
    # Get list of files
    flist=get_filelist(ftype,**input_kw)
    # Loop over groups
    for grp,infolist in flist.iteritems():
        # Loop over files
        for info in infolist:
            idir=os.path.dirname(info['path'])
            # Remove dir if empty
            if os.path.isdir(idir) and len(glob.glob(os.path.join(idir,'*')))==0:
                print '{} directory empty. Removing: {}'.format(grp,idir)
                os.removedirs(idir)
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
                   inclmktag =None,nprocmk =None,
                   inclevtag =None,nprocev =None,
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
                'inclmktag'  :False,'nprocmk' :0,
                'inclevtag'  :False,'nprocev' :0,
                'inclscftag' :False,'nprocscf':0}
    memcomp=mmlpars.mml_pars(memcomp,default=simstr['memcomp'],type=str)
    runtag=mmlpars.mml_pars(runtag,default=simstr['runtag'],type=str)
    runtyp=mmlpars.mml_pars(runtyp,default=simstr['runtyp'],type=str)
    subtag=mmlpars.mml_pars(subtag,default=simstr['subtag'],type=str)
    inclmktag =mmlpars.mml_pars(inclmktag ,default=simstr['inclmktag' ],type=bool)
    inclevtag =mmlpars.mml_pars(inclevtag ,default=simstr['inclevtag' ],type=bool)
    inclscftag=mmlpars.mml_pars(inclscftag,default=simstr['inclscftag'],type=bool)
    nprocmk =mmlpars.mml_pars(nprocmk ,default=simstr['nprocmk' ],type=int,min=0)
    nprocev =mmlpars.mml_pars(nprocev ,default=simstr['nprocev' ],type=int,min=0)
    nprocscf=mmlpars.mml_pars(nprocscf,default=simstr['nprocscf'],type=int,min=0)
    # Get directories
    topdirDEF=get_topdir(memcomp=memcomp,checkaccess=checkaccess)
    topdir=mmlpars.mml_pars(topdir,default=topdirDEF,type=str)
    simdirDEF=get_simdir(runtyp=runtyp,topdir=topdir)
    simdir=mmlpars.mml_pars(simdir,default=simdirDEF,type=str)
    rundirDEF=get_rundir(runtag=runtag,simdir=simdir)
    rundir=mmlpars.mml_pars(rundir,default=rundirDEF,type=str)
    # Get added tag
    subtagproc=get_subtag(subtag=subtag,
                          inclmktag =inclmktag ,nprocmk =nprocmk ,
                          inclevtag =inclevtag ,nprocev =nprocev ,
                          inclscftag=inclscftag,nprocscf=nprocscf)
    # Get prefixes
    pfix,pfix_shrt=get_prefix(runtag,subtagproc)
    # Get other file stuff
    mkinfo=os.path.join(rundir,runtag+'.mkinfo')
    archiv=os.path.join(simdir,runtag+'.tar.gz')
    rapsht=mmlparam.par2fname('mmlyt.simlist','mmlsim',tagstr=runtag)
    # Create and return dictionary
    input_kw=dict(rundir=rundir,runtag=runtag,subtag=subtagproc,memcomp=memcomp,
                  topdir=topdir,simdir=simdir,
                  inclmktag =inclmktag ,nprocmk =nprocmk ,
                  inclevtag =inclevtag ,nprocev =nprocev ,
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
    # Create subrunid including # of processors and subtag
    subtag=""
    if simstr['inclmktag' ]: subtag+="_MK%(nprocmk)d"   % simstr
    if simstr['inclevtag' ]: subtag+="_EV%(nprocev)d"   % simstr
    if simstr['inclscftag']: subtag+="_SCF%(nprocscf)d" % simstr
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
    tagkw=mmlparam.parspar('mmlyt.simlist',simtype,inpar=tagkw,askuser=True)
    # Proceed based on simtype
    runtag=mmlparam.par2tag('mmlyt.simlist',simtype,tagkw)
    # Return runtag
    return runtag
