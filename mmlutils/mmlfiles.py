####################################################################################################################################
#
# MEAGAN LANG'S FILE METHODS
#
####################################################################################################################################
import mmlpars,mmlclass,mmlio,mmlinfo
import os,shutil,subprocess,types,copy,itertools,fnmatch,stat,glob


def pathtype_local(path):
    """
    Returns the type of local paths
    """
    if   os.path.isdir(path) : pathtype='dir'
    elif os.path.isfile(path): pathtype='file'
    else                     : pathtype='none'
    return pathtype
def pathtype_remote(path,sftpclient):
    """
    Returns the type of remote paths
    """
    try:
        pathstat=sftpclient.stat(path)
        if stat.S_ISDIR(pathstat.st_mode): pathtype='dir'
        else                             : pathtype='file'
    except IOError: pathtype='none'
    return pathtype

def chkhost(host):
    if not host: return host
    else:
        host=mmlpars.mml_pars(host,type=str)
        memlist=mmlinfo.computers()['memlist']
        if host in memlist: return None
        else              : return host

####################################################################################################################################
def ls(temp,host=None):
    """
    Returns a list of files matching the template
    """
    host=chkhost(host)
    if   isinstance(temp,str ): temp=[temp]
    elif isinstance(temp,list): pass
    else: raise Exception('template ({}) is incorrect format ({})'.format(temp,type(temp)))
    if host:
        try   : flist=subprocess.check_output(['ssh',host,'ls']+temp).split('\n')
        except: flist=[]
    else:
        flist=[]
        for itemp in temp: flist+=glob.glob(itemp)
    return sorted(flist)

def rm(temps,host=None):
    """
    Removes files
    """
    host=chkhost(host)
    if isinstance(temps,str): temps=[temps]
    for t in temps:
        for f in glob.glob(t):
            # flist=ls(ftemp,host=host)
            if host: subprocess.call(['ssh',host,'rm',f])
            else   : os.remove(f)

####################################################################################################################################
# METHOD TO SEARCH A DIRECTORY FOR A FILE
def search(topdir,extlist,exact=False,nlvlup=0):
    """
    Searchs directory tree for a file
    """
    fname=False
    if not topdir: return fname
    if os.path.isfile(topdir): topdir=os.path.dirname(topdir)
    for ilvl in range(nlvlup): topdir=os.path.dirname(topdir)
    for root,dirs,files in os.walk(topdir):
        for ext in extlist:
            if exact: mkcand=fnmatch.filter(files,ext)
            else    : mkcand=fnmatch.filter(files,'*'+ext+'*')
            if len(mkcand)>0: fname=os.path.join(root,mkcand[0])
            if isinstance(fname,str): break
        if isinstance(fname,str): break
    return fname

####################################################################################################################################
# METHOD TO PREPARE FOR FILE WRITE
def prep_write(fname,overwrite=None,askuser=None):
    if fname is None: return False
    fname=mmlpars.mml_pars(fname,type=str)
    askuser=mmlpars.mml_pars(askuser,default=False,type=bool)
    if fname in ['stdout','stderr']: return True
    if not isinstance(fname,str): return False
    if len(fname)==0: return False
    if os.path.isfile(fname):
        if not isinstance(overwrite,bool):
            if askuser:
                print fname
                overwrite=mmlio.yorn('File already exists. Overwrite?')
            else      : overwrite=mmlpars.mml_pars(overwrite,default=False,type=bool)
        if overwrite:
            os.remove(fname)
        else:
            print 'File already exists: {}'.format(fname)
            return False
    if len(os.path.dirname(fname)) > 0: mkdirs(os.path.dirname(fname))
    return True

####################################################################################################################################
# METHOD TO MAKE DIRECTORIES
####################################################################################################################################
def mkdir(*args,**kwargs): return mkdirs(*args,**kwargs)
def mkdirs(dirlist0,host=None):
    '''
    NAME:
        mmlfiles.mkdirs
    PURPOSE:
        To create directories if they do not exist.
    CALLING:
        mmlfiles.mkdirs(dirlist)
    ARGUMENTS:
        dirlist: List of directories to create.
    '''
    host=chkhost(host)
    # Ensure that even singular directory is iterable
    if not isinstance(dirlist0,list):
        dirlist=[dirlist0]
    else:
        dirlist=dirlist0
    # Iterate over directories, creating them if they don't exist
    for idir in dirlist:
        if host: subprocess.call(['ssh',host,'mkdir','-p',idir])
        else:
            if not os.path.exists(idir):
                os.makedirs(idir)


####################################################################################################################################
# METHOD TO COPY FILES/DIRECTORIES
def cp(*args,**kwargs): return cpfiles(*args,**kwargs)
def mv(*args,**kwargs):
    kwargs['move']=True
    return cpfiles(*args,**kwargs)
def cpfiles(srclist0,dstlist0,overwrite=False,move=False,cpmeta=True,verbose=False,
            srchost=None,dsthost=None):
    '''
    NAME:
        mmlfiles.cpfiles
    PURPOSE:
        To copy a list of files.
    CALLING:
        cpfiles(srclist,dstlist[,overwrite=,move=,cpmeta=,srchost=,dsthost=])
    ARGUMENTS:
        srclist:   List of source file paths.
        dstlist:   List of destination file paths.
    KEYWORDS:
        overwrite: Boolean determining if existing destination files are overwritten.
        move:      Boolean determining if existing source files are deleted.
        cpmeta:    Boolean determining if meta data should/shouldn"t be copied.
        srchost:   Str address of computer hosting source files
        dsthost:   Str address of computer hosting destination files
    '''
    # Pars input
    srchost=chkhost(srchost)
    dsthost=chkhost(dsthost)
    # Ensure input is in list
    if not isinstance(srclist0,list):
        srclist=[srclist0]
    else:
        srclist=srclist0
    if not isinstance(dstlist0,list):
        dstlist=[dstlist0]
    else:
        dstlist=dstlist0
    # Make sure lists are same lengths
    nsrc=len(srclist)
    ndst=len(dstlist)
    if nsrc != ndst:
        raise Exception('Number of destination paths ({}) must match number of source paths ({}).'.format(ndst,nsrc))
    # Open sftp if source/destination remote
    if   srchost and dsthost:
        raise Exception('Remote source to remote destination not supported.')
    elif srchost or  dsthost:
        if move  : mmlio.verbose('Source files will not be removed during remote transfer.')
        if cpmeta: mmlio.verbose('Meta data will not be copied during remote transfers')
        import paramiko
        # Load info from SSH Config file
        sshconf=paramiko.SSHConfig()
        sshconf.parse(open(os.path.expanduser('~/.ssh/config')))
        # Get info on the remote host
        if   srchost: sshinfo=sshconf.lookup(srchost)
        elif dsthost: sshinfo=sshconf.lookup(dsthost)
        # Get SSH Client
        sshclnt=paramiko.SSHClient()
        sshclnt.load_system_host_keys()
        sshclnt.connect(sshinfo['hostname'],username=sshinfo['user'],
                        key_filename=sshinfo['identityfile'])
        # Get SFTP Client
        sftpclnt=sshclnt.open_sftp()
    # Define functions
    if srchost: srctype=lambda path: pathtype_remote(path,sftpclnt)
    else      : srctype=lambda path: pathtype_local(path)
    if dsthost: dsttype=lambda path: pathtype_remote(path,sftpclnt)
    else      : dsttype=lambda path: pathtype_local(path)
    # Loop over files
    for isrc,idst in zip(srclist,dstlist):
        if isrc == idst: continue
        # Ensure files do not get overwritten unless desired
        ifilecopy=False ; idircopy=False
        # Identify file types
        isrctype=srctype(isrc)
        idsttype=dsttype(idst)
        # Handle errors
        if isrctype=='none': raise Exception('Source is not a valid path: {}'.format(isrc))
        elif isrctype=='file' and idsttype=='dir' : 
            idst=os.path.join(idst,os.path.basename(isrc))
            idsttype=dsttype(idst)
        elif isrctype=='dir' and idsttype=='file': raise Exception('Cannot copy a directory ({}) to a file ({}).'.format(isrc,idst))
        # Prevent overwrite
        if idsttype!='none' and not overwrite: continue
        # Create directory
        idstdir=os.path.dirname(idst)
        if dsttype(idstdir)!='dir':
            # if dsthost: dirout=subprocess.call(['ssh',dsthost,'mkdir','-p',idstdir])
            if dsthost: dirout=sshclnt.exec_command('mkdir -p {}'.format(idstdir))
            else      : mkdirs(idstdir)
        # Move files
        if verbose: print '{} ---> {}'.format(isrc,idst)
        if srchost or dsthost:
            spwnarg=['rsync','-a']
            if srchost: spwnarg.append(srchost+':'+isrc)
            else      : spwnarg.append(isrc)
            if dsthost: spwnarg.append(dsthost+':'+idst)
            else      : spwnarg.append(idst)
            subprocess.call(spwnarg)
        # if   srchost: sftpclnt.get(isrc,idst)
        # elif dsthost: sftpclnt.put(isrc,idst)
        else:
            if move: shutil.move(isrc,idst)
            else:
                if cpmeta: shutil.copy2(isrc,idst)
                else     : shutil.copy(isrc,idst)

        # # Handle both files being on remote host
        # if len(srchost)>0 and len(dsthost)>0:
        #     # Initialize command
        #     spwnarg=['scp','-r']
        #     # Handle source
        #     if len(srchost)==0: spwnarg.append(isrc)
        #     else              : spwnarg.append(srchost+':'+os.path.basename(isrc))
        #     # Handle destination
        #     if len(dsthost)==0: spwnarg.append(idst)
        #     else              : spwnarg.append(dsthost+':'+os.path.basename(idst))
        #     # Call copy command
        #     subprocess.call(spwnarg)
        # # Handle source file
        # elif os.path.isfile(isrc):
        #     if os.path.isdir(idst):
        #         idst=os.path.join(idst,os.path.basename(isrc))
        #     elif os.path.isfile(idst):
        #         if overwrite: os.remove(idst)
        #         else: continue
        #     if move: shutil.move(isrc,idst)
        #     else:
        #         if cpmeta: shutil.copy2(isrc,idst)
        #         else: shutil.copy(isrc,idst)
        # # Handle source directory
        # elif os.path.isdir(isrc):
        #     if os.path.isdir(idst):
        #         if overwrite: shutil.rmtree(idst)
        #         else: continue
        #     elif os.path.isfile(idst):
        #         raise Exception('Source is a directory ({}) and destination is a file({}).'.format(isrc,idst))
        #     shutil.copytree(isrc,idst)
        #     if move: shutil.rmtree(isrc)
        # # Handle invalide source
        # else:
        #     if not os.path.exists(isrc):
        #        raise Exception('Source is not a valid path: {}'.format(isrc))
        #     else:
        #         raise Exception('Source is a path, but not file or directory. What?: {}'.format(isrc))

    # Close sftp if source/destination remote
    if   srchost and dsthost: pass
    elif srchost or  dsthost:
        sftpclnt.close()
        sshclnt.close()
    # Return control
    return

####################################################################################################################################
# METHOD TO MAKE A DICTIONARY FULL OF FILE NAMES
####################################################################################################################################
def filedict(topdir,filekeys=[],filelist=[],dirkeys=[],dirlist=[],pfix='',fdict={}):
    '''
    NAME:
        mmlfiles.filedict
    PURPOSE:
        To create a dictionary of full file paths.
    CALLING:
        fdict=filedict(topdir,filekeys=[],filelist=[],dirkeys=[],dirlist=[],pfix="",fdict={})
    ARGUMENTS:
        topdir:   Full path directory in which files/directories should be placed
    KEYWORDS:
        filekeys: List of keywords used for generated file paths in created dict
        filelist: List of file base names to be used in created dict
        dirkeys:  List of keywords used for generated dir paths in created dict
        dirlist:  List of directory base names to be used in created dict
        pfix:     String prefix to append to beginning of file names.
        fdict:    Dictionary to added file/dir key/value pairs to.
    OUTPUT:
        fdict:    Dictionary with file/dir key/vale pairs in it.
    '''
    nfile=len(filekeys)
    ndir=len(dirkeys)
    if len(filelist) != nfile:
        raise Exception('Length of filekeys list must be the same as filelist.')
    if len(dirlist) != ndir:
        raise Exception('Length of dirkeys list must be the same as dirlist.')

    for ifile in range(nfile):
        fdict[filekeys[ifile]]=os.path.join(topdir,pfix+filelist[ifile])
    for idir in range(ndir):
        fdict[dirkeys[idir]]=os.path.join(topdir,dirlist[idir])

    return fdict
