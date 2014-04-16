#!/usr/bin/python
import mmlclass,mmlpars,mmlio
import os

####################################################################################################################################
#
# MMLINFO METHODS
#
####################################################################################################################################
COMPDICT_PPN={'accre':2,'stampede':1}
COMPDICT_MPI={'accre'   :'/usr/local/mpiexec/latest/x86_64/gcc46/nonet/bin/mpiexec',
              'stampede':'ibrun'}
LIST_COMPUTERS=['bender','flexo','vpac36','vpac38','accre','stampede']
LIST_COSMOSIM=['SINHA1024','SINHA512']

####################################################################################################################################
def dir2compid(dir):
    """
    Finds a computers ID string given a directory
    """
#    raise Exception('Do not use this anymore!!')
    if dir==None: compid=hostname2compid()
    else:
        # Pars input
        dir=mmlpars.mml_pars(dir,type=str)
        # Proceed based on substring
        if   'vpac' in dir:
            idxvpac=dir.find('vpac')
            compid=dir[idxvpac:idxvpac+6]
        elif 'scratch' in dir:
            compid='accre'
        elif 'fs0' in dir:
            compid='bender'
        elif 'data1' in dir:
            compid='flexo'
        else:
            raise Exception('Cannot pars directory into computer ID: {}'.format(dir))
    # Return computer ID
    return compid

####################################################################################################################################
def hostname2compid(hostid=None):
    """
    Converts the hostname returned by os.environ['HOSTNAME'] to a short computer ID
    """
    # Pars input
    hostid=mmlpars.mml_pars(hostid,type=str,default=os.environ['HOSTNAME'])
    # Check for ACCRE hosted machines
    if '.vampire' in hostid:
        if   'bender' in hostid: compid='bender'
        elif 'flexo'  in hostid: compid='bender'
        else: compid='accre'
    # Otherwise take it as is
    else:
        compid=hostid
    # Return computer ID
    return compid

####################################################################################################################################
def complist(update=False):
    """
    Returns dictionary of computer properties
    """
    # Set constants
    compfile='/home/langmm/code/python/mine/files/complist.txt'
    proplist=['compid','owner','userid','host','simdir','memlist']
    # Read complist from file
    complist=mmlio.rwtable('R',compfile)
    vpacLIST_astro=[ivpac for ivpac in complist['compid'] if 'vpac' in ivpac]
    # Add missing info
    for iprop in proplist:
        if iprop not in complist: complist[iprop]=['?' for icomp in complist['compid']]
    for idx in range(len(complist['compid'])):
        # User ID
        if '?' == complist['userid'][idx]: complist['userid'][idx]=os.getlogin()
        # Host address
        if '?' == complist['host'  ][idx]:
            if 'vpac' in complist['compid'][idx]:
                host='{}.phy.vanderbilt.edu'.format(complist['compid'][idx])
            elif complist['compid'][idx] in ['bender','flexo']:
                host='{}.phy.vanderbilt.edu'.format(complist['compid'][idx])
            elif complist['compid'][idx] == 'accre':
                host='vmplogin.accre.vanderbilt.edu'
            elif complist['compid'][idx] == 'cfhmac':
                host='68.53.59.122'
            elif complist['compid'][idx] == 'local':
                host=os.environ['SSH_CLIENT'].split()[0]
            complist['host'][idx]=host
        # Memory list
        if '?' == complist['memlist'][idx]:
            if 'vpac' in complist['compid'][idx]: 
                memlist=vpacLIST_astro
            elif complist['compid'][idx] in ['bender','flexo']:
                memlist=['vpac36','vpac38','accre']
            elif complist['compid'][idx] == 'accre':
                memlist=['accre']
            elif complist['compid'][idx] == 'cfhmac':
                memlist=['cfhmac']
            elif complist['compid'][idx] == 'local':
                memlist=['local']
            complist['memlist'][idx]=memlist
        # Simulation directory
        if '?' == complist['simdir'][idx]:
            if 'vpac' in complist['compid'][idx]:
                simdir='/net/{}/astro1/{}'.format(complist['compid'][idx],complist['userid'][idx])
            elif complist['compid'][idx] == 'bender':
                simdir='/fs0/{}'.format(complist['userid'][idx])
            elif complist['compid'][idx] == 'flexo':
                simdir='/data1/{}'.format(complist['userid'][idx])
                print '[mmlinfo.complist] Check where to store files on FLEXO'
            elif complist['compid'][idx] == 'accre':
                simdir='/scratch/{}'.format(complist['userid'][idx])
            elif complist['compid'][idx] == 'cfhmac':
                simdir='/Users/{}'.format(complist['userid'][idx])
            elif complist['compid'][idx] == 'local':
                simdir='~/'
            complist['simdir'][idx]=simdir
    # Update file
    if update: mmlio.rwtable('W',compfile,tabdict=complist,overwrite=True,keylist=proplist)
    # Return list
    return complist

####################################################################################################################################
def computers(usrcomp=None,memcomp=None,compusr=None,compmem=None,checkaccess=None):
    '''
    NAME:
        mmlinfo.computers
    CALLING:
        compdict=computers(usrcomp,memcomp=None)
    ARGUMENTS:
        usrcomp:   String specifying what computer to return info on.
    KEYWORDS:
        memcomp:  String specifying what VPAC machine to use for storage on vpac network
    OUTPUT:
        compdict: Dictionary of information on a computer.
    '''
    # Pars input
    if compusr != None: usrcomp=compusr
    if compmem != None: memcomp=compmem
    usrcomp=mmlpars.mml_pars(usrcomp,default=hostname2compid(),type=str)
    memcomp=mmlpars.mml_pars(memcomp,default=usrcomp,type=str)
    checkaccess=mmlpars.mml_pars(checkaccess,default=False,type=bool)
    # Force usrcomp to uppercase to ensure matching
    usrcomp=usrcomp.lower()
    memcomp=memcomp.lower()
    # Read compinfo from list
    compinfo=complist()
    vpacLIST_astro=[ivpac for ivpac in compinfo['compid'] if 'vpac' in ivpac]
    if usrcomp not in compinfo['compid']: raise Exception('Unrecognized computer ID for user: {}'.format(usrcomp))
    if memcomp not in compinfo['compid']: raise Exception('Unrecognized computer ID for mem : {}'.format(usrcomp))
    idxusr=compinfo['compid'].index(usrcomp)
    idxmem=compinfo['compid'].index(memcomp)
    # Set user comp dependent stuff
    owner=compinfo['owner'][idxusr]
    userid_usr=compinfo['userid'][idxusr]
    host=compinfo['host'][idxusr]
    memlist=compinfo['memlist'][idxusr]
    # Set mem comp dependent stuff
    if checkaccess and memcomp not in memlist:
        raise Exception('Memory on {} is not accessible by {}.'.format(memcomp,usrcomp))
    userid_mem=compinfo['userid'][idxmem]
    if usrcomp in ['bender','flexo'] and memcomp in vpacLIST_astro:
        simdir='/net/{}.phy.vanderbilt.edu/astro1/{}'.format(memcomp,userid_mem)
    else:
        simdir=compinfo['simdir'][idxmem]
    # Create dictionary
    userhost=userid_usr+'@'+host
    compdict=mmlclass.mydict(compid=usrcomp,owner=owner,host=userhost,simdir=simdir,memlist=memlist)
    return compdict

####################################################################################################################################
# METHOD TO RETURN INFO ON MANUDEEPS COSMOLOGICAL SIMULATIONS
def cosmosim(simid=None,askuser=None):
    """
    Returns a dictionary with information on a cosmological simulation
    """
    # Pars input
    askuser=mmlpars.mml_pars(askuser,default=False)
    if simid not in LIST_COSMOSIM and askuser:
        simid=mmlio.askselect('Select a valid cosmological simulation:',LIST_COSMOSIM)
    simid=mmlpars.mml_pars(simid,list=LIST_COSMOSIM,default=LIST_COSMOSIM[0])
    # Add info based on id string
    if   simid=='SINHA1024': simstr={'simid'   :simid,
                                     'boxsize' :50.0 ,
                                     'npart3'  :1024 ,
                                     'zstart'  :249. ,
                                     'groupdir':'/data1/sinham/isolatedhalos/50mpc/1024/new/groups_links'}
    elif simid=='SINHA512' : simstr={'simid'   :simid,
                                     'boxsize' :50.0 ,
                                     'npart3'  :512  ,
                                     'zstart'  :249. ,
                                     'groupdir':'/data1/sinham/isolatedhalos/50mpc/512/z249/groups'}
    else: raise Exception('Invalid simulation ID: {}'.format(simid))
    # Return output
    return simstr



