#!/usr/bin/env python
####################################################################################################################################
#
# MEAGAN LANG'S SIMFILE METHODS
#
####################################################################################################################################
import sys,os,shutil,glob,copy,pprint,scipy,math,itertools,re,importlib,fnmatch,inspect
import numpy as np
LIST_METHODS=['getfdict','getflist','mvrun','mvtree','ziprun','unziprun','procnum','cleantree']
from mmlutils import *
from .. import simlist


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
        outdict=mmlparam.parspar('pysim.files.simfile','mklog_'+method,inpar=outdict)
        return outdict
    # Write
    elif rwid=='W': 
        method='{ftype}_{method}'.format(**indict)
        outdict=mmlparam.parspar('pysim.files.simfile','mklog_'+method,inpar=indict)
        outdict['keylist']=['timestamp','makecomp']+mmlparam.listpar('pysim.files.simfile','mklog_'+method)
        outdict.update(timestamp=strftime("%Y%b%d-%H:%M:%S",gmtime()),
                       makecomp =mmlinfo.hostname2compid())
        mmlio.rwdict('W',fname,outdict,overwrite=overwrite)
        return
    # Error
    else: raise Exception('Invalid rwid: {}'.format(rwid))
