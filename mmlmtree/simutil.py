#!/usr/bin/python    
###################################################################################################################################
#
# MEAGAN LANG'S SIMSTAT METHODS
#
####################################################################################################################################
import matplotlib
#matplotlib.use('Agg')
import matplotlib.pyplot as plt
import sys,os,shutil,glob,copy,pprint,scipy,math,cPickle
import numpy as np
from mmlutils import *
import main as mmlmtree
import simfile

LIST_METHODS=['haloid2groupnum','haloid2lastint']

def main():
    """
    Provides command line access to STAT methods
    """
    out = mmlnbody.walk(mtype='simutil')
    return out

####################################################################################################################################
####################################################################################################################################
# HIGH LEVEL METHODS REQUIRING MMLSIM

####################################################################################################################################
# METHOD FOR RUNNING DIFFERENT STAT METHODS
def run(simstr,method=None,**method_kw):
    """
    Provides interface for running different methods
    """
    if method not in LIST_METHODS: method=mmlio.askselect('Select a valid simutil method:',LIST_METHODS)
    # Proceed based on method
    if   method=='haloid2groupnum': haloid2groupnum(simstr,**method_kw)
    elif method=='haloid2lastint' : haloid2lastint(simstr,**method_kw)
    else: raise Exception('Invalid method: {}'.format(method))
    # Return output
    return out

###################################################################################################################################
# METHOD FOR CONVERTING HALO IDS TO INFO ON THE LAST INTERACTION
def haloid2lastint(simstr,snapnum=None,haloid=None,askuser=None,verbose=None,outflags=None,
                   intlist=None,cutpar=None,cuttag=None,owcut=None,**exkw):
    """
    Method for converting halo IDs to group numbers
    """
    singval=False
    # Pars input
    askuser=mmlpars.mml_pars(askuser,default=False,type=bool)
    verbose=mmlpars.mml_pars(verbose,default=False,type=bool)
    outflags=mmlpars.mml_pars(outflags,default=True,type=bool)
    # Get snapshot
    if not isinstance(snapnum,int) and askuser:
        snapnum=mmlio.askquest('Entre a snapshot number:',default=0,dtype='int')
    else:
        snapnum=mmlpars.mml_pars(snapnum,type=int)
    # Pars input list/variable
    if askuser and haloid==None:
        fin=False ; haloid=[]
        while not fin:
            iid=mmlio.askquest('Enter a haloid to get the last interaction for (-1 to stop):',default=-1L,dtype='long')
            if iid==-1: 
                if len(haloid)==0: print 'At least one haloid must be selected...'
                else: fin=True
            else     : haloid.append(iid)
        if len(haloid)==1: singval=True
    if not isinstance(haloid,list): 
        haloid=[haloid]
        singval=True
    # Load data
    if intlist==None: intlist=simstat.mkcuts(simstr,cutpar=cutpar,cuttag=cuttag,overwrite=owcut,askuser=askuser)
    # Sort interactions by ending redshift
    idxkeep=(np.array(intlist['isnap'])<=snapnum)
    for ikey in intlist['keylist']: intlist[ikey]=list(np.array(intlist[ikey])[idxkeep])
    idxsort=np.argsort(np.array(intlist['fsnap']))[::-1]
    for ikey in intlist['keylist']: intlist[ikey]=list(np.array(intlist[ikey])[idxsort])
    hidarr=np.array(intlist['haloid1'])
    # ID halos
    outvar={ikey:[] for ikey in intlist['keylist']}
    outvar['keylist']=intlist['keylist']
    flaglist=[]
    for ihaloid in haloid:
        if ihaloid in intlist['haloid1']:
            ihidx=np.arange(len(hidarr))[hidarr==ihaloid]
            iisnap=np.array(intlist['isnap'])[ihidx]
            ifsnap=np.array(intlist['fsnap'])[ihidx]
            idxcur=ihidx[np.logical_and(iisnap<=snapnum,ifsnap>=snapnum)] ; ncur=len(idxcur)
            if ncur==0:
                flaglist.append(0)
                lastint=ihidx[ifsnap==ifsnap[0]]
                lastidx=lastint[np.argmin(np.array(intlist['q'])[lastint])]
                for ikey in intlist['keylist']: outvar[ikey].append(intlist[ikey][lastidx])
            elif ncur==1:
                flaglist.append(-2)
                if verbose: mmlio.verbose('Haloid {} undergoing an interaction at snapshot {}'.format(ihaloid,snapnum))
                lastidx=idxcur[0]
                for ikey in intlist['keylist']: outvar[ikey].append(intlist[ikey][lastidx])
            else:
                flaglist.append(-3)
                if verbose: mmlio.verbose('Haloid {} undergoing multiple interactions at snapshot {}'.format(ihaloid,snapnum))
                lastidx=idxcur[np.argmin(np.array(intlist['q'])[idxcur])]
                for ikey in intlist['keylist']: outvar[ikey].append(intlist[ikey][lastidx])
        else:
            flaglist.append(-1)
            if verbose: mmlio.verbose('Haloid {} could not be found.'.format(ihaloid))
            for ikey in intlist['keylist']: outvar[ikey].append(-1)
            outvar['haloid1'][-1]=ihaloid
    # Return output
    if singval: 
        for ikey in intlist['keylist']: outvar[ikey]=outvar[ikey][0]
    if outflags: return outvar,flaglist
    else       : return outvar

###################################################################################################################################
# METHOD FOR CONVERTING HALO IDS TO GROUP NUMBERS
def haloid2groupnum(simstr,snapnum=None,haloid=None,groupnum=None,askuser=None,**exkw):
    """
    Method for converting halo IDs to group numbers
    """
    singval=False
    kyvar={'haloid':'Groupid','groupnum':'GroupNum'}
    # Pars input
    askuser=mmlpars.mml_pars(askuser,default=False,type=bool)
    # Get snapshot
    if not isinstance(snapnum,int) and askuser:
        snapnum=mmlio.askquest('Entre a snapshot number:',default=0,dtype='int')
    else:
        snapnum=mmlpars.mml_pars(snapnum,type=int)
    # Pars variable that should be converted
    if haloid==None and groupnum==None and askuser:
        var_i=mmlio.askselect('What variable do you want to start with?',['haloid','groupnum'])
    elif haloid  !=None: var_i='haloid'
    elif groupnum!=None: var_i='groupnum'
    else: raise Exception('Neither haloid or groupnum provided.')
    if   var_i=='haloid'  : invar=haloid   ; var_o='groupnum'
    elif var_i=='groupnum': invar=groupnum ; var_o='haloid'
    # Pars input list/variable
    if askuser and invar==None:
        fin=False ; invar=[]
        while not fin:
            iid=mmlio.askquest('Enter a {} to convert (-1 to stop):'.format(var_i),default=-1L,dtype='long')
            if iid==-1: 
                if len(invar)==0: print 'At least one {} must be selected...'.format(var_i)
                else: fin=True
            else     : invar.append(iid)
        if len(invar)==1: singval=True
    if not isinstance(invar,list): 
        invar=[invar]
        singval=True
    # If its the last snapshot haloid=groupnum
    if snapnum==(len(simstr.fdict['redlist'])-1): outvar=invar
    # Otherwise use info from table
    else:
        # Load data
        ftable=simstr.fdict['tabbase'].format(snapnum)
        dtable=simfile.rwtable('R',ftable,'snaptab')
        # ID halos
        outvar=[]
        for iinvar in invar:
            if iinvar in dtable[kyvar[var_i]]:
                iidx=dtable[kyvar[var_i]].index(iinvar)
                outvar.append(dtable[kyvar[var_o]][iidx])
            else:
                mmlio.verbose('{} {} could not be found.'.format(var_i,iinvar))
                outvar.append(-1)
    # Return output
    if singval: return outvar[0]
    else      : return outvar


###################################################################################################################################
###################################################################################################################################
# PROVIDE COMMAND LINE ACCESS
if __name__ == '__main__': main()
