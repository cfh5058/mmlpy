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
import pNbody
LIST_STATS=['general','profile','bar','bulk','com','var']
LIST_METHODS_RUN=LIST_STATS
LIST_METHODS_SNP=[]
LIST_METHODS=LIST_METHODS_RUN+LIST_METHODS_SNP
DICT_VARSTAT={'angmom':['Ltot','Lx','Ly','Lz','spin'],
              'spin'  :['spintot']}
from mmlutils import *
import main as mmlnbody
import simlist,simlyze,simplot,simfile

def main():
    """
    Provides command line access to STAT methods
    """
    out = mmlnbody.walk(mtype='stat')
    return out

####################################################################################################################################
####################################################################################################################################
# METHODS FOR GETTING LISTS
def list_varstat(method):
    """
    Returns a list of variables grouped for a method
    """
    method=mmlpars.mml_pars(method,type=str)
    if method in DICT_VARSTAT: varlist=DICT_VARSTAT[method]
    else                     : varlist=[method]
    return varlist
def list_ppar(method):
    """
    Returns a list of supported plot parameters
    """
    method=mmlpars.mml_pars(method,list=LIST_STATS)
    oppar={}
    ppar=dict(oppar,**ippar)
    return ppar
def list_fpar(method,plot=None):
    """
    Returns a list of supported parameters
    """
    method=mmlpars.mml_pars(method,list=LIST_STATS)
    plot=mmlpars.mml_pars(plot,default=False,type=bool)
    ofpar={}
    if   method=='general': ifpar={'ptyp'  :simlist.LIST_PTYPS      ,
                                   'cmeth' :simlist.LIST_COMMETHS   ,
                                   'ctyp'  :simlist.LIST_COMPTYPS   ,
                                   'cgal'  :simlist.LIST_COMPGALS   }
    elif method=='profile': ifpar={'ptyp'  :simlist.LIST_PTYPS      ,
                                   'cmeth' :simlist.LIST_COMMETHS   ,
                                   'ctyp'  :simlist.LIST_COMPTYPS   ,
                                   'cgal'  :simlist.LIST_COMPGALS   }
    elif method=='bar'    : ifpar={'ptyp'  :simlist.LIST_PTYPS      ,
                                   'pgal'  :simlist.LIST_PGALS      ,
                                   'smeth' :['dens','scfm2','ellip']}
    elif method=='bulk'   : ifpar={'ptyp'  :simlist.LIST_PTYPS      ,
                                   'pgal'  :simlist.LIST_PGALS      ,
                                   'cmeth' :simlist.LIST_COMMETHS   ,
                                   'rmax'  :float                   }
    elif method=='com'    : ifpar={'cmeth' :simlist.LIST_COMMETHS   ,
                                   'cmeth0':simlist.LIST_COMMETHS   ,
                                   'cgal0' :simlist.LIST_COMPGALS   ,
                                   'ctyp0' :simlist.LIST_COMPTYPS   ,
                                   'rmax'  :float                   }
    elif method=='var'    : ifpar={'smeth' :DICT_VARSTAT.keys()     ,
                                   'cmeth' :simlist.LIST_COMMETHS   ,
                                   'ctyp'  :simlist.LIST_COMPTYPS   ,
                                   'cgal'  :simlist.LIST_COMPGALS   ,
                                   'rmax'  :float                   }
    else: raise Exception('Invalid method: {}'.format(method))
    if plot:
        if   method=='var': ippar={'rowvar':['gal','typ']}
        else              : ippar={}
        ifpar=dict(ifpar,**ippar)
    fpar=dict(ofpar,**ifpar)
    return fpar
def list_fparDEF(method,plot=None):
    """
    Returns a dictionary with default parameter values
    """
    method=mmlpars.mml_pars(method,list=LIST_STATS)
    plot=mmlpars.mml_pars(plot,default=False,type=bool)
    ofpar={}
    if   method=='general': ifpar={'ptyp'  :'all'   ,
                                   'cmeth' :'mass'  ,
                                   'ctyp'  :'all'   ,
                                   'cgal'  :1       }
    elif method=='profile': ifpar={'ptyp'  :'all'   ,
                                   'cmeth' :'mass'  ,
                                   'ctyp'  :'all'   ,
                                   'cgal'  :1       }
    elif method=='bar'    : ifpar={'ptyp'  :'all'   ,
                                   'pgal'  :1       ,
                                   'smeth' :'ellip' }
    elif method=='bulk'   : ifpar={'ptyp'  :'all'   ,
                                   'pgal'  :1       ,
                                   'cmeth' :'mass'  ,
                                   'rmax'  :100.0   }
    elif method=='com'    : ifpar={'cmeth' :'mass'  ,
                                   'cmeth0':'mass'  ,
                                   'cgal0' :1       ,
                                   'ctyp0' :'all'   ,
                                   'rmax'  :100.0   }
    elif method=='var'    : ifpar={'smeth' :'angmom',
                                   'cmeth' :'mass'  ,
                                   'ctyp'  :'all'   ,
                                   'cgal'  :1       ,
                                   'rmax'  :100.0   }
    else: raise Exception('Invalid method: {}'.format(method))
    if plot:
        if   method=='var': ippar={'rowvar':'gal'}
        else              : ippar={}
        ifpar=dict(ifpar,**ippar)
    fpar=dict(ofpar,**ifpar)
    return fpar
def fpar2tag(method,inpar,plot=None):
    """
    Returns a file tag given the set of parameters
    """
    method=mmlpars.mml_pars(method,list=LIST_STATS)
    plot=mmlpars.mml_pars(plot,default=False,type=bool)
    outpar=simfile.parsfpar('stat',method,inpar,plot=plot)
    if 'rmax' in outpar: outpar['rmaxstr']=mmlstring.dec2str(outpar['rmax'])
    if   method=='general': tagstr='{ptyp}_c{ctyp}{cgal}{cmeth}'.format(**outpar)
    elif method=='profile': tagstr='{ptyp}_c{ctyp}{cgal}{cmeth}'.format(**outpar)
    elif method=='bar'    : tagstr='{smeth}_{ptyp}{pgal}'.format(**outpar)
    elif method=='bulk'   : tagstr='{ptyp}{pgal}_{cmeth}to{rmaxstr}'.format(**outpar)
    elif method=='com'    : tagstr='{cmeth}to{rmaxstr}_rel2{ctyp0}{cgal0}{cmeth0}'.format(**outpar)
    elif method=='var'    : tagstr='{smeth}_c{ctyp}{cgal}{cmeth}'.format(**outpar)
    else: raise Exception('Invalid method: {}'.format(method))
    if plot:
        if   method=='var': pltstr='_{rowvar}'.format(**outpar)
        else              : pltstr=''
        tagstr+=pltstr
    return tagstr
def dict_ftype():
    """
    """
    fdict={'shortout':[iclass+'stats' for iclass in LIST_STATS],
           'longout' :LIST_METHODS_SNP}
    return fdict
def dict_fgroup():
    """
    """
    fdict={'clean' :['shortout','longout'],
           'input' :['shortin' ,'longin' ],
           'output':['shortout','longout']}
    return fdict

####################################################################################################################################
####################################################################################################################################
# HIGH LEVEL METHODS REQUIRING MMLSIM

####################################################################################################################################
# METHOD FOR RUNNING DIFFERENT STAT METHODS
def run(simstr,method,**method_kw):
    """
    Provides interface for running different STAT methods
    """
    # Set constants
    methLIST=LIST_METHODS
    # Pars input
    method=mmlpars.mml_pars(method,list=methLIST)
    # Initialize default output
    out=None
    # Proceed based on method
    if method in methLIST:
        simlyze.analyze(simstr,method='stat_'+method,**method_kw)
    else:
        raise Exception('Option for method {} needs to be added to simcalc.run.'.format(method))
    # Return output
    return out

####################################################################################################################################
####################################################################################################################################
# FILE CLASSES AND METHODS
                                                                    
####################################################################################################################################
# METHOD TO RETURN DICTIONARY OF SIMCALC FILES FOR CALCULATING SNAPSHOT STATISTICS 
def snapfiles(simstr,method,fpar,plot=None,**exkw):
    """
    Returns a dictionary containing a snapshot hist directory and file base 
    """
    snfdict=simfile.fpar2file(simstr,'stat',method,fpar,sing=True,plot=plot,**exkw)
    return snfdict

####################################################################################################################################
# METHOD TO RETURN CALC FILES
def files(fdict):
    """
    Returns a dictionary of STAT files & directories
    """
    # Set constants
    statlist=LIST_METHODS_RUN
    # Pars input
    fdict=mmlpars.mml_pars(fdict,type=dict)
    pfix=fdict['pfix']
    # Initialize dictionary with directory & default/static files
    files=fdict['stat']
    # Add input files
    fkeys=['stat_'+istat for istat in statlist]
    flist=['stat_'+istat for istat in statlist]
    # Directories
    dkeys=[]
    dlist=[]
    files['subdirs']=dkeys
    files=mmlfiles.filedict(files['dir'],dirkeys=dkeys,dirlist=dlist,filekeys=fkeys,filelist=flist,pfix=pfix,fdict=files)
    # Return file dictionary
    fdict['stat']=files
    return fdict

####################################################################################################################################
# METHOD TO RETURN LIST OF FILES
def get_filelist(fdict,keylist):
    """
    Returns a list of relavent CALC files
    """
    # Pars input
    fdict=mmlpars.mml_pars(fdict,type=dict)
    keylist=mmlpars.mml_pars(keylist,type=list)
    # Add files
    filelist={}
    for iclass in LIST_STATS:
        iclassstr=iclass+'stats'
        istatlist=list_stattypes(iclass)
        if iclassstr in keylist:
            filelist[iclassstr]=[fdict['stat_'+istat] for istat in istatlist]
    # Return output
    return filelist

####################################################################################################################################
# METHOD TO PERFORM STATISTICS
def rwstats(method,stattype=None,substattype=None,statfile=None,statdict=None,addstatdict=None,owentry=None):
    # Pars input
    if isinstance(method,str): method=method.upper()
    method=mmlpars.mml_pars(method,type=str,list=['R','W','A'])
    if isinstance(stattype,str): stattype=stattype.lower()
    stattype=mmlpars.mml_pars(stattype,type=str,list=LIST_STATS)
    owentry=mmlpars.mml_pars(owentry,type=bool,default=False)
    # Get list of keys
    statkeys=['time']
    if stattype=='general':
        statkeys+=['comx','comy','comz','vcomx','vcomy','vcomz','clause','sigma']
    elif stattype=='profile':
        statkeys+=['disk_A','disk_r','disk_z','halo_A','halo_r','bulge_A','bulge_r']
    elif stattype=='bar':
        statkeys+=['A','Aalt','pa','omega','r','x','y','z']
    elif stattype=='bulk':
        statkeys+=['Ep','Ek','vsig','Lx','Ly','Lz','Lx_circ','Ly_circ','Lz_circ']
    elif stattype=='com':
        posstr=['x','y','z']
        for ipgal in simlist.LIST_PGALS:
            for iptyp in simlist.LIST_PTYPS:
                for ipos in posstr:
                    statkeys+=['{}{}_{}'.format(iptyp,ipgal,ipos)]
    elif stattype=='var':
        varlist=list_varstat(substattype)
        for ivar in varlist:
            for ipgal in simlist.LIST_PGALS:
                for iptyp in simlist.LIST_PTYPS:
                    statkeys+=['{}_{}{}'.format(ivar,iptyp,ipgal)]
    statkeys+=['fname']
    # Read
    if method=='R':
        statfile=mmlpars.mml_pars(statfile,type=str)
        if os.path.isfile(statfile):
            statdict=mmlio.rwtable('R',statfile)
        else:
            statdict={ikey:[] for ikey in statkeys}
    # Append
    elif method=='A':
        addstatdict=mmlpars.mml_pars(addstatdict,type=dict)
        if statdict==None: statdict=rwstats('R',stattype=stattype,substattype=substattype,statfile=statfile)
        if addstatdict['fname'] not in statdict['fname'] or owentry:
            if addstatdict['fname'] in statdict['fname']:
                idxadd=statdict['fname'].index(addstatdict['fname'])
                for ikey in statkeys:
                    statdict[ikey][idxadd]=addstatdict[ikey]
            else:
                for ikey in statkeys:
                    statdict[ikey].append(addstatdict[ikey])
            statdict=rwstats('W',stattype=stattype,substattype=substattype,statfile=statfile,statdict=statdict)
    # Write
    elif method=='W':
        statfile=mmlpars.mml_pars(statfile,type=str)
        statdict=mmlpars.mml_pars(statdict,type=dict)
        if not os.path.isfile(statfile):
            mmlfiles.mkdirs(os.path.dirname(statfile))
        statdict['keylist']=statkeys
        mmlio.rwtable('W',statfile,statdict,overwrite=True)
    # Error
    else: raise Exception('Invalid method: {}'.format(method))
    # Return output
    return statdict

###################################################################################################################################
###################################################################################################################################
# PLOTTING METHODS

###################################################################################################################################
# METHOD TO PLOT VARIABLE STATS
def plot_varstat(method,statDict,fname=None,overwrite=None,rowvar=None,**plot_kw):
    """
    Plots variable statistics
    """
    # Set constants
    fnameTAG='plot_{}stats'.format(method)
    # Pars input
    options=mmlplot.parsopt(plot_kw)
    statDict=mmlpars.mml_pars(statDict,type=dict)
    fname=mmlpars.mml_pars(fname,default=os.path.join(options['plotdir'],fnameTAG+options['plotext']),type=str)
    if not mmlfiles.prep_write(fname,overwrite): return None
    rowvar=mmlpars.mml_pars(rowvar,list=['typ','gal'],default='gal')
    # Count variables
    galList=[] ; typList=[] ; varList=list_varstat(method)
    for igal in simlist.LIST_PGALS:
        for ityp in simlist.LIST_PTYPS:
            nonzero=False
            for ivar in varList:
                if nonzero: break
                iy=np.array(statDict['{}_{}{}'.format(ivar,ityp,igal)])
                if iy.min()!=iy.max(): nonzero=True
            if nonzero: 
                if igal not in galList: galList.append(igal)
                if ityp not in typList: typList.append(ityp)
    ntyp=len(typList) ; ngal=len(galList) ; nvar=len(varList)
    # Assign variable for rows, cols, & colors
    ncol=nvar+1
    if   rowvar=='gal': 
        nrow=ngal ; clrvar='typ' ; clrDict=mmlplot.get_typcolor(outstr=True)
    elif rowvar=='typ': 
        nrow=ntyp ; clrvar='gal' ; clrDict=mmlplot.get_galcolor()
    # Initialize figure & axes
    options['pad']=2.0
    options['xlabpos']=(0.5,-0.06)
    options['ylabpos']=(-0.1,0.5)
    options['figsize']=(5.*ncol,5.*nrow)
    fig=plt.figure()
    fig.suptitle(os.path.basename(statDict['fname'][0]).split('.')[0])
    # Find index to sort by time
    timeidx=np.argsort(np.array(statDict['time']))
    xlab='time'
    time=np.array(statDict[xlab])[timeidx]
    xlim=(time.min(),time.max())
    # Loop over galaxies
    axDict={} ; mnDict={} ; mxDict={}
    for galidx in range(ngal):
        igal=galList[galidx]
        for typidx in range(ntyp):
            ityp=typList[typidx]
            for varidx in range(nvar):
                ivar=varList[varidx]
                icol=varidx
                # Set row and color
                if   rowvar=='gal': 
                    irow=galidx ; iclr=mmlplot.get_mmlcolor(clrDict[ityp]) ; ilab=ityp ; itit=simlist.DICT_PGALS[igal]
                elif rowvar=='typ': 
                    irow=typidx ; iclr=mmlplot.get_mmlcolor(clrDict[igal]) ; ilab=simlist.DICT_PGALS[igal] ; itit=ityp
                # Create legend
                if irow not in axDict: 
                    axDict[irow]={} ; mnDict[irow]={} ; mxDict[irow]={}
                    axDict[irow]['leg']=plt.subplot(nrow,ncol,ncol*(irow+1))
                    axDict[irow]['leg'].set_xlim((-1,1))
                    axDict[irow]['leg'].set_ylim((-1,1))
                # Create axes for variable
                if icol not in axDict[irow]:
                    axDict[irow][icol]=plt.subplot(nrow,ncol,ncol*irow+icol+1)
                    mnDict[irow][icol]=float('inf')
                    mxDict[irow][icol]=float('-inf')
                # Get y data
                iykey='{}_{}{}'.format(ivar,ityp,igal)
                iy=np.array(statDict[iykey])[timeidx]
                if iy.min()==iy.max(): continue
                # Plot
                iplot=axDict[irow][icol].plot(time,iy,label=ilab,color=iclr)
                if icol==0: ilegn=axDict[irow]['leg'].plot(time[0],iy[0],label=ilab,color=iclr)
                # Limits
                mnDict[irow][icol]=min(mnDict[irow][icol],iy.min())
                mxDict[irow][icol]=max(mxDict[irow][icol],iy.max())
                axDict[irow][icol].set_xlim(xlim)
                axDict[irow][icol].set_ylim((mnDict[irow][icol],mxDict[irow][icol]))
                # Labels
                axDict[irow][icol].set_xlabel(xlab)
                axDict[irow][icol].set_ylabel(ivar)
                axDict[irow][icol].set_title(itit)
    # Create legends
    for irow in range(nrow): axDict[irow]['leg'].legend()
    # Set axes properties
    mmlplot.set_figprop(fig,options)
    fig.tight_layout(pad=options['pad'])
    # Display plot
    if options['showplot']: fig.show()
    # Save/output plot
    if options['outplot']: return fig
    else:
        fig.savefig(fname,facecolor=options['bgclr'],edgecolor=options['bgclr'])
        return None

###################################################################################################################################
def varstattrash():
    axList=[]
    for galidx in range(ngal):
        igal=galList[galidx]
        iax={}
        iax['leg']=plt.subplot(nrow,ncol,galidx*ncol+nvar+1)
        iax['leg'].set_xlim((-1,1))
        iax['leg'].set_ylim((-1,1))
        # Loop over variables
        for varidx in range(nvar):
            imax=None ; imin=None
            ivar=varList[varidx]
            iax[ivar]=plt.subplot(nrow,ncol,galidx*ncol+varidx+1)
            # Loop over particle types
            for typidx in range(ntyp):
                ityp=typList[typidx]
                iclr=mmlplot.get_mmlcolor(clrDict[ityp])
                iykey='{}_{}{}'.format(ivar,ityp,igal)
                iy=np.array(statDict[iykey])[timeidx]
                # Plot
                if iy.min()==iy.max(): continue
                iplot=iax[ivar].plot(time,iy,label=ityp,color=iclr)
                if varidx==0: ilegn=iax['leg'].plot(time[0],iy[0],label=ityp,color=iclr)
                # Get limits
                if imax==None: imax=iy.max()
                else         : imax=max(imax,iy.max())
                if imin==None: imin=iy.min()
                else         : imin=min(imin,iy.min())
            # Set labels & limits
            iax[ivar].set_xlabel(xlab)
            iax[ivar].set_ylabel(ivar)
            if imax==None or imin==None: ilim=(-1,1)
            else                       : ilim=(imin,imax)
            iax[ivar].set_xlim(xlim)
            iax[ivar].set_ylim(ilim)
        # Add limits, title & legend
        iax['leg'].legend()
        # Add axes to dictionary
        axList.append(iax)
    # Set axes properties

###################################################################################################################################
# METHOD TO PLOT COM STATS
def plot_comstat(statDict,fname=None,overwrite=None,**plot_kw):
    """
    Plots center of mass statistics
    """
    # Set constants 
    fnameTAG='plot_comstats'
    typList=simlist.LIST_PTYPS
    galList=simlist.LIST_PGALS
    posList=['xy','xz','yz']
    ntyp=len(typList)
    ngal=len(galList)
    npos=len(posList)
    ncol=npos+1
    nrow=ngal
    clrDict=mmlplot.get_typcolor(outstr=True)
    # Pars input
    options=mmlplot.parsopt(plot_kw)
    statDict=mmlpars.mml_pars(statDict,type=dict)
    fname=mmlpars.mml_pars(fname,default=os.path.join(options['plotdir'],fnameTAG+options['plotext']),type=str)
    if not mmlfiles.prep_write(fname,overwrite): return None
    # Set figure options
    options['figsize']=(15.,8.)
    # Initialize figure & axes
    fig=plt.figure()
    fig.suptitle(os.path.basename(statDict['fname'][0]).split('.')[0])
    # Find index to sort by time
    timeidx=np.argsort(np.array(statDict['time']))
    # Loop over galaxies
    axList=[]
    for galidx in range(ngal):
        igal=galList[galidx]
        iax={}
        imax=0.
        iax['leg']=plt.subplot(nrow,ncol,galidx*ncol+npos+1)
        iax['leg'].set_xlim((-1,1))
        iax['leg'].set_ylim((-1,1))
        # Loop over positions
        for posidx in range(npos):
            ipos=posList[posidx]
            iax[ipos]=plt.subplot(nrow,ncol,galidx*ncol+posidx+1)
            iax[ipos].set_xlim((-1,1))
            iax[ipos].set_ylim((-1,1))
            ixlab=ipos[0]
            iylab=ipos[1]
            # Loop over particle types
            for typidx in range(ntyp):
                ityp=typList[typidx]
                iclr=mmlplot.get_mmlcolor(clrDict[ityp])
                ix=np.array(statDict['{}{}_{}'.format(ityp,igal,ixlab)])[timeidx]
                iy=np.array(statDict['{}{}_{}'.format(ityp,igal,iylab)])[timeidx]
                # Plot
                if ix.min()==ix.max(): continue
                if iy.min()==iy.max(): continue
                imax=max(imax,abs(ix).max(),abs(iy).max())
                iplot=iax[ipos].plot(ix,iy,label=ityp,color=iclr)
                if posidx==0: ilegn=iax['leg'].plot(ix[0],iy[0],label=ityp,color=iclr)
                # Plot first point
                icirc=iax[ipos].plot(ix[0],iy[0],'o',color=iclr)
            # Set labels
            iax[ipos].set_xlabel(ixlab)
            iax[ipos].set_ylabel(iylab)
        # Add limits, title & legend
        if imax!=0:
            for ipos in posList:
                iax[ipos].set_xlim((-imax,imax))
                iax[ipos].set_ylim((-imax,imax))
            iax[posList[0]].set_title('Galaxy {}'.format(igal))
            iax['leg'].legend()
        # Add axes to dictionary
        axList.append(iax)
    # Set axes properties
    mmlplot.set_figprop(fig,options)
    fig.tight_layout()
    # Display plot
    if options['showplot']: fig.show()
    # Save/output plot
    if options['outplot']: return fig
    else:
        fig.savefig(fname,facecolor=options['bgclr'],edgecolor=options['bgclr'])
        return None

###################################################################################################################################
###################################################################################################################################
# PROVIDE COMMAND LINE ACCESS
if __name__ == '__main__': main()
