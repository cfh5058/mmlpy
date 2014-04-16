#!/usr/bin/python    
###################################################################################################################################
#
# MEAGAN LANG'S SIMLYZE METHODS
#
####################################################################################################################################
#import matplotlib
#matplotlib.use('Agg')
#import matplotlib.pyplot as plt
import sys,os,shutil,glob,copy,pprint,scipy,math
import numpy as np
import pNbody
import cPickle
LIST_METHODS=['all','stat','calc','hist','plot','plotsnap','plotrun']
from mmlutils import *
import main as mmlnbody
import simstat,simcalc,simhist,simprof,simplot,simfile,mmlbuildgal,mmlgadget,mmlscf,mmlgalfit,mml_pNbody

def main():
    """
    Provides command line access to CALC methods
    """
    out = mmlnbody.walk(mtype='analyze')
    return out

####################################################################################################################################
####################################################################################################################################
# METHODS FOR GETTING LISTS
def dict_methods():
    """
    Returns a dictionary of supported ANALYZE methods for analyze.run
    """
    methDICT={}
    methDICT['stat']=['all','stat']+['stat_'+imeth for imeth in simstat.LIST_METHODS]
    methDICT['calc']=['all','calc']+['calc_'+imeth for imeth in simcalc.LIST_METHODS]
    methDICT['hist']=['all','hist']+['hist_'+imeth for imeth in simhist.LIST_METHODS]
    methDICT['plotsnap']=['all','plot','plotsnap']+['plotsnap_'+imeth for imeth in simplot.LIST_METHODS_SNP]
    methDICT['plotrun' ]=['all','plot','plotrun' ]+['plotrun_' +imeth for imeth in simplot.LIST_METHODS_RUN]
    return methDICT
def list_fullmeths():
    """
    Returns a complete list of supported ANALYZE methods and sumbethods
    """
    methLIST=copy.deepcopy(LIST_METHODS)
    methDICT=dict_methods()
    for ikey in methDICT.keys(): methLIST+=methDICT[ikey]
    return methLIST

####################################################################################################################################
####################################################################################################################################
# HIGH LEVEL METHODS REQUIRING MMLSIM

####################################################################################################################################
# METHOD FOR RUNNING DIFFERENT CALC METHODS
def run(simstr,method,**method_kw):
    """
    Provides interface for running different CALC methods
    """
    # Set constants
    methLIST=LIST_METHODS
    # Pars input
    method=mmlpars.mml_pars(method,list=methLIST)
    # Initialize default output
    out=None
    # Proceed based on method
    if method in methLIST:
        analyze(simstr,method=method,**method_kw)
    else:
        raise Exception('Option for method {} needs to be added to simcalc.run.'.format(method))
    # Return output
    return out

####################################################################################################################################
####################################################################################################################################
# RUN METHODS
                                                                    
####################################################################################################################################
# METHOD TO PERFORM SIMULATION CALCULATIONS
def analyze(simstr,method=None,ftype=None,ptype=None,overwrite=None,verbose=None,shortlist=None,
            calcdict0=None,idgal2=None,plot_kw={},anim_kw={},hist_fparfile=None,stat_fparfile=None,**extra_kw):
    """
    NAME:
        mmlnbody.main.analyze
    PURPOSE:
        To run calculations on a Nbody simulation.
    CALLING:
        analyze(simstr[,method=,ftype=,ptype=,overwrite=,verbose=,shortlist=,calcdict0=,
                idgal2=,scielist=,displist=,ellilist=])
    ARGUMENTS:
        simstr:    mmlsim object describing basic parameters of the simulation
    KEYWORDS:
        method:    String specifying a specific part of the analysis to perform.
        ftype:     String specifying what type of simulation it was.
        ptype:     String specifying what type of particle snapshot to load
                   [Only valide for ftype='scf']
        overwrite: Bool specifying whether or not to overwrite existing files
        verbose:   Bool specifying if additional output should be printed to the screen
        shortlist: Bool specifying if analysis should be limited to IC file + most recent snapshot
        calcdict0: Dictionary containing info on calcdict defaults
        idgal2:    Long minimum ID of second galaxy
    """
    # Set constants
    verbwid=150
    methodLIST=list_fullmeths()
    # Pars input
    method=mmlpars.mml_pars(method,default='all',list=methodLIST)
    ftype=mmlpars.mml_pars(ftype,default='gadget',type=str)
    ptype=mmlpars.mml_pars(ptype,default='disk',type=str)
    verbose=mmlpars.mml_pars(verbose,default=True,type=bool)
    method=method.lower()
    ftype=ftype.lower()
    ptype=ptype.lower()
    # Initialize stuff
    statdict={}
    histdict=None ; histlim=None
    flagdict=dict(stat={},calc={},hist={},plotsnap={},plotrun={})
    # Get flag
    if not isinstance(shortlist,bool):
        shortlist=not mmlio.yorn('Perform method {} for all snapshots?'.format(method))
    # Get file names
    files=simstr.fdict
    # Get list of snapshots & icfile
    if shortlist:
        validfile=False
        while not validfile:
            fext=mmlio.askquest('What snapshot extension should be loaded?',default='_000ic',dtype='str')
            fsnap,flag_ic=simstr.get_snapname(fext=fext,ftype=ftype,ptype=ptype)
            validfile=os.path.isfile(fsnap)
        snaplist=[fsnap]
    else:
        if   ftype=='gadget':
            snaplist=sorted(glob.glob(files[ftype]['output']['snapbase']+'*'),reverse=True)
            icfile=files[ftype]['input']['ic']
        elif ftype=='scf':
            snaplist=sorted(glob.glob(files[ftype][ptype]['output']['snapbase']+'*'),reverse=True)
            icfile=snaplist.pop()
        elif ftype=='buildgal':
            raise Exception('Buildgal analysis is a work in progress.')
            snaplist=[]
            icfile=files[ftype]['output']['ic']
        else: raise Exception('Invalid file type: {}'.format(ftype))
        # Determine if IC file should be included
        if os.path.isfile(icfile): snaplist.insert(0,icfile)
    # Raise error if no files found
    if len(snaplist)==0: raise Exception('There are not any IC/snapshot files for the specified run.')
    # Ask user what to do
    if isinstance(overwrite,bool):
        if overwrite: overwrite=mmlio.yorn('Are you sure you want to mass overwrite all {} files?'.format(method.upper()),width=verbwid)
    # Get keys
    methdict=dict_methods()
    # Get flags for each method
    loadkeys=dict(fname=snaplist[0],ftype=ftype,ptype=ptype,idgal2=idgal2)
    deftkeys=dict(setflags=True,overwrite=overwrite,method=method)
    if method in methdict['stat'    ]: flagdict['stat'    ]=snapstat(simstr,loadkeys,flagdict=flagdict['stat'    ],fparfile=stat_fparfile,**deftkeys)
    if method in methdict['calc'    ]: flagdict['calc'    ]=snapcalc(simstr,loadkeys,flagdict=flagdict['calc'    ],**deftkeys)
    if method in methdict['hist'    ]: flagdict['hist'    ]=snaphist(simstr,loadkeys,flagdict=flagdict['hist'    ],fparfile=hist_fparfile,**deftkeys)
    if method in methdict['plotsnap']: flagdict['plotsnap']=snapplot(simstr,loadkeys,flagdict=flagdict['plotsnap'],**deftkeys)
    if method in methdict['plotrun' ]: flagdict['plotrun' ]= runplot(simstr,loadkeys,flagdict=flagdict['plotrun' ],**deftkeys)
    # Loop over snapshots
    if method in methdict['stat']+methdict['calc']+methdict['hist']+methdict['plotsnap']:
        for isnapfile in snaplist:
            # Stop after 2 snapshots if shortlist option selected
#            if shortlist and isnapfile==snaplist[2]: break
            # Print intro output to screen
            if verbose: mmlio.verbose('Beginning to analyze: ',addval=os.path.basename(isnapfile),
                                      border=True,width=verbwid,valpad=50)
            # Make ID string
            fext=simstr.get_fext(isnapfile,ftype=ftype,ptype=ptype)
            if verbose: mmlio.verbose(fext)
            # Initialize nbody objects
            inb=None ; calcdict=None
            loadkeys=dict(fname=isnapfile,ftype=ftype,ptype=ptype,idgal2=idgal2)
            # Statistics
            if method in methdict['stat']:
                flagdict['stat'],inb,statdict = snapstat(simstr,loadkeys,flagdict=flagdict['stat'],inb=inb,method=method,
                                                         statdict=statdict,verbose=verbose,verbwid=verbwid)
            # Calculations
            if method in methdict['calc']:
                flagdict['calc'],inb,calcdict,calcdict0 = snapcalc(simstr,loadkeys,flagdict=flagdict['calc'],inb=inb,
                                                                   calcdict=calcdict,calcdict0=calcdict0,method=method,
                                                                   verbose=verbose,verbwid=verbwid)
            # Histograms
            if method in methdict['hist']:
                flagdict['hist'],inb,histdict,histlim = snaphist(simstr,loadkeys,flagdict=flagdict['hist'],inb=inb,histdict=histdict,
                                                                 limdict=histlim,method=method,verbose=verbose,verbwid=verbwid)
            # Plots
            if method in methdict['plotsnap']:
                flagdict['plotsnap'],inb,calcdict,calcdict0 = snapplot(simstr,loadkeys,flagdict=flagdict['plotsnap'],inb=inb,
                                                                       calcdict=calcdict,calcdict0=calcdict0,method=method,
                                                                       verbose=verbose,verbwid=verbwid,**plot_kw)
    # Create plots
    if method in methdict['plotrun']:
        # Plot
        flagdict['plotrun'],statdict=runplot(simstr,loadkeys,flagdict=flagdict['plotrun'],statdict=statdict,histdict=histdict,
                                             verbose=verbose,verbwid=verbwid,stat_fparfile=stat_fparfile,**plot_kw)
    return

####################################################################################################################################
# METHOD TO HANDLE CREATING CALC FILES
def snapcalc(simstr,loadkeys,flagdict=None,inb=None,calcdict=None,calcdict0=None,
             setflags=None,overwrite=None,verbose=None,verbwid=None,method=None,
             Rscale=None,profmode=None,profspace=None,**calc_kw):
    # Set constants
    methlistVAL=dict_methods()['calc']
    ptyplistDEF=['visi']
    viewlistDEF=simplot.list_plotviews()
    pgallistDEF=[1]
    RscaleDEF='log'
    profmodeDEF='spin'
    profspaceDEF='R'
    # Pars input
    loadkeys=simstr.get_loadkeys(**loadkeys)
    method=mmlpars.mml_pars(method,default='calc_all',list=methlistVAL)
    flagdict=mmlpars.mml_pars(flagdict,default={},type=dict)
    setflags=mmlpars.mml_pars(setflags,default=False,type=bool)
    verbose=mmlpars.mml_pars(verbose,default=False,type=bool)
    ptyplist=mmlpars.mml_pars(ptyplist,default=ptyplistDEF,type=list)
    pgallist=mmlpars.mml_pars(pgallist,default=pgallistDEF,type=list)
    viewlist=mmlpars.mml_pars(viewlist,default=viewlistDEF,type=list)
    Rscale=mmlpars.mml_pars(Rscale,default=RscaleDEF)
    profmode=mmlpars.mml_pars(profmode,default=profmodeDEF)
    profspace=mmlpars.mml_pars(profspace,default=profspaceDEF)
    if setflags: verbose=False
    # Find calc method
    if method.startswith('calc_'): calcmethod=method.split('calc_')[-1]
    else                         : calcmethod='all'
    # Get file dictionary
    filedict=simstr.fdict
    # Set flags if not provided
    if isinstance(overwrite,bool):
        for imeth in simcalc.LIST_METHODS_SNP: flagdict[imeth]=overwrite
    if calcmethod in ['all','calcsnap']:
        if 'stat' not in flagdict: flagdict['stat']=mmlio.yorn('Overwrite existing static calc file?',width=verbwid)
    for imeth in simcalc.LIST_METHODS_SNP:
        if calcmethod in ['all',imeth]:
            if imeth not in flagdict: flagdict[imeth]=mmlio.yorn('Overwrite existing {} files?'.format(imeth),width=verbwid)
    if not setflags:
        # Calcsnap
        if calcmethod in ['all','calcsnap']:
            # Only continue if overwrite selected or file does not exist
            calcfile=simcalc.snapfiles(simstr,'calcsnap',fext=loadkeys['fext'])['file']
            if flagdict['calcsnap'] or not os.path.isfile(calcfile):
                # Get info
                calcdict,calcdict0,inb=simcalc.get_calcsnap(simstr,snapdict=inb,calcdict0=calcdict0,
                                                            overwrite=flagdict['calcsnap'],owstat=flagdict['stat'],
                                                            fulloutput=True,**loadkeys)
                if verbose: mmlio.verbose('Created calc file:')
                if verbose: print '    '+calcfile
        # Ellipse fit
        if calcmethod in ['all','ellipfit']:
            for iptyp in ptyplist:
                for ipgal in pgallist:
                    for iview in viewlist:
                        # Only continue if overwrite selected or file does not exist
                        if 'mode' not in calc_kw: calc_kw['mode']='pos_m'
                        efitfile=simcalc.snapfiles(simstr,'ellipfit',iptyp,calc_kw['mode'],ipgal,iview,fext=loadkeys['fext'])['file'] #+'.txt'
                        if flagdict['ellipfit'] or not os.path.isfile(efitfile):
                            icalckw=dict(calc_kw,type=iptyp,galaxy=ipgal,view=iview,shape=100,**loadkeys)
                            calcdict,calcdict0,inb=simcalc.get_ellipfit(simstr,snapdict=inb,calcdict0=calcdict0,
                                                                        overwrite=flagdict['ellipfit'],**icalckw)
                            if verbose: mmlio.verbose('Create ellipfit file:')
                            if verbose: print '    '+efitfile
        # Profile
        if calcmethod in ['all','profile']:
            for ipgal in pgallist:
                # Only continue if overwrite selected or file does not exist
                fpar=dict(space=profspace,scale=Rscale,mode=profmode,galaxy=ipgal)
                proffile=simfile.fpar2file(simstr,'prof','profile',fext=loadkeys['fext'],**fpar)['file']
                if flagdict['profile'] or not os.path.isfile(proffile):
                    # Get info
                    icalckw=dict(fpar,fulloutput=True,**loadkeys)
                    iprofdict,inb=simprof.get_profile(simstr,overwrite=flagdict['profile'],**icalckw)
                    if verbose: mmlio.verbose('Created profile file:')
                    if verbose: print '    '+proffile
    # Return info
    if setflags:
        return flagdict
    else:
        return flagdict,inb,calcdict,calcdict0
    
####################################################################################################################################
# METHOD TO HANDLE CREATING CALC FILES
def snaphist(simstr,loadkeys,flagdict=None,inb=None,method=None,histdict=None,limdict=None,
             setflags=None,overwrite=None,verbose=None,verbwid=None,fparfile=None):
    """
    Creates a histogram for one snapshot
    """
    # Set constants
    methlistVAL=dict_methods()['hist']
    # Pars input
    loadkeys=simstr.get_loadkeys(**loadkeys)
    method=mmlpars.mml_pars(method,default='hist_all',list=methlistVAL)
    flagdict=mmlpars.mml_pars(flagdict,default={},type=dict)
    limdict=mmlpars.mml_pars(limdict,default={},type=dict)
    setflags=mmlpars.mml_pars(setflags,default=False,type=bool)
    verbose=mmlpars.mml_pars(verbose,default=False,type=bool)
    if setflags: verbose=False
    # Find hist method
    if method.startswith('hist_'): histmethod=method.split('hist_')[-1]
    else                         : histmethod='all'
    # Select plot type
    if histmethod in simhist.LIST_METHODS_PLT:
        hmeth=histmethod[-2:]
        if 'plotmeth' not in flagdict: flagdict['plotmeth']=simhist.askplot(hmeth)
        pmeth=flagdict['plotmeth']
    else:
        hmeth=histmethod
        pmeth=None
    if pmeth==None: plot=False
    else          : plot=True
    # Set flags if not provided
    if isinstance(overwrite,bool):
        for imeth in simhist.LIST_METHODS_SNP: flagdict[imeth]=overwrite
    for imeth in simhist.LIST_METHODS_SNP:
        if histmethod in ['all',imeth]:
            if imeth not in flagdict: flagdict[imeth]=mmlio.yorn('Overwrite existing {} files?'.format(imeth),width=verbwid)
            if plot and hmeth not in flagdict: flagdict[hmeth]=mmlio.yorn('Overwrite existing {} files?'.format(hmeth),width=verbwid)
    # Initialize histdict if not provided
    if histdict==None: 
        if 'histdict' in flagdict: histdict=flagdict['histdict']
        else                     : histdict=simfile.parsfpar('hist',hmeth,askuser=True,tagstr=fparfile,plot=plot)
    # Create histogram/plot
    if not setflags:
        # Histogram
        if not plot:
            # Only continue if overwrite selected or file does not exist
            histfile=simhist.snapfiles(simstr,histmethod,histdict,fext=loadkeys['fext'])['file']
            if flagdict[histmethod] or not os.path.isfile(histfile):
                # Get info
                histdict,inb=simhist.get_histsnap(simstr,histmethod,histdict=histdict,
                                                  overwrite=flagdict[histmethod],fulloutput=True,askuser=False,
                                                  snapdict=inb,**loadkeys)
                if verbose: mmlio.verbose('Created hist file:')
                if verbose: print '    '+histfile
        # Plot
        elif 'plot' in histmethod:
            # Only continue if overwrite selected or file does not exist
            plotfile=simhist.snapfiles(simstr,hmeth,histdict,fext=loadkeys['fext'],plot=True,plotmeth=pmeth)['file']
            if flagdict[histmethod] or not os.path.isfile(plotfile):
                # Get plot
                histdict,limdict,inb=simhist.plothist(simstr,hmeth,pmeth,histdict=histdict,limdict=limdict,
                                                      owplot=flagdict[histmethod],owhist=flagdict[hmeth],askuser=False,
                                                      snapdict=inb,loadkeys=loadkeys)
                if verbose: mmlio.verbose('Created hist plot:')
                if verbose: print '    '+plotfile
    # Return info
    if setflags:
        flagdict['histdict']=histdict
        return flagdict
    else:
        return flagdict,inb,histdict,limdict
    
####################################################################################################################################
# METHOD TO HANDLE PLOTTING FOR EACH SNAPSHOT
def snapplot(simstr,loadkeys,flagdict=None,inb=None,calcdict=None,calcdict0=None,
             setflags=None,overwrite=None,verbose=None,verbwid=None,method=None,
             parlist=None,**plot_kw):
    # Set constants
    methlistVAL=dict_methods()['plotsnap']
    plotmethVAL=simplot.LIST_METHODS_SNP
    parlistDEF={}
    for imeth in plotmethVAL:
        if   imeth=='science': parlistDEF[imeth]={'type':['mult']}
        elif imeth=='profile': parlistDEF[imeth]={'space':['R'],'scale':['log'],'mode':['spin'],'galaxy':[1]}
        elif imeth=='contour': parlistDEF[imeth]={'type':['visi'],'mode':['pos_scfpotm2'],'view':['xy','xz'],'galaxy':[1]}
        else                 : parlistDEF[imeth]={'type':['visi'],'mode':['pos_m'],'view':['xy','xz'],'galaxy':[1]}
    flagdict_calc={'snap':False,'stat':False}
    # Pars input
    flagdict=mmlpars.mml_pars(flagdict,default={},type=dict)
    setflags=mmlpars.mml_pars(setflags,default=False,type=bool)
    verbose=mmlpars.mml_pars(verbose,default=False,type=bool)
    method=mmlpars.mml_pars(method,default='plotsnap_all',list=methlistVAL)
    if setflags: verbose=False
    # Create file dictionary, file extension, & plotting windows
    filedict=simstr.fdict
    fext=simstr.get_fext(loadkeys['fname'],ftype=loadkeys['ftype'],ptype=loadkeys['ptype'])
    calcfile=simcalc.snapfiles(simstr,'calcsnap',fext=fext)['file']
    calcfile0=filedict['calc']['calcstat']
    # Handle flags for isolating plots
    if method.startswith('plotsnap_'): plotmethod=method.split('plotsnap_')[-1]
    else                             : plotmethod='all'
    for imeth in plotmethVAL:
        if imeth=='profile': pass
        else               : 
            if plotmethod.split(imeth)[-1] in simfile.list_fpar('plots',imeth)['type']:
                parlistDEF[imeth]['type']=[plotmethod.split(imeth)[-1]]
    # Set flags
    flagdict=simplot.plotsnap(simstr,method=plotmethod,flagdict=flagdict,setflags=True,overwrite=overwrite,
                              askuser=True,allsnap=False,**loadkeys)
    # Plots
    if not setflags:
        for isnplot in plotmethVAL:
            # Don't plot if flag not set
            if not flagdict[isnplot]: continue
            # Fill in parameters
            fpar=parlistDEF[isnplot]
            for ikey in fpar.keys():
                if ikey in plot_kw: fpar[ikey]=[plot_kw[ikey]]
            # Loop over list of plots
            iterpar=simfile.iter_fpar('plots',isnplot,fpar=fpar)
            for ifpar in iterpar:
                iplotkw=dict(plot_kw,loadcalc=False,**loadkeys)
                for ikey in ifpar.keys(): iplotkw[ikey]=ifpar[ikey]
                if isnplot=='science': iplotkw['loadcalc']=True
                else                 : iplotkw['shape']=512
                out=simplot.plotsnap(simstr,method=imeth,snapdict=inb,calcdict=calcdict,calcdict0=calcdict0,
                                     flagdict=flagdict,askuser=False,verbose=verbose,verbwid=verbwid,**iplotkw)
                flagdict,inb,calcdict,calcdict0=out
    # Return info
    if setflags:
        return flagdict
    else:
        return flagdict,inb,calcdict,calcdict0

####################################################################################################################################
# METHOD TO HANDLE STATISTICS FOR EACH SNAPSHOT
def snapstat(simstr,loadkeys,flagdict=None,inb=None,statdict=None,fparfile=None,
             setflags=None,overwrite=None,verbose=None,verbwid=None,method=None):
    # Set constants
    methlistVAL=dict_methods()['stat']
    # Pars input
    method=mmlpars.mml_pars(method,default='stat_all',list=methlistVAL)
    statdict=mmlpars.mml_pars(statdict,default={},type=dict)
    flagdict=mmlpars.mml_pars(flagdict,default={},type=dict)
    setflags=mmlpars.mml_pars(setflags,default=False,type=bool)
    verbose=mmlpars.mml_pars(verbose,default=False,type=bool)
    if setflags: verbose=False
    # Create file dictionary
    filedict=simstr.fdict
    # Create directory
    mmlfiles.mkdirs(filedict['stat']['dir'])
    # Isolate stat method
    if method.startswith('stat_'): statmethod=method.split('stat_')[-1]
    else                         : statmethod='all'
    if statmethod=='all': statlist=simstat.LIST_STATS
    else                : statlist=[statmethod]
    # Loop over stat files
    for istat in statlist:
        # Get parameters
        ifpar_deft=simstat.list_fparDEF(istat)
        ifpar_deft['rmax']=simstr.get_Rvir(1)
        if istat+'_fpar' not in flagdict: flagdict[istat+'_fpar']=simfile.parsfpar('stat',istat,askuser=True,fpar_deft=ifpar_deft,tagstr=fparfile)
        istatfpar=flagdict[istat+'_fpar']
        if istat=='var': isubstat=istatfpar['smeth']
        else           : isubstat=None
        # Select file
        statfile=simstat.snapfiles(simstr,istat,istatfpar)['file']
#        statfile=filedict['stat']['stat_'+istat]
        # Add flag
        if isinstance(overwrite,bool): flagdict[istat]=overwrite
        if not istat in flagdict: flagdict[istat]=mmlio.yorn('Overwrite existing {} statistics entries?'.format(istat),width=verbwid)
        if not setflags:
            # Load statistics dictionary if not already
            if not istat in statdict: statdict[istat]=simstat.rwstats('R',stattype=istat,substattype=isubstat,statfile=statfile)
            # Do statistics as required
            if loadkeys['fname'] not in statdict[istat]['fname'] or flagdict[istat]:
                # Load pNbody object if not provided
                if inb==None: inb=simstr.get_pnbody(**loadkeys)
                # Prevent scfm2 stats if file not there
                if istat=='bar':
                    if istatfpar['smeth']=='scfm2' and inb.scfpotm2.min()==inb.scfpotm2.max(): continue
                # Prevent bulk stats if particle type not there
                if istat=='bulk':
                    if istatfpar['ptyp'] not in inb.get_typlist(): continue
                # Perform statistics
                addstatdict=inb.get_statdict(method=istat,**istatfpar)
                # Add statistics to file
                statdict[istat]=simstat.rwstats('A',stattype=istat,substattype=isubstat,statfile=statfile,statdict=statdict[istat],
                                                addstatdict=addstatdict,owentry=flagdict[istat])
                if verbose: mmlio.verbose('Added {} statistics to file: {}'.format(istat,os.path.basename(statfile)))
    # Return output
    if setflags: return flagdict
    else       : return flagdict,inb,statdict

####################################################################################################################################
# METHOD FOR RUN PLOTS
def runplot(simstr,loadkeys,flagdict=None,statdict=None,histdict=None,
            setflags=None,overwrite=None,verbose=None,verbwid=None,
            method=None,stat_fparfile=None,**plot_kw):
    # Set constants
    statlist=list_stats()
    methlistVAL=dict_methods()['plotrun']
    # Pars input
    method=mmlpars.mml_pars(method,default='plotrun_all',list=methlistVAL)
    flagdict=mmlpars.mml_pars(flagdict,default={},type=dict)
    statdict=mmlpars.mml_pars(statdict,default={},type=dict)
    setflags=mmlpars.mml_pars(setflags,default=False,type=bool)
    verbose=mmlpars.mml_pars(verbose,default=False,type=bool)
    if setflags: verbose=False
    # Create file dictionary
    filedict=simstr.fdict
    # Isolate plot method
    if method.startswith('plotrun_'): plotmethod=method.split('plotrun_')[-1]
    else                            : plotmethod='all'
    # Print some info to the screen
    if verbose: mmlio.verbose('Beginning total run plots',border=True,width=verbwid)
    # Gadget specific plots
    if plotmethod in ['all','energy']:
        if loadkeys['ftype'] == 'gadget':
            flagdict=simplot.plotrun_energy(simstr,flagdict=flagdict,setflags=setflags,overwrite=overwrite,
                                            verbose=verbose,verbwid=verbwid,**plot_kw)
    # Plot statistics
    if plotmethod in ['all','stats']:
        out=simplot.plotrun_stats(simstr,statdict=statdict,flagdict=flagdict,fparfile=stat_fparfile,
                                  setflags=setflags,overwrite=overwrite,
                                  verbose=verbose,verbwid=verbwid,**plot_kw)
        if setflags: flagdict=out
        else       : flagdict,statdict=out
    # Create animation
    if plotmethod in ['all','anim']:
        flagdict['anim']=simplot.plotrun_anim(simstr,loadkeys,flagdict=flagdict['anim'],
                                              verbose=verbose,verbwid=verbwid,overwrite=overwrite,
                                              **plot_kw)
    # Return info
    if setflags: return flagdict
    else       : return flagdict,statdict

###################################################################################################################################
###################################################################################################################################
# PROVIDE COMMAND LINE ACCESS
if __name__ == '__main__': main()
