####################################################################################################################################
#
# MEAGAN LANG'S SIMPLOT METHODS
#
####################################################################################################################################
import matplotlib
import matplotlib.pyplot as plt
import sys,os,shutil,glob,copy,pprint,scipy,math,itertools
import numpy as np
import pNbody
from mpl_toolkits.axes_grid1 import make_axes_locatable
import ImageDraw,ImageFont,Image
import scipy.interpolate as interpolate
LIST_METHODS_RUN=['anim','energy','stats']
LIST_METHODS_SNP=['display','science','ellipse','imagpil','imagshw','contour','snapdif','profile']
LIST_METHODS_OTR=['retrieve']
LIST_METHODS_TST=['test'+imeth for imeth in LIST_METHODS_SNP]
LIST_METHODS=LIST_METHODS_RUN+LIST_METHODS_SNP+LIST_METHODS_OTR+LIST_METHODS_TST
from mmlutils import *
from mmlastro import mmlellipse
import main as mmlnbody
import simlist,simlyze,simstat,simcalc,simhist,simprof,simfile,mmlgadget

def main():
    """
    Provides command line access to PLOT methods
    """
    out = mmlnbody.walk(mtype='plot')
    return out

####################################################################################################################################
####################################################################################################################################
# METHODS FOR GETTING LISTS
def dict_ftype():
    """
    Returns a dictionary of files listed by type 
    """
    fdict={'shortout': LIST_METHODS_RUN+['short'+isnap for isnap in LIST_METHODS_SNP],
           'longout' : LIST_METHODS_SNP}
    return fdict
def dict_fgroup():
    """
    Returns a dictionary of file types listed by group
    """
    fdict={'input'  :['shortin','longin'],
           'output' :['shortout','longout']}
    return fdict
def list_plotmodes(method=None):
    """
    Returns a list of supported PLOT modes for 2D histograms
    """
    if method=='profile':
        modeLIST=['spin']
    else:
        modeLIST=['pos_m','pos_vr','pos_gal','pos_scfpot','pos_scfpotm2',
                  'EkLz_Lz','EpLz_Lz','EtotLz_Lz','EtotLtot_Ltot']
    return modeLIST
def list_plotptyps():
    """
    Returns a list of supported PLOT particle types for 2D histograms
    """
    typLIST=simlist.LIST_PTYPS
    return typLIST
def list_plotpgals():
    """
    Returns a list of supported PLOT galaxies for 2D histograms
    """
    galLIST=[0,1,2]
    return galLIST
def list_plotviews():
    """
    Returns a list of supported PLOT views for 2D histograms
    """
    viewLIST=['xy','xz']#,'yz']
    return viewLIST
def list_snaptypes(mtype):
    """
    Returns a list of supported snap plot types
    """
    mtype=mmlpars.mml_pars(mtype,list=LIST_METHODS_SNP)
    if   mtype=='science': typeLIST=['mult','hist1d','hist2d']
    elif mtype=='display': typeLIST=simlist.LIST_PTYPS
    elif mtype=='ellipse': typeLIST=simlist.LIST_PTYPS
    elif mtype=='imagpil': typeLIST=simlist.LIST_PTYPS
    elif mtype=='imagshw': typeLIST=simlist.LIST_PTYPS
    elif mtype=='contour': typeLIST=simlist.LIST_PTYPS
    elif mtype=='snapdif': typeLIST=simlist.LIST_PTYPS
    elif mtype=='profile': typeLIST=['log','lin']
    else: raise Exception('Invalid method class: {}'.format(mtype))
    return typeLIST
def list_snapplots(mtype):
    """
    Returns a list of supprted snap plots
    """
    typeLIST=list_snaptypes(mtype)
    return [mtype+itype for itype in typeLIST]
def list_fpar(method,anim=None):
    """
    Returns a list of supported parameters
    """
    quadlist=['display','ellipse','imagpil','imagshw','contour','snapdif']
    method=mmlpars.mml_pars(method,list=LIST_METHODS_SNP)
    if   method in quadlist: fpar=simcalc.list_fpar('ellipfit')
    elif method=='science' : fpar={'type':['mult','hist1d','hist2d']}
    elif method=='profile' : fpar=simprof.list_fpar('profile')
    else: raise Exception('Invalid method: {}'.format(method))
    return fpar
def list_fparDEF(method,anim=None):
    """
    Returns a dictionary of supported parameter defaults
    """
    quadlist=['display','ellipse','imagpil','imagshw','contour','snapdif']
    method=mmlpars.mml_pars(method,list=LIST_METHODS_SNP)
    if   method in quadlist: fpar=simcalc.list_fparDEF('ellipfit')
    elif method=='science' : fpar={'type':'mult'}
    elif method=='profile' : fpar=simprof.list_fparDEF('profile')
    else: raise Exception('Invalid method: {}'.format(method))
    return fpar
def fpar2tag(method,inpar):
    """
    Returns a file tag given the set of parameters 
    """
    quadlist=['display','ellipse','imagpil','imagshw','contour','snapdif']
    method=mmlpars.mml_pars(method,list=LIST_METHODS_SNP)
    outpar=simfile.parsfpar('plots',method,fpar=inpar)
    if   method in quadlist: tagstr='{type}_{mode}_{galaxy}{view}'.format(**outpar)
    elif method=='science' : tagstr='{type}'.format(**outpar)
    elif method=='profile' : tagstr=simprof.fpar2tag('profile',outpar)
    else: raise Exception('Invalid method: {}'.format(method))
    return tagstr
def testkw(simstr,method,**method_kw):
    """
    Sets keywords necessary for testing plotting methods
    """
    method_kw['loadcalc']=False
    method_kw['allsnap']=False
    method_kw['askuser']=False
    method_kw['fext']=simstr.get_lastsnap(ftype='gadget',outfext=True)
    if   method=='contour':
        method_kw['mode']='pos_scfpotm2'
    elif method=='imagshw':
        method_kw['mode']='pos_scfpotm2'
        method_kw['palette']='pnbody_isophot'
    elif method=='snapdif':
        method_kw['mode']='EkLz_Lz'
    else:
        method_kw['mode']='pos_m'
    method_kw['type']='visi'
    method_kw['galaxy']=1
    method_kw['view']='xy'
#    method_kw['size']=30.
    method_kw['shape']=100
    method_kw['center']=True
    method_kw['centermeth']=simlist.LIST_COMMETHS
    method_kw['centertyp']='disk'
    method_kw['centergal']=1
    method_kw['outplot']=False
    method_kw['showplot']=True
    method_kw['overwrite']=True
    method_kw['owdata']=True
    method_kw['test']=True
    return method_kw

####################################################################################################################################
####################################################################################################################################
# HIGH LEVEL METHODS REQUIRING MMLSIM

####################################################################################################################################
# METHOD FOR RUNNING DIFFERENT PLOT METHODS
def run(simstr,method,verbose=None,askuser=None,**method_kw):
    """
    Provides interface for running different PLOT methods
    """
    # Set constants
    methLIST=LIST_METHODS
    snapmethLIST=LIST_METHODS_SNP
    snapplotLIST=[]
    for imeth in snapmethLIST: snapplotLIST+=list_snapplots(imeth)
    # Pars input
    method=mmlpars.mml_pars(method,list=methLIST+snapplotLIST)
    verbose=mmlpars.mml_pars(verbose,type=bool,default=True)
    askuser=mmlpars.mml_pars(askuser,type=bool,default=True)
    method_kw['verbose']=verbose
    method_kw['askuser']=askuser
    # Initialize default output
    out=None
    # Proceed based on method
    if   method=='retrieve':
        simstr.mvfiles_cfhmac('TO','plots','all')
    elif method in LIST_METHODS_RUN:
        plotrun(simstr,method=method,**method_kw)
    elif method in LIST_METHODS_SNP+snapplotLIST:
        plotsnap(simstr,method=method,**method_kw)
    elif method.startswith('test'):
        method_kw=testkw(simstr,method.split('test')[-1],**method_kw)
        plotsnap(simstr,method.split('test')[-1],**method_kw)
##     elif method in methLIST:
##         shortlist=mmlio.yorn('Restrict analysis to a short list of snapshots?')
##         method_kw['shortlist']=shortlist
##         method_kw['overwrite']=None
##         simlyze.analyze(simstr,method=method,**method_kw)
    else:
        raise Exception('Option for method {} needs to be added to simplot.run.'.format(method))
    # Return output
    return out

####################################################################################################################################
####################################################################################################################################
# FILE CLASSES AND METHODS

####################################################################################################################################
# METHOD FOR PARSING SIMPLOT FILE OPTIONS
def pars_fileopt(fopt):
    """
    Returns parsed PLOT file option dictionary with keys:
        plotext: Str extension to use for plot files
        animext: Str extension to use for animation files
    """
    # Define keys
    plotopt=mmlplot.parsopt()
    form={
        'plotext': mmlpars.parsdict(default=plotopt['plotext'],type=str),
        'animext': mmlpars.parsdict(default=plotopt['animext'],type=str)
        }
    # Pars input
    fopt=mmlpars.mml_formpars(fopt,form)
    # Return parsed input
    return fopt

####################################################################################################################################
# METHOD TO RETURN DICTIONARY OF SIMPLOT FILES FOR SNAPSHOT PLOTS
def snapfiles(simstr,method,inpar,**exkw):
    """
    Returns a dictionary containing a snapshot plot directory and file base
    """
    if method=='profile':
        snfdict=simfile.fpar2file(simstr,'prof',method,inpar,plot=True,**exkw)
    else:
        snfdict=simfile.fpar2file(simstr,'plots',method,inpar,plot=True,**exkw)
    return snfdict
    
def olsnapfiles(simstr,method,type=None,mode=None,galaxy=None,view=None,fext=None,test=None):
    # Pars input
    method=mmlpars.mml_pars(method,list=LIST_METHODS_SNP)
    type=mmlpars.mml_pars(type,list=list_snaptypes(method))
    if method!='science':
        mode=mmlpars.mml_pars(mode,list=list_plotmodes(method))
        galaxy=mmlpars.mml_pars(galaxy,list=list_plotpgals())
        if method != 'profile':
            view=mmlpars.mml_pars(view,list=list_plotviews())
    fext=mmlpars.mml_pars(fext,type=str,default='_*')
    test=mmlpars.mml_pars(test,type=bool,default=False)
    if test: teststr='test_'
    else   : teststr=''
    # Get file info
    fdict=simstr.fdict
    plotdir=fdict['plots']['dir']
    animdir=fdict['plots']['animdir']
    plotext=fdict['plots']['plotext']
    animext=fdict['plots']['animext']
    pfix=fdict['pfix']
    # Initialize dictionary
    snfdict={'plotdir':plotdir,'plotext':plotext,'animdir':animdir,'animext':animext}
    # Directory
    if   method=='science': dir='{}_{}'.format(method,type)
    else                  : dir='{}_{}_{}_{}_{}'.format(method,type,mode,galaxy,view)
    snfdict['dir']=os.path.join(plotdir,method,dir)
    # File base
    if   method=='science': base='{}{}{}_{}plot'.format(teststr,pfix,method,type)
    else                  : base='{}{}{}_{}plot_{}_{}_{}'.format(teststr,pfix,method,type,mode,galaxy,view)
    snfdict['base']=os.path.join(snfdict['dir'],base)
    # File name
    snfdict['file']='{}{}.{}'.format(snfdict['base'],fext,plotext)
    # Animation
    if   method=='science': anim='{}{}{}_{}anim.{}'.format(teststr,pfix,method,type,animext)
    else                  : anim='{}{}{}_{}anim_{}_{}_{}.{}'.format(teststr,pfix,method,type,mode,galaxy,view,animext)
    snfdict['anim']=os.path.join(animdir,anim)
    # Return dictionary
    return snfdict

####################################################################################################################################
# METHOD TO RETURN DICTIONARY OF SIMPLOT FILES AND DICTIONARY
def files(fdict):
    '''
    Returns a dictionary of PLOT files & directories with keys:
      dir:        Directory containing all simulation calc files
      energy:     File containing plot of energy stats vs. time.
    '''
    # Set constants
    statlist=simstat.LIST_STATS
    # Pars input
    fdict=mmlpars.mml_pars(fdict,type=dict)
    options=pars_fileopt({})
    rundir=fdict['rundir']
    pfix=fdict['pfix']
    pfix_shrt=fdict['pfix_shrt']
    plotext=options['plotext']
    animext=options['animext']
    # Initialize output dictionary
    files=dict(fdict['plots'],
               plotext=plotext,animext=animext,
               plotstat=os.path.join(rundir,pfix_shrt+'static_plot')+'.'+plotext)
    # Singular plot key/value pairs
    fkeys=['energy']+['stat_'+istat for istat in statlist]
    flist=[ifkey+'.'+plotext for ifkey in fkeys]
    # Snapshot plot directory key/value pairs
    dkeys=['animdir','compdir']
    dlist=['anim','compare']
    # Add input files
    files['subdirs']=dkeys
    files=mmlfiles.filedict(files['dir'],dirkeys=dkeys,dirlist=dlist,filekeys=fkeys,filelist=flist,pfix=pfix,fdict=files)
    # Return file dictionary
    fdict['plots']=files
    return fdict

####################################################################################################################################
# METHOD TO RETURN LIST OF RELAVENT FILES
def get_filelist(fdict,keylist):
    """
    Returns a list of relevant PLOT files.
    """
    # Set constants
    statlist=simstat.LIST_STATS
    # Pars input
    fdict=mmlpars.mml_pars(fdict,type=dict)
    keylist=mmlpars.mml_pars(keylist,type=list)
    # Add files
    filelist={}
    if 'static'    in keylist:
        filelist['static'  ]=[fdict['plotstat']]
    if 'energy'    in keylist:
        filelist['energy'  ]=[fdict['energy']]
    if 'stats'     in keylist:
        filelist['stats'   ]=[fdict['stat_'+istat] for istat in statlist]
    if 'anim' in keylist: filelist['anim']=[]
    for ilong in dict_ftyp()['longout']:
        ishrt='short'+ilong
        if 'anim' in keylist or ilong in keylist or ishrt in keylist:
            if ilong in keylist: filelist[ilong]=[]
            if ishrt in keylist: filelist[ishrt]=[]
            fpar=list_fpar(ilong)
            if len(fpar.keys())==0: fpar={'none':[None]}
            iterpar=simfile.iter_fpar('plots',ilong,fpar=fpar)
            for ifpar in iterpar:
                isnfdict=snapfiles(fdict,ilong,ifpar)
                if 'anim' in keylist: filelist['anim'].append(isnfdict['anim'])
                if ilong  in keylist: filelist[ilong ].append(isnfdict['file'])
                if ishrt  in keylist: 
                    ilongfiles=glob.glob(isnfdict['file'])
                    if len(ilongfiles)>0: filelist[ishrt].append(sorted(ilongfiles)[-1])
    # Return output
    return filelist

####################################################################################################################################
####################################################################################################################################
# HIGH LEVEL PLOTTING METHODS
def plotrun(simstr,method=None,**method_kw):
    """
    Creates plots for entire run
    """
    # Set constants
    methLIST=LIST_METHODS_RUN
    # Pars input
    method=mmlpars.mml_pars(method,default='all',list=['all']+methLIST)
    # Proceed based on method
    if   method=='all'   :
        for imeth in methLIST: plotrun(simstr,method=imeth,**method_kw)
    elif method=='anim'  : plotrun_anim(simstr,**method_kw)
    elif method=='energy': plotrun_energy(simstr,**method_kw)
    elif method=='stats' : plotrun_stats(simstr,**method_kw)
    else: raise Exception('[simplot.plotrun] Invalid method: {}'.format(method))
    # Return control
    return

####################################################################################################################################
# METHOD TO GENERATE ANIMATIONS
def plotrun_anim(simstr,loadkeys=None,flagdict=None,setflags=None,
                 verbose=None,verbwid=None,**anim_kw):
    """
    Creates animations for run
    """
    # Set constants
    methlist=LIST_METHODS_SNP+['hist_'+ihst for ihst in simhist.LIST_METHODS_HST]
    # Pars input
    flagdict=mmlpars.mml_pars(flagdict,default={},type=dict)
    setflags=mmlpars.mml_pars(setflags,default=False,type=bool)
    verbose=mmlpars.mml_pars(verbose,default=False,type=bool)
    if setflags: verbose=False
#    anim_kw['loglevel']='debug'
    # Get file dictionary
    filedict=simstr.fdict
    # Create directory
    mmlfiles.mkdirs(filedict['plots']['animdir'])
    # Set flags and print info
    if 'anim' not in flagdict: flagdict['anim']=False
    if verbose: mmlio.verbose('[{}] Beginning animations'.format(simstr['runtag']),border=True,width=verbwid)
    # Create list of frames
    baselist=[] ; framDict={} ; animDict={}
    ometh=mmlio.askselect('What type of animation should be looked for?',['all']+methlist)
    if ometh!='all': methlist=[ometh]
    for imeth in methlist:
        if 'hist' in imeth:
            if ometh=='all':
                if not mmlio.yorn('Look for {} animations?'.format(imeth)): continue
            hmeth=imeth.split('hist_')[-1]
            movestr=simhist.animhist(simstr,hmeth,verbose=verbose,**anim_kw)
        else:
            parlist=simfile.iter_fpar('plots',imeth)
            for ipar in parlist:
                ibase='{}_'.format(imeth)
                for ikey in ipar: ibase+=str(ipar[ikey]) 
                ifdict=snapfiles(simstr,imeth,ipar)
                if len(glob.glob(ifdict['file']))==0: continue
                framDict[ibase]=ifdict['file']
                animDict[ibase]=ifdict['anim']
                baselist.append(ibase)
                # Ask if animation should be overwritten
                if ibase not in flagdict:
                    if flagdict['anim']: flagdict[ibase]=True
                    else: flagdict[ibase]=mmlio.yorn('Overwrite {} animation?'.format(ibase),width=verbwid)
                # Create animation
                if not setflags:
                    movestr=mmlplot.ffmpeg(framDict[ibase],animDict[ibase],verbose=verbose,
                                           overwrite=flagdict[ibase],rmtemp=True,**anim_kw)
    # Return info
    if setflags:
        return flagdict
    else:
        return flagdict

####################################################################################################################################
# METHOD TO PLOT ENERGY FOR RUN
def plotrun_energy(simstr,loadkeys=None,flagdict=None,setflags=None,verbose=None,verbwid=None,**plot_kw):
    """
    Plots energy as a function of time
    """
    # Pars input
    flagdict=mmlpars.mml_pars(flagdict,default={},type=dict)
    setflags=mmlpars.mml_pars(setflags,default=False,type=bool)
    verbose=mmlpars.mml_pars(verbose,default=False,type=bool)
    if setflags: verbose=False
    # Get file dictionary
    filedict=simstr.fdict
    # Set flag and print info
    if 'energy' not in flagdict: flagdict['energy']=mmlio.yorn('Overwrite existing energy plot?',width=verbwid)
    if verbose: mmlio.verbose('[{}] Beginning energy plot'.format(simstr['runtag']),border=True,width=verbwid)
    # Continue if not setting flags
    if not setflags:
        # Get file names
        energyfile=filedict['gadget']['output']['energy']
        energyplot=filedict['plots']['energy']
        # Only continue if overwrite set or plot dosn't exist
        if not os.path.isfile(energyplot) or flagdict['energy']:
            # Make directory
            mmlfiles.mkdirs(os.path.dirname(energyplot))
            # Plot
            mmlgadget.plot_energy(energyfile,plotfile=energyplot,overwrite=flagdict['energy'])
            if verbose: mmlio.verbose('Energy plot:')
            if verbose: print '    {}'.format(energyplot)
    # Return info
    if setflags:
        return flagdict
    else:
        return flagdict

####################################################################################################################################
# PLOT STATISTICS
def plotrun_stats(simstr,method=None,fpar=None,fparfile=None,loadkeys=None,statdict=None,flagdict=None,
                  setflags=None,verbose=None,verbwid=None,**plot_kw):
    """
    Plots statistics
    """
    # Pars input
    if method not in simstat.LIST_STATS:
        method=mmlio.askselect('Select a stat type to plot:',simstat.LIST_STATS)
    method=mmlpars.mml_pars(method,list=simstat.LIST_STATS)
    statdict=mmlpars.mml_pars(statdict,default={},type=dict)
    flagdict=mmlpars.mml_pars(flagdict,default={},type=dict)
    setflags=mmlpars.mml_pars(setflags,default=False,type=bool)
    verbose=mmlpars.mml_pars(verbose,default=False,type=bool)
    if setflags: verbose=False
    # Get file dictionary
    filedict=simstr.fdict
    # Print info
    if verbose: mmlio.verbose('[{}] Beginning stat plots'.format(simstr['runtag']),border=True,width=verbwid)
    # Get file parameters
    if fpar==None: fpar=simfile.parsfpar('stat',method,askuser=True,tagstr=fparfile,plot=True)
    # Select files
    statfile=simstat.snapfiles(simstr,method,fpar)['file']
    statplot=simstat.snapfiles(simstr,method,fpar,plot=True)['file']
    # Continue if file dosn't exist
    if not os.path.isfile(statfile): raise Exception('{} statistics file does not exist.'.format(method))
    # Set flag if not provided
    if method not in flagdict: flagdict[method]=mmlio.yorn('Overwrite existing {} statistics plot?'.format(method),width=verbwid)
    # Continue if not setting flags
    if not setflags:
        # Only continue if overwrite set or plot does not exist
        if not os.path.isfile(statplot) or flagdict[method]:
            # Load missing statistics dictionaries
            if method=='var': submeth=fpar['smeth']
            else            : submeth=None
            if not method in statdict: statdict[method]=simstat.rwstats('R',stattype=method,substattype=submeth,statfile=statfile)
            # Plot
            if   method=='com' : out=simstat.plot_comstat(statdict[method],fname=statplot,overwrite=flagdict[method],**plot_kw)
            elif method=='var' : out=simstat.plot_varstat(submeth,statdict[method],fname=statplot,overwrite=flagdict[method],rowvar=fpar['rowvar'],**plot_kw)
            else               : out=plot_stats(statdict[method],method=method,fname=statplot,overwrite=flagdict[method],**plot_kw)
            if verbose: mmlio.verbose('Statistics plot [{}]:'.format(method))
            if verbose: print '     {}'.format(statplot)
    # Return info
    if setflags:
        return flagdict
    else:
        return flagdict,statdict
        
####################################################################################################################################
####################################################################################################################################
# PLOTTING CLASSES AND METHODS FOR SNAPSHOTS
def plotsnap(simstr,method=None,snapdict=None,calcdict=None,calcdict0=None,profdict=None,loadcalc=None,
             flagdict=None,setflags=None,allsnap=None,overwrite=None,askuser=None,**input_kw):
    """
    Creates plots for single snapshots
    """
    # Set constants
    snapmethLIST=LIST_METHODS_SNP
    methlistDEF=copy.deepcopy(LIST_METHODS_SNP)
    for isnapmeth in snapmethLIST: methlistDEF+=list_snapplots(isnapmeth)
    # Pars input
    method=mmlpars.mml_pars(method,default='all',list=['all']+methlistDEF)
    askuser=mmlpars.mml_pars(askuser,default=False,type=bool)
    if not isinstance(allsnap,bool):
        if askuser:
            allsnap=mmlio.yorn('Plot all snapshots for this run?')
        else:
            allsnap=False
    flagdict=mmlpars.mml_pars(flagdict,default={},type=dict)
    setflags=mmlpars.mml_pars(setflags,default=False,type=bool)
    plot_kw=mmlplot.parsopt(input_kw)
    # Call calc to handle looping if allsnap selected
    if allsnap:
        input_kw['shortlist']=False
        simlyze.analyze(simstr,method='plotsnap_'+method,**input_kw)
        out=None
    # Proceed based on method
    else:
        # Loadkeys
        loadkeys,method_kw=simstr.get_loadkeys(snapdict=snapdict,noidgal2=True,askuser=askuser,outextra=True,**input_kw)
        # Set flags
        flagdict=plotsnap_setflags(method,flagdict=flagdict,askuser=askuser,overwrite=overwrite)
        # Do things that require loading data
        if not setflags:
            # Loadcalc flag
            if not isinstance(loadcalc,bool):
                if askuser:
                    loadcalc=mmlio.yorn('Load histogram info from file?')
                else:
                    loadcalc=False
            # Assemble keywords
            method_kw['loadcalc']=loadcalc
            method_kw['loadkeys']=loadkeys
            method_kw['askuser' ]=askuser
            # Plot
            outtot=[]
            for ikey in snapmethLIST:
                if not flagdict[ikey]: continue
                method_kw['overwrite']=flagdict['ow'+ikey]
                if method in list_snapplots(ikey): method_kw['type']=method.split(ikey)[-1]
                if   ikey == 'science': out=plotsnap_science(simstr,snapdict=snapdict,calcdict=calcdict,calcdict0=calcdict0,**method_kw)
                elif ikey == 'profile': out=plotsnap_profile(simstr,fext=loadkeys['fext'],**method_kw)
                else                  : out=plotsnap_wrapper(ikey,simstr,snapdict=snapdict,calcdict=calcdict,calcdict0=calcdict0,**method_kw)
                if plot_kw['outplot']: outtot.append(out)
                else:
                    if   ikey == 'profile': pass
                    else                  : snapdict,calcdict,calcdict0=out
    # Return output
    if setflags: return flagdict
    if not plot_kw['outplot']: outtot=flagdict,snapdict,calcdict,calcdict0
    return outtot

####################################################################################################################################
# METHOD FOR SETTING PLOTSNAP FLAGS
def plotsnap_setflags(method,flagdict=None,askuser=None,overwrite=None):
    """
    Sets flags for a type of plot
    """
    # Set constants
    mtypes=LIST_METHODS_SNP
    # Pars input
    flagdict=mmlpars.mml_pars(flagdict,type=dict,default={})
    askuser=mmlpars.mml_pars(askuser,default=False,type=bool)
    # Loop over method types
    for imtype in mtypes:
        # Overall flag
        if imtype not in flagdict: flagdict[imtype]=(method in ['all',imtype]+list_snaptypes(imtype))
        # Overwrite flag
        if flagdict[imtype] and 'ow'+imtype not in flagdict:
            if isinstance(overwrite,bool):
                flagdict['ow'+imtype]=overwrite
            else:
                if askuser:
                    flagdict['ow'+imtype]=mmlio.yorn('Overwrite existing {} plots for each snapshot?'.format(imtype))
                else:
                    flagdict['ow'+imtype]=False
    # Return flag dictionary
    return flagdict
            
####################################################################################################################################
# METHOD FOR PLOTTING DISPLAY PLOTS
def plotsnap_wrapper(method,simstr,snapdict=None,calcdict=None,calcdict0=None,loadcalc=None,loadkeys=None,
                     plotfile=None,overwrite=None,askuser=None,verbose=None,verbwid=None,test=None,
                     **input_kw):
    """
    Provides command line interface for creating display plots
    """
    # Set constants
    snapmeths=LIST_METHODS_SNP
    fpar=simfile.list_fpar('plots',method)
    # Pars input
    askuser=mmlpars.mml_pars(askuser,default=False,type=bool)
    verbose=mmlpars.mml_pars(verbose,default=False,type=bool)
    loadcalc=mmlpars.mml_pars(loadcalc,default=False,type=bool)
    loadkeys=mmlpars.mml_pars(loadkeys,default={},type=dict)
    loadkeys=simstr.get_loadkeys(snapdict=snapdict,noidgal2=True,outextra=False,askuser=askuser,**loadkeys)
    test=mmlpars.mml_pars(test,default=False,type=bool)
    if test:
        overwrite=True
        input_kw['owdata']=True
    plot_kw,input_kw=mmlplot.parsopt(input_kw,outextra=True)
    # Ask for user input
    if askuser and method not in snapmeths:
        method=mmlio.askselect('What method should be used to create plot?',snapmeths)
    else:
        method=mmlpars.mml_pars(method,default='display',list=snapmeths)
#    for ipar in fpar['keylist']:
#        if ipar not in input_kw: input_kw[ipar]=None
#        if askuser and input_kw[ipar] not in fpar[ipar]:
#            input_kw[ipar]=mmlio.askselect('What plot {} should be created?'.format(ipar),fpar[ipar])
#        else:
#            input_kw[ipar]=mmlpars.mml_pars(input_kw[ipar],default=fpar[ipar][0],list=fpar[ipar])
    inpar=simfile.parsfpar('plots',method,fpar=input_kw,askuser=askuser)
    # Create lists of things to calc/load
    if   method=='science':
        singflag=False
        if   inpar['type']=='mult': pass
        elif inpar['type']=='hist1d': input_kw['h2DList']=[]
        elif inpar['type']=='hist2d': input_kw['h1DList']=[]
    else:
        singflag=True
        input_kw['h1DList']=[]
        input_kw['h2DList']=[inpar['mode']]
        input_kw['typList']=[inpar['type']]
        input_kw['galList']=[inpar['galaxy']]
        input_kw['rayList']=[inpar['view']]
    if   method=='contour': 
        if 'pos_m' not in input_kw['h2DList']: input_kw['h2DList'].append('pos_m')
    # Set histogram package
    if   method=='display': histpackage='pNbody'
    else                  : histpackage='numpy'
    # Get required info
    snapdict,calcdict,calckw=simcalc.askcalc(simstr,askuser=askuser,snapdict=snapdict,calcdict=calcdict,loadcalc=loadcalc,
                                             singflag=singflag,**input_kw)
    plotkw=dict(plot_kw,**inpar)
    # Files
    fdict=simstr.fdict
    if not plot_kw['outplot']:
        snfdict=snapfiles(simstr,method,inpar,fext=loadkeys['fext'],test=test)
        plotfile=mmlpars.mml_pars(plotfile,default=snfdict['file'],type=str)
        if not mmlfiles.prep_write(plotfile,overwrite=overwrite,askuser=askuser): return snapdict,calcdict,calcdict0
    # Load required objects
    snapdict,calcdict,calcdict0,calckw=simcalc.getcalc(simstr,askuser=askuser,
                                                       snapdict=snapdict,calcdict=calcdict,calcdict0=calcdict0,
                                                       loadcalc=loadcalc,loadkeys=loadkeys,
                                                       histpackage=histpackage,**calckw)
    # Plot
    if   method=='science':
        pass
    elif method=='display': out=plotcalc_impil(calcdict,plotfile=plotfile,overwrite=overwrite,verbose=verbose,**plotkw)
    elif method=='ellipse':
        ellipfileDEF=simcalc.snapfiles(simstr,'ellipfit',inpar,fext=loadkeys['fext'],test=test)['file'] #+'.txt'
        if 'datafile' not in plotkw: plotkw['datafile']=ellipfileDEF
        out=plotcalc_ellip(calcdict,plotfile=plotfile,overwrite=overwrite,verbose=verbose,**plotkw)
    elif method=='contour':
        implotkw=dict(plot_kw,**inpar) ; implotkw['mode']='pos_m' ; implotkw['outplot']=True ; implotkw['showplot']=False
        figobj=plotcalc_imshw(calcdict,verbose=verbose,**implotkw)
        out=plotcalc_contr(calcdict,plotfile=plotfile,overwrite=overwrite,verbose=verbose,figobj=figobj,**plotkw)
    elif method=='imagpil': out=plotcalc_impil(calcdict,plotfile=plotfile,overwrite=overwrite,verbose=verbose,**plotkw)
    elif method=='imagshw': out=plotcalc_imshw(calcdict,plotfile=plotfile,overwrite=overwrite,verbose=verbose,**plotkw)
    elif method=='snapdif':
        snaplist=[fdict['gadget']['input']['ic']]+sorted(glob.glob(fdict['gadget']['output']['snapbase']+'*'))
        snapfile1,icflag1=simstr.get_snapname(**loadkeys)
#        mmlio.verbose(snapfile1)
#        mmlio.yorn(fdict['gadget']['output']['snapbase']+'*')
        idxsnap1=snaplist.index(snapfile1)
        if idxsnap1==0: calcdict2=calcdict
        else:
            idxsnap2=0 # idxsnap1-1
            snapfile2=snaplist[idxsnap2]
            loadkeys2=copy.deepcopy(loadkeys) ; loadkeys2['fname']=snapfile2 ; del loadkeys2['fext']
            loadkeys2['fext']=simstr.get_fext(**loadkeys)
            snapdict2,calcdict2,calcdict20,calckw=simcalc.getcalc(simstr,askuser=askuser,
                                                                  loadcalc=loadcalc,loadkeys=loadkeys2,
                                                                  histpackage=histpackage,**calckw)
        out=plotcalc_sndif(calcdict,calcdict2,plotfile=plotfile,overwrite=overwrite,verbose=verbose,
                           residflag=True,**plotkw)
    else: raise Exception('Invalid plotting method: {}'.format(method))
    # Print info
    if not plot_kw['outplot'] and verbose:
        pltstr=''
        for ipar in inpar.keys(): pltstr+='{} '.format(inpar[ipar])
        mmlio.verbose('Created {}{} plot:'.format(pltstr,method))
        print '    '+plotfile
    # Return output
    if plot_kw['outplot']: return out
    else                 : return snapdict,calcdict,calcdict0

####################################################################################################################################
# METHOD FOR PLOTTING PROFILE PLOTS
def plotsnap_profile(simstr,**prof_kw):
    """
    Provides command line interface for creating profile plots
    """
    # Plot
    out=simprof.plot_profile(simstr,**prof_kw)
    # Return output
    return out

####################################################################################################################################
# METHOD FOR PLOTTING SCIENCE PLOTS
def plotsnap_science(simstr,type=None,snapdict=None,calcdict=None,calcdict0=None,loadcalc=None,loadkeys=None,
                     plotfile=None,overwrite=None,askuser=None,verbose=None,verbwid=None,**input_kw):
    """
    Provides command line interface for creating science plots
    """
    # Pars input
    type=mmlpars.mml_pars(type,default='mult',list=list_snaptypes('science'))
    askuser=mmlpars.mml_pars(askuser,default=False,type=bool)
    verbose=mmlpars.mml_pars(verbose,default=False,type=bool)
    loadcalc=mmlpars.mml_pars(loadcalc,default=False,type=bool)
    loadkeys=mmlpars.mml_pars(loadkeys,default={},type=dict)
    loadkeys=simstr.get_loadkeys(snapdict=snapdict,noidgal2=True,outextra=False,askuser=askuser,**loadkeys)
    plot_kw=mmlplot.parsopt(input_kw)
    # Files
    fdict=simstr.fdict
    if not plot_kw['outplot']:
        plotbase=fdict['plots']['science_'+type+'base']
        plotfileDEF=plotbase+loadkeys['fext']+'.'+fdict['plots']['plotext']
        plotfile=mmlpars.mml_pars(plotfile,default=plotfileDEF,type=str)
        if os.path.isfile(plotfile) and not isinstance(overwrite,bool) and askuser:
            overwrite=mmlio.yorn('That file already exists. Overwrite?')
        overwrite=mmlpars.mml_pars(overwrite,default=False,type=bool)
        if os.path.isfile(plotfile) and not overwrite: return snapdict,calcdict,calcdict0
    # Load required objects
    snapdict,calcdict,calcdict0,refkw,extrakw=simcalc.multcalc(simstr,method=type,askuser=askuser,
                                                               snapdict=snapdict,calcdict=calcdict,calcdict0=calcdict0,
                                                               loadcalc=loadcalc,loadkeys=loadkeys,outextra=True,
                                                               histpackage='numpy',**input_kw)
    # Plot
    plotcalckw=dict(refkw,**extrakw)
    if type=='mult':
        out=plotcalc_multi(calcdict,calcdict0=calcdict0,plotfile=plotfile,overwrite=overwrite,**plotcalckw)
    else: raise Exception('Invalid science plot type: {}'.format(type))
    # Print info
    if not plot_kw['outplot'] and verbose:
        mmlio.verbose('Created {} science plot: '.format(type))
        print '    '+plotfile
    # Return output
    if plot_kw['outplot']: return out
    else                 : return snapdict,calcdict,calcdict0

####################################################################################################################################
# METHOD FOR PLOTTING ELLIPSES OVER IMAGE
def plot_ellipse(nbObj,plotfile=None,overwrite=None,datafile=None,outimage=None,showimage=None,
                 semimajorlist=None,nsemimajor=None,params0=None,
                 interpmeth=None,nphi=None,errtol=None,miniter=None,maxiter=None,
                 **plot_kw):
    """
    Plots simple images of individual particle types
    """
    # Pars input
    saveflag=(plotfile != None)
    plotfile=mmlpars.mml_pars(plotfile,default='ellipse.png',type=str)
    datafile=mmlpars.mml_pars(datafile,default='ellipse.txt',type=str)
    overwrite=mmlpars.mml_pars(overwrite,default=False,type=bool)
    outimage=mmlpars.mml_pars(outimage,default=False,type=bool)
    showimage=mmlpars.mml_pars(showimage,default=False,type=bool)
    # Get image data
    plot_kw['outimage']=True
    imgpil,view_window=plot_sing(nbObj,**plot_kw)
    imgdat=np.array(imgpil.convert('L')).astype(float)
    imgrgb=np.array(imgdat)
    imgdim=[[-view_window[0]/2,view_window[0]/2],[-view_window[1]/2,view_window[1]/2]]
    # Get ellipse data
    if os.path.isfile(datafile) and not overwrite:
        elliplist=mmlellipse.readellipse(datafile)
    else:
        elliplist=mmlellipse.fitmultellip(semimajorlist=semimajorlist,nsemimajor=nsemimajor,params0=params0,
                                          imgdat=imgdat,imgdim=imgdim,interpmeth=interpmeth,
                                          nphi=nphi,errtol=errtol,miniter=miniter,maxiter=maxiter)
        mmlellipse.saveellipse(datafile,elliplist,overwrite=overwrite)
    # Plot
    out=mmlellipse.plotellipse(elliplist,imgdat=imgrgb,imgdim=imgdim,nphi=nphi)
    # Display image
    if showimage: out.show()
    # Return
    if outimage:
        return out
    else:
        if overwrite or not os.path.isfile(plotfile):
            out.savefig(plotfile)
    # Return
    return None

####################################################################################################################################
# METHOD TO PLOT SIMULATION
def plot_mult(plotfile=None,calcfile=None,calcdict=None,calcfile0=None,calcdict0=None,
              h1DList=None,h2DList_imag=None,h2DList_cont=None,typList=None,galList=None,rayList=None,
              overwrite=None,plotobj=None,**plot_kw):
    '''
    Plots info on a snapshot
    '''
    # Set constants
    fnameTAG='calcplot'
    h1DListDEF=['rxyz_rho','z_n','vz_n','phi_n','rxyz_rho_resid']
    h2DList_imagDEF=['pos_m','pos_vr']#,'pos_vr','pos_gal']
    h2DList_contDEF=['pos_scfpotm2',None]
    typListDEF=['gas','halo','disk','bulge','stars','bndry']
    galListDEF=0
    rayListDEF=['xy','yz','xz']
    # Pars file input
    if not isinstance(calcdict ,dict) and not isinstance(calcfile ,str): raise Exceptions('File name not provided and dictionary is empty.')
    if calcdict==None: calcdict=loadcalc(calcfile)
    if not isinstance(calcdict0,dict) and     isinstance(calcfile0,str): calcdict0=loadcalc(calcfile0)
    overwrite=mmlpars.mml_pars(overwrite,default=False,type=bool)
    if not mmlfiles.prep_write(plotfile,overwrite): return
    # Pars plotting input
    options=mmlplot.parsopt(plot_kw)
    options['pad']=5. ; options['wpad']=2.2 ; options['tickfontsize']=12. ; options['textfontsize']=14.
    options['xlabpos']=(0.5,-0.12) ; options['ylabpos']=(-0.12,0.5)
    h1DList=mmlpars.mml_pars(h1DList,default=h1DListDEF,type=list) ; nh1D=len(h1DList)
    h2DList_imag=mmlpars.mml_pars(h2DList_imag,default=h2DList_imagDEF,type=list,nelements=2) ; nh2D_imag=len(h2DList_imag)
    h2DList_cont=mmlpars.mml_pars(h2DList_cont,default=h2DList_contDEF,type=list,nelements=2) ; nh2D_cont=len(h2DList_cont)
    typList=mmlpars.mml_pars(typList,default=typListDEF,type=list) ; ntyp=len(typList)
#    galList=mmlpars.mml_pars(galList,default=galListDEF,type=list) ; ngal=len(galList)
    rayList=mmlpars.mml_pars(rayList,default=rayListDEF,type=list) ; nray=len(rayList)
    # Count rows and columns
    ntyp_plot=0
    for ityp in typList:
        if ityp in calcdict['typList']: ntyp_plot+=1
    nray_plot=0
    for iray in rayList:
        if iray in calcdict['rayList']: nray_plot+=1
    nh1D_plot=0
    for ih1D in h1DList:
        if ih1D in calcdict['h1DList']: nh1D_plot+=1
        else:
            if ih1D.split('_')[-1] == 'resid': nh1D_plot+=1
    nrow=max([ntyp_plot,nh1D_plot])
    ncol=nray_plot
    # Determine size of figure
    fontsize_inc=max(options['tickfontsize'],options['textfontsize'])/72.
    opad_inc=options['pad']*fontsize_inc ; ipad_inc=options['wpad']*fontsize_inc
    options['figsize']=(4*opad_inc+2*(ncol+1)*options['axsize'][0]+2*(ncol-1)*ipad_inc,
                        2*opad_inc+nrow*options['axsize'][1]+(nrow-1)*ipad_inc)
    # Set dimensions relative to figure
    opad_relx=opad_inc/options['figsize'][0] ; opad_rely=opad_inc/options['figsize'][1]
    ipad_relx=ipad_inc/options['figsize'][0] ; ipad_rely=ipad_inc/options['figsize'][1]
    rect_rho=np.array([0.,1.-opad_rely,options['axsize'][0]/options['figsize'][0],options['axsize'][1]/options['figsize'][1]])
    rect_rad=np.array([0.,1.-opad_rely,2.*options['axsize'][0]/options['figsize'][0],options['axsize'][1]/options['figsize'][1]])
    rect_vel=np.array([0.,1.-opad_rely,options['axsize'][0]/options['figsize'][0],options['axsize'][1]/options['figsize'][1]])
    rect_rho[0]=opad_relx
    rect_rad[0]=rect_rho[0]+ncol*rect_rho[2]+(ncol-1)*ipad_relx+opad_relx
    rect_vel[0]=rect_rad[0]+rect_rad[2]+opad_relx
    rect_rho_cb=np.array([rect_rho[0],rect_rho[1]-ntyp_plot*rect_rho[3]-ntyp_plot*ipad_rely,ncol*rect_rho[2]+(ncol-1)*ipad_relx,ipad_rely])
    rect_vel_cb=np.array([rect_vel[0],rect_vel[1]-ntyp_plot*rect_vel[3]-ntyp_plot*ipad_rely,ncol*rect_vel[2]+(ncol-1)*ipad_relx,ipad_rely])
    rect_txt=np.array([rect_vel_cb[0],rect_vel_cb[1]-ntyp_plot*rect_vel_cb[3]-(ntyp_plot-1)*ipad_rely-ipad_rely,rect_vel[2],rect_vel[3]])
    # Create plot objects
    if plotobj==None:
        newfig=True
        plotobj=dict(fig=plt.figure(),rho=None,rad=None,vel=None,rho_cb=None,vel_cb=None,txt=None)
    else:
        newfig=False
    # Plots of mass vs. position
    imag_mode1=h2DList_imag[0] ; cont_mode1=h2DList_cont[0]
#    imag_mode1='pos_m'
#    cont_mode1='pos_scfpotm2'
    plotobj['rho']=panels_hist2D(plotobj['fig'],typList,rayList,0,calcdict['hist2D'],imag_mode=imag_mode1,cont_mode=cont_mode1,
                                 unitdict=calcdict['unitDict'],rect=rect_rho,wpad=ipad_relx,hpad=ipad_rely)
    plotobj['rho_cb']=panel_colorbar(plotobj['fig'],[plotobj['rho']['imagList'],plotobj['rho']['contList']],
                                     rect=rect_rho_cb,wpad=ipad_relx,hpad=ipad_rely,label=[imag_mode1.split('_')[1],cont_mode1.split('_')[1]])
    # 1D histograms
    if calcdict0!=None:
        plotobj['rad']=panels_hist1D(plotobj['fig'],typList,h1DList,1,calcdict['hist1D'],calcdict0['hist1D'],
                                     unitdict=calcdict['unitDict'],rect=rect_rad,wpad=ipad_relx,hpad=ipad_rely)
    else:
        plotobj['rad']=panels_hist1D(plotobj['fig'],typList,h1DList,1,calcdict['hist1D'],
                                     unitdict=calcdict['unitDict'],rect=rect_rad,wpad=ipad_relx,hpad=ipad_rely)
    # Plots of velocity vs. position
    imag_mode2=h2DList_imag[1] ; cont_mode2=h2DList_cont[1]
#    imag_mode2='pos_vr'
    plotobj['vel']=panels_hist2D(plotobj['fig'],typList,rayList,0,calcdict['hist2D'],imag_mode=imag_mode2,cont_mode=cont_mode2,
                                 unitdict=calcdict['unitDict'],rect=rect_vel,wpad=ipad_relx,hpad=ipad_rely)
    plotobj['vel_cb']=panel_colorbar(plotobj['fig'],[plotobj['vel']['imagList'],plotobj['vel']['contList']],
                                     rect=rect_vel_cb,wpad=ipad_relx,hpad=ipad_rely,label=[imag_mode2.split('_')[1]])
    # Text panel
    plotobj['txt']=panel_text(plotobj['fig'],calcdict['statDict'],calcdict['unitDict'],rect=rect_txt,wpad=ipad_relx,hpad=ipad_rely)
    # Save figure
    if newfig: mmlplot.set_figprop(plotobj['fig'],options)
    plotobj['fig'].savefig(plotfile,edgecolor=options['fgclr'],facecolor=options['bgclr'])
    # Return plot objects
    return plotobj

####################################################################################################################################
# METHOD FOR PLOTTING ELLIPSE OVER IMAGE
def plotcalc_ellip(calcdict,plotfile=None,overwrite=None,verbose=None,
                   datafile=None,owdata=None,semimajorlist=None,nsemimajor=None,params0=None,
                   interpmeth=None,nphi=None,errtol=None,dertol=None,miniter=None,maxiter=None,
                   palette=None,**imagekw):
    """
    Plot ellipses over images for a snapshot
    """
    # Pars input
    plotfile=mmlpars.mml_pars(plotfile,default='ellipse.png',type=str)
    overwrite=mmlpars.mml_pars(overwrite,default=False,type=bool)
    datafile=mmlpars.mml_pars(datafile,default='ellipse.txt',type=str)
    owdata=mmlpars.mml_pars(owdata,default=False,type=bool)
    palette=mmlpars.mml_pars(palette,default='pnbody_light',type=str)
    plot_kw=mmlplot.parsopt(imagekw)
    if not plot_kw['outplot']:
        if not mmlfiles.prep_write(plotfile,overwrite): return None
    # Get image data
    imagekw['outdata']=True
    imgdat,imgscl,view_window=plotcalc_impil(calcdict,verbose=verbose,**imagekw)
    imgdim=[[-view_window[0]/2.,view_window[0]/2.],[-view_window[1]/2.,view_window[1]/2.]]
    # Get ellipse data
    if os.path.isfile(datafile) and not owdata:
        elliplist=mmlellipse.readellipse(datafile)
    else:
        elliplist=mmlellipse.fitmultellip(semimajorlist=semimajorlist,nsemimajor=nsemimajor,params0=params0,
                                          imgdat=imgdat,imgdim=imgdim,interpmeth=interpmeth,
                                          nphi=nphi,errtol=errtol,dertol=dertol,miniter=miniter,maxiter=maxiter,verbose=verbose)
        mmlellipse.saveellipse(datafile,elliplist,overwrite=owdata,verbose=verbose)
    # Plot
    out=mmlellipse.plotellipse(elliplist,imgdat=imgdat,imgdim=imgdim,nphi=nphi,palette=palette,verbose=verbose)
    # Display image
    if plot_kw['showplot']: out.show()
    # Return
    if plot_kw['outplot']: return out
    else                 : out.savefig(plotfile)
    return None

####################################################################################################################################
# METHOD TO PLOT DIFFERENCE BETWEEN 2D HISTOGRAMS FOR SIMULATION
def plotcalc_sndif(calcdict,calcdict2,plotmeth=None,**plotkw):
    """
    Plot imag of difference of histograms for a snapshot
    """
    # Pars input
    plotmeth=mmlpars.mml_pars(plotmeth,default='imshow',list=['pil','imshow'])
    # Plot
    if   plotmeth=='pil':
        out=plotcalc_impil(calcdict,calcdict2=calcdict2,**plotkw)
    elif plotmeth=='imshow':
        out=plotcalc_imshw(calcdict,calcdict2=calcdict2,**plotkw)
    else: raise Exception('Invalid plotting method: {}'.format(plotmeth))
    # Return output
    return out

####################################################################################################################################
# METHOD TO PLOT 2D HISTOGRAMS FOR SIMULATION
def plotcalc_impil(calcdict,calcdict2=None,mode=None,type=None,view=None,galaxy=None,
                   plotfile=None,overwrite=None,verbose=None,annotate=None,residflag=None,
                   palette=None,outdata=None,**plot_kw):
    """
    Plots images for a snapshot
    """
    # Set constants
    fnameTAG='calcimage'
    imag_modeDEF='pos_m'
    cont_modeDEF='pos_scfpotm2'
    typListDEF=simlist.LIST_PTYPS
    rayListDEF=list_plotviews()
    # Pars general input
    calcdict=mmlpars.mml_pars(calcdict,type=dict)
    if calcdict2==None: residflag=False
    annotate=mmlpars.mml_pars(annotate,default=True,type=bool)
    residflag=mmlpars.mml_pars(residflag,default=False,type=bool)
    outdata=mmlpars.mml_pars(outdata,default=False,type=bool)
    verbose=mmlpars.mml_pars(verbose,default=False,type=bool)
    # Pars plotting info
    mode=mmlpars.mml_pars(mode,default='pos_m',list=calcdict['h2DList'])
    type=mmlpars.mml_pars(type,default='visi' ,list=calcdict['typList'])
    view=mmlpars.mml_pars(view,default='xy'   ,list=calcdict['rayList'])
    galaxy=mmlpars.mml_pars(galaxy,default=0  ,list=range(3))
    palette=mmlpars.mml_pars(palette,default='pnbody_light',type=str)
    plot_kw=mmlplot.parsopt(plot_kw)
    # Pars file input
    plotfile=mmlpars.mml_pars(plotfile,default='pNbodyplot_{}.png'.format(fnameTAG),type=str)
    overwrite=mmlpars.mml_pars(overwrite,default=False,type=bool)
    if not outdata and not plot_kw['outplot']:
        if not mmlfiles.prep_write(plotfile,overwrite): return None
    # Isolate calcdict options
    if residflag:
        histopt=resid_hist2D(calcdict['hist2D'][mode][type][galaxy][view],calcdict2['hist2D'][mode][type][galaxy][view])
    else:
        histopt=copy.deepcopy(calcdict['hist2D'][mode][type][galaxy][view])
        if galaxy==0 and calcdict['ngal']>1:
            histopt['zlim']=(histopt['zlim'][0],histopt['zlim'][1]/2.)
        histopt['zlim']=(histopt['zlim'][0],histopt['zlim'][1]/(10.**2))
        if histopt['zscale']=='log':
            histopt['zmean']=10.**((np.log10(histopt['zlim'][1])+np.log10(histopt['zlim'][0]))/2.-1.)
    view_window=(histopt['xlim'][1]-histopt['xlim'][0],histopt['ylim'][1]-histopt['ylim'][0])
    # Make image
    mat=np.ma.masked_array(histopt['z'],mask=histopt['zmask'])
    inkw_mat={'mn':float(mat.min())  ,'mx':float(mat.max())  ,'cd':float(mat.mean())}
    inkw_def={'mn':histopt['zlim'][0],'mx':histopt['zlim'][1],'cd':histopt['zmean']}
    sclinput=dict(scale=histopt['zscale'],clrmin=0,clip=True,**inkw_def)
#    print 'lim mat: {mn:.3e},{mx:.3e},{cd:.3e}'.format(**inkw_mat)
#    print 'lim def: {mn:.3e},{mx:.3e},{cd:.3e}'.format(**inkw_def)
    matint=mmlplot.map_array(mat,**sclinput).filled(0)
    if outdata: return mat.filled(0),matint,view_window
#    if outdata: return matint,view_window
    matint_T=np.transpose(matint)
    imgpil=Image.fromstring("P",(matint_T.shape[1],matint_T.shape[0]),matint_T.tostring())
    # Add colormap
    cmap=mmlplot.get_cmap(palette)
    imgpil.putpalette(mmlplot.cmap2palette(cmap))
#    imgpil.convert('RGB')
    # Annotate
    if annotate:
      imgsiz = imgpil.size
      fntclr = 255
      draw   = ImageDraw.Draw(imgpil)
      font   = ImageFont.load_default()
      fntwid, fnthgt = font.getsize('test')
      locfct = 20
      xlocpad=imgsiz[0]/locfct
      ylocpad=imgsiz[1]/locfct
#      time2Gyr=nbObj.localsystem_of_units.get_UnitTime_in_s()*3.16888e-8/1.0e9
#      leng2kpc=nbObj.localsystem_of_units.get_UnitLength_in_cm()/3.0856780e+21
      # Label specifying time
      tloc = (xlocpad,ylocpad)
      if residflag:
          tstr = "t1 = {:5.2f} {}, t2 = {:5.2f} {}".format(calcdict['statDict']['time'],calcdict['unitDict']['UnitTime'].symbol,
                                                          calcdict2['statDict']['time'],calcdict2['unitDict']['UnitTime'].symbol)
      else:
          tstr = "t = {:5.2f} {}".format(calcdict['statDict']['time'],calcdict['unitDict']['UnitTime'].symbol)
      draw.text(tloc,tstr,fill=fntclr,font=font)
      # Label specifying size
      dstr = "{:3.0f} x {:3.0f} {}".format(view_window[0],view_window[1],calcdict['unitDict']['UnitLength'].symbol)
      dwid,dhgt=font.getsize(dstr)
      dloc = (imgsiz[0]-xlocpad-dwid,ylocpad)
      draw.text(dloc,dstr,fill=fntclr,font=font)
    # Display image
    if plot_kw['showplot']: imgpil.show()
    # Return output
    if plot_kw['outplot']: return imgpil
    else                 : imgpil.save(plotfile)
#      libutil.mplot(matint,palette=view_palette,save=fname)
    return None

####################################################################################################################################
# METHOD TO PLOT 2D HISTOGRAMS FOR SIMULATION (IMSHOW)
def plotcalc_imshw(calcdict,calcdict2=None,mode=None,type=None,view=None,galaxy=None,
                   plotfile=None,overwrite=None,verbose=None,annotate=None,residflag=None,
                   palette=None,outdata=None,figobj=None,**plot_kw):
    """
    Plots images for a snapshot
    """
    # Set constants
    fnameTAG='calcimshw'
    typListDEF=simlist.LIST_PTYPS
    rayListDEF=list_plotviews()
    # Pars general input
    calcdict=mmlpars.mml_pars(calcdict,type=dict)
    if calcdict2==None: residflag=False
    annotate=mmlpars.mml_pars(annotate,default=True,type=bool)
    residflag=mmlpars.mml_pars(residflag,default=False,type=bool)
    outdata=mmlpars.mml_pars(outdata,default=False,type=bool)
    verbose=mmlpars.mml_pars(verbose,default=False,type=bool)
    # Pars plotting info
    mode=mmlpars.mml_pars(mode,default='pos_m',list=calcdict['h2DList'])
    type=mmlpars.mml_pars(type,default='visi' ,list=calcdict['typList'])
    view=mmlpars.mml_pars(view,default='xy'   ,list=calcdict['rayList'])
    galaxy=mmlpars.mml_pars(galaxy,default=0  ,list=range(3))
    palette=mmlpars.mml_pars(palette,default='pnbody_light',type=str)
    plot_kw=mmlplot.parsopt(plot_kw)
    # Pars file input
    plotfile=mmlpars.mml_pars(plotfile,default='pNbodyplot_{}.png'.format(fnameTAG),type=str)
    overwrite=mmlpars.mml_pars(overwrite,default=False,type=bool)
    if not outdata and not plot_kw['outplot']:
        if not mmlfiles.prep_write(plotfile,overwrite): return None
    # Isolate calcdict options
    if residflag:
        histopt=resid_hist2D(calcdict['hist2D'][mode][type][galaxy][view],calcdict2['hist2D'][mode][type][galaxy][view])
        print 'xlim={}'.format(histopt['xlim'])
        print 'ylim={}'.format(histopt['ylim'])
        print 'zlim={}'.format(histopt['zlim'])
    else:
        histopt=copy.deepcopy(calcdict['hist2D'][mode][type][galaxy][view])
        if galaxy==0 and calcdict['ngal']>1:
            histopt['zlim']=(histopt['zlim'][0],histopt['zlim'][1]/2.)
        histopt['zlim']=(histopt['zlim'][0],histopt['zlim'][1]/(10.**2))
        if histopt['zscale']=='log':
            histopt['zmean']=10.**((np.log10(histopt['zlim'][1])+np.log10(histopt['zlim'][0]))/2.-1.)
    view_window=(histopt['xlim'][1]-histopt['xlim'][0],histopt['ylim'][1]-histopt['ylim'][0])
    # Make image
    mat=np.ma.masked_array(histopt['z'],mask=histopt['zmask'])
#    inkw_mat={'mn':float(mat.min())  ,'mx':float(mat.max())  ,'cd':float(mat.mean())}
    inkw_def={'mn':histopt['zlim'][0],'mx':histopt['zlim'][1],'cd':histopt['zmean']}
    sclinput=dict(scale=histopt['zscale'],clrmin=0,clip=True,**inkw_def)
#    print 'lim mat: {mn:.3e},{mx:.3e},{cd:.3e}'.format(**inkw_mat)
#    print 'lim def: {mn:.3e},{mx:.3e},{cd:.3e}'.format(**inkw_def)
    matint=mmlplot.map_array(mat,**sclinput).filled(0)
    if outdata: return mat.filled(0),matint,view_window
    # Create figure
    if figobj==None:
        figobj=plt.figure()
        ax1=plt.subplot(1,1,1)
    else:
        ax1=figobj.gca()
    extent=(histopt['xlim'][0],histopt['xlim'][1],histopt['ylim'][0],histopt['ylim'][1])
    mmlplot.hist2imshow(matint,cmapname=palette,extent=extent)
    # Annotate
    if annotate:
      fntclr = 'w'
      locfct = 20.
      xlocpad=view_window[0]/locfct
      ylocpad=view_window[1]/locfct
      # Label specifying time
      tloc = (histopt['xlim'][0]+xlocpad,histopt['ylim'][1]-ylocpad)
      if residflag:
          tstr = "t1 = {:5.2f} {}, t2 = {:5.2f} {}".format(calcdict['statDict']['time'],calcdict['unitDict']['UnitTime'].symbol,
                                                           calcdict2['statDict']['time'],calcdict2['unitDict']['UnitTime'].symbol)
      else:
          tstr = "t = {:5.2f} {}".format(calcdict['statDict']['time'],calcdict['unitDict']['UnitTime'].symbol)
      ax1.text(tloc[0],tloc[1],tstr,color=fntclr,ha='left' )
      # Label specifying size
      dstr = "{:3.0f} x {:3.0f} {}".format(view_window[0],view_window[1],calcdict['unitDict']['UnitLength'].symbol)
      dloc = (histopt['xlim'][1]-xlocpad,histopt['ylim'][1]-ylocpad)
      ax1.text(dloc[0],dloc[1],dstr,color=fntclr,ha='right')
    # Display image
    if plot_kw['showplot']: plt.show()
    # Return output
    if plot_kw['outplot']: return figobj
    else                 : figobj.savefig(plotfile)
    return None

####################################################################################################################################
# METHOD TO PLOT 2D HISTOGRAMS FOR SIMULATION (IMSHOW)
def plotcalc_contr(calcdict,calcdict2=None,mode=None,type=None,view=None,galaxy=None,
                   plotfile=None,overwrite=None,verbose=None,annotate=None,residflag=None,
                   palette=None,outdata=None,figobj=None,ncont=None,**plot_kw):
    """
    Plots images for a snapshot
    """
    # Set constants
    fnameTAG='calcimshw'
    typListDEF=simlist.LIST_PTYPS
    rayListDEF=list_plotviews()
    # Pars general input
    calcdict=mmlpars.mml_pars(calcdict,type=dict)
    if calcdict2==None: residflag=False
    annotate=mmlpars.mml_pars(annotate,default=True,type=bool)
    residflag=mmlpars.mml_pars(residflag,default=False,type=bool)
    outdata=mmlpars.mml_pars(outdata,default=False,type=bool)
    verbose=mmlpars.mml_pars(verbose,default=False,type=bool)
    ncont=mmlpars.mml_pars(ncont,default=11,type=int,min=1)
    # Pars plotting info
    mode=mmlpars.mml_pars(mode,default='pos_m',list=calcdict['h2DList'])
    type=mmlpars.mml_pars(type,default='visi' ,list=calcdict['typList'])
    view=mmlpars.mml_pars(view,default='xy'   ,list=calcdict['rayList'])
    galaxy=mmlpars.mml_pars(galaxy,default=0  ,list=range(3))
    palette=mmlpars.mml_pars(palette,default='pnbody_rainbow4',type=str)
    plot_kw=mmlplot.parsopt(plot_kw)
    # Pars file input
    plotfile=mmlpars.mml_pars(plotfile,default='pNbodyplot_{}.png'.format(fnameTAG),type=str)
    overwrite=mmlpars.mml_pars(overwrite,default=False,type=bool)
    if not outdata and not plot_kw['outplot']:
        if not mmlfiles.prep_write(plotfile,overwrite): return None
    # Isolate calcdict options
    if residflag:
        histopt=resid_hist2D(calcdict['hist2D'][mode][type][galaxy][view],calcdict2['hist2D'][mode][type][galaxy][view])
    else:
        histopt=copy.deepcopy(calcdict['hist2D'][mode][type][galaxy][view])
        if histopt['zmode'] not in ['scfpot','scfpotm2']:
            if galaxy==0 and calcdict['ngal']>1:
                histopt['zlim']=(histopt['zlim'][0],histopt['zlim'][1]/2.)
            histopt['zlim']=(histopt['zlim'][0],histopt['zlim'][1]/(10.**2))
        if histopt['zscale']=='log':
            histopt['zmean']=10.**((np.log10(histopt['zlim'][1])+np.log10(histopt['zlim'][0]))/2.-1.)
    view_window=(histopt['xlim'][1]-histopt['xlim'][0],histopt['ylim'][1]-histopt['ylim'][0])
    # Make image
#    mmlio.yorn(str((histopt['zmode'],histopt['z'].min(),histopt['z'].max())))
    mat=np.ma.masked_array(histopt['z'],mask=histopt['zmask'])
    inkw_mat={'mn':float(mat.min())  ,'mx':float(mat.max())  ,'cd':float(mat.mean())}
    inkw_def={'mn':histopt['zlim'][0],'mx':histopt['zlim'][1],'cd':histopt['zmean']}
    sclinput=dict(scale=histopt['zscale'],clrmin=0,clip=True,**inkw_def)
#    sclinput['scale']='linear'
#    mmlio.yorn('mode={}, lim={}, scale={}'.format(histopt['zmode'],histopt['zlim'],histopt['zscale']))
#    print 'lim mat: {mn:.3e},{mx:.3e},{cd:.3e}'.format(**inkw_mat)
#    print 'lim def: {mn:.3e},{mx:.3e},{cd:.3e}'.format(**inkw_def)
#    matint=mmlplot.map_array(mat,**sclinput).filled(0)
    matint=mat
    if outdata: return mat.filled(0),matint,view_window
    # Create figure
    if figobj==None:
        figobj=plt.figure()
        ax1=plt.subplot(1,1,1)
    else:
        ax1=figobj.gca()
    extent=(histopt['xlim'][0],histopt['xlim'][1],histopt['ylim'][0],histopt['ylim'][1])
    mmlplot.hist2contour(matint,cmapname=palette,extent=extent,ncont=ncont,**sclinput)
    # Annotate
    if annotate:
      fntclr = 'w'
      locfct = 20.
      xlocpad=view_window[0]/locfct
      ylocpad=view_window[1]/locfct
      # Label specifying time
      tloc = (histopt['xlim'][0]+xlocpad,histopt['ylim'][1]-ylocpad)
      if residflag:
          tstr = "t1 = {:5.2f} {},t2 = {:5.2f} {}".format(calcdict['statDict']['time'],calcdict['unitDict']['UnitTime'].symbol,
                                                          calcdict2['statDict']['time'],calcdict2['unitDict']['UnitTime'].symbol)
      else:
          tstr = "t = {:5.2f} {}".format(calcdict['statDict']['time'],calcdict['unitDict']['UnitTime'].symbol)
      ax1.text(tloc[0],tloc[1],tstr,color=fntclr,ha='left' )
      # Label specifying size
      dstr = "{:3.0f} x {:3.0f} {}".format(view_window[0],view_window[1],calcdict['unitDict']['UnitLength'].symbol)
      dloc = (histopt['xlim'][1]-xlocpad,histopt['ylim'][1]-ylocpad)
      ax1.text(dloc[0],dloc[1],dstr,color=fntclr,ha='right')
    # Display image
    if plot_kw['showplot']: plt.show()
    # Return output
    if plot_kw['outplot']: return figobj
    else                 : figobj.savefig(plotfile)
    return None

####################################################################################################################################
# METHOD TO CREATE MULTIPANEL PLOT FOR SIMULATION SNAPSHOT
def plotcalc_multi(calcdict,calcdict0=None,
                   h1DList=None,h2DList_imag=None,h2DList_cont=None,typList=None,galList=None,rayList=None,
                   plotfile=None,overwrite=None,plotobj=None,**plot_kw):
    """
    Plot multiple panels of info on a snapshot
    """
    # Set constants
    fnameTAG='calcplot'
    h1DListDEF=['rxyz_rho','z_n','vz_n','phi_n','rxyz_rho_resid']
    h2DList_imagDEF=['pos_m','pos_vr']#,'pos_vr','pos_gal']
    h2DList_contDEF=['pos_scfpotm2',None]
    typListDEF=simlist.LIST_PTYPBASE
    galListDEF=0
    rayListDEF=['xy','yz','xz']
    # Pars plotting input
    options=mmlplot.parsopt(plot_kw)
    options['pad']=5. ; options['wpad']=2.2 ; options['tickfontsize']=12. ; options['textfontsize']=14.
    options['xlabpos']=(0.5,-0.12) ; options['ylabpos']=(-0.12,0.5)
    h1DList=mmlpars.mml_pars(h1DList,default=h1DListDEF,type=list) ; nh1D=len(h1DList)
    h2DList_imag=mmlpars.mml_pars(h2DList_imag,default=h2DList_imagDEF,type=list,nelements=2) ; nh2D_imag=len(h2DList_imag)
    h2DList_cont=mmlpars.mml_pars(h2DList_cont,default=h2DList_contDEF,type=list,nelements=2) ; nh2D_cont=len(h2DList_cont)
    typList=mmlpars.mml_pars(typList,default=typListDEF,type=list) ; ntyp=len(typList)
#    galList=mmlpars.mml_pars(galList,default=galListDEF,type=list) ; ngal=len(galList)
    rayList=mmlpars.mml_pars(rayList,default=rayListDEF,type=list) ; nray=len(rayList)
    # Pars file input
    overwrite=mmlpars.mml_pars(overwrite,default=False,type=bool)
    if not options['outplot']:
        if not mmlfiles.prep_write(plotfile,overwrite): return
    # Count rows and columns
    ntyp_plot=0
    for ityp in typList:
        if ityp in calcdict['typList']: ntyp_plot+=1
    nray_plot=0
    for iray in rayList:
        if iray in calcdict['rayList']: nray_plot+=1
    nh1D_plot=0
    for ih1D in h1DList:
        if ih1D in calcdict['h1DList']: nh1D_plot+=1
        else:
            if ih1D.split('_')[-1] == 'resid': nh1D_plot+=1
    nrow=max([ntyp_plot,nh1D_plot])
    ncol=nray_plot
    # Determine size of figure
    fontsize_inc=max(options['tickfontsize'],options['textfontsize'])/72.
    opad_inc=options['pad']*fontsize_inc ; ipad_inc=options['wpad']*fontsize_inc
    options['figsize']=(4*opad_inc+2*(ncol+1)*options['axsize'][0]+2*(ncol-1)*ipad_inc,
                        2*opad_inc+nrow*options['axsize'][1]+(nrow-1)*ipad_inc)
    # Set dimensions relative to figure
    opad_relx=opad_inc/options['figsize'][0] ; opad_rely=opad_inc/options['figsize'][1]
    ipad_relx=ipad_inc/options['figsize'][0] ; ipad_rely=ipad_inc/options['figsize'][1]
    rect_rho=np.array([0.,1.-opad_rely,options['axsize'][0]/options['figsize'][0],options['axsize'][1]/options['figsize'][1]])
    rect_rad=np.array([0.,1.-opad_rely,2.*options['axsize'][0]/options['figsize'][0],options['axsize'][1]/options['figsize'][1]])
    rect_vel=np.array([0.,1.-opad_rely,options['axsize'][0]/options['figsize'][0],options['axsize'][1]/options['figsize'][1]])
    rect_rho[0]=opad_relx
    rect_rad[0]=rect_rho[0]+ncol*rect_rho[2]+(ncol-1)*ipad_relx+opad_relx
    rect_vel[0]=rect_rad[0]+rect_rad[2]+opad_relx
    rect_rho_cb=np.array([rect_rho[0],rect_rho[1]-ntyp_plot*rect_rho[3]-ntyp_plot*ipad_rely,ncol*rect_rho[2]+(ncol-1)*ipad_relx,ipad_rely])
    rect_vel_cb=np.array([rect_vel[0],rect_vel[1]-ntyp_plot*rect_vel[3]-ntyp_plot*ipad_rely,ncol*rect_vel[2]+(ncol-1)*ipad_relx,ipad_rely])
    rect_txt=np.array([rect_vel_cb[0],rect_vel_cb[1]-ntyp_plot*rect_vel_cb[3]-(ntyp_plot-1)*ipad_rely-ipad_rely,rect_vel[2],rect_vel[3]])
    # Create plot objects
    if plotobj==None:
        newfig=True
        plotobj=dict(fig=plt.figure(),rho=None,rad=None,vel=None,rho_cb=None,vel_cb=None,txt=None)
    else:
        newfig=False
    # Plots of mass vs. position
    modList1=[h2DList_imag[0],h2DList_cont[0]]
#    imag_mode1='pos_m'
#    cont_mode1='pos_scfpotm2'
    plotobj['rho']=panels_hist2D(plotobj['fig'],typList,rayList,modList1,0,calcdict['hist2D'],
                                 unitdict=calcdict['unitDict'],rect=rect_rho,wpad=ipad_relx,hpad=ipad_rely)
    plotobj['rho_cb']=panel_colorbar(plotobj['fig'],[plotobj['rho']['imagList'],plotobj['rho']['contList']],
                                     rect=rect_rho_cb,wpad=ipad_relx,hpad=ipad_rely,label=[modList1[0].split('_')[1],modList1[1].split('_')[1]])
    # 1D histograms
    if calcdict0!=None:
        plotobj['rad']=panels_hist1D(plotobj['fig'],typList,h1DList,1,calcdict['hist1D'],calcdict0['hist1D'],
                                     unitdict=calcdict['unitDict'],rect=rect_rad,wpad=ipad_relx,hpad=ipad_rely)
    else:
        plotobj['rad']=panels_hist1D(plotobj['fig'],typList,h1DList,1,calcdict['hist1D'],
                                     unitdict=calcdict['unitDict'],rect=rect_rad,wpad=ipad_relx,hpad=ipad_rely)
    # Plots of velocity vs. position
    modList2=[h2DList_imag[1],h2DList_cont[1]]
#    imag_mode2='pos_vr'
    plotobj['vel']=panels_hist2D(plotobj['fig'],typList,rayList,modList2,0,calcdict['hist2D'],
                                 unitdict=calcdict['unitDict'],rect=rect_vel,wpad=ipad_relx,hpad=ipad_rely)
    plotobj['vel_cb']=panel_colorbar(plotobj['fig'],[plotobj['vel']['imagList'],plotobj['vel']['contList']],
                                     rect=rect_vel_cb,wpad=ipad_relx,hpad=ipad_rely,label=[modList2[0].split('_')[1]])
    # Text panel
    plotobj['txt']=panel_text(plotobj['fig'],calcdict['statDict'],calcdict['unitDict'],rect=rect_txt,wpad=ipad_relx,hpad=ipad_rely)
    # Set figure options
    if newfig: mmlplot.set_figprop(plotobj['fig'],options)
    # Display image
    if options['showplot']: plt.show()
    # Return output
    if options['outplot']: return plotobj['fig']
    else                 : plotobj['fig'].savefig(plotfile,edgecolor=options['fgclr'],facecolor=options['bgclr'])
    return plotobj

####################################################################################################################################
# METHOD TO PLOT 1D HISTOGRAMS FOR SIMULATION
def plotcalc_hist1D(calcdict,calcdict2=None,
                    typList=None,modList=None,galList=None,residualflag=None,
                    plotfile=None,overwrite=None,verbose=None,
                    labelList=None,flag_legend=None,
                    rowvar=None,colvar=None,plotobj=None,**plot_kw):
##     plotfile=None,calcfile=None,calcdict=None,calcfile2=None,calcdict2=None,
##                     typList=None,pgal=None,hist_mode=None,residualflag=None,
##                     overwrite=None,verbose=None,flag_legend=None,labelList=None,
##                     plotobj=None,**plot_kw):
    """
    Plots 1D histograms for a snapshot
    """
    # Set constants
    varLIST=['type','mode','galaxy']
    fnameTAG='snaphist1D'
    typListDEF=simlist.LIST_PTYPBASE
    modListDEF=['rxyz_rho']
    galListDEF=[1]
    # Pars general input
    verbose=mmlpars.mml_pars(verbose,default=False,type=bool)
    if verbose: mmlio.verbose('Parsing input...')
    # Pars file input
    overwrite=mmlpars.mml_pars(overwrite,default=False,type=bool)
    residualflag=mmlpars.mml_pars(residualflag,default=False,type=bool)
    if not mmlfiles.prep_write(plotfile,overwrite): return None
    if not isinstance(calcdict,dict): raise Exception('Calcdict must be a dictionary.')
    if calcdict2==None:
        residualflag=False
        calcdict2={'hist1D':None}
    else:
        if not isinstance(calcdict2,dict): raise Exception('Calcdict2 must be a dictionary')
    # Pars variable options
    if verbose: mmlio.verbose('Setting list options...')
    modList=mmlpars.mml_pars(modList,default=modListDEF,type=list) ; nmod=len(modList)
    typList=mmlpars.mml_pars(typList,default=typListDEF,type=list) ; ntyp=len(typList)
    galList=mmlpars.mml_pars(galList,default=galListDEF,type=list) ; ngal=len(galList)
    varDict={'mode':modList,'type':typList,'galaxy':galList}
    if verbose: mmlio.verbose('modList={}'.format(modList))
    if verbose: mmlio.verbose('typList={}'.format(typList))
    if verbose: mmlio.verbose('galList={}'.format(galList))
    # Pars row/col options
    rowvar=mmlpars.mml_pars(rowvar,default=varLIST[0],list=varLIST)
    colvar=mmlpars.mml_pars(colvar,default=varLIST[1],list=varLIST)
    for ivar in varLIST:
        if ivar!=rowvar and ivar!=colvar: otrvar=ivar
    if rowvar==colvar:
        raise Exception('Both rows and columns cannot vary {}'.format(rowvar))
    if len(varDict[otrvar])>1:
        raise Exception('{} cannot vary across axes because both rows & columns are taken.')
    # Pars plotting input
    if verbose: mmlio.verbose('Setting plot options...')
    options=mmlplot.parsopt(plot_kw)
    options['axsize']=(5.,2.)
    options['pad']=5. ; options['wpad']=4. ; options['tickfontsize']=12. ; options['textfontsize']=14.
    options['xlabpos']=(0.5,-0.18) ; options['ylabpos']=(-0.10,0.5)
    # Count rows and columns
    nvar_plot={'mode':0,'type':0,'galaxy':0}
    varDict_calc={'mode'  :calcdict['h1DList'],
                  'type'  :calcdict['typList'],
                  'galaxy':calcdict['galList']}
    for ivar in varLIST:
        for iele in varDict[ivar]:
            if iele in varDict_calc[ivar]: nvar_plot[ivar]+=1
    nrow=nvar_plot[rowvar]
    ncol=nvar_plot[colvar]
    if residualflag: ncol*=2
    # Add residual flags to modes
    if residualflag:
        modList_his=[[imod,imod+'_compr'] for imod in modList]
        modList_res=[imod+'_resid' for imod in modList]
    else           :
        modList_his=modList
        modList_res=modList
    # Determine size of figure
    fontsize_inc=max(options['tickfontsize'],options['textfontsize'])/72.
    opad_inc=options['pad']*fontsize_inc ; ipad_inc=options['wpad']*fontsize_inc
    options['figsize']=(2*opad_inc+ncol*options['axsize'][0]+(ncol-1)*ipad_inc,
                        2*opad_inc+nrow*options['axsize'][1]+(nrow-1)*ipad_inc)
    # Set dimensions relative to figure
    opad_relx=opad_inc/options['figsize'][0] ; opad_rely=opad_inc/options['figsize'][1]
    ipad_relx=ipad_inc/options['figsize'][0] ; ipad_rely=ipad_inc/options['figsize'][1]
    rect_his=np.array([1.2*opad_relx,1.-opad_rely,
                       options['axsize'][0]/options['figsize'][0],options['axsize'][1]/options['figsize'][1]])
    if residualflag:
        rect_res=copy.deepcopy(rect_his)
        rect_res[0]+=1.25*ipad_relx+options['axsize'][0]/options['figsize'][0]
    # Create plot objects
    if plotobj==None:
        newfig=True
        plotobj=dict(fig=plt.figure(),hist=None,resid=None)
    else:
        newfig=False
    # 1D histograms
    if verbose: mmlio.verbose('Plotting histograms...')
    plotobj['hist' ]=panels_hist1D(plotobj['fig'],typList,modList_his,galList,
                                   calcdict=calcdict['hist1D'],calcdict0=calcdict2['hist1D'],unitdict=calcdict['unitDict'],
                                   rect=rect_his,wpad=ipad_relx,hpad=ipad_rely,rowvar=rowvar,colvar=colvar,
                                   flag_legend=flag_legend,labelList=labelList,verbose=verbose)
    if residualflag:
        if verbose: mmlio.verbose('Plotting residuals...')
        plotobj['resid']=panels_hist1D(plotobj['fig'],typList,modList_res,galList,
                                       calcdict=calcdict['hist1D'],calcdict0=calcdict2['hist1D'],unitdict=calcdict['unitDict'],
                                       rect=rect_res,wpad=ipad_relx,hpad=ipad_rely,rowvar=rowvar,colvar=colvar,
                                       flag_legend=flag_legend,verbose=verbose)
    # Set properties
    if newfig:
        mmlplot.set_figprop(plotobj['fig'],options)
        plotobj['fig'].suptitle(options['title'],fontsize=options['textfontsize'],color=options['fgclr'])
    # Display plot
    if options['showplot']: plotobj['fig'].show()
    # Save figure/return plot objects
    if options['outplot']:
        return plotobj
    else:
        plotobj['fig'].savefig(plotfile,edgecolor=options['fgclr'],facecolor=options['bgclr'])
        return None

####################################################################################################################################
# METHOD TO PLOT 2D HISTOGRAMS FOR SIMULATION
def plotcalc_hist2D(calcdict,calcdict2=None,
                    typList=None,rayList=None,modList=None,galList=None,residualflag=None,
                    plotfile=None,overwrite=None,verbose=None,
                    rowvar=None,colvar=None,plotobj=None,**plot_kw):
##     plotfile=None,calcfile=None,calcdict=None,calcdict2=None,
##                         typList=None,rayList=None,pgal=None,imag_mode=None,cont_mode=None,
##                         overwrite=None,verbose=None,plotobj=None,
##                         imag_residflag=None,cont_residflag=None,skipimag=None,skipcont=None,**plot_kw):
    """
    Plots 2D histograms for a snapshot
    """
    # Set constants
    varLIST=['type','view','mode']
    fnameTAG='snaphist2D'
    imag_modeDEF='pos_m'
    cont_modeDEF='pos_scfpotm2'
    typListDEF=simlist.LIST_PTYPBASE
    rayListDEF=['xy','yz','xz']
    modListDEF=['pos_m','pos_scfpotm2']
    galListDEF=[1]
    # Pars general input
    verbose=mmlpars.mml_pars(verbose,default=False,type=bool)
    if verbose: mmlio.verbose('Parsing input...')
    # Pars file input
    overwrite=mmlpars.mml_pars(overwrite,default=False,type=bool)
    residualflag=mmlpars.mml_pars(residualflag,default=False,type=bool)
    if not mmlfiles.prep_write(plotfile,overwrite): return
    if not isinstance(calcdict,dict): raise Exception('Calcdict must be a dictionary.')
    if calcdict2==None:
        residualflag=False
        calcdict2={'hist2D':None}
    else:
        if not isinstance(calcdict2,dict): raise Exception('Calcdict2 must be a dictionary')
    # Pars variable options
    if verbose: mmlio.verbose('Setting list options...')
    modList=mmlpars.mml_pars(modList,default=modListDEF,type=list) ; nmod=len(modList)
    typList=mmlpars.mml_pars(typList,default=typListDEF,type=list) ; ntyp=len(typList)
    galList=mmlpars.mml_pars(galList,default=galListDEF,type=list) ; ngal=len(galList)
    rayList=mmlpars.mml_pars(rayList,default=rayListDEF,type=list) ; nray=len(rayList)
    varDict={'mode':modList,'type':typList,'galaxy':galList,'view':rayList}
    if verbose: mmlio.verbose('modList={}'.format(modList))
    if verbose: mmlio.verbose('typList={}'.format(typList))
    if verbose: mmlio.verbose('galList={}'.format(galList))
    if verbose: mmlio.verbose('rayList={}'.format(rayList))
    # Pars row/col options
    rowvar=mmlpars.mml_pars(rowvar,default=varLIST[0],list=varLIST)
    colvar=mmlpars.mml_pars(colvar,default=varLIST[1],list=varLIST)
    for ivar in varLIST:
        if ivar!=rowvar and ivar!=colvar: otrvar=ivar
    if rowvar==colvar:
        raise Exception('Both rows and columns cannot vary {}'.format(rowvar))
    if len(varDict[otrvar])>1:
        raise Exception('{} cannot vary across axes because both rows & columns are taken.')
    # Pars plotting input
    if verbose: mmlio.verbose('Setting plot options...')
    options=mmlplot.parsopt(plot_kw)
    options['pad']=5. ; options['wpad']=2.2 ; options['tickfontsize']=12. ; options['textfontsize']=14.
    options['xlabpos']=(0.5,-0.12) ; options['ylabpos']=(-0.12,0.5)
    # Count rows and columns
    nvar_plot={'mode':0,'type':0,'galaxy':0,'view':0}
    varDict_calc={'mode'  :calcdict['h1DList'],
                  'type'  :calcdict['typList'],
                  'galaxy':calcdict['galList'],
                  'view'  :calcdict['rayList']}
    for ivar in varLIST:
        for iele in varDict[ivar]:
            if iele in varDict_calc[ivar]: nvar_plot[ivar]+=1
    nrow=nvar_plot[rowvar]
    ncol=nvar_plot[colvar]
    notr=nvar_plot[otrvar]
    if notr>2: raise Exception('Only 2 methods are supported for z direction plotting.')
    # Add residual info to modes
    if residualflag: modList=[imod+'_resid' for imod in modList]
    # Determine size of figure
    fontsize_inc=max(options['tickfontsize'],options['textfontsize'])/72.
    opad_inc=options['pad']*fontsize_inc ; ipad_inc=options['wpad']*fontsize_inc
    options['figsize']=(2*opad_inc+ncol*options['axsize'][0]+(ncol-1)*ipad_inc,
                        2*opad_inc+nrow*options['axsize'][1]+(nrow-1)*ipad_inc+2*nrow*ipad_inc)
    # Set dimensions relative to figure
    opad_relx=opad_inc/options['figsize'][0] ; opad_rely=opad_inc/options['figsize'][1]
    ipad_relx=ipad_inc/options['figsize'][0] ; ipad_rely=ipad_inc/options['figsize'][1]
    rect=np.array([opad_relx,1.-opad_rely,options['axsize'][0]/options['figsize'][0],options['axsize'][1]/options['figsize'][1]])
    rect_cb=np.array([rect[0],rect[1]-ntyp_plot*rect[3]-ntyp_plot*ipad_rely,ncol*rect[2]+(ncol-1)*ipad_relx,ipad_rely])
    # Create plot objects
    if plotobj==None:
        newfig=True
        plotobj=dict(fig=plt.figure(),hist=None,cbar=None)
    else:
        newfig=False
    # Plots of mass vs. position
    if verbose: mmlio.verbose('Plotting histograms...')
    plotobj['hist']=panels_hist2D(plotobj['fig'],typList,rayList,modList,galaxy,rowvar=rowvar,colvar=colvar,
                                  calcdict=calcdict['hist2D'],calcdict0=calcdict2['hist2D'],unitdict=calcdict['unitDict'],
                                  rect=rect,wpad=ipad_relx,hpad=ipad_rely,verbose=verbose)
    if verbose: mmlio.verbose('Adding colorbar...')
    labList=[ivar for ivar in varDict[otrvar]]
    plotobj['cbar']=panel_colorbar(plotobj['fig'],[plotobj['hist']['imagList'],plotobj['hist']['contList']],
                                   rect=rect_cb,wpad=ipad_relx,hpad=ipad_rely,label=labList)
    # Set figure options
    if newfig:
        mmlplot.set_figprop(plotobj['fig'],options)
        plotobj['fig'].suptitle(options['title'],fontsize=options['textfontsize'],color=options['fgclr'])
    # Display figure
    if options['showplot']: plotobj['fig'].show()
    # Output/save figure
    if options['outplot']:
        return plotobj
    else:
        plotobj['fig'].savefig(plotfile,edgecolor=options['fgclr'],facecolor=options['bgclr'])
        return None

####################################################################################################################################
####################################################################################################################################
# METHODS TO MAKE HISTOGRAM PANELS

####################################################################################################################################
# METHOD TO MAKE COLORBAR PANEL
def panel_colorbar(figObj,mapObjArray=None,rect=[1.,1.,1.,1.],wpad=0.,hpad=0.,label='',ntick=3):
    """
    Creates a panel for a colorbar
    """
    # Pars input
    if not isinstance(mapObjArray,list): mapObjArray=[mapObjArray]
    ncol=0
    nrow_list=[]
    for icol in range(len(mapObjArray)):
        if not np.all(np.array(mapObjArray[icol])==None): ncol+=1
        nrow_list.append(np.ones(np.array(mapObjArray).shape)[np.array(mapObjArray)!=None].sum())
    plotobj=dict.fromkeys(['axes','cbar'],[])
    for icol in range(ncol):
        plotobj['axes'].append(max(nrow_list)*[None])
        plotobj['cbar'].append(max(nrow_list)*[None])
    plotobj['text']=ncol*[None]
    # Create colorbars
    rect_end=ncol*[None]
    icol=-1
    for mapObjList in mapObjArray:
        irow=-1
        for mapObj in mapObjList:
            if mapObj==None: continue
            irow+=1
            if irow==0:
                icol+=1
            # Set position
            irect=copy.deepcopy(rect)
            irect[2]=(rect[2]-(ncol-1)*wpad)/ncol
            irect[0]+=icol*(irect[2]+wpad)
            irect[1]+=(-(irow+1)*irect[3]-irow*hpad)
            # mmlio.verbose('{}: {} rect={}'.format(label[icol],(irow,icol),irect))
            # Create colorbar axes
            if plotobj['axes'][icol][irow]==None:
                plotobj['axes'][icol][irow]=figObj.add_axes(list(irect))
            # Create tick locator
#            cbar_tickLocator.set_params(nbins=6)
            vmin=mapObj.norm.vmin
            vmax=mapObj.norm.vmax
            if mapObj.norm.__class__ == matplotlib.colors.LogNorm().__class__:
                tick_formatter=matplotlib.ticker.LogFormatter(labelOnlyBase=False)
                vpad=(np.log10(vmax)-np.log10(vmin))/(2.*ntick)
                ticks=np.logspace(np.log10(vmin)+vpad,np.log10(vmax)-vpad,ntick)
                # mmlio.verbose('LogNorm: {}'.format(ticks))
            else:
                tick_formatter=matplotlib.ticker.ScalarFormatter()
                vpad=(vmax-vmin)/(2.*ntick)
                ticks=np.linspace(vmin+vpad,vmax-vpad,ntick)
                tick_formatter.set_scientific(True)
                tick_formatter.set_powerlimits((-1,1))
                #mmlio.verbose('Normalize: {}'.format(ticks))
#                tick_formatter=matplotlib.ticker.LogFormatter(labelOnlyBase=False)
            # Create colorbar object
#            if plotobj['cbar'][icol][irow]==None:
#                plotobj['cbar'][icol][irow]=
            figObj.colorbar(mapObj,plotobj['axes'][icol][irow],orientation='horizontal',
                            ticks=list(ticks),format='%1.1g')#tick_formatter)
#                            ticks=cbar_tickLocator)
            # Set axes properties
            plotobj['axes'][icol][irow].tick_params(labelleft=False,labelright=False,
                                                    labeltop=False,labelbottom=True)
        # Create text object
        rect_end[icol]=[irect[0]+irect[2],irect[1]-irect[3],irect[2],irect[3]]
        plotobj['text'][icol]=figObj.text(irect[0]-wpad/2.5,rect_end[icol][1]+(rect[1]-rect_end[icol][1])/2.+hpad/2.,label[icol],
                                          rotation='vertical',ha='center',va='baseline')
    # Return output
    return plotobj


####################################################################################################################################
# METHOD TO MAKE TEXT PANEL
def panel_text(figObj,sdict,udict,rect=[1.,1.,1.,1.],wpad=0.,hpad=0.):
    """
    Creates a panel of text info on a simulation
    """
    figObj.suptitle('%(fname)s' % sdict)
    linekeys=['time','m2','com','vcom','halo','disk','bulge','bar']
    linelist=dict(
        time =r'%(time)#6.2f '+udict['UnitTime'].symbol,
        m2   =r'%(bar_m2)#6.2f',
        com  =r'(%(comx)#5.2f,%(comy)#5.2f,%(comz)#5.2f) '+udict['UnitLength'].symbol,
        vcom =r'(%(vcomx)#8.2f,%(vcomy)#8.2f,%(vcomz)#8.2f) '+udict['UnitVelocity'].symbol,
        halo =r'$A$=%(halo_A)#8.2f ; r=%(halo_r)#8.2f',
        disk =r'$A$=%(disk_A)#8.2f ; r=%(disk_r)#8.2f ; z=%(disk_z)#8.2f',
        bulge=r'$A$=%(bulge_A)#8.2f ; r=%(bulge_r)#8.2f',
        bar  =r'$A$=%(bar_amp)#8.2f ; r=%(bar_x)#8.2f ; z=%(bar_z)#8.2f ; $\phi$=%(bar_phi)#8.2f')
#        dim=(%(bar_x)#8.2f,%(bar_y)#8.2f,%(bar_z)86.2f) ; m2=%(bar_m2)#5.2f')
    iline=''
    for ilinekey in linekeys: iline+='{:8s}: '.format(ilinekey)+linelist[ilinekey]+' \n'
    textObj=figObj.text(rect[0],rect[1],iline % sdict,va='top',
                        fontproperties=matplotlib.font_manager.FontProperties('monospace'))
    return textObj


####################################################################################################################################
# METHOD TO MAKE 1D HISTOGRAM PANEL
def panels_hist1D(figObj,typList,modList,galList,calcdict={},calcdict0={},unitdict=None,
                  rect=[1.,1.,0.,0.],wpad=0.,hpad=0.,flag_legend=None,labelList=None,
                  rowvar=None,colvar=None,clrvar=None,styvar=None,verbose=None):
    """
    Create 1D histogram panels
    """
    # Set constants
    varLIST=['mode','type','galaxy']
    typListDEF=simlist.LIST_PTYPBASE
    # Pars input
    verbose=mmlpars.mml_pars(verbose,default=False,type=bool)
    if verbose: mmlio.verbose('Parsing input...')
    typList=mmlpars.mml_pars(typList,type=list)
    modList=mmlpars.mml_pars(modList,type=list)
    galList=mmlpars.mml_pars(galList,type=list)
    varDict={'mode':modList,'type':typList,'galaxy':galList}
    flag_legend=mmlpars.mml_pars(flag_legend,default=False,type=bool)
    # Check row/col variables
    rowvar=mmlpars.mml_pars(rowvar,default=varLIST[0],list=varLIST)
    colvar=mmlpars.mml_pars(colvar,default=varLIST[1],list=varLIST)
    if rowvar==colvar: raise Excpetion('Both rows and columns cannot vary {}'.format(rowvar))
    for ivar in varLIST:
        if ivar!=rowvar and ivar!=colvar: otrvar=ivar
    if len(varDict[otrvar])>1:
        raise Exception('Variable {} not allowed for rowvar={} and colvar={}'.format(otrvar,rowvar,colvar))
    clrvarDEF='line'
    clrvar=mmlpars.mml_pars(clrvar,default=clrvarDEF,list=varLIST+['line','none'])
    if clrvar=='line': styvarDEF='none'
    else             : styvarDEF='line'
    styvar=mmlpars.mml_pars(styvar,default=styvarDEF,list=varLIST+['line','none'])
    # Initialize plot object dictionary
    modList_plotobj=[] ; typList_plotobj=[] ; galList_plotobj=[]
    for imod0 in modList:
        if isinstance(imod0,list):
            for imod in imod0: modList_plotobj.append(imod)
        else: modList_plotobj.append(imod0)
    for ityp0 in typList:
        if isinstance(ityp0,list):
            for ityp in ityp0: typList_plotobj.append(ityp)
        else: typList_plotobj.append(ityp0)
    for igal0 in galList:
        if isinstance(igal0,list):
            for igal in igal0: galList_plotobj.append(igal)
        else: galList_plotobj.append(igal0)
    plotobj=mmlclass.mydict_nd({},modList_plotobj,typList_plotobj,galList_plotobj)
    # Initalize variables
    plotobj['lineList']=[]
    modList_plot=[] ; modLine_plot={}
    typList_plot=[] ; typLine_plot={}
    galList_plot=[] ; galLine_plot={}
    if verbose: mmlio.verbose('Beginning loop over modes & types...')
    # Loop over modList
    for imod0 in modList:
        # Handle multiple mode lines
        imod0_idx=modList.index(imod0)
        if imod0_idx not in modLine_plot: modLine_plot[imod0_idx]=[]
        if isinstance(imod0,list):
            imodList=imod0   ; imodLines=True
        else:
            imodList=[imod0] ; imodLines=False
        for imod in imodList:
            # Handle option for residual/comparison plotting
            flag_resid=False ; flag_compr=False ; imod_short=imod
            if 'resid' in imod.split('_')[-1]:
                flag_resid=True
                imod_short='_'.join(imod.split('_')[:-1])
            if 'compr' in imod.split('_')[-1]:
                flag_compr=True
                imod_short='_'.join(imod.split('_')[:-1])
            # Skip if data dosn't exist for a mode
            if imod_short not in calcdict.keys():
                if verbose: mmlio.verbose('Mode {} not in calcdict. Skipping...'.format(imod_short))
                continue
            if (flag_resid or flag_compr) and imod_short not in calcdict0.keys():
                if verbose: mmlio.verbose('Mode {} not in calcdict0. Skipping...'.format(imod_short))
                continue
            if imod0 not in modList_plot: modList_plot.append(imod0)
            if imod not in modLine_plot[imod0_idx]: modLine_plot[imod0_idx].append(imod)
            # Loop over typList
            for ityp0 in typList:
                # Handle multiple type lines
                ityp0_idx=typList.index(ityp0)
                if ityp0_idx not in typLine_plot: typLine_plot[ityp0_idx]=[]
                if isinstance(ityp0,list):
                    itypList=ityp0   ; itypLines=True
                else:
                    itypList=[ityp0] ; itypLines=False
                for ityp in itypList:
                    # Skip if data dosn't exist for a type
                    if ityp not in calcdict[imod_short]:
                        if verbose: mmlio.verbose('Type {} not in calcdict. Skipping...'.format(ityp))
                        continue
                    if (flag_resid or flag_compr) and ityp not in calcdict0[imod_short]:
                        if verbose: mmlio.verbose('Type {} not in calcdict0. Skipping...'.format(ityp))
                        continue
                    if ityp0 not in typList_plot: typList_plot.append(ityp0)
                    if ityp not in typLine_plot[ityp0_idx]: typLine_plot[ityp0_idx].append(ityp)
                    # Loop over galaxies
                    for igal0 in galList:
                        # Handle multiple galaxy lines
                        igal0_idx=galList.index(igal0)
                        if igal0_idx not in galLine_plot: galLine_plot[igal0_idx]=[]
                        if isinstance(igal0,list):
                            igalList=igal0   ; igalLines=True
                        else:
                            igalList=[igal0] ; igalLines=False
                        for igal in igalList:
                            # Skip if data dosn't exist for a galaxy
                            if igal not in calcdict[imod_short][ityp]:
                                if verbose: mmlio.verbose('Galaxy {} not in calcdict. Skipping...'.format(igal))
                                continue
                            if (flag_resid or flag_compr) and igal not in calcdict0[imod_short][ityp]:
                                if verbose: mmlio.verbose('Galaxy {} not in calcdict0. Skipping...'.format(igal))
                                continue
                            if igal0 not in galList_plot: galList_plot.append(igal0)
                            if igal not in galLine_plot[igal0_idx]: galLine_plot[igal0_idx].append(igal)
                            # Get info on mode & type
                            if imodLines+itypLines+igalLines>1:
                                raise Exception('Multiple mode/type/galaxy lines not supported.')
                            varList_plot={'mode':modList_plot,'type':typList_plot,'galaxy':galList_plot}
                            ivar ={'mode':imod ,'type':ityp ,'galaxy':igal }
                            ivar0={'mode':imod0,'type':ityp0,'galaxy':igal0}
                            if verbose: mmlio.verbose('mode={} type={} galaxy={}'.format(imod,ityp,igal))
                            # Determine location of plot
                            if verbose: mmlio.verbose('Getting plot location...')
                            irow=varList_plot[rowvar].index(ivar0[rowvar])
                            icol=varList_plot[colvar].index(ivar0[colvar])
                            irect=copy.deepcopy(rect)
                            irect[1]+=-(irow+1)*rect[3]-irow*hpad
                            # Determine line number
                            ilin=0
                            if imodLines: ilin=modLine_plot[imod0_idx].index(imod)
                            if itypLines: ilin=typLine_plot[ityp0_idx].index(ityp)
                            if igalLines: ilin=galLine_plot[igal0_idx].index(igal)
                            if   clrvar=='line': iclr=ilin
                            elif clrvar=='none': iclr=0
                            elif clrvar=='type': iclr=typListDEF.index(ityp)+igal*len(typListDEF)
                            else: iclr=varList_plot[clrvar].index(ivar0[clrvar])
                            if   styvar=='line': isty=ilin
                            elif styvar=='none': isty=0
                            else: isty=varList_plot[styvar].index(ivar0[styvar])
                            # Determine axes
                            imodax=modLine_plot[imod0_idx][0]
                            itypax=typLine_plot[ityp0_idx][0]
                            igalax=galLine_plot[igal0_idx][0]
                            # Create residual if necessary
                            if (flag_resid or flag_compr) and calcdict0.has_key(imod_short):
                                if flag_resid: iopt=resid_hist1D(calcdict[imod_short][ityp][igal],calcdict0[imod_short][ityp][igal])
                                if flag_compr:
                                    iopt0=calcdict0[imod_short][ityp][igal]
                                    iopt=copy.deepcopy(calcdict[imod_short][ityp][igal])
                                    iopt['y']=iopt0['y']
                                    iopt['ymask']=iopt0['ymask']
                                    iopt['pgal']+=1
                            else: iopt=copy.deepcopy(calcdict[imod_short][ityp][igal])
                            # Add correct tag for row
                            if imodLines: iopt['label']=imod
                            if itypLines: iopt['label']=ityp
                            if igalLines: iopt['label']=str(igal)
                            if isinstance(labelList,list): iopt['label']=labelList[ilin]
                            if   rowvar=='type':
                                iopt['ymode']+=' ({})'.format(ityp)
                                iopt['pgal' ]+=4
                            # Add color & style index
                            iopt['clridx']=iclr
                            iopt['styidx']=isty
                            # Plot
                            if verbose: mmlio.verbose('Plotting...')
                            if ilin!=0: plotobj[imod][ityp][igal]['axes']=plotobj[imodax][itypax][igalax]['axes']
                            plotobj[imod][ityp][igal]=plot_hist1D(figObj,iopt,unitdict,rect=irect,objdict=plotobj[imod][ityp][igal])
                            if ilin!=0 and flag_legend and irow==icol==0: plotobj[imod][ityp][igal]['axes'].legend(loc=3)
                            # Add plot object to list
                            if irow==icol==0:
                                plotobj['lineList'].append(plotobj[imod][ityp][igal]['line'])
    # Return plot objects
    return plotobj
        

####################################################################################################################################
# METHOD TO MAKE 2D HISTOGRAM PANEL
def panels_hist2D(figObj,typList,rayList,modList,galList,
                  calcdict={},calcdict0={},unitdict=None,
                  rect=[1.,1.,1.,1.],wpad=0.,hpad=0.,hcb=0.,
                  rowvar=None,colvar=None,clrvar=None,verbose=None):
    """
    Creates 2D histogram panel
    """
    # Set constants
    varLIST=['type','view','mode','galaxy']
    # Pars input
    verbose=mmlpars.mml_pars(verbose,default=False,type=bool)
    rowvar=mmlpars.mml_pars(rowvar,default=varLIST[0],list=varLIST)
    colvar=mmlpars.mml_pars(colvar,default=varLIST[1],list=varLIST)
    clrvar=mmlpars.mml_pars(clrvar,default='type',list=varLIST)
    if verbose: mmlio.verbose('Parsing input...')
    # Initialize plot object dictionary
    plotobj=mmlclass.mydict_nd(None,modList,typList,galList,rayList)
    # Intialize stuff
    plotobj['imagList']=[]
    plotobj['contList']=[]
    modList_plot=[] ; modLine_plot={}
    typList_plot=[] ; typLine_plot={}
    galList_plot=[] ; galLine_plot={}
    rayList_plot=[] ; rayLine_plot={}
    if verbose: mmlio.verbose('Beginning loop over modes, types & views...')
    # Loop over modList
    for imod in modList:
        # Handle multiple mode plots per axes
        imod0_idx=modList.index(imod0)
        if imod0_idx not in modLine_plot: modLine_plot[imod0_idx]=[]
        if isinstance(imod0,list):
            imodList=imod0   ; imodLines=True
        else:
            imodList=[imod0] ; imodLines=False
        for imod in imodList:
            # Handle option for residual/comparison plotting
            flag_resid=False ; flag_compr=False ; imod_short=imod
            if 'resid' in imod.split('_')[-1]:
                flag_resid=True
                imod_short='_'.join(imod.split('_')[:-1])
            if 'compr' in imod.split('_')[-1]:
                flag_compr=True
                imod_short='_'.join(imod.split('_')[:-1])
            # Skip if data dosn't exist for a mode
            if imod_short not in calcdict:
                if verbose: mmlio.verbose('Mode {} not in calcdict. Skipping...'.format(imod_short))
                continue
            if (flag_resid or flag_compr) and imod_short not in calcdict0:
                if verbose: mmlio.verbose('Mode {} not in calcdict0. Skipping...'.format(imod_short))
                continue
            if imod0 not in modList_plot: modList_plot.append(imod0)
            if imod not in modLine_plot[imod0_idx]: modLine_plot[imod0_idx].append(imod)
            # Loop over typList
            for ityp0 in typList:
                # Handle multiple type plots per axes
                ityp0_idx=typList.index(ityp0)
                if ityp0_idx not in typLine_plot: typLine_plot[ityp0_idx]=[]
                if isinstance(ityp0,list):
                    itypList=ityp0   ; itypLines=True
                else:
                    itypList=[ityp0] ; itypLines=False
                for ityp in itypList:
                    # Skip if data dosn't exist for a type
                    if ityp not in calcdict[imod_short]:
                        if verbose: mmlio.verbose('Type {} not in calcdict. Skipping...'.format(ityp))
                        continue
                    if (flag_resid or flag_compr) and ityp not in calcdict0[imod_short]:
                        if verbose: mmlio.verbose('Type {} not in calcdict0. Skipping...'.format(ityp))
                        continue
                    if ityp0 not in typList_plot: typList_plot.append(ityp0)
                    if ityp not in typLine_plot[ityp0_idx]: typLine_plot[ityp0_idx].append(ityp)
                    # Loop over galList
                    for igal0 in galList:
                        # Handle multiple galaxy plots per axes
                        igal0_idx=galList.index(igal0)
                        if igal0_idx not in galLine_plot: galLine_plot[igal0_idx]=[]
                        if isinstance(igal0,list):
                            igalList=igal0   ; igalLines=True
                        else:
                            igalList=[igal0] ; igalLines=False
                        for igal in igalList:
                            # Skip if data dosn't exist for a galaxy
                            if igal not in calcdict[imod_short][ityp]:
                                if verbose: mmlio.verbose('Galaxy {} not in calcdict. Skipping...'.format(igal))
                                continue
                            if (flag_resid or flag_compr) and igal not in calcdict0[imod_short][ityp]:
                                if verbose: mmlio.verbose('Galaxy {} not in calcdict0. Skipping...'.format(igal))
                                continue
                            if igal0 not in galList_plot: galList_plot.append(igal0)
                            if igal not in galLine_plot[igal0_idx]: galLine_plot[igal0_idx].append(igal)
                            # Loop over rayList
                            for iray0 in rayList:
                                # Handle multiple ray plots per axes
                                iray0_idx=rayList.index(iray0)
                                if iray0_idx not in rayLine_plot: rayLine_plot[iray0_idx]=[]
                                if isinstance(iray0,list):
                                    irayList=iray0   ; irayLines=True
                                else:
                                    irayList=[iray0] ; irayLines=False
                                for iray in irayList:
                                    # Skip if data dosn't exist for a ray
                                    if iray not in calcdict[imod_short][ityp][igal]:
                                        if verbose: mmlio.verbose('View {} not in calcdict. Skipping...'.format(iray))
                                        continue
                                    if (flag_resid or flag_compr) and iray not in calcdict0[imod_short][ityp][igal]:
                                        if verbose: mmlio.verbose('View {} not in calcdict0. Skipping...'.format(iray))
                                        continue
                                    if iray0 not in rayList_plot: rayList_plot.append(iray0)
                                    if iray not in rayLine_plot[iray0_idx]: rayLine_plot[iray0_idx].append(iray)
                                    
                                    # Get some info
                                    if imodLines+itypLines+igalLines+irayLines > 1: raise Exception('Only one variable can have multiple plots in one axes.')
                                    varList_plot={'mode':modList_plot,'type':typList_plot,'galaxy':galList_plot,'view':rayList_plot}
                                    ivar ={'mode':imod ,'type':ityp ,'galaxy':igal ,'view':iray }
                                    ivar0={'mode':imod0,'type':ityp0,'galaxy':igal0,'view':iray0}
                                    if verbose: mmlio.verbose('mode={} type={} galaxy={} view={}'.format(imod_short,ityp,igal,iray))
                                    # Determine location of plot
                                    if verbose: mmlio.verbose('Getting plot location...')
                                    irow=varList_plot[rowvar].index(ivar0[rowvar])
                                    icol=varList_plot[colvar].index(ivar0[colvar])
                                    if clrvar=='type':
                                        iclr=simlist.LIST_PTYPBASE.index(ityp)
                                    else:
                                        iclr=varList_plot[clrvar].index(ivar0[clrvar])
                                    irect=copy.deepcopy(rect)
                                    irect[0]+=icol*rect[2]+icol*wpad
                                    irect[1]+=-(irow+1)*rect[3]-irow*hpad
                                    # Determine line number
                                    ilin=0
                                    if imodLines: ilin=modLine_plot[imod0_idx].index(imod)
                                    if itypLines: ilin=typLine_plot[ityp0_idx].index(ityp)
                                    if igalLines: ilin=galLine_plot[igal0_idx].index(igal)
                                    if irayLines: ilin=rayLine_plot[ityp0_idx].index(iray)
                                    # Determine axes
                                    imodax=modLine_plot[imod0_idx][0]
                                    itypax=typLine_plot[ityp0_idx][0]
                                    igalax=galLine_plot[igal0_idx][0]
                                    irayax=rayLine_plot[iray0_idx][0]
                                    # Create residual if necessary
                                    if flag_resid: iopt=resid_hist2D(calcdict[imod_short][ityp][igal][iray],calcdict0[imod_short][ityp][igal][iray])
                                    else         : iopt=copy.deepcopy(calcdict[imod_short][ityp][igal][iray])
                                    # Set options
                                    iopt_imag=None ; iopt_cont=None
                                    if ilin==0: iopt_imag=iopt
                                    else      : iopt_cont=iopt
                                    # Set default axes
                                    if ilin!=0: plotobj[imod][ityp][igal][iray]['axes']=plotobj[imodax][itypax][igalax][irayax]['axes']
                                    # Plot
                                    plotobj[imod][ityp][igal][iray]=plot_hist2D(figObj,iopt_imag,iopt_cont,unitdict,irect,
                                                                               firstrow=(irow==0),firstcol=(icol==0),
                                                                               objdict=plotobj[imod][ityp][igal][iray])
                                    # Add plot objects to list
                                    if icol==0:
                                        plotobj['imagList'].append(plotobj[imod][ityp][igal][iray]['imag_smap'])
                                        plotobj['contList'].append(plotobj[imod][ityp][igal][iray]['cont_smap'])
    # Tally number of rows and columns
    nrow=len(varList_plot[rowvar])
    ncol=len(varList_plot[colvar])
    rect_end=copy.deepcopy(rect)
    rect_end[0]+=(+ncol*rect[2]+(ncol-1)*wpad)
    rect_end[1]+=(-nrow*rect[3]-(nrow-1)*hpad)
    # Return plot objects
    return plotobj
            
####################################################################################################################################
####################################################################################################################################
# METHOD TO RETURN RESIDUAL INFO

####################################################################################################################################
# METHOD TO RETURN PARTICLE RESIDUAL INFO
def resid_scatter(inopt1,inopt2):
    """
    Creates residual object for pNbody
    """
    # Check that lengths are equal
    if inopt1.nbody != inopt2.nbody:
        raise Exception('Nbody residual requires the same number of particles in each pNbody object.')
    # Initialize output dictionary
    outopt=copy.deepcopy(inopt1)
    # Get variable list
    varlist=inopt1.get_list_of_array()
    # Do residual of each
    for ivar in varlist:
        iarr1=getattr(inopt1,ivar)
        iarr2=getattr(inopt2,ivar)
        setattr(outopt,ivar,(iarr1-iarr2)/iarr1)
    # Return residual dictionary
    return outopt

####################################################################################################################################
# METHOD TO RETURN 1D RESIDUAL INFO
def resid_hist1D(inopt1,inopt2):
    """
    Creates options for 1D histogram residual
    """
    # Initialize output dictionary
    outopt=copy.deepcopy(inopt1)
    # Adjust general info
    outopt['ymode' ]+=' residual'
    outopt['yscale'] = 'linear'
    outopt['ylim'  ] = (0.,1.)
    # Get data from primary options dict
    x1=inopt1['x'][:-1]+(inopt1['x'][1:]-inopt1['x'][:-1])/2.
    y1=inopt1['y']
    # Get data from secondary options dict
    x2=inopt2['x'][:-1]+(inopt2['x'][1:]-inopt2['x'][:-1])/2.
    y2=inopt2['y']
    # Interpolate
    res=np.abs(np.ma.masked_invalid((y1-interpolate.griddata(x2,y2,x1))/y1))
    outopt['y'    ]=res.data
    outopt['ymask']=res.mask
    # Return residual dictionary
    return outopt

####################################################################################################################################
# METHOD TO RETURN 2D RESIDUAL INFO
def resid_hist2D(inopt1,inopt2):
    """
    Creates options for 2D histogram residual
    """
    verbose=False
    # Initialize output dictionary
    outopt=copy.deepcopy(inopt1)
    # Addjust general info
#    outopt['zmode' ]+= ' residual'
    outopt['zscale'] = 'linear'
    outopt['zlim'  ] = (0.,1.)
    outopt['zmean' ] = 0.5
    # Get data from primary options dict
    if verbose:
        print 'opt 1'
        print 'xlim={}'.format(inopt1['xlim'])
        print 'ylim={}'.format(inopt1['ylim'])
        print 'zlim={}'.format(inopt1['zlim'])
    x1=inopt1['x'][:-1]+(inopt1['x'][1:]-inopt1['x'][:-1])/2.
    y1=inopt1['y'][:-1]+(inopt1['y'][1:]-inopt1['y'][:-1])/2.
    dim1=(len(x1),len(y1)) ; len1=dim1[0]*dim1[1]
    yy1,xx1=np.meshgrid(y1,x1)
#    xx1,yy1=np.meshgrid(x1,y1)
    xx1=np.reshape(xx1,len1)
    yy1=np.reshape(yy1,len1)
    z1=np.reshape(inopt1['z'],len1)
    # Get data from secondary options dict
    if verbose:
        print 'opt 2'
        print 'xlim={}'.format(inopt2['xlim'])
        print 'ylim={}'.format(inopt2['ylim'])
        print 'zlim={}'.format(inopt2['zlim'])
    x2=inopt2['x'][:-1]+(inopt2['x'][1:]-inopt2['x'][:-1])/2.
    y2=inopt2['y'][:-1]+(inopt2['y'][1:]-inopt2['y'][:-1])/2.
    dim2=(len(x2),len(y2)) ; len2=dim2[0]*dim2[1]
    yy2,xx2=np.meshgrid(y2,x2)
#    xx2,yy2=np.meshgrid(x2,y2)
    xx2=np.reshape(xx2,len2)
    yy2=np.reshape(yy2,len2)
    z2=np.reshape(inopt2['z'],len2)
    # Interpolate
    res=np.abs(np.ma.masked_invalid((z1-interpolate.griddata((xx2,yy2),z2,(xx1,yy1)))/z1))
    if verbose:
        print 'res lim={}'.format([res.min(),np.mean(res),res.max()])
    outopt['z'    ]=np.reshape(res.data,dim1)
    outopt['zmask']=np.reshape(res.mask,dim1)
    # Return residual dictionary
    return outopt

####################################################################################################################################
####################################################################################################################################
# METHODS TO PLOT HISTOGRAMS

####################################################################################################################################
# METHOD TO PLOT 1D HISTOGRAMS
def plot_hist1D(figObj,opt=None,unitdict=None,rect=[0.,0.,1.,1.],objdict=None):
    '''
    Plots 1D histograms
    '''
    # Set constants
    plot_kw=['axes','line','line_cmap']
    unitlist=['UnitDensity','UnitEnergy','UnitLength','UnitMass',
              'UnitSpecEnergy','UnitSurfaceDensity','UnitTime','UnitVelocity']
    # Pars input
    if not isinstance(unitdict,dict): unitdict=dict.fromkeys(unitlist,None)
    # Object dictionary
    if not isinstance(objdict,dict): objdict=dict.fromkeys(plot_kw,None)
    for iplot_kw in plot_kw:
        if not objdict.has_key(iplot_kw): objdict[iplot_kw]=None
    # Axes
    if objdict['axes']==None:
        objdict['axes']=figObj.add_axes(list(rect))
        objdict['axes'].set_xlim(opt['xlim'])
        objdict['axes'].set_ylim(opt['ylim'])
        objdict['axes'].set_xscale(opt['xscale'])
        objdict['axes'].set_yscale(opt['yscale'])
        objdict['axes'].autoscale(enable=False)
        objdict['axes'].set_xlabel(opt['xmode'])
        objdict['axes'].set_ylabel(opt['ymode'])
    # Colormap
    if objdict['line_cmap']==None: objdict['line_cmap']=mmlplot.get_cmap(list_cmaps()[opt['clridx']])
    color=list_colors(clrsim=opt['clrsim'])[opt['clridx']]
    # Line styles
    style=list_linestyles()[opt['styidx']]
    # Eensure that dimensions match
    if len(opt['x']) == len(opt['y']):
        x=opt['x']
        y=opt['y']
    else:
        x=opt['x'][1:]-(opt['x'][1:]-opt['x'][:-1])/2.
        y=opt['y']
    # Apply mask
    if 'xmask' in opt: x=np.ma.masked_array(x,mask=opt['xmask'])
    if 'ymask' in opt: y=np.ma.masked_array(y,mask=opt['ymask'])
    # Plot
    if objdict['line']==None: objdict['line']=objdict['axes'].plot(x,y,color=color,ls=style,label=opt['label'])
    else:                     objdict['line'].set_data(x,y)
    # Return output
    return objdict

####################################################################################################################################
# METHOD TO PLOT 2D HISTOGRAMS
def plot_hist2D(figObj,imagopt=None,contopt=None,unitdict=None,rect=[0.,0.,1.,1.],
                firstrow=False,firstcol=False,objdict=None,skipimag=None,skipcont=None):
    '''
    Plots 2D histograms
    '''
    # Set constants
    plot_kw=['axes',
             'imag','imag_cmap','imag_norm','imag_cbax','imag_cbar','imag_smap',
             'cont','cont_cmap','cont_norm','cont_cbax','cont_cbar','cont_smap']
    cmapdict=dict(m='DEF',vr='matplotlib_jet',pot='matplotlib_PuOr',scfpotm2='matplotlib_PuOr',scfpot='matplotlib_PuOr',nan='DEF')
    cmapdict_dir=dict(m=False,vr=False,pot=True,scfpotm2=True,scfpot=True,nan=False)
    ncont=10
    unitlist=['UnitDensity','UnitEnergy','UnitLength','UnitMass',
              'UnitSpecEnergy','UnitSurfaceDensity','UnitTime','UnitVelocity']
    # Pars input
    if not isinstance(unitdict,dict): unitdict=dict.fromkeys(unitlist,None)
    skipimag=mmlpars.mml_pars(skipimag,default=False,type=bool)
    skipcont=mmlpars.mml_pars(skipcont,default=False,type=bool)
    # Determine what to plot
    flag_imag=False
    if not skipimag and imagopt!=None:
        if imagopt['zmode']!='nan': flag_imag=True
    flag_cont=False
    if not skipcont and contopt!=None:
        if contopt['zmode']!='nan': flag_cont=True
    # Raise error if neither
    if not flag_imag and not flag_cont:
        raise Exception('Either image or contour data must be provided.')
    # Set axes option dict
    if flag_imag: opt=imagopt
    else:         opt=contopt
    # Check that imag and cont are in same space if both provided
    if flag_imag and flag_cont:
#        if imagopt['ptyp'  ]!=contopt['ptyp'  ]: raise Exception('Image and contour particle types do not match.')
#        if imagopt['pgal'  ]!=contopt['pgal'  ]: raise Exception('Image and contour particle galaxies do not match.')
        if imagopt['view'  ]!=contopt['view'  ]: raise Exception('Image and contour views do not match.')
#        if imagopt['xymode']!=contopt['xymode']: raise Exception('Image and contour xymodes do not match.')
    # Object dictionary
    if objdict==None: objdict=dict.fromkeys(plot_kw,None)
    # Axes
    if objdict['axes']==None:
        objdict['axes']=figObj.add_axes(list(rect))
        objdict['axes'].set_xlim(opt['xlim'])
        objdict['axes'].set_ylim(opt['ylim'])
        objdict['axes'].set_xscale(opt['xscale'])
        objdict['axes'].set_yscale(opt['yscale'])
        objdict['axes'].autoscale(enable=False)
        objdict['axes'].tick_params(labelbottom=False,labeltop=False,labelright=False,labelleft=False)
        if opt['xymode']=='pos':
            objdict['axes'].set_xlabel(opt['view'][0])
            objdict['axes'].set_ylabel(opt['view'][1])
        elif opt['xymode']=='vel':
            objdict['axes'].set_xlabel('v'+opt['view'][0])
            objdict['axes'].set_ylabel('v'+opt['view'][1])
        else: raise Exception('Unsupported xymode of {}'.format(opt['xymode']))
        if firstrow and firstcol:
            objdict['axes'].set_xlabel('')
            if   opt['xymode']=='pos': objdict['axes'].set_ylabel(unitdict['UnitLength'].symbol)
            elif opt['xymode']=='vel': objdict['axes'].set_ylabel(unitdict['UnitVelocity'].symbol)
        else:
            objdict['axes'].set_xlabel('')
            objdict['axes'].set_ylabel('')
        if firstcol:
            objdict['axes'].tick_params(labelleft=True)
            objdict['axes'].text(0.1,0.9,opt['ptyp'].title(),ha='left',va='top',transform=objdict['axes'].transAxes)
        if firstrow:
            if   isinstance(opt['view'],str  ): iraytitle=opt['view']
            elif isinstance(opt['view'],tuple): iraytitle='({},{})'.format(opt['view'][0],opt['view'][1])
            else: raise Exception('Unsuported view type of {}'.format(type(opt['view'])))
            objdict['axes'].text(0.9,0.9,iraytitle,ha='right',va='top',transform=objdict['axes'].transAxes)
    # Colormap
    if flag_imag:
        if objdict['imag_cmap']==None:
            if cmapdict[imagopt['zmode']]=='DEF':
                objdict['imag_cmap']=mmlplot.get_cmap(list_cmaps()[imagopt['clridx']],reverse=cmapdict_dir[imagopt['zmode']])
            else:
                objdict['imag_cmap']=mmlplot.get_cmap(cmapdict[imagopt['zmode']])
##         if objdict['imag_cmap']==None: objdict['imag_cmap']=get_cmaps(cmapdict[imagopt['zmode']],imagopt['ptyp'],imagopt['pgal'],
##                                                                       reverse=cmapdict_dir[imagopt['zmode']])
    if flag_cont:
        if objdict['cont_cmap']==None:
            if cmapdict[contopt['zmode']]=='DEF':
                objdict['cont_cmap']=mmlplot.get_cmap(list_cmaps()[contopt['clridx']],reverse=cmapdict_dir[contopt['zmode']])
            else:
                objdict['cont_cmap']=mmlplot.get_cmap(cmapdict[contopt['zmode']])
##         if objdict['cont_cmap']==None: objdict['cont_cmap']=get_cmaps(cmapdict[contopt['zmode']],contopt['ptyp'],contopt['pgal'],
##                                                                       reverse=cmapdict_dir[contopt['zmode']])
    # Normalize object
    if flag_imag:
        if objdict['imag_norm']==None:
            objdict['imag_norm']=matplotlib.colors.Normalize(vmin=0,vmax=255,clip=True)
##             if imagopt['zscale']=='linear':
##                 objdict['imag_norm']=matplotlib.colors.Normalize(vmin=imagopt['zlim'][0],vmax=imagopt['zlim'][1],clip=True)
##             elif imagopt['zscale']=='log':
##                 if imagopt['zlim'][0] <= 0 or imagopt['zlim'][1] <= 0:
##                     objdict['imag_norm']=matplotlib.colors.Normalize(vmin=imagopt['zlim'][0],vmax=imagopt['zlim'][1],clip=True)
##                 else:
##                     objdict['imag_norm']=matplotlib.colors.LogNorm(vmin=imagopt['zlim'][0],vmax=imagopt['zlim'][1],clip=True)
##             else:
##                 raise Exception('Unsupported zscale option for image ({}).'.format(imagopt['zscale']))
    if flag_cont:
        if objdict['cont_norm']==None:
            objdict['cont_norm']=matplotlib.colors.Normalize(vmin=0,vmax=255,clip=False)
## #            contopt['zscale']='linear'
##             if contopt['zscale']=='linear':
##                 objdict['cont_norm']=matplotlib.colors.Normalize(vmin=contopt['zlim'][0],vmax=contopt['zlim'][1],clip=False)
##             elif contopt['zscale']=='log':
##                 if contopt['zlim'][0] <= 0 or contopt['zlim'][1] <= 0:
##                     objdict['cont_norm']=matplotlib.colors.Normalize(vmin=contopt['zlim'][0],vmax=contopt['zlim'][1],clip=False)
##                 else:
##                     objdict['cont_norm']=matplotlib.colors.LogNorm(vmin=contopt['zlim'][0],vmax=contopt['zlim'][1],clip=False)
##             else:
##                 raise Exception('Unsupported zscale option for contour ({}).'.format(contopt['zscale']))
    # Scalar mappable object
#    if flag_imag:
#        if objdict['imag_smap']==None: objdict['imag_smap']=matplotlib.cm.ScalarMappable(norm=objdict['imag_norm'],cmap=objdict['imag_cmap'])
#    if flag_cont:
#        if objdict['cont_smap']==None: objdict['cont_smap']=matplotlib.cm.ScalarMappable(norm=objdict['cont_norm'],cmap=objdict['cont_cmap'])
    # Plot (ensure that dimensions match)
    if flag_imag:
        if imagopt['zlim'][0]==imagopt['zlim'][1]:
            mmlio.verbose('Image data is singular. Skipping plot.')
        else:
            # Select data and add mask
            ximag=imagopt['x']
            yimag=imagopt['y']
            zimag=imagopt['z']
            if 'xmask' in imagopt: ximag=np.ma.masked_array(ximag,mask=imagopt['xmask'])
            if 'ymask' in imagopt: yimag=np.ma.masked_array(yimag,mask=imagopt['ymask'])
            if 'zmask' in imagopt: zimag=np.ma.masked_array(zimag,mask=imagopt['zmask'])
            # Scale image data
            zimag=map_array(zimag,scale=imagopt['zscale'],cd=imagopt['zmean'],
                            mn=imagopt['zlim'][0],mx=imagopt['zlim'][1])
            # Set image extent based on x & y data
            extent=(ximag.min(),ximag.max(),yimag.min(),yimag.max())
            # Create object/adjust data
            if objdict['imag']==None:
                objdict['imag']=objdict['axes'].imshow(np.transpose(zimag),aspect='auto',extent=extent,interpolation='bicubic',
                                                       cmap=objdict['imag_cmap'],norm=objdict['imag_norm'],origin='lower')
            else:
                objdict['imag'].set_data(np.transpose(zimag))
                objdict['imag'].set_extent(extent)
            if objdict['imag_smap']==None: objdict['imag_smap']=objdict['imag']
    if flag_cont:
        if contopt['zlim'][0]==contopt['zlim'][1]:
            mmlio.verbose('contour data is singular. Skipping plot.')
        else:
#            if contopt['zscale']=='log': zcont=np.ma.masked_less_equal(contopt['z'],0.)
#            else:                        zcont=contopt['z']
            xcont=contopt['x']
            ycont=contopt['y']
            zcont=contopt['z']
            if 'xmask' in contopt: xcont=np.ma.masked_array(xcont,mask=contopt['xmask'])
            if 'ymask' in contopt: ycont=np.ma.masked_array(ycont,mask=contopt['ymask'])
            if 'zmask' in contopt: zcont=np.ma.masked_array(zcont,mask=contopt['zmask'])
            # Scale contour data
            zcont=map_array(zcont,scale=contopt['zscale'],cd=contopt['zmean'],
                            mn=contopt['zlim'][0],mx=contopt['zlim'][1])
            # Set contour extent based on x & y data
            extent=(xcont.min(),xcont.max(),ycont.min(),ycont.max())
            # Create contour object/adjust data
            if objdict['cont']==None:
                contLvls=np.linspace(0,255,ncont)
##                 if   contopt['zscale']=='linear': contLvls=np.linspace(contopt['zlim'][0],contopt['zlim'][1],ncont)
##                 elif contopt['zscale']=='log'   : contLvls=mmlmath.logspace(contopt['zlim'][0],contopt['zlim'][1],ncont,addzero=True)
##                 contLvls=map_array(zcont,scale=contopt['zscale'],cd=contopt['zcd'],
##                                    mn=contopt['zlim'][0],mx=contopt['zlim'][1])
                objdict['cont']=objdict['axes'].contour(np.transpose(zcont),extent=extent,levels=contLvls,
                                                        cmap=objdict['cont_cmap'],norm=objdict['cont_norm'],origin='lower')
            else:
                objdict['cont'].set_data(zcont)
                objdict['cont'].set_extent(extent)
            # Create scalar mappable for contours
            if objdict['cont_smap']==None:
                objdict['cont_smap']=matplotlib.cm.ScalarMappable(norm=objdict['cont_norm'],cmap=objdict['cont_cmap'])
                objdict['cont_smap']._A=contLvls
#    print imagopt['xymode']+'_'+imagopt['zmode']
#    print 'zlim={}'.format(imagopt['zlim'])
#    print 'dlim={}'.format((imagopt['z'].min(),imagopt['z'].max()))
    # Return output
    return objdict

####################################################################################################################################
####################################################################################################################################
# METHODS TO PARS HISTOGRAM OPTIONS

####################################################################################################################################
# METHOD TO COUNT DIMENSIONS
def dict_mode1d():
    mdict={
        'rxy'   : ('rxy'   ,'z'    )
        }
def dict_mode2d():
    mdict={
        'xy'    : ('x'     ,'y'    ,'z'  ),
        'xz'    : ('x'     ,'z'    ,'y'  ),
        'yz'    : ('y'     ,'z'    ,'x'  ),
        'Rz'    : ('rxy'   ,'z'    ),
        'vxy'   : ('vx'    ,'vy'   ),
        'vxz'   : ('vx'    ,'vz'   ),
        'vyz'   : ('vy'    ,'vz'   ),
        'scfxy' : ('scfx'  ,'scfy' ),
        'scfxz' : ('scfx'  ,'scfz' ),
        'scfyz' : ('scfy'  ,'scfz' ),
        'scfvxy': ('scfvx' ,'scfvy'),
        'scfvxz': ('scfvx' ,'scfvz'),
        'scfvyz': ('scfvy' ,'scfvz'),
        'scfRz' : ('scfrxy','scfz' )
        }
    return mdict
def dict_mode3d():
    mdict={
        'xyz'      : ('x'      ,'y'     ,'z'       ),
        'pos'      : ('x'      ,'y'     ,'z'       ),
        'cylpos'   : ('rxy'    ,'phi'   ,'z'       ),
        'sphpos'   : ('rxyz'   ,'phi'   ,'theta'   ),
        'vxyz'     : ('vx'     ,'vy'    ,'vz'      ),
        'vel'      : ('vx'     ,'vy'    ,'vz'      ),
        'scfxyz'   : ('scfx'   ,'scfy'  ,'scfz'    ),
        'scfpos'   : ('scfx'   ,'scfy'  ,'scfz'    ),
        'scfcylpos': ('scfrxy' ,'scfphi','scfz'    ),
        'scfsphpos': ('scfrxyz','scfphi','scftheta'),
        'scfvxyz'  : ('scfvx'  ,'scfvy' ,'scfvz'   ),
        'scfvel'   : ('scfvx'  ,'scfvy' ,'scfvz'   )
        }
def list_mode2d():
    return dict_mode2d.keys()
def list_mode3d():
    return dict_mode3d.keys()
def list_modedim(mode):
    if   mode in dict_mode1d(): dims=dict_mode1d()[mode]
    elif mode in dict_mode2d(): dims=dict_mode2d()[mode]
    elif mode in dict_mode3d(): dims=dict_mode3d()[mode]
    else: dims=(mode)
    return dims
def nmodedim(mode):
    if   mode in dict_mode1d(): ndim=1 
    elif mode in dict_mode2d(): ndim=2
    elif mode in dict_mode3d(): ndim=3
    else: ndim=1
    return ndim

####################################################################################################################################
# METHOD TO COMPUTE ND HISTOGRAM
def histnd(nbobj,outdict=None,totidx=None,limdict=None):
    """
    Returns ND histogram
    """
    # Pars input
    opt=parsopt_histND(hist_kw)
    outdict=mmlpars.mml_pars(outdict,default=False,type=bool)
    ndat=len(nbobj.pos[:,0])
    # Get space modes
    smodes=list_modedim(opt['smode']) ; ndim=nmodedim(opt['smode']) ; dlist=[]
    for ismod in smodes: dlist.append(nbobj.get_var(ismod))
    d=np.vstack(dlist).T
    # Get weight mode
    w=nbobj.get_var(opt['wmode'])
    # Select correct particles
    if totidx==None: totidx=nbobj.get_index(opt['ptyp'],opt['pgal'])
    d=d[totidx,:]
    w=w[totidx,:]
    # Set limits
    for idim in range(len(smodes)):
        if opt['slim'][idim][0]==opt['slim'][idim][1]:
            opt['slim'][idim]=nbobj.get_limits(smodes[idim],opt['slimpad'][idim],opt['slimmir'][idim],totindex=totidx,lim=limdict)
    if opt['wlim'][0]==opt['wlim'][1]:
        opt['wlim']=nbobj.get_limits(opt['wmode'],opt['wlimpad'],opt['wlimmir'],totindex=totidx,lim=limdict)
    # Create bins
    bins=[]
    for idim in range(len(smodes)):
        if   opt['sscale'][idim]=='linear': bins.append(np.linspace(opt['slim'][idim][0],opt['slim'][idim][1],opt['nbins'][idim]))
        elif opt['sscale'][idim]=='log'   : bins.append(mmlmath.logspace(opt['slim'][idim][0],opt['slim'][idim][1],opt['nbins'][idim],addzero=True))
        else: raise Exception('Unsupported scale for mode {}: {}'.format(smodes[idim],opt['sscale'][idim]))
    # Create histogram
    hist,bins=np.histogramdd(d,weights=w,bins=bins)
    if opt['wmode']=='n' or opt['wmode']=='m':
        nhist=copy.deepcopy(hist)
    else:
        nhist,nbins=np.histogramdd(d,bins=bins)
    # Mask empty bins
    hist =np.ma.masked_where(nhist==0,hist )
    nhist=np.ma.masked_where(nhist==0,nhist)
    # Adjust histogram based on method
    if   opt['method']=='hist': pass
    elif opt['method']=='mean': hist/=nhist
    elif opt['method']=='dens':
        binsmean=[]
        for ibin in bins: binsmean.append(ibin[:-1]+np.diff(ibin)/2.)
        if   opt['smode']=='rxyz'  : # Spherical shells
            vols=(4.0/3.0)*math.pi*(np.power(bins[0][1:],3)-np.power(bins[0][:-1],3))
        elif opt['smode']=='rxy'   : # Circular rings
            vols=math.pi*(np.power(bins[0][1:],2)-np.power(bins[0][:-1],2))*(bins[1][1:]-bins[1][:-1])
        elif opt['smode']=='Rz'    : # Circular rings
            vols=math.pi*(np.power(bins[0][1:],2)-np.power(bins[0][:-1],2))*(bins[1][1:]-bins[1][:-1])
        elif opt['smode']=='cylpos': # Cylindrical volume elements
            pass
        elif opt['smode']=='sphpos': # Spherical volume elements
            pass
        else:                        # Linear density
            vols=1.
            for ibin in bins: vols*=(ibin[1:]-ibin[:-1])
        raise Exception('Density histogram not currently supported.')
    # Set limits if singular
    for idim in range(len(smodes)):
        if opt['slim'][idim][0]==opt['slim'][idim][1]: opt['slim'][idim]=(bins[idim].min(),bins[idim].max())
    if opt['wlim'][0]==opt['wlim'][1]: opt['wlim']=(hist.data.flatten()[hist.argmin()],hist.data.flatten()[hist.argmax()])
    # Set mean
    for idim in range(ndim):
        opt['smean'][idim]=bins[idim].mean()
    if np.all(hist.mask): opt['wmean']=hist.data.mean()
    else                : opt['wmean']=hist.mean()
    # Add output to dictionary
    opt['bins']=bins
    opt['hist']=hist.data ; opt['mask']=hist.mask
    # Return output
    if outdict: return opt
    else      : return hist,bins
    
####################################################################################################################################
# METHOD TO RETURN DEFAULT HISTOGRAM OPTIONS FOR DIFFERENT VARIABLES                                                                                                               
def get_histdef(mode):
    # Create category lists
    mirrlist=['phi','x','y','z','pos','vx','vy','vz','vr','vt','vcir','vrad','pos','vel','pot','scfpot','scfpotm2','scfpotm2rel']
    loglist=['n','m','rho','Rz','rxyz','pot','scfpot','scfpotm2','scfpotm2rel']
    linlist=['pos','vel','phi']
    histlist=['n','m']
    denslist=['rho']
    # Pars mode for scale
    if mode in loglist:    scale='log'
    else:                  scale='linear'
    # Pars mode for method
    if mode in histlist:   method='hist'
    elif mode in denslist: method='dens'
    else:                  method='mean'
    # Pars mode for mirror
    if mode in mirrlist:   mirror=True
    else:                  mirror=False
    # Return output
    return method,scale,mirror

####################################################################################################################################
# METHOD TO PARS ND HISTOGRAM OPTIONS
def parsopt_histND(inopt=None,ndim=3):
    inopt=mmlpars.mml_pars(inopt,default={},type=dict)
    optform={
      'ptyp'     : mmlpars.parsdict(default='all'         ,type=str                ),
      'pgal'     : mmlpars.parsdict(default=0             ,type=int, min=0         ),
      'method'   : mmlpars.parsdict(default='DEF'         ,type=str                ),
      'nbins'    : mmlpars.parsdict(default=ndim*[512]    ,type=list               ),
      'smode'    : mmlpars.parsdict(default='xyz'         ,type=str                ),
      'wmode'    : mmlpars.parsdict(default='n'           ,type=str                ),
      'slim'     : mmlpars.parsdict(default=ndim*[[0.,0.]],type=list               ),
      'wlim'     : mmlpars.parsdict(default=[0.,0.]       ,type=list, nelements=2  ),
      'slimpad'  : mmlpars.parsdict(default=ndim*[0.]     ,type=list               ),
      'wlimpad'  : mmlpars.parsdict(default=0.0           ,type=float, min=0.      ),
      'slimmir'  : mmlpars.parsdict(default=ndim*['DEF']  ,type=list               ),
      'wlimmir'  : mmlpars.parsdict(default='DEF'                                  ),
      'sscale'   : mmlpars.parsdict(default=ndim*['DEF']  ,type=list               ),
      'wscale'   : mmlpars.parsdict(default='DEF'         ,type=str                ),
      'smean'    : mmlpars.parsdict(default=ndim*[0.]     ,type=list               ),
      'wmean'    : mmlpars.parsdict(default=0.0           ,type=float              ),
      'package'  : mmlpars.parsdict(default='numpy'       ,list=['numpy','pNbody'] ),
      'clridx'   : mmlpars.parsdict(default=0             ,type=int,min=0          ),
      'styidx'   : mmlpars.parsdict(default=0             ,type=int,min=0          ),
      'label'    : mmlpars.parsdict(default=''            ,type=str                )
      }
    outopt=mmlpars.mml_formpars(inopt,optform)
    # Check dimensions
    ndim=nmodedim(outopt['smode'])
    if nmodedim(outopt['wmode']) != 1: raise Exception('Mode weighting histogram must be 1D.')
    smodes=list_modedim(outopt['smode'])
    # Set space defaults
    for idim in range(ndim):
        ismethDEF,issclDEF,ismirDEF=get_histdef(smodes[idim])
        if outopt['sscale' ][idim]=='DEF': outopt['sscale' ][idim]=issclDEF
        if outopt['slimmir'][idim]=='DEF': outopt['slimmir'][idim]=ismirDEF
    # Set mode defaults
    mmethodDEF,mscaleDEF,mmirrorDEF=get_histdef(outopt['wmode'])
    if outopt['wscale' ]=='DEF': outopt['wscale' ]=mscaleDEF
    if outopt['method' ]=='DEF': outopt['method' ]=mmethodDEF
    if outopt['wlimmir']=='DEF': outopt['wlimmir']=mmirrorDEF
    return outopt

####################################################################################################################################
# METHOD TO PARS 1D HISTOGRAM OPTIONS
def parsopt_hist1D(inopt=None):
    inopt=mmlpars.mml_pars(inopt,default={},type=dict)
    optform={
      'clrsim'   : mmlpars.parsdict(default=False   ,type=bool               ),
      'clridx'   : mmlpars.parsdict(default=0       ,type=int,min=0          ),
      'styidx'   : mmlpars.parsdict(default=0       ,type=int,min=0          ),
      'ptyp'     : mmlpars.parsdict(default='all'   ,type=str                ),
      'pgal'     : mmlpars.parsdict(default=0       ,type=int, min=0         ),
      'view'     : mmlpars.parsdict(default='xz'    ,type=str                ),
      'method'   : mmlpars.parsdict(default='DEF'   ,type=str                ),
      'xnbins'   : mmlpars.parsdict(default=512     ,type=int,min=1          ),
      'xmode'    : mmlpars.parsdict(default='rxyz'  ,type=str                ),
      'ymode'    : mmlpars.parsdict(default='n'     ,type=str                ),
      'xlim'     : mmlpars.parsdict(default=(0.,0.) ,type=tuple,nelements=2  ),
      'ylim'     : mmlpars.parsdict(default=(0.,0.) ,type=tuple,nelements=2  ),
      'xlimpad'  : mmlpars.parsdict(default=0.0     ,type=float,min=0.       ),
      'ylimpad'  : mmlpars.parsdict(default=0.0     ,type=float,min=0.       ),
      'xlimmir'  : mmlpars.parsdict(default='DEF'                            ), 
      'ylimmir'  : mmlpars.parsdict(default='DEF'                            ), 
      'xscale'   : mmlpars.parsdict(default='DEF'   ,type=str                ),
      'yscale'   : mmlpars.parsdict(default='DEF'   ,type=str                ),
      'xmean'    : mmlpars.parsdict(default=0.0     ,type=float              ),
      'ymean'    : mmlpars.parsdict(default=0.0     ,type=float              ),
      'label'    : mmlpars.parsdict(default=''      ,type=str                )
      }
    outopt=mmlpars.mml_formpars(inopt,optform)
    xmethodDEF,xscaleDEF,xmirrorDEF=get_histdef(outopt['xmode'])
    ymethodDEF,yscaleDEF,ymirrorDEF=get_histdef(outopt['ymode'])
    if outopt['xscale']=='DEF': outopt['xscale']=xscaleDEF
    if outopt['yscale']=='DEF': outopt['yscale']=yscaleDEF
    if outopt['method']=='DEF': outopt['method']=ymethodDEF
    if not isinstance(outopt['xlimmir'],bool): outopt['xlimmir']=xmirrorDEF
    if not isinstance(outopt['ylimmir'],bool): outopt['ylimmir']=ymirrorDEF
    return outopt

####################################################################################################################################
# METHOD TO PARS 2D HISTOGRAM OPTIONS
def parsopt_hist2D(inopt=None):
    inopt=mmlpars.mml_pars(inopt,default={},type=dict)
    optform={
      'package'  : mmlpars.parsdict(default='numpy' ,list=['numpy','pNbody'] ),
      'clridx'   : mmlpars.parsdict(default=0       ,type=int,min=0          ),
      'styidx'   : mmlpars.parsdict(default=0       ,type=int,min=0          ),
      'ptyp'     : mmlpars.parsdict(default='all'   ,type=str                ),
      'pgal'     : mmlpars.parsdict(default=0       ,type=int, min=0         ),
      'view'     : mmlpars.parsdict(default='xy'    ,type=str                ),
      'method'   : mmlpars.parsdict(default='DEF'   ,type=str                ),
      'xnbins'   : mmlpars.parsdict(default=512     ,type=int,min=1          ),
      'ynbins'   : mmlpars.parsdict(default=512     ,type=int,min=1          ),
      'xymode'   : mmlpars.parsdict(default='pos'   ,type=str                ),
      'zmode'    : mmlpars.parsdict(default='n'     ,type=str                ),
      'xlim'     : mmlpars.parsdict(default=(0.,0.) ,type=tuple,nelements=2  ),
      'ylim'     : mmlpars.parsdict(default=(0.,0.) ,type=tuple,nelements=2  ),
      'zlim'     : mmlpars.parsdict(default=(0.,0.) ,type=tuple,nelements=2  ),
      'xlimpad'  : mmlpars.parsdict(default=0.0     ,type=float,min=0.       ),
      'ylimpad'  : mmlpars.parsdict(default=0.0     ,type=float,min=0.       ),
      'zlimpad'  : mmlpars.parsdict(default=0.0     ,type=float,min=0.       ),
      'xlimmir'  : mmlpars.parsdict(default='DEF'                            ), 
      'ylimmir'  : mmlpars.parsdict(default='DEF'                            ), 
      'zlimmir'  : mmlpars.parsdict(default='DEF'                            ), 
      'xscale'   : mmlpars.parsdict(default='DEF'   ,type=str                ),
      'yscale'   : mmlpars.parsdict(default='DEF'   ,type=str                ),
      'zscale'   : mmlpars.parsdict(default='DEF'   ,type=str                ),
      'xmean'    : mmlpars.parsdict(default=0.0     ,type=float              ),
      'ymean'    : mmlpars.parsdict(default=0.0     ,type=float              ),
      'zmean'    : mmlpars.parsdict(default=0.0     ,type=float              ),
      'label'    : mmlpars.parsdict(default=''      ,type=str                )
      }
    outopt=mmlpars.mml_formpars(inopt,optform)
    xymethodDEF,xyscaleDEF,xymirrorDEF=get_histdef(outopt['xymode'])
    zmethodDEF,zscaleDEF,zmirrorDEF=get_histdef(outopt['zmode'])
    if outopt['xscale']=='DEF': outopt['xscale']=xyscaleDEF
    if outopt['yscale']=='DEF': outopt['yscale']=xyscaleDEF
    if outopt['zscale']=='DEF': outopt['zscale']=zscaleDEF
    if outopt['method']=='DEF': outopt['method']=zmethodDEF
    if not isinstance(outopt['xlimmir'],bool): outopt['xlimmir']=xymirrorDEF
    if not isinstance(outopt['ylimmir'],bool): outopt['ylimmir']=xymirrorDEF
    if not isinstance(outopt['zlimmir'],bool): outopt['zlimmir']=zmirrorDEF
    return outopt

####################################################################################################################################
####################################################################################################################################
# COLORMAP METHODS

####################################################################################################################################
# METHOD TO MAP ARRAY TO A COLORMAP
def map_array(mat,scale=None,mn=None,mx=None,cd=None):
    """
    Maps an array onto an integer array corresponding to a colormap
    """
    # Pars input
    rm=np.ravel(mat)
    scale=mmlpars.mml_pars(scale,default='log',list=['linear','log'])
    mn=mmlpars.mml_pars(mn,default=min(rm),type=float)
    mx=mmlpars.mml_pars(mx,default=max(rm),type=float,min=mn)
    if mn==mx:
        mn=min(rm)
        mx=max(rm)
    cd=mmlpars.mml_pars(cd,default=rm.mean(),type=float)
    if cd==0: cd=rm.mean()
    if scale=='linear': cd==0
    # Scale array
    if scale=='log':
        matscl = 255.*np.log(1.+(mat-mn)/(cd)) / np.log(1.+(mx-mn)/(cd))
    elif scale=='linear':
        matscl = 255.*(mat-mn)/(mx-mn)
    else: raise Exception('[simplot.map_array] Invalid scale: {}'.format(scale))
    # Return scaled array
    return matscl.astype(int)

####################################################################################################################################
# METHOD TO LIST COLORMAPS
def list_cmaps():
    """
    Lists default order for selecting colormaps
    """
    cmapLIST=['mmlparttype_{}'.format(imap) for imap in mmlplot.list_cmaps('mmlparttype',inclrev=False)]
    return cmapLIST

####################################################################################################################################
# METHOD TO LIST COLORS
def list_colors(clrsim=False):
    """
    Lists default order for selecting colors
    """
    if clrsim:
        cmapLIST=list_cmaps()
        inorm=matplotlib.colors.Normalize(vmin=0.,vmax=1.)
        colorLIST=[]
        for icmapname in cmapLIST:
            icmap=mmlplot.get_cmap(icmapname)
            ismapObj=matplotlib.cm.ScalarMappable(cmap=icmap,norm=inorm)
            colorLIST.append(ismapObj.to_rgba(0.5))
    else:
        cnameLIST=['magenta','cyan','green','blue','red2','yellow']
        colorLIST=[mmlplot.get_mmlcolor(icname) for icname in cnameLIST]
    return colorLIST

####################################################################################################################################
# METHOD TO LIST LINE STYLES
def list_linestyles():
    """
    Lists default order for selecting line styles
    """
    styleLIST=['-','--','-.',':']
    return styleLIST

####################################################################################################################################
# METHOD TO CREATE COLORMAPS
def get_cmaps(cmapColors=None,ptyps=[0,1,2,3,4,5],pgal=0,reverse=False,
              iColor=(0,0,0),fColor=(1,1,1),mColor=(0,0,0),
              iAlpha=1.0    ,fAlpha=1.0    ,mAlpha=0.0):
    '''
    get_cmaps(imag_cmap,ityp,igal)
    '''
    # Set constants
    typList=['gas','halo','disk','bulge','stars','bndry']
    typDict={typList[ityp]:ityp for ityp in range(len(typList))}
    cmapColors0DEF=[(0.,1.,1.),(0.,1.,0.),(0.,0.,1.), # cyan, green, blue                                                              
                    (1.,0.,1.),(1.,0.,0.),(1.,1.,0.)] # magenta, red, yellow                                                  
    cmapColors1DEF=[(0.,0.6,1.),(0.6,1.,0.),(0.,0.6,1.),
                    (1.,0.,0.6),(1.,0.6,0.),(0.8,1.,0.)]
    cmapColors2DEF=6*[(1.,0.84,0.)]
    cmapColors3DEF=6*[(0.64,0.,1.)]
    cmapColors4DEF=6*[(0.,0.,1.)]
    cmapColors5DEF=6*[(1.,0.,0.)]
    cmapColors6DEF=6*[(0.,1.,0.)]
    cmapColors7DEF=6*[(0.,1.,1.)]
    cmapColors8DEF=6*[(1.,0.,1.)]
    cmapColors9DEF=6*[(1.,1.,0.)]
    cmapColorsDEF=[cmapColors0DEF,cmapColors0DEF,cmapColors1DEF,cmapColors2DEF,
                   cmapColors3DEF,cmapColors4DEF,cmapColors5DEF,cmapColors6DEF,
                   cmapColors7DEF,cmapColors8DEF,cmapColors9DEF]
    cmapColorsMML1=['cyan' ,'green' ,'blue' ,'magenta' ,'red' ,'yellow' ]
    cmapColorsMML2=['cyan2','green2','blue2','magenta2','red2','yellow2']
    cmapColorsLGT1=['pnbody_mmllight5','pnbody_mmllight2','pnbody_mmllight0',
                    'pnbody_mmllight1','pnbody_mmllight3','pnbody_mmllight4']
    cmapColorsLGT2=['pnbody_greenlut','pnbody_bluelut','pnbody_heat']
    # Set defaults
    if not isinstance(ptyps,list): ptyps=[ptyps]
    if cmapColors==None or cmapColors=='DEF': cmapColors=cmapColorsDEF[pgal]
    if cmapColors==None or cmapColors=='LGT': cmapColors=cmapColorsLGT
    if not isinstance(cmapColors,list): cmapColors=[cmapColors]*(len(typList)+1)
    # Create color maps
    cmaps=[]
    for itypCount in range(len(ptyps)):
        if    isinstance(ptyps[itypCount],int): ityp=ptyps[itypCount]
        elif  isinstance(ptyps[itypCount],str): ityp=typDict[ptyps[itypCount].lower()]
        else: raise Exception('Particle types must be strings or integers not {}'.format(type(ptyps[itypCount])))
        # Primary colormap
        if isinstance(cmapColors[ityp],str):
            cmaps.append(mmlplot.get_cmap(cmapColors[ityp]))
#            cmaps.append(matplotlib.cm.get_cmap(name=cmapColors[ityp]))
        else:
            if reverse:
                cend=tuple(np.array(cmapColors[ityp])/2.)
                cbeg=(1,1,1)
            else:
                cbeg=tuple(np.array(cmapColors[ityp])/2.)
                cend=(1,1,1)
            cmaps.append(matplotlib.colors.LinearSegmentedColormap.from_list('defGal',[cbeg,cend])) # Fad to white
        cmaps[-1].set_under(iColor,iAlpha)
        cmaps[-1].set_over(fColor,fAlpha)
        cmaps[-1].set_bad(mColor,mAlpha)
#        if iColor != None: cmaps[-1].set_under(iColor)
#        if fColor != None: cmaps[-1].set_over(fColor)
#        if mColor != None: cmaps[-1].set_bad(mColor)
    # Return output
    if len(cmaps)==1: cmaps=cmaps[0]
    return cmaps


####################################################################################################################################
# METHOD TO PLOT OTHER STATISTICS
####################################################################################################################################
def plot_stats(statDict,statDict2=None,method=None,fname=None,overwrite=None,residflag=None,
               colors=None,labels=None,**plot_kw):
    """
    NAME:
        simplot.plot_stats
    PURPOSE:
        To plot various simulation statistics.
    CALLING:
        plot_stats(statDict[,fname=None,overwrite=False,**plot_kw])
    ARGUMENTS:
        statDict:  Dictionary of simulation statistics (See pNbody.gadget.list_stats)
    KEYWORDS:
        fname:     Str absolute path to file where energy plot should be saved.
        overwrite: Bool specifying if existing file should be overwritten 
        Additional keywords are assumed to be mmlplot options (See mmlplot.parsopt)
    """
    # Set constants
    fnameTAG='plot_stats'
    # Pars input
    method=mmlpars.mml_pars(method,default='calclist',type=str)
    options=mmlplot.parsopt(plot_kw)
#    options['ylabpos']=(-0.2,0.5)
    options['pad']=2.0
    options['ylabpos']=(0.06,0.5)
    options['xlabpos']=(0.5,-0.2)
    options['figsize']=(15.,8.)
    colorsDEF=2*[options['fgclr']]
    labelsDEF=['1','2']
    fname=mmlpars.mml_pars(fname,default=os.path.join(options['plotdir'],fnameTAG+options['plotext']),type=str)
    overwrite=mmlpars.mml_pars(overwrite,default=False,type=bool)
    residflag=mmlpars.mml_pars(residflag,default=False,type=bool)
    colors=mmlpars.mml_pars(colors,default=colorsDEF,type=list)
    labels=mmlpars.mml_pars(labels,default=labelsDEF,type=list)
    if   len(colors)==0: colors=colorsDEF
    elif len(colors)==1: colors=colors+colorsDEF[1]
    if   len(labels)==0: labels=labelsDEF
    elif len(labels)==1: labels=labels+labelsDEF[1]
    if not mmlfiles.prep_write(fname,overwrite): return
    if statDict2==None: residflag=False
    # Determine order of plots
    method=method.lower()
    if method=='calclist':
        xTags=[['comx'   ,'comx'   ,'comy'   ,'time'   ],
               ['vcomx'  ,'vcomx'  ,'vcomy'  ,'time'   ],
               ['time'   ,'time'   ,'time'   ,'time'   ],
               ['time'   ,'time'   ,'time'   ,'time'   ],
               ['time'   ,'time'   ,'time'   ,'time'   ],
               ['time'   ,'time'   ,'time'   ,'time'   ]]
        yTags=[['comy'   ,'comz'   ,'comz'   ,'clause' ],
               ['vcomy'  ,'vcomz'  ,'vcomz'  ,'sigma'  ],
               ['halo_A' ,'disk_A' ,'bar_amp','bulge_A'],
               ['halo_r' ,'disk_r' ,'bar_x'  ,'bulge_r'],
               [None     ,'disk_z' ,'bar_phi',None     ],
               [None     ,None     ,'bar_m2' ,None     ]]
    elif method=='general':
        xTags=[['comx'   ,'comx'   ,'comy'   ],
               ['vcomx'  ,'vcomx'  ,'vcomy'  ],
               ['time'   ,'time'   ,'time'   ]]
        yTags=[['comy'   ,'comz'   ,'comz'   ],
               ['vcomy'  ,'vcomz'  ,'vcomz'  ],
               ['sigma'  ,'clause' ,None     ]]
    elif method=='profile':
        xTags=[['time'   ,'time'   ,'time'   ],
               ['time'   ,'time'   ,'time'   ],
               ['time'   ,'time'   ,'time'   ]]
        yTags=[['halo_A' ,'disk_A' ,'bulge_A'],
               ['halo_r' ,'disk_r' ,'bulge_r'],
               [None     ,'disk_z' ,None     ]]
    elif method == 'spin':
        xTags=[['time'   ,'time'   ,'time'   ,'time'   ]]
        yTags=[['all'    ,'halo'   ,'disk'   ,'bulge'  ]]
    elif 'bar' in method:
        xTags=[['time'   ,'time'   ,'time'   ,'time'   ],
               ['time'   ,'time'   ,'time'   ,'time'   ]]
        yTags=[['r'      ,'x'      ,'y'      ,'z'      ],
               ['A'      ,'Aalt'   ,'pa'     ,'omega'  ]]
    elif 'bulk' in method:
        timeidx=np.argsort(np.array(statDict['time']))
        for ikey in statDict['keylist']: statDict[ikey]=np.array(statDict[ikey])[timeidx]
        for ikey in statDict['keylist']: statDict[ikey]=statDict[ikey][1:] # zero potential in IC
        statDict['dtime']=statDict['time'][1:]
        statDict['dEp/Ep']=(statDict['Ep'][1:]-statDict['Ep'][0])/statDict['Ep'][0]
        statDict['dEk/Ek']=(statDict['Ek'][1:]-statDict['Ek'][0])/statDict['Ek'][0]
        statDict['dvsig/vsig']=(statDict['vsig'][1:]-statDict['vsig'][0])/statDict['vsig'][0]
        statDict['dLz/Lz']=(statDict['Lz'][1:]-statDict['Lz'][0])/statDict['Lz'][0]
        for ikey in statDict['keylist']: statDict[ikey]=list(statDict[ikey])
        print (statDict2!=None)
        if statDict2!=None:
            timeidx2=np.argsort(np.array(statDict2['time']))
            for ikey in statDict2['keylist']: statDict2[ikey]=np.array(statDict2[ikey])[timeidx2]
            for ikey in statDict2['keylist']: statDict2[ikey]=statDict2[ikey][1:] # zero potential in IC
            statDict2['dtime']=statDict2['time'][1:]
            statDict2['dEp/Ep']=(statDict2['Ep'][1:]-statDict2['Ep'][0])/statDict2['Ep'][0]
            statDict2['dEk/Ek']=(statDict2['Ek'][1:]-statDict2['Ek'][0])/statDict2['Ek'][0]
            statDict2['dvsig/vsig']=(statDict2['vsig'][1:]-statDict2['vsig'][0])/statDict2['vsig'][0]
            statDict2['dLz/Lz']=(statDict2['Lz'][1:]-statDict2['Lz'][0])/statDict2['Lz'][0]
            for ikey in statDict2['keylist']: statDict2[ikey]=list(statDict2[ikey])
            
        xTags=[['time'   ,'time'   ,'time'      ],
               ['dtime'  ,'dtime'  ,'dtime'     ],
               ['time'   ,'time'   ,'dtime'     ]]
        yTags=[['Ep'     ,'Ek'     ,'vsig'      ],
               ['dEp/Ep' ,'dEk/Ek' ,'dvsig/vsig'],
               ['Lz_circ','Lz'     ,'dLz/Lz'    ]]
    else: raise Exception('Invalid method: {}'.format(method))
    nrow,ncol=np.array(xTags).shape
    # Get limits
    if method in ['calclist','general']:
        rlim=np.absolute(np.array([statDict['comx' ],statDict['comy' ],statDict['comz' ]])).max()*(1.+options['limpad'])
        vlim=np.absolute(np.array([statDict['vcomx'],statDict['vcomy'],statDict['vcomz']])).max()*(1.+options['limpad'])
    # Initialize figure
    figObj=plt.figure()
    figObj.suptitle(os.path.basename(statDict['fname'][0]).split('.')[0])
    # Plot
#    mmlio.verbose('File range:')
#    print '    '+statDict['fname'][0]
#    print '    '+statDict['fname'][1]
    timeidx=np.argsort(np.array(statDict['time']))
    if statDict2!=None: timeidx2=np.argsort(np.array(statDict2['time']))
    axList=[]
    iplt=0
    for irow in range(nrow):
        for icol in range(ncol):
            iplt+=1
            if yTags[irow][icol]==None:
                pass
            else:
                # Create axes
                axList.append(plt.subplot(nrow,ncol,iplt))
                # Get time
                if xTags[irow][icol]=='dtime':
                    itimeidx=np.argsort(np.array(statDict['dtime']))
                    if statDict2!=None: itimeidx2=np.argsort(np.array(statDict2['dtime']))
                else:
                    itimeidx=timeidx
                    if statDict2!=None: itimeidx2=timeidx2
                # Select data
                ix=np.array(statDict[xTags[irow][icol]])[itimeidx]
                iy=np.array(statDict[yTags[irow][icol]])[itimeidx]
                ixlab=xTags[irow][icol] ; iylab=yTags[irow][icol]
                if statDict2!=None:
                    ix2=np.array(statDict2[xTags[irow][icol]])[itimeidx2]
                    iy2=np.array(statDict2[yTags[irow][icol]])[itimeidx2]
                    if residflag:
                        iylab+=' resid' ; iy1=iy
                        iy=(iy1-interpolate.griddata(ix2,iy2,ix))/iy1
                # Plot
                iplot=axList[-1].plot(ix,iy,color=colors[0],label=labels[0])
                if statDict2!=None and not residflag:
                    iplot2=axList[-1].plot(ix2,iy2,color=colors[1],label=labels[1])
                # Add labels
                axList[-1].set_xlabel(ixlab) ; axList[-1].set_ylabel(iylab)
                # Legend
                if statDict2!=None and not residflag and iplt==1: plt.legend()
                # Plot first point for center of mass stuff
                if iylab.startswith('com') or iylab.startswith('vcom'):
                    data=iplot[0].get_data()
                    axList[-1].plot(data[0][0],data[1][0],'o',color=colors[0])
                    if statDict2!=None and not residflag:
                        data2=iplot2[0].get_data()
                        axList[-1].plot(data2[0][0],data2[1][0],'o',color=colors[1])
                    if iylab.startswith('com' ) and not residflag:
                        axList[-1].set_xlim(-rlim,rlim) ; axList[-1].set_ylim(-rlim,rlim)
                    if iylab.startswith('vcom') and not residflag:
                        axList[-1].set_xlim(-vlim,vlim) ; axList[-1].set_ylim(-vlim,vlim)
    # Set axes properties
    mmlplot.set_figprop(figObj,options)
    figObj.tight_layout(pad=options['pad'],h_pad=options['hpad'],w_pad=options['wpad'])
    # Display plot
    if options['showplot']: figObj.show()
    # Save plot
    if options['outplot']:
        return figObj
    else:
        figObj.savefig(fname,facecolor=options['bgclr'],edgecolor=options['bgclr'])
        return None

####################################################################################################################################
# METHOD TO PLOT ENERGY STATISTICS
def plot_energy(t,ek,ep,fname=None,overwrite=False,labels=None,colors=None,
                deemethod=None,residualflag=None,**plot_kw):
    '''
    NAME:
        simplot.plot_energy
    PURPOSE:
        To plot energy statistics (2*T/W and Delta E/E)
    CALLING:
        plot_energy(t,ek,ep[,fname=None,overwrite=False,**plot_kw])
    INPUT:
        t:         Array of times (Gyr)
        ek:        Array of kinetic energies
        ep:        Array of potential energies
    KEYWORDS:
        fname:     Str absolute path to file where energy plot should be saved.
        overwrite: Bool specifying if existing file should be overwritten
        Additional keywords are assumed to be mmlplot options (See mmlplot.parsopt)
    '''
    # Set constants
    fnameTAG='plot_energy'
    clrLIST=['w','b','r','g','m','y']
    # Pars input
    if not isinstance(t,tuple): t=(t) ; ek=(ek) ; ep=(ep)
    for it,iek,iep in zip(t,ek,ep):
        inin=len(it)
        if len(iek) != inin:
            raise Exception('List of kinetic energies has a different length ({}) than list of times ({}).'.format(len(iek),inin))
        if len(iep) != inin:
            raise Exception('List of potential energies has a different length ({}) than list of times ({}).'.format(len(iep),inin))
    options=mmlplot.parsopt(plot_kw)
    options['ylabpos']=(1.1,0.5)
    options['xlabpos']=(0.5,-0.1)
    options['pad']=2.0
    clrLIST[0]=options['fgclr']
    labels=mmlpars.mml_pars(labels,type=list,default=[str(idx) for idx in range(len(t))],nelements=len(t))
    colors=mmlpars.mml_pars(colors,type=list,default=clrLIST[:len(t)])
    fname=mmlpars.mml_pars(fname,default=os.path.join(options['plotdir'],fnameTAG+options['plotext']),type=str)
    overwrite=mmlpars.mml_pars(overwrite,type=bool,default=False)
    deemethod=mmlpars.mml_pars(deemethod,type=str,default='previous',list=['previous','first'])
    if   deemethod == 'previous': errdeelim=0.00005
    elif deemethod == 'first'   : errdeelim=0.01
    else: raise Exception('[simplot.plot_energy] Invalid dee method: {}'.format(deemethod))
    residualflag=mmlpars.mml_pars(residualflag,type=bool,default=False)
    if len(t)!=2: residualflag=False
    if not mmlfiles.prep_write(fname,overwrite): return
    # Create new arrays
    ttw=[] ; etot=[] ; dee=[] ; max_tlist=[] ; max_ttwlist=[] ; max_deelist=[]
    for it,iek,iep in zip(t,ek,ep):
        ittw=abs(2.*iek/iep)
        ietot=iek+iep
        if   deemethod == 'previous': idee=(ietot[1:]-ietot[:-1])/ietot[:-1]
        elif deemethod == 'first'   : idee=(ietot[1:]-ietot[0])/ietot[0]
        else: raise Exception('[simplot.plot_energy] Invalid dee method: {}'.format(deemethod))
        idee=np.ma.masked_outside(idee[1:],-errdeelim,errdeelim)
        max_tlist.append(it.max())
        max_ttwlist.append(abs(ittw-1.).max()+1.)
        max_deelist.append(abs(idee).max())
        ttw.append(ittw)
        etot.append(ietot)
        dee.append(idee)
    ttw=tuple(ttw)
    etot=tuple(etot)
    dee=tuple(dee)
    if residualflag:
        ttwresid=np.abs(np.ma.masked_invalid((ttw[0]-interpolate.griddata(t[1],ttw[1],t[0]))/ttw[0]))
        deeresid=np.abs(np.ma.masked_invalid((dee[0]-interpolate.griddata(t[1][2:],dee[1],t[0][2:]))/dee[0]))
    # Get limits
    max_t=min(max_tlist) ; min_t=0.
    max_ttw=max(max_ttwlist) ; min_ttw=1.-(max_ttw-1.)
    max_dee=max(max_deelist) ; min_dee=-max_dee
#    max_ttw=abs(ttw-1.).max()+1. ; min_ttw=1.-(max_ttw-1.)
#    max_dee=abs(dee).max() ; min_dee=-max_dee
    if residualflag:
        max_ttwres=np.abs(ttwresid).max() ; min_ttwres=0.#-max_ttwres
        max_deeres=np.abs(deeresid).max() ; min_deeres=0.#-max_deeres
    ypad_ttw=options['limpad']*(max_ttw-min_ttw)
    ypad_dee=options['limpad']*(max_dee-min_dee)
    if residualflag:
        ypad_ttwres=options['limpad']*(max_ttwres-min_ttwres)
        ypad_deeres=options['limpad']*(max_deeres-min_deeres)
    xlim=(min_t,max_t)
    ylim_ttw=(min_ttw-ypad_ttw,max_ttw+ypad_ttw)
    ylim_dee=(min_dee-ypad_dee,max_dee+ypad_dee)
    if residualflag:
        ylim_ttwres=(min_ttwres-ypad_ttwres,max_ttwres+ypad_ttwres)
        ylim_deeres=(min_deeres-ypad_deeres,max_deeres+ypad_deeres)
    # Create labels
    xlab='Time (Gyr)'
    ylab_ttw='2*T/W'
    ylab_dee='[E(t)-E(t-1)]/E(t)'
    # Initialize figure
    fig=plt.figure()
    if residualflag: plotdim=(2,2)
    else           : plotdim=(2,1)
    iplot=1
    # Plot 2*T/W
    ax_ttw=plt.subplot(plotdim[0],plotdim[1],iplot) ; iplot+=1
    for idx in range(len(t)):
        ax_ttw.plot(t[idx],ttw[idx],'-',color=colors[idx],label=labels[idx])
    ax_ttw.axhline(1.0,linestyle='--',color=options['fgclr'])
    ax_ttw.set_title(ylab_ttw,axes=ax_ttw,figure=fig)
#    ax_ttw.set_ylabel(ylab_ttw,axes=ax_ttw,figure=fig)
    ax_ttw.set_ylim(ylim_ttw)
    ax_ttw.set_xlim(xlim)
    ax_ttw.tick_params(labelbottom=False)
    if len(t) > 1: plt.legend()
    if residualflag:
        ax_ttwres=plt.subplot(plotdim[0],plotdim[1],iplot) ; iplot+=1
        ax_ttwres.plot(t[0],ttwresid,'-',color=options['fgclr'])
        ax_ttwres.axhline(0.0,linestyle='--',color=options['fgclr'])
        ax_ttwres.set_title(ylab_ttw+' Residual',axes=ax_ttwres,figure=fig)
#        ax_ttwres.set_ylabel(ylab_ttw+' Residual',axes=ax_ttwres,figure=fig)
        ax_ttwres.set_ylim(ylim_ttwres)
        ax_ttwres.set_xlim(xlim)
        ax_ttwres.tick_params(labelbottom=False)
    # Plot (E2-E1)/E1
    ax_dee=plt.subplot(plotdim[0],plotdim[1],iplot,sharex=ax_ttw) ; iplot+=1
    for idx in range(len(t)):
        ax_dee.plot(t[idx][2:],dee[idx],'-',color=colors[idx],label=labels[idx])
    ax_dee.axhline(0.0,linestyle='--',color=options['fgclr'])
    ax_dee.set_xlabel(xlab)
    ax_dee.set_title(ylab_dee)
#    ax_dee.set_ylabel(ylab_dee)
    ax_dee.set_ylim(ylim_dee)
    ax_dee.set_xlim(xlim)
    if residualflag:
        ax_deeres=plt.subplot(plotdim[0],plotdim[1],iplot,sharex=ax_ttwres) ; iplot+=1
        ax_deeres.plot(t[0][2:],deeresid,'-',color=options['fgclr'])
        ax_deeres.axhline(0.0,linestyle='--',color=options['fgclr'])
        ax_deeres.set_title(ylab_dee+' Residual')
#        ax_deeres.set_ylabel(ylab_dee+' Residual')
        ax_deeres.set_ylim(ylim_deeres)
        ax_deeres.set_xlim(xlim)
    # Set axes properties
    mmlplot.set_figprop(fig,options)
    # Display plot
    if options['showplot']: fig.show()
    # Save plot
    if options['outplot']:
        return fig
    else:
        fig.savefig(fname,facecolor=options['bgclr'],edgecolor=options['bgclr'])
        return None

###################################################################################################################################
###################################################################################################################################
# PROVIDE COMMAND LINE ACCESS
if __name__ == '__main__': main()
