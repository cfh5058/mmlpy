#!/usr/bin/python
import sys,os,shutil,glob,copy,pprint,scipy,math
import numpy as np
import matplotlib as mplib
LIST_METHODS_RUN=['energy','profparam','anim','stat']
LIST_METHODS_SNP0=['hist1d','hist2d']
LIST_METHODS_SNP=LIST_METHODS_SNP0+[ilong+'_soft' for ilong in LIST_METHODS_SNP0]
LIST_METHODS=LIST_METHODS_RUN+LIST_METHODS_SNP
from mmlutils import *
import main as mmlnbody
import simlist,simstat,simcalc,simplot,simfile,mmlgadget

def main():
    """
    Provides command line access to COMPARE methods
    """
    out = mmlnbody.walk(mtype='compare')
    return out

####################################################################################################################################
####################################################################################################################################
# METHODS FOR GETTING LISTS
def dict_ftype():
    """
    Returns a dictionary of files listed by type
    """
    fdict={'shortout':LIST_METHODS_RUN+['short'+ilong for ilong in LIST_METHODS_SNP],
           'shortlon':['short'+ilong for ilong in LIST_METHODS_SNP],
           'longout' :LIST_METHODS_SNP}
    return fdict
def dict_fgroup():
    """
    Returns a dictionary of file types listed by group  
    """
    fdict={'input' :['shortin' ,'longin' ],
           'output':['shortout','longout']}
    return fdict

####################################################################################################################################
####################################################################################################################################
# FILE CLASSES AND METHODS

####################################################################################################################################
# METHOD FOR PARSING COMPARE FILE OPTIONS
def pars_fileopt(fopt):
    """
    Returns parsed COMPARE file option dictionary with keys:
        plotext: Str extension to use for plot files
        animext: Str extension to use for animation files
    """
    # Define keys
    form={
        'plotext': mmlpars.parsdict(default='png',type=str),
        'animext': mmlpars.parsdict(default='mp4',type=str),
        'runtag2': mmlpars.parsdict(default='None',type=str)
        }
    # Pars input
    fopt=mmlpars.mml_formpars(fopt,form)
    # Return parsed input
    return fopt

####################################################################################################################################
# METHOD TO CREATE FILE DICTIONARY
def files(fdict,options=None,**input_kw):
    """
    Returns a dictionary of COMPARE files and directories
    """
    # Set constants
    shrtLIST=LIST_METHODS_RUN
    longLIST=LIST_METHODS_SNP
    statLIST=simstat.LIST_STATS
    # Pars input
    options=pars_fileopt(options)
    comptag=fdict['runtag']+'_VS_'+options['runtag2']
    pfix=comptag+'.'
    plotext=options['plotext']
    animext=options['animext']
    # Initialize directory
    files=dict(fdict['compare'],plotext=plotext,animext=animext)
    files['pfix']=pfix
    # Singular plots
    fkeys=[] ; flist=[]
    for ishrt in shrtLIST:
        if   ishrt=='anim': continue
        elif ishrt=='stat':
            for istat in statLIST:
                fkeys.append('stat_{}'.format(istat))
                flist.append('stat_{}.{}'.format(istat,plotext))
                fkeys.append('stat_{}_resid'.format(istat))
                flist.append('stat_{}_resid.{}'.format(istat,plotext))
        else:
            fkeys.append(ishrt)
            flist.append('{}.{}'.format(ishrt,plotext))
    # Multiple plot directories
    dkeys=[] ; dlist=[]
    for ilong in longLIST+['anim']:
        dkeys.append(ilong+'dir')
        dlist.append(ilong)
    # Add plots and directories
    files['subdirs']=dkeys
    files=mmlfiles.filedict(files['dir'],dirkeys=dkeys,dirlist=dlist,filekeys=fkeys,filelist=flist,pfix=pfix,fdict=files)
    # Multiple plot basenames & animations
    for ilong in longLIST:
        files[ilong+'base']=os.path.join(files[ilong+'dir'],pfix+ilong)
        files[ilong+'anim']=os.path.join(files['animdir'],pfix+ilong+'anim.'+animext)
    # Return file dictionary
    fdict['compare']=files
    return fdict

####################################################################################################################################
# METHOD TO RETURN RELEVANT LIST OF FILES
def get_filelist(fdict,keylist):
    """
    Returns a list of relevant COMPARE files.
    """
    # Set constants
    shrtLIST=LIST_METHODS_RUN
    longLIST=LIST_METHODS_SNP
    # Pars input
    fdict=mmlpars.mml_pars(fdict,type=dict)
    keylist=mmlpars.mml_pars(keylist,type=list)
    # Add files
    filelist={}
    shrtLIST=['energy','profparam','anim','stat']
    for ishrt in shrtLIST:
        if ishrt in keylist:
            if ishrt=='anim': 
                filelist[ishrt]=[]
            else:
                filelist[ishrt]=[fdict[ishrt]]
    for ilong in longLIST:
        ishrt='short'+ilong
        if ilong in keylist or ishrt in keylist or 'anim' in keylist:
            ibase=fdict[ilong+'base']+'*'
            if ilong in keylist: filelist[ilong]=ibase
            if ishrt in keylist:
                ilongfiles=glob.glob(ibase)
                if len(ilongfiles)>0: filelist[ishrt]=[sorted(ilongfiles)[-1]]
            if 'anim' in keylist: filelist['anim'].append(fdict[ilong+'anim'])
    # Return output
    return filelist

####################################################################################################################################
####################################################################################################################################
# HIGH LEVEL METHODS REQUIRING MMLSIM 

####################################################################################################################################
# METHOD FOR RUNNING DIFFERENT COMPARISON METHODS
def run(simstr1,method,verbose=None,simstr2=None,runtag2=None,overwrite=None,**method_kw):
    """
    Provides interface for running different COMPARE methods
    """
    # Set constants
    methLIST=LIST_METHODS
    # Pars input
    method=mmlpars.mml_pars(method,list=methLIST)
    verbose=mmlpars.mml_pars(verbose,type=bool,default=True)
    method_kw['verbose']=verbose
    # Initialize default output
    out=None
    # Get second run for comparison
    if simstr2 == None:
        if runtag2==None:
            simtyp2=mmlio.askselect('What type of simulation should {} be compared with?'.format(simstr1['runtag']),simlist.LIST_RUNTYPS)
            runtag2=mmlio.askselect('Which {} simulation should {} be compared with?'.format(simtyp2,simstr1['runtag']),simfile.fpar_taglist('icgen',simtyp2))
        simstr2=mmlnbody.runtag2simobj(runtag2)
    if not isinstance(overwrite,bool): overwrite=mmlio.yorn('Overwrite existing files?')
    method_kw['overwrite']=overwrite
    # Proceed based on method
    if   method in LIST_METHODS_RUN: out=comprun( simstr1,simstr2,method=method,**method_kw)
    elif method in LIST_METHODS_SNP: out=compsnap(simstr1,simstr2,method=method,**method_kw)
    else: raise Exception('Option for method {} needs to be added to compare.run.'.format(method))
    # Return output
    return out

####################################################################################################################################
# METHOD FOR COMPARING THINGS FOR AN ENTIRE RUN
def comprun(simstr1,simstr2,method=None,**method_kw):
    """
    Compares statistics for entire runs
    """
    # Set constants
    methLIST=LIST_METHODS_RUN
    # Pars input
    method=mmlpars.mml_pars(method,default='all',list=['all']+methLIST)
    # Proceed based on method
    if method=='all':
        out=[]
        for imeth in methLIST: out.append(comprun(simstr1,simstr2,method=imeth,**method_kw))
    elif method=='energy'   : out=comprun_energy(simstr1,simstr2,**method_kw)
    elif method=='profparam': out=comprun_profparam(simstr1,simstr2,**method_kw)
    elif method=='anim'     : out=comprun_anim(simstr1,simstr2,**method_kw)
    elif method=='stat'     : out=comprun_stat(simstr1,simstr2,**method_kw)
    else: raise Exception('[compare.comprun] Invalid method: {}'.format(method))
    # Return output
    return out

####################################################################################################################################
# METHOD FOR COMPARING RUN STATISTICS
def comprun_stat(simstr1,simstr2,verbose=None,**plot_kw):
    """
    Compares statistics between runs
    """
    # Get info
    statlist=simstat.LIST_STATS
    labels=[simstr1['runtag'],simstr2['runtag']]
    colors=['b','r']
    # Pars input
    plot_kw=mmlplot.parsopt(plot_kw)
    verbose=mmlpars.mml_pars(verbose,default=True,type=bool)
    if verbose: mmlio.verbose('Beginning comparison of energies...')
    # Get filenames
    fdict1=simstr1.mkfiledict(options={'runtag2':simstr2['runtag']})
    fdict2=simstr2.mkfiledict()
    # Loop over stats
    out={}
    for istat in statlist:
        # Get file names
        plotfile_compr=fdict1['compare']['stat_'+istat]
        plotfile_resid=fdict1['compare']['stat_'+istat+'_resid']
        fname1=fdict1['stat']['stat_'+istat]
        fname2=fdict1['stat']['stat_'+istat]
        if not os.path.isfile(fname1) or not os.path.isfile(fname2):
            if verbose: mmlio.verbose('One or both of the {} statistics files do not exist.'.format(istat))
            if not os.path.isfile(fname1) and verbose: print '    '+fname1
            if not os.path.isfile(fname2) and verbose: print '    '+fname2
            continue
        # Load dictionaries
        statdict1=simstat.rwstats('R',stattype=istat,statfile=fname1)
        statdict2=simstat.rwstats('R',stattype=istat,statfile=fname2)
        # Plot comparison
        out[istat+'_compr']=simplot.plot_stats(statdict1,statDict2=statdict2,method=istat,
                                               fname=plotfile_compr,residflag=False,
                                               colors=colors,labels=labels,**plot_kw)
        if verbose: print '    '+plotfile_compr
        # Plot residual
        out[istat+'_resid']=simplot.plot_stats(statdict1,statDict2=statdict2,method=istat,
                                               fname=plotfile_resid,residflag=True,
                                               colors=colors,labels=labels,**plot_kw)
        if verbose: print '    '+plotfile_resid
    # Return output
    return out
        
####################################################################################################################################
# METHOD FOR COMPARING ENERGY EVOLUTION BETWEEN RUNS
def comprun_energy(simstr1,simstr2,verbose=None,**plot_kw):
    """
    Compares energy statistics between runs
    """
    # Pars input
    plot_kw=mmlplot.parsopt(plot_kw)
    verbose=mmlpars.mml_pars(verbose,default=True,type=bool)
    if verbose: mmlio.verbose('Beginning comparison of energies...')
    # Get filenames
    fdict1=simstr1.mkfiledict(options={'runtag2':simstr2['runtag']})
    fdict2=simstr2.mkfiledict()
    fname1=fdict1['gadget']['output']['energy']
    fname2=fdict2['gadget']['output']['energy']
    plotfile=fdict1['compare']['energy']
    # Load data
    edict1=mmlgadget.read_energy(fname1)
    edict2=mmlgadget.read_energy(fname2)
    t=(np.array(edict1['t'],float),np.array(edict2['t'],float))
    ek=(np.array(edict1['ek'],float),np.array(edict2['ek'],float))
    ep=(np.array(edict1['ep'],float),np.array(edict2['ep'],float))
    labels=[simstr1['runtag'],simstr2['runtag']]
    colors=['b','r']
    plot_kw['figsize']=(20.,10.)
    out=simplot.plot_energy(t,ek,ep,fname=plotfile,residualflag=True,
                            labels=labels,colors=colors,**plot_kw)
    if verbose: print '    '+plotfile
    # Return
    return out

####################################################################################################################################
# METHOD FOR COMPARING PROFILE PARAMETERS FOR A RUN
def comprun_profparam(simstr1,simstr2,**method_kw):
    """
    Compares profile parameters between runs
    """
    raise Exception('comprun_profparam is a work in progress.')

####################################################################################################################################
# METHOD FOR CREATING COMPARISON ANIMATIONS FOR TWO RUNS
def comprun_anim(simstr1,simstr2,flagdict=None,setflags=None,verbose=None,overwrite=None,**anim_kw):
    """
    Creates comparison animations for two runs
    """
    # Set constants
    longLIST=LIST_METHODS_SNP
    # Pars input
    flagdict=mmlpars.mml_pars(flagdict,default={},type=dict)
    setflags=mmlpars.mml_pars(setflags,default=False,type=bool)
    verbose=mmlpars.mml_pars(verbose,default=False,type=bool)
    overwrite=mmlpars.mml_pars(overwrite,default=False,type=bool)
    if setflags: verbose=False
    # Get filenames
    fdict1=simstr1.mkfiledict(options={'runtag2':simstr2['runtag']})
    fdict2=simstr2.mkfiledict()
    # Create directory
    mmlfiles.mkdirs(fdict1['compare']['animdir'])
    # Set flags and print info
    if 'anim' not in flagdict: flagdict['anim']=overwrite
    if verbose: mmlio.verbose('[{}_{}] Beginning animations'.format(simstr1['runtag'],simstr2['runtag']))
    # Create list of frames
    baselist=[] ; framDict={} ; animDict={}
    for ilong in longLIST:
        ibase=ilong
        baselist.append(ibase)
        framDict[ibase]=fdict1['compare'][ilong+'base']+'_*.'+fdict1['compare']['plotext']
        animDict[ibase]=fdict1['compare'][ilong+'anim']
    # Loop over list of animations
    for ibase in baselist:
        # Continue if frames do not exist
        if len(glob.glob(framDict[ibase])) == 0: continue
        # Ask in animation should be overwritten
        if ibase not in flagdict:
            if flagdict['anim']: flagdict[ibase]=True
            else: flagdict[ibase]=mmlio.yorn('Overwrite {} animation?'.format(ibase))
        # Create animation
        if not setflags:
            if not os.path.isfile(animDict[ibase]) or flagdict[ibase]:
                movestr=mmlplot.ffmpeg(framDict[ibase],animDict[ibase],
                                       overwrite=flagdict[ibase],rmtemp=True,**anim_kw)
                if verbose: print 'Created animation:'
                if verbose: print '    '+animDict[ibase]
    # Return info
    if setflags: return flagdict
    else       : return flagdict

####################################################################################################################################
# METHOD FOR COMPARING THINGS BETWEEN SNAPSHOTS
def compsnap(simstr1,simstr2,method=None,fextlist1=None,fextlist2=None,allsnap=None,**method_kw):
    """
    Compares statistics between individual snapshots
    """
    # Set constants
    methLIST=LIST_METHODS_SNP
    # Pars input
    method=mmlpars.mml_pars(method,default='all',list=['all']+methLIST)
    if not isinstance(fextlist1,list):
        if not isinstance(allsnap,bool): allsnap=mmlio.yorn('Compare all snapshots?')
        if allsnap:
            fdict1=simstr1.fdict
            fdict2=simstr2.fdict
            snapbase1=fdict1['gadget']['output']['snapbase']
            snapbase2=fdict2['gadget']['output']['snapbase']
            snaplist1=sorted(glob.glob(snapbase1+'*'))
            snaplist2=sorted(glob.glob(snapbase2+'*'))
            nfext=min(len(snaplist1),len(snaplist2))
            fextlist1=sorted([isnap.split(snapbase1)[-1] for isnap in snaplist1],reverse=True)
            mmlio.verbose('WARNING: allsnap comparisons will only be meaningfull if simulations have the same sampling.')
        else:
            fext1=mmlio.askquest('What snapshot extension should be compared from {}?'.format(simstr1['runtag']),dtype='str',default='_000ic')
            fextlist1=[fext1]
    if not isinstance(fextlist2,list):
        if not isinstance(allsnap,bool): allsnap=mmlio.yorn('Compare all snapshots?')
        if allsnap:
            fextlist2=fextlist1
        else:
            fext2=mmlio.askquest('What snapshot extension should be compared from {}?'.format(simstr2['runtag']),dtype='str',default=fextlist1[0])
            fextlist2=[fext2]
    if len(fextlist1) != len(fextlist2):
        raise Exception('[compare.compsnap] List of extensions for each run must be the same length.')
    # Proceed based on method
    for ifext1,ifext2 in zip(fextlist1,fextlist2):
        method_kw['fext1']=ifext1
        method_kw['fext2']=ifext2
        if method=='all':
            out=[]
            for imeth in methLIST: out.append(compsnap(simstr1,simstr2,method=imeth,**method_kw))
        elif method=='hist1d': out=compsnap_hist1D(simstr1,simstr2,**method_kw)
        elif method=='hist2d': out=compsnap_hist2D(simstr1,simstr2,**method_kw)
        elif method=='hist1d_soft': out=compsnap_hist1D(simstr1,simstr2,softflag=True,**method_kw)
        elif method=='hist2d_soft': out=compsnap_hist2D(simstr1,simstr2,softflag=True,**method_kw)
        else: raise Exception('[compare.compsnap] Invalid method: {}'.format(method))
    # Return output
    return out

####################################################################################################################################
# METHOD FOR COMPARING 1D HISTOGRAMS BETWEEN SNAPSHOTS
def compsnap_hist1D(simstr1,simstr2,fext1=None,fext2=None,verbose=None,softflag=None,
                    modList=None,typList=None,galList=None,**plot_kw):
    """
    Compares 1D histograms for particle types between snapshots
    """
    # Set constants
    typListDEF=simlist.LIST_PTYPBASE
    modListDEF=['rxyz_rho']
    galListDEF=[1]
    # Pars input
    verbose=mmlpars.mml_pars(verbose,default=True,type=bool)
    if verbose: mmlio.verbose('Beginning 1D histogram comparison of snapshots {}{} and {}{}...'.format(simstr1['runtag'],fext1,simstr2['runtag'],fext2))
    fext1=mmlpars.mml_pars(fext1,type=str,default='_000ic')
    fext2=mmlpars.mml_pars(fext2,type=str,default='_000ic')
    softflag=mmlpars.mml_pars(softflag,default=False,type=bool)
    modList=mmlpars.mml_pars(modList,default=modListDEF,type=list)
    typList=mmlpars.mml_pars(typList,default=typListDEF,type=list)
    galList=mmlpars.mml_pars(galList,default=galListDEF,type=list)
    # Get filenames
    fdict1=simstr1.mkfiledict(options={'runtag2':simstr2['runtag']})
    fdict2=simstr2.fdict
    if softflag:
        plotfile=fdict1['compare']['hist1d_softbase']+fext1+fext2+'.'+fdict1['compare']['plotext']
    else:
        plotfile=fdict1['compare']['hist1dbase']+fext1+fext2+'.'+fdict1['compare']['plotext']
    # Get load keys
    loadkeys1=simstr1.get_loadkeys(fext=fext1)
    loadkeys2=simstr2.get_loadkeys(fext=fext2)
    # Load data
    if softflag:
        loadcalc=False
        soft1=simstr1.get_softenings()
        soft2=simstr2.get_softenings()
        softmin1=min(soft1.values()) ; softmax1=max(soft1.values())
        softmin2=min(soft2.values()) ; softmax2=max(soft2.values())
        softmax=max(softmax1,softmax2)
        calcdict1=simcalc.initcalc(h1DList=modList,typList=typList,galList=galList)
        calcdict2=simcalc.initcalc(h1DList=modList,typList=typList,galList=galList)
        logkeys=['rxy','rxyz']
        linkeys=['x','y','z']
        limdict={}
        for ilogkey in logkeys: limdict[ilogkey]=(softmax/10.,softmax*10.)
        for ilinkey in linkeys: limdict[ilinkey]=(-softmax*10.,softmax*10.)
        for imod in modList:
            if imod.split('_')[0] not in logkeys+linkeys:
                mmlio.verbose('WARNING: softflag meaningless for mode {}'.format(imod))
    else:
        loadcalc=True
        limdict=None
        calcdict1=None
        calcdict2=None
    out1=simcalc.multcalc(simstr1,calcdict=calcdict1,loadkeys=loadkeys1,limdict=limdict,loadcalc=loadcalc,
                          h1DList=modList,h2DList=[],typList=typList,rayList=[],galList=galList)
    out2=simcalc.multcalc(simstr2,calcdict=calcdict2,loadkeys=loadkeys2,limdict=limdict,loadcalc=loadcalc,
                          h1DList=modList,h2DList=[],typList=typList,rayList=[],galList=galList)
    snapdict1,calcdict1,calcdict01,refkw1=out1
    snapdict2,calcdict2,calcdict02,refkw2=out2
    # Plot
    plot_kw['residualflag']=True
    plot_kw['flag_legend']=True
    plot_kw['labelList']=[simstr1['runtag'],simstr2['runtag']]
    plot_kw['title']=""
    plot_kw['title']+="t1 = {:5.2f} {}".format(calcdict1['statDict']['time'],calcdict1['unitDict']['UnitTime'].symbol)+', '
    plot_kw['title']+="t2 = {:5.2f} {}".format(calcdict2['statDict']['time'],calcdict2['unitDict']['UnitTime'].symbol)
    out=simplot.plotcalc_hist1D(calcdict1,calcdict2=calcdict2,
                                modList=modList,typList=typList,galList=galList,
                                plotfile=plotfile,verbose=False,**plot_kw)
    if verbose: print '    '+plotfile
    # Return
    return out

####################################################################################################################################
# METHOD FOR COMPARING 2D HISTOGRAMS BETWEEN SNAPSHOTS
def compsnap_hist2D(simstr1,simstr2,fext1=None,fext2=None,verbose=None,**plot_kw):
    """
    Compares 2D histograms for particle types between snapshots
    """
    # Pars input
    fext1=mmlpars.mml_pars(fext1,type=str,default='_000ic')
    fext2=mmlpars.mml_pars(fext2,type=str,default='_000ic')
    verbose=mmlpars.mml_pars(verbose,default=True,type=bool)
    if verbose: print 'Beginning 2D histogram comparison of snapshots {}{} and {}{}...'.format(simstr1['runtag'],fext1,simstr2['runtag'],fext2)
    # Get filenames
    fdict1=simstr1.mkfiledict(options={'runtag2':simstr2['runtag']})
    fdict2=simstr2.mkfiledict()
    plotfile=fdict1['compare']['hist2dbase']+fext1+fext2+'.'+fdict1['compare']['plotext']
    # Load data
    calcdict1=simcalc.get_calcsnap(simstr1,fext=fext1,ftype='gadget')
    calcdict2=simcalc.get_calcsnap(simstr2,fext=fext2,ftype='gadget')
    # Plot
    plot_kw['imag_mode']='pos_m'
    plot_kw['cont_mode']='pos_scfpotm2'
    plot_kw['imag_residflag']=True
    plot_kw['cont_residflag']=False
    plot_kw['title']=""
    plot_kw['title']+="t1 = {:5.2f} {}".format(calcdict1['statDict']['time'],calcdict1['unitDict']['UnitTime'].symbol)+', '
    plot_kw['title']+="t2 = {:5.2f} {}".format(calcdict2['statDict']['time'],calcdict2['unitDict']['UnitTime'].symbol)
    plot_kw['title']+=" residual"
    out=simplot.plot_snaphist2D(plotfile=plotfile,calcdict=calcdict1,calcdict2=calcdict2,verbose=False,**plot_kw)
    if verbose: print '    '+plotfile
    # Return
    return out

if __name__ == '__main__':
    main()
