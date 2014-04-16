#!/usr/bin/python 
####################################################################################################################################
#
# METHODS FOR REPORTING ON THE STATUS OF RUNS
#
####################################################################################################################################
import sys,os,shutil,glob,copy,pprint,scipy,math,re,subprocess
import numpy as np
from mmlutils import *
from pysim import files

####################################################################################################################################
# WRITE STATS
def writestat(runlist=None,fname0=None):
    """
    Formated output of run stats
    """
    progsym='='
    nprog0=10
    # Get filename
    if isinstance(runlist,(str,list)): fname0='stdout'
    if not isinstance(fname0,str):
        fname0=os.path.join(os.path.expanduser('~'),'runstatus_')
    # Get stats
    stats=allruns(runlist=runlist)
    keylist=copy.deepcopy(stats['keylist'])
    keylist.remove('runtyp')
    # Make lines
    lines={}
    for i in range(len(stats['runtag'])):
        if stats['runtyp'][i] not in lines: lines[stats['runtyp'][i]]={key:[] for key in keylist}
        for key in keylist:
            if isinstance(stats[key][i],dict):
                if stats[key][i]['frac']>=1: 
                    nprog=nprog0
                else:
                    nprog=int(stats[key][i]['frac']*float(nprog0))
                lines[stats['runtyp'][i]][key].append('['+(progsym*nprog).ljust(nprog0)+']')
            else: lines[stats['runtyp'][i]][key].append(stats[key][i])
    # Write dictionaries
    for runtyp,fdict in lines.iteritems():
        if fname0 in ['stdout','stderr']: fname=fname0
        elif isinstnace(fname0,str)     : fname=fname0+'{}.txt'.format(runtyp)
        fdict['keylist']=keylist
        mmlio.rwtable('W',fname,fdict,overwrite=True)
        if fname not in ['stdout','stderr']:
            print '    '+fname
    return

####################################################################################################################################
# STATUS OF ALL RUNS
def allruns(runlist=None):
    """
    Return list of status dictionaries for all runs
    """
    # Get runlist
    if isinstance(runlist,str): 
        runlist=[runlist]
    if not isinstance(runlist,list):
        runlist=mmlparam.par_taglist('pysim.simlist','mmlsim')
    # Loop over runlist adding stats
    stats={}
    for runtag in runlist: 
        istat=runstat(runtag)
        if 'keylist' not in stats: stats['keylist']=istat['keylist']
        for key in istat['keylist']:
            if key in stats: stats[key].append(istat[key])
            else           : stats[key]=[istat[key]]
    # Return stats
    return stats

####################################################################################################################################
# STATUS OF A SINGLE RUN
def runstat(runtag,makemeth='buildgal',evolmeth='gadget',alzymethlist=None):
    """
    Return dictionary describing status of run
    """
    import pysim
    simstr=pysim.loadsim(runtag)
    stats={'runtag':simstr['runtag'],'runtyp':simstr['runtyp']}
    keylist=['runtag','runtyp']
    # MAKE
    if os.path.isfile(simstr.fdict[evolmeth]['icgen']['stat_ic']):
        stats['make']={'code':'done','frac':1.,'method':makemeth}
    else:
        stats['make']=stat_make(simstr,makemeth)
    keylist.append('make')
    # EVOLUTION
    evolmod=globals()['stat_'+evolmeth]
    stats['evolve']=evolmod(simstr)
    stats['evolve']['method']=evolmeth
    lastext=stats['evolve']['ext']
    stats['ext' ]=stats['evolve']['ext' ]
    stats['time']='%3.2f' % stats['evolve']['time']
    keylist+=['ext','time','evolve']
    # ANALYSIS
    if alzymethlist is None: alzymethlist=['pilimg','scf','ellipse','barprop_Am2']
    for a in alzymethlist: 
        alzymod=globals()['stat_'+a]
        stats[a]=alzymod(simstr,evolmeth,lastext)
        keylist.append(a)
    # Return stats
    stats['keylist']=keylist
    return stats
    
####################################################################################################################################
# BUILDGAL STATUS
def stat_make(simstr,makemeth):
    """
    Return make status
    """
    if   len(mmlfiles.ls(simstr.finfo[makemeth]['stat_ic']['path']))>0: code='done'
    elif len(mmlfiles.ls(os.path.join(simstr.finfo[makemeth]['output']['path'],'*')))>0: code='endrun'
    elif len(mmlfiles.ls(os.path.join(simstr.finfo[makemeth]['input' ]['path'],'*')))>0: code='setup'
    else: code=''
    stat={'code':code,'frac':float(code=='done'),'method':makemeth}
    return stat

####################################################################################################################################
# GADGET STATUS
def stat_gadget(simstr,partime=False):
    """
    Return GADGET status
    """
    # Handle no progress
    if not os.path.isfile(simstr.finfo['gadget']['param']['path']): stat={'code':'','frac':0.,'ext':'','time':0.}
    # Handle partial progress
    else:
        from pysim.files import gadget
        par=gadget.rw_param('R',simstr.finfo['gadget']['param']['path'])
        time2Gyr=(par['UnitLength_in_cm']/par['UnitVelocity_in_cm_per_s'])/((3.15569e7)*(1.0e9))
        if not partime:
            par['TimeMax']=par['TimeBegin']+gadget.runtyp2runtime(simstr['runtyp'])/time2Gyr
        # Get last snapshot
        snapkey,snaptemp,snaplist=files.listsnap(simstr,'gadget')
        if len(snaplist)==0:
            ext=''
            time=par['TimeBegin']
        else:
            ext=files.get_ext(snaplist[-1],snaptemp)
            time=par['TimeOfFirstSnapshot']+par['TimeBetSnapshot']*float(ext)
        frac=float(time-par['TimeBegin'])/float(par['TimeMax']-par['TimeBegin'])
        time_Gyr=time*time2Gyr
        # Status based on fraction
        if   frac==0: code='setup'
        elif frac <1: code='running'
        elif frac>=1: code='done'
        stat={'code':code,'frac':frac,'ext':ext,'time':time_Gyr}
    # Return 
    return stat
        
####################################################################################################################################
# SCF STATUS
def stat_scf(simstr,evolmeth,lastext=None):
    """
    Returns SCF status
    """
    # Get last snap of run if not provided
    if not lastext:
        module=globals()['stat_'+evolmeth]
        lastext=module(simstr)['ext']
    # Get list of io snapshots
    otemp=simstr.finfo['scf']['osnap']['path']
    olist=mmlfiles.ls(files.add_ext(otemp,'*'))
    itemp=simstr.finfo['scf']['isnap']['path']
    ilist=mmlfiles.ls(files.add_ext(itemp,'*'))
    # Find last io snapshot
    frameno=lambda x: float(max(re.findall(r'\d+',x),key=len))
    iext=None ; i=0
    while True:
        i+=1
        if i>len(ilist): iext='000'
        else           : iext=files.get_ext(ilist[-i],itemp)
        if iext: break
    oext=None ; i=0
    while True:
        i+=1
        if i>len(olist): oext='000'
        else           : oext=files.get_ext(olist[-i],otemp)
        if oext: break
    iextval=frameno(iext)
    oextval=frameno(oext)
    try              : lastextval=float(lastext)
    except ValueError: lastextval=0.
    # Get fractions
    if lastextval==0:
        ifrac=0.
        ofrac=0.
    else:
        ifrac=iextval/lastextval
        ofrac=oextval/lastextval
    # Set code based on frac
    if lastextval==0: code=''
    else:
        if   ofrac>=1   : code='done'
        elif ifrac>=1   : code='run'
        elif ifrac< 1   : code='update'
    # Return stat
    stat={'code':code,'frac':ofrac}
    return stat

####################################################################################################################################
# PLOT STATUS
def stat_plot(simstr,plotmeth,tagstr,evolmeth,lastext=None,plotmodule=None):
    """
    Returns PLOT status
    """
    fext0='{test}'
    # Get parameters
    parfile=mmlparam.par2fname('pysim.simlist',plotmeth,tagstr=tagstr)
    if   plotmeth=='pilimg' : pNbody=True
    elif 'pNbody' in parfile: pNbody=True
    else                    : pNbody=False
    if plotmodule is None:
        plotmodule=files.ftype2snaptyp(evolmeth,pNbody=pNbody)
    if pNbody: parfile=parfile.replace('pysim/simlist','pNbody/'+plotmodule)
    else     : parfile=parfile.replace('pysim/simlist','pysim/'+plotmodule)
    inpar=mmlio.rwdict('R',parfile)
    inpar['tagstr']=tagstr
    # Get last snap of run if not provided
    if not lastext:
        module=globals()['stat_'+evolmeth]
        lastext=module(simstr)['ext']
    # Get list of plots
    fdict=mmlparam.par2file(simstr,'pysim.simlist',plotmeth,plot=True,inpar=inpar,fext=fext0)
    ftemp=fdict['file'].replace('simlist',plotmodule)
    fanim=fdict['anim'].replace('simlist',plotmodule)
    flist=mmlfiles.ls(ftemp.replace(fext0,'*'))
    # Get extension of last snapshot
    ext=None ; i=0
    while True:
        i+=1
        if i>len(flist): ext='000'
        else           : ext=files.get_ext(flist[-i],ftemp)
        if ext: break
    try              : extval=float(ext.split('_')[0])
    except ValueError: extval=0.
    try              : lastextval=float(lastext)
    except ValueError: lastextval=0.
    # Get fraction
    if lastextval==0: frac=0.
    else            : frac=extval/lastextval
    # Set code based on frac
    if lastextval==0: code=''
    else:
        if   frac< 1: code='update'
        elif frac>=1: 
            if os.path.isfile(fanim):
                tsnap=os.path.getmtime(files.add_ext(ftemp,ext))
                tanim=os.path.getmtime(fanim)
                if tanim>tsnap: code='done'
                else          : code='animate'
            else:
                code='animate'
    # Return stat
    stat={'code':code,'frac':frac}
    return stat

####################################################################################################################################
# PILIMG STATUS
def stat_pilimg(simstr,evolmeth,lastext=None):
    """
    Returns PILIMG status
    """
    tagstr='disk_galaxy1_x_y_logm_potCenter'
    return stat_plot(simstr,'pilimg',tagstr,evolmeth,lastext=lastext)

####################################################################################################################################
# ELLIPSE STATUS
def stat_ellipse(simstr,evolmeth,lastext=None):
    """
    Returns ELLIPSE status
    """
    tagstr='numpy_disk_galaxy1_x_y_logmass_hybCenter'
    return stat_plot(simstr,'ellipse',tagstr,evolmeth,lastext=lastext)

####################################################################################################################################
# BARPROP STATUS
def stat_barprop_Am2(simstr,evolmeth,lastext=None):
    """
    Returns Am2 BARPROP status
    """
    tagstr='disk_galaxy1_logr0p01to30_Center_phi_phi0p1_100bins'
    return stat_plot(simstr,'barprop_Am2',tagstr,evolmeth,lastext=lastext,
                     plotmodule='simlist')
