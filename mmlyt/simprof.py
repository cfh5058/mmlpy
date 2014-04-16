#!/usr/bin/python
import sys,os,shutil,glob,copy,pprint,scipy,math
import numpy as np
import matplotlib as mplib
import pNbody
LIST_METHODS=['profile']
from mmlutils import *
import main as mmlnbody
import simlyze,simstat,simplot,simfile,mmlbuildgal,mmlgadget,mmlscf,mmlgalfit,mml_pNbody

####################################################################################################################################
####################################################################################################################################
# METHODS FOR GETTING LISTS
def list_fpar(method):
    """
    Returns a list of supported parameters
    """
    method=mmlpars.mml_pars(method,list=LIST_METHODS)
    if method=='profile': pardict={'space' :['r','R'],
                                   'scale' :['lin','log'],
                                   'mode'  :['M','Menc','spin','omega','Vcirc'],
                                   'galaxy':[1,2]}
    else: raise Exception('Invalid method: {}'.format(method))
    return pardict
def fpar2tag(method,inpar):
    """
    Returns a file tag given the set of parameters
    """
    method=mmlpars.mml_pars(method,list=LIST_METHODS)
    outpar=simfile.parsfpar('prof',method,fpar=inpar)
    if method=='profile': tagstr='{space}{scale}_{mode}_{galaxy}'.format(**outpar)
    else: raise Exception('Invalid method: {}'.format(method))
    return tagstr

####################################################################################################################################
####################################################################################################################################
# METHOD FOR RUNNING DIFFERENT PROFILES

####################################################################################################################################
# METHOD FOR RUNNING DIFFERENT CALC METHODS
def run(simstr,method,**method_kw):
    """
    Provides interface for running different PROF methods
    """
    # Pars input & initialize output
    method=mmlpars.mml_pars(method,list=LIST_METHODS)
    out=None
    # Proceed based on method
    if   method in LIST_METHODS_RUN: mmlio.verbose('Run methods not available.')
    elif method in LIST_METHODS_SNP: simlyze.analyze(simstr,method='prof_'+method,**method_kw)
    else: raise Exception('Invalid method: {}'.format(method))
    # Return output
    return out

####################################################################################################################################
# METHOD TO READ/WRITE PROFILE FILES
def rwprofs(method,profspace,proffile=None,profdict=None,addprofdict=None,overwrite=None):
    """
    Read/writes/appends entries to profile
    """
    # Pars input
    if isinstance(method,str): method=method.upper()
    method=mmlpars.mml_pars(method,list=['R','W','A'])
    profspace=mmlpars.mml_pars(profspace,list=list_fpar('profile')['space'])
    overwrite=mmlpars.mml_pars(overwrite,default=False,type=bool)
    # Set keys
    keylist=[profspace]+mmlgadget.LIST_PTYPS
    # Read
    if method=='R':
        proffile=mmlpars.mml_pars(proffile,type=str)
        if os.path.isfile(proffile):
            profdict=mmlio.rwtable('R',proffile)
        else:
            profdict={ikey:[] for ikey in keylist}
    # Append
    elif method=='A':
        addprofdict=mmlpars.mml_pars(addprofdict,type=dict)
        nadd=len(addprofdict[profspace])
        if profdict==None: profdict=rwprofs('R',profspace,proffile=proffile)
        added=False
        for idx in range(nadd):
            if addprofdict[profspace][idx] not in profdict[profspace] or overwrite:
                added=True
                if addprofdict[profspace][idx] in profdict[profspace]:
                    idx0=profdict[profspace].index(addprofdict[profspace][idx])
                    for ikey in keylist:
                        profdict[ikey][idx0]=addprofdict[ikey][idx]
                else:
                    for ikey in keylist:
                        profdict[ikey].append(addprofdict[ikey][idx])
        if added: profdict=rwprofs('W',profspace,proffile=proffile,profdict=profdict)
    # Write
    elif method=='W':
        proffile=mmlpars.mml_pars(proffile,type=str)
        profdict=mmlpars.mml_pars(profdict,type=dict)
        if not os.path.isfile(proffile):
            mmlfiles.mkdirs(os.path.dirname(proffile))
        profdict['keylist']=keylist
        mmlio.rwtable('W',proffile,profdict,overwrite=True)
    # Error
    else: raise Exception('Invalid method: {}'.format(method))
    # Return output
    return profdict

####################################################################################################################################
# METHOD TO RETRUN PROFILE LIMITS
def get_proflim(simstr,galaxy=None):
    """
    Returns profile limits
    """
    # Set constants
    fparlist=list_fpar('profile')
    modlist=fparlist['mode' ]
    spalist=fparlist['space']
    # Pars input
    galaxy=mmlpars.mml_pars(galaxy,default=fparlist['galaxy'][0],list=fparlist['galaxy'])
    # Initalize dictionary
    limdict={}
    # Radial limits
    softdict=simstr.get_softenings()
    sft=np.array(softdict.values())
    Smin=sft[sft>0].min()
    Rvir=simstr.get_Rvir(galaxy=galaxy)
    Rmin=10.*Smin
    Rmax=Rvir/10.
    limdict['r'    ]=(Rmin,Rmax)
    limdict['R'    ]=(Rmin,Rmax)
    # Mass limits
    Mvir=simstr.get_Rvir(galaxy=galaxy)
    Mmin=Mvir/100.
    Mmax=Mvir
    Mscl='log'
    limdict['M'    ]=(Mmin,Mmax,Mscl)
    limdict['Menc' ]=(Mmin,Mmax,Mscl)
    # Spin
    limdict['spin' ]=(0.0,10.0,'lin')
    # Omega
    limdict['omega']=(0.0,5.0,'lin')
    # Velocity limits
    limdict['Vcirc']=(0.0,200.0,'lin')
    # Check that all modes and spaces are accounted for
    modmiss=[]
    for imod in modlist:
        if imod not in limdict: modmiss.append(imod)
    spamiss=[]
    for ispa in spalist:
        if ispa not in limdict: spamiss.append(ispa)
    if len(modmiss)+len(spamiss)>0:
        raise Exception('Missing mode(s): {} and space(s): {}'.format(modmiss,spamiss))
    # Return dictionary
    return limdict

####################################################################################################################################
# METHOD TO COMPUTE PROFILE
def get_profile(simstr,space=None,scale=None,mode=None,galaxy=None,
                fext=None,fname=None,idgal2=None,ftype=None,ptype=None,snapdict=None,
                centermeth=None,nbin=None,
                overwrite=None,fulloutput=None,test=None,**calc_kw):
    """
    Generic wrapper for returning profiles 
    """
    # Set constats
    fpar=list_fpar('profile')
    proflim=get_proflim(simstr,galaxy=galaxy)
    # Pars input
    if not isinstance(fext,str) and not isinstance(fname,str):
        raise Exception('Must provide either fext or fname')
    if not isinstance(fext ,str): fext=simstr.get_fext(fname=fname,ftype=ftype,ptype=ptype)
    if not isinstance(fname,str): fname=simstr.get_snapname(fext,ftype=ftype,ptype=ptype)[0]
    loadkeys=dict(fname=fname,fext=fext,ftype=ftype,ptype=ptype,idgal2=idgal2)
    space=mmlpars.mml_pars(space,default=fpar['space'][0],list=fpar['space'])
    scale=mmlpars.mml_pars(scale,default=fpar['scale'][0],list=fpar['scale'])
    mode=mmlpars.mml_pars(mode,default=fpar['mode'][0],list=fpar['mode'])
    galaxy=mmlpars.mml_pars(galaxy,default=fpar['galaxy'][0],list=fpar['galaxy'])
    nbin=mmlpars.mml_pars(nbin,default=100,type=int,min=2)
    if simstr['runtyp']=='galaxy': galaxy=1
    overwrite=mmlpars.mml_pars(overwrite,default=False,type=bool)
    # Select files
    snapfile=loadkeys['fname']
    proffile=simfile.fpar2file(simstr,'prof','profile',fext=loadkeys['fext'],test=test,
                               space=space,scale=scale,mode=mode,galaxy=galaxy)['file']
    if not os.path.isfile(snapfile): raise Exception('Invalid file extension: {}'.format(loadkeys['fext']))
    # Load profile if it exists
    if os.path.isfile(proffile) and not overwrite:
        profdict=rwprofs('R',space,proffile)
    # Otherwise create it
    else:
        # Load pNbody object if not provided
        if snapdict==None: snapdict=simstr.get_pnbody(**loadkeys)
        # Create profile
        Rmin,Rmax=proflim[space]
#        sft=np.array(snapdict.get_softenings())
#        Rmin=10.*sft[sft>0].min()
#        Rmax=simstr.get_Rvir(galaxy=galaxy)/10.
        if   scale=='lin': Rlist=np.linspace(Rmin,Rmax,nbin)
        elif scale=='log': Rlist=np.logspace(np.log10(Rmin),np.log10(Rmax),nbin)
        else: raise Exception('Invalid scale: {}'.format(scale))
        Rlist=list(Rlist)
        profdict=snapdict.get_profile(mode,Rlist,space=space,galaxy=galaxy,centermeth=centermeth)
        # Save the dictionary to file
        profdict=rwprofs('W',space,proffile,profdict,overwrite=overwrite)
    # Return dictionaries
    if fulloutput:
        return profdict,snapdict
    else:
        return profdict

####################################################################################################################################
# METHOD TO PLOT PROFILES
def plot_profile(simstr,space=None,scale=None,mode=None,galaxy=None,fext=None,
                 overwrite=None,verbose=None,askuser=None,test=None,**prof_kw):
    """
    Plots profile
    """
    # Set constats
    fpar=list_fpar('profile')
    proflim=get_proflim(simstr,galaxy=galaxy)
    # Pars input
    space=mmlpars.mml_pars(space,default=fpar['space'][0],list=fpar['space'])
    scale=mmlpars.mml_pars(scale,default=fpar['scale'][0],list=fpar['scale'])
    mode=mmlpars.mml_pars(mode,default=fpar['mode'][0],list=fpar['mode'])
    galaxy=mmlpars.mml_pars(galaxy,default=fpar['galaxy'][0],list=fpar['galaxy'])
    fext=mmlpars.mml_pars(fext,default='_000ic',type=str)
    verbose=mmlpars.mml_pars(verbose,default=False,type=bool)
    askuser=mmlpars.mml_pars(askuser,default=False,type=bool)
    options=mmlplot.parsopt(prof_kw)
    # Get plot file name
    plotfile=simfile.fpar2file(simstr,'prof','profile',fext=fext,test=test,plot=True,
                               space=space,scale=scale,mode=mode,galaxy=galaxy)['file']
    if not mmlfiles.prep_write(plotfile,overwrite=overwrite,askuser=askuser): return
    # Get profile
    prof_kw['fulloutput']=False
    profdict=get_profile(simstr,space=space,scale=scale,mode=mode,galaxy=galaxy,fext=fext,**prof_kw)
    # Get sorted index
    idxsrt=np.array(sorted(range(len(profdict[space])),key=lambda k: profdict[space][k]))
    # Set limits
    xlim=proflim[space]
    ylim=proflim[mode ][0:2]
    xscl=scale
    yscl=proflim[mode ][2]
    if xscl=='lin': xscl='linear'
    if yscl=='lin': yscl='linear'
    # Create figure and plot
    fig=mplib.pyplot.figure()
    axs=mplib.pyplot.subplot(1,1,1)
    # Add plots
    typelist=mmlgadget.LIST_PTYPS
    typecolr=mmlplot.get_typcolor(outstr=True)
    x=np.array(profdict[space])[idxsrt]
    for itype in typelist:
        iy=np.array(profdict[itype])[idxsrt]
        if iy.min()==iy.max(): continue
        iclr=mmlplot.get_mmlcolor(typecolr[itype])
        axs.plot(x,iy,color=iclr,label=itype)
    # Add legend
    mplib.pyplot.legend()
    # Set labels and limits
    axs.set_xlabel(space)
    axs.set_ylabel(mode )
    axs.set_xlim(xlim)
    axs.set_ylim(ylim)
    axs.set_xscale(xscl)
    axs.set_yscale(yscl)
    # Set axes properties
    mmlplot.set_figprop(fig,options)
    # Display plot
    if options['showplot']: fig.show()
    # Save plot
    if options['outplot']: return fig
    else:
        fig.savefig(plotfile,facecolor=options['bgclr'],edgecolor=options['bgclr'])
        mplib.pyplot.close(fig)
        if verbose:
            mmlio.verbose('Created profile plot:')
            print '    '+plotfile
        return None

if __name__ == '__main__':
    main()
