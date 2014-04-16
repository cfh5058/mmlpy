#!/usr/bin/python    
###################################################################################################################################
#
# MEAGAN LANG'S SIMHIST METHODS
#
####################################################################################################################################
#import matplotlib
#matplotlib.use('Agg')
import matplotlib.pyplot as plt
import sys,os,shutil,glob,copy,pprint,scipy,math,itertools
import ImageDraw,ImageFont,Image
import numpy as np
import cPickle as pickle
import pNbody
LIST_MODES_ENC=['m_enc','Ltot_enc','spin','vdyn']
LIST_MODES_CLC=['Boort_calc','epicycle_calc','toomreQ_calc','dvcirdr_calc']
LIST_MODES_ND=['x','y','z','rxyz','rxy','phi','theta',
               'n','m','gal','vr','pot','vomega','vphi','vcir',
               'scfpot','scfpotm2','scfpotm2rel',
               'Ek','Ep','Etot','Lz','Ltot']+LIST_MODES_ENC+LIST_MODES_CLC
DICT_PLOTMETHODS={'1d':['plot1d'],
                  '2d':['ellipse','imagpil','imagshw','contour']}
LIST_METHODS_HST=['1d','2d']
LIST_METHODS_PLT=['plot'+ihst for ihst in LIST_METHODS_HST]
LIST_METHODS_RUN=[]
LIST_METHODS_SNP=LIST_METHODS_HST+LIST_METHODS_PLT
LIST_METHODS=LIST_METHODS_RUN+LIST_METHODS_SNP
DIR_FPAR='/home/langmm/analysis/histparams'
from mmlutils import *
import main as mmlnbody
import simlist,simlyze,simstat,simplot,simfile

def main():
    """
    Provides command line access to HIST methods
    """
    out = mmlnbody.walk(mtype='hist')
    return out

####################################################################################################################################
####################################################################################################################################
# METHODS FOR GETTING LISTS
def dict_ftype():
    """
    Returns a dictionary of files listed by type
    """
    fdict={'shortout':LIST_METHODS_RUN,
           'longout' :LIST_METHODS_SNP}
    return fdict
def dict_fgroup():
    """
    Returns a dictionary of file types listed by group
    """
    fdict={'clean'  :['shortout','longout'],
           'output' :['shortout','longout']}
    return fdict
def list_fpar(method,plot=None):
    """
    Returns a list of supported parameters
    """
    method=mmlpars.mml_pars(method,list=LIST_METHODS)
    plot=mmlpars.mml_pars(plot,default=False,type=bool)
    opardict={'keylist':['ptyp','pgal','cmeth','ctyp','cgal'],
              'ptyp'   :simlist.LIST_PTYPS    ,
              'pgal'   :simlist.LIST_PGALS    ,
              'cmeth'  :simlist.LIST_COMMETHS ,
              'ctyp'   :simlist.LIST_COMPTYPS ,
              'cgal'   :simlist.LIST_COMPGALS ,
              'y'      :str                   }#LIST_MODES_ND         }
    ndim=int(float(method[0]))
    ipardict={}
    for idim in range(1,ndim+1):
        opardict['keylist']+=['x{}'.format(idim),'xscl{}'.format(idim),'xlim{}'.format(idim),'xbin{}'.format(idim)]
        ipardict[   'x{}'.format(idim)]=str #LIST_MODES_ND
        ipardict['xscl{}'.format(idim)]=simlist.LIST_SCALES
        ipardict['xlim{}'.format(idim)]=tuple
        ipardict['xbin{}'.format(idim)]=int
    if plot:
        opardict['keylist']+=['y','yscl','ylim','clrvar','hmeth']
        ipardict['yscl'  ]=simlist.LIST_SCALES
        ipardict['ylim'  ]=tuple
        ipardict['clrvar']=simlist.LIST_DISPVARS
        ipardict['hmeth' ]=simlist.LIST_HISTMETHS
    else: opardict['keylist'].append('y')
    if method=='2d': ipardict['histmodule']=['numpy','pNbody']
    pardict=dict(opardict,**ipardict)
    return pardict
def list_fparDEF(method,plot=None):
    """
    Returns a dictionary with default parameter values
    """
    method=mmlpars.mml_pars(method,list=LIST_METHODS)
    plot=mmlpars.mml_pars(plot,default=False,type=bool)
    opardict={'ptyp' :'visi',
              'pgal' :1,
              'cmeth':'mass',
              'ctyp' :'visi',
              'cgal' :1     ,
              'y'    :'m'   }
    if   method=='1d': ipardict={'x1'   :'r'}
    elif method=='2d': ipardict={'x1'   :'x',
                                 'x2'   :'y'}
    else: raise Exception('Invalid method: {}'.format(method))
    ndim=int(float(method[0]))
    for idim in range(1,ndim+1):
        ipardict['xscl{}'.format(idim)]='lin'
        ipardict['xlim{}'.format(idim)]=(0.,0.)
        ipardict['xbin{}'.format(idim)]=100
    if plot:
        ipardict['yscl'  ]='lin'
        ipardict['ylim'  ]=(0.,0.)
        ipardict['clrvar']='otr'
        ipardict['hmeth' ]='dir'
    if method=='2d': ipardict['histmodule']='numpy'
    pardict=dict(opardict,**ipardict)
    return pardict
def fpar2tag(method,inpar,plot=None):
    """
    Returns a file tag given the set of parameters 
    """
    method=mmlpars.mml_pars(method,list=LIST_METHODS)
    plot=mmlpars.mml_pars(plot,default=False,type=bool)
    outpar=simfile.parsfpar('hist',method,fpar=inpar,plot=plot)
    ndim=int(float(method[0])) ; tagstr=''
    for idim in range(1,ndim+1):
        ilim='xlim{}'.format(idim)
        if ilim in outpar: outpar[ilim+'str']='{}to{}'.format(mmlstring.dec2str(outpar[ilim][0]),mmlstring.dec2str(outpar[ilim][1]))
        tagstr+=str('{x'+str(idim)+'}{xlim'+str(idim)+'str}{xscl'+str(idim)+'}{xbin'+str(idim)+'}_').format(**outpar)
    if plot: 
        if 'ylim' in outpar: outpar['ylimstr']='{}to{}'.format(mmlstring.dec2str(outpar['ylim'][0]),mmlstring.dec2str(outpar['ylim'][1]))
        tagstr+='{y}{ylimstr}{yscl}_'.format(**outpar)
    else:
        tagstr+='{y}_'.format(**outpar)
    tagstr+='{ptyp}{pgal}_{cmeth}{ctyp}{cgal}'.format(**outpar)
    if method=='2d' and outpar['histmodule']!='numpy': tagstr+='_{histmodule}'.format(**outpar)
#    if   method=='1d': tagstr='{x1}{xbin1}{xscl1}{xlim1str}_{y}_{ptyp}{pgal}_{cmeth}{ctyp}{cgal}'.format(**outpar)
#    elif method=='2d': tagstr='{x1}{xbin1}{xscl1}{xlim1str}_{x2}{xbin2}{xscl2}{xlim2str}_{y}_{ptyp}{pgal}_{cmeth}{ctyp}{cgal}'.format(**outpar)
#    else: raise Exception('Invalid method: {}'.format(method))
#    if plot: tagstr+='_{clrvar}{hmeth}'.format(**outpar)
    return tagstr

####################################################################################################################################
####################################################################################################################################
# HIGH LEVEL METHODS REQUIRING MMLSIM

####################################################################################################################################
# METHOD FOR RUNNING DIFFERENT CALC METHODS
def run(simstr,method,**method_kw):
    """
    Provides interface for running different CALC methods
    """
    # Pars input & initialize output
    method=mmlpars.mml_pars(method,list=LIST_METHODS)
    out=None
    # Proceed based on method
    if   method in LIST_METHODS_RUN: mmlio.verbose('Run methods not available.')
    elif method in LIST_METHODS_SNP: simlyze.analyze(simstr,method='hist_'+method,**method_kw)
    else: raise Exception('Invalid method: {}'.format(method))
    # Return output
    return out

####################################################################################################################################
####################################################################################################################################
# FILE CLASSES AND METHODS
                                                                    
####################################################################################################################################
# METHOD TO RETURN DICTIONARY OF SIMCALC FILES FOR CALCULATING SNAPSHOT STATISTICS
def snapfiles(simstr,method,fpar,plot=None,plotmeth=None,**exkw):
    """
    Returns a dictionary containing a snapshot hist directory and file base
    """
    plot=mmlpars.mml_pars(plot,default=False,type=bool)
    if plot:
        plotmeth=mmlpars.mml_pars(plotmeth,type=str)
        exkw['plot']=plot
        exkw['plotstr0']='{}_{}{}'.format(plotmeth,fpar['clrvar'],fpar['hmeth'])
        exkw['animstr0']='anim'+plotmeth
    snfdict=simfile.fpar2file(simstr,'hist',method,fpar,**exkw)
    return snfdict
        
####################################################################################################################################
# METHOD TO RETURN LIST OF FILES
def get_filelist(fdict,keylist):
    """
    Returns a list of relavent HIST files
    """
    # Set constants
    ftypDICT=simfile.dict_ftype(fmeth='hist')
    fgrpDICT=simfile.dict_fgroup(fmeth='hist')
    ftypLIST=ftypDICT.keys()
    fgrpLIST=fgrpDICT.keys()
    # Pars input
    fdict=mmlpars.mml_pars(fdict,type=dict)
    keylist=mmlpars.mml_pars(keylist,type=list)
    # Add files
    filelist={}
    for ilong in ftypLIST['longout']:
        if ilong  in keylist:
            filelist[ilong]=[]
            fpar=list_fpar(ilong)
            if len(fpar.keys())==0: fpar={'none':[None]}
            iterpar=simfile.iter_fpar('hist',ilong,fpar=fpar)
            for ifpar in iterpar:
                filelist[ilong].append(snapfiles(fdict,ilong,ifpar)['file'])
    # Return output
    return filelist

####################################################################################################################################
# METHOD TO RETURN A HISTSNAP DICTIONARY FOR A SNAPSHOT
def get_histsnap(simstr,method,histdict=None,overwrite=None,fulloutput=None,askuser=None,
                 fext=None,fname=None,idgal2=None,ftype=None,ptype=None,snapdict=None,**hist_kw):
    """
    Returns a histsnap dictionary for a snapshot
    """
    # Pars input
    method=mmlpars.mml_pars(method,list=LIST_METHODS_SNP)
    overwrite=mmlpars.mml_pars(overwrite,default=False,type=bool)
    fulloutput=mmlpars.mml_pars(fulloutput,default=False,type=bool)
    askuser=mmlpars.mml_pars(askuser,default=False,type=bool)
    loadkeys=dict(fname=fname,fext=fext,ftype=ftype,ptype=ptype,idgal2=idgal2)
    loadkeys=simstr.get_loadkeys(snapdict=snapdict,noidgal2=True,askuser=askuser,**loadkeys)
    # Get parameters
    fpar=simfile.parsfpar('hist',method,fpar=histdict,askuser=askuser)
    # Select files
    snfdict=snapfiles(simstr,method,fpar,fext=loadkeys['fext'])
    snapfile=loadkeys['fname']
    histfile=snfdict['file']
    if not os.path.isfile(snapfile): raise Exception('Invalid file extension: {}'.format(loadkeys['fext']))
    # Load histsnap if file exists
    if os.path.isfile(histfile) and not overwrite: histdict=loadhist(histfile)
    # Otherwise create it
    else: histdict,snapdict=gethist(simstr,method,fpar,fname=histfile,overwrite=overwrite,
                                    loadkeys=loadkeys,snapdict=snapdict,
                                    fulloutput=True,**hist_kw)
    # Return dictionaries
#    print 'get_histsnap: {}'.format(histdict['y'])
#    print histdict['ydat'].min(),histdict['ydat'].max()
    if fulloutput: return histdict,snapdict
    else         : return histdict
            
####################################################################################################################################
# METHOD TO UNMASK AN ARRAY
def unmask(invar,verbose=False):
    """
    Unmasks elements in a variable
    """
    unmasked=False
    outvar=invar
#    outvar=copy.deepcopy(invar)
    if verbose: mmlio.verbose('Beginning unmasking process...')
    if issubclass(outvar.__class__,dict):
#    if isinstance(outvar,dict) or issubclass(outvar,dict) or isinstance(outvar,mmlclass.mydict_nd):
        for ikey in outvar:
            if isinstance(outvar[ikey],np.ma.masked_array):
                mmlio.verbose('Unmasking key {} in variable'.format(ikey))
                outvar[ikey]=outvar[ikey].data
                unmasked=True
            else:
                outvar[ikey],iunmasked=unmask(outvar[ikey])
                if iunmasked:
                    mmlio.verbose('Unmasked part of key {} in variable'.format(ikey))
                    unmasked=True
    elif issubclass(outvar.__class__,list) or issubclass(outvar.__class__,tuple):
#    elif isinstance(outvar,list) or isinstance(outvar,tuple) or issubclass(outvar,list) or issubclass(outvar,tuple):
        newlist=[]
        for idx in range(len(outvar)):
            if isinstance(outvar[idx],np.ma.masked_array):
                mmlio.verbose('Unmasking element {} of variable'.format(idx))
                ivar=outvar[idx].data
                unmasked=True
            else:
                ivar,iunmasked=unmask(outvar[idx])
                if iunmasked:
                    mmlio.verbose('Unmasked part of element {} in variable'.format(idx))
                    unmasked=True
            newlist.append(ivar)
        if issubclass(outvar.__class__,tuple): outvar=tuple(newlist)
        else                                 : outvar=newlist
    elif isinstance(outvar,np.ndarray): pass
    elif issubclass(outvar.__class__,float     ): pass
    elif issubclass(outvar.__class__,long      ): pass
    elif issubclass(outvar.__class__,int       ): pass
    elif issubclass(outvar.__class__,str       ): pass
    elif issubclass(outvar.__class__,type      ): pass
    elif issubclass(outvar.__class__,bool      ): pass
    elif issubclass(outvar.__class__,pNbody.units.Units): pass
    elif issubclass(outvar.__class__,np.bool_  ): pass
    elif issubclass(outvar.__class__,np.float32): pass
    else:
        mmlio.verbose('Unrecognized object type/class: {}/{}'.format(type(outvar),outvar.__class__))
    if not unmasked:
        if verbose: mmlio.verbose('No masked variables found')
    return outvar,unmasked
    
####################################################################################################################################
# METHOD TO ASK FOR PLOTTING PARAMETERS
def askplot(method):
    """
    Asks user for infor required to create histogram plot
    """
    # Pars input
    method=mmlpars.mml_pars(method,list=LIST_METHODS_HST)
    # Get plot base
    list_pmeths=DICT_PLOTMETHODS[method]
    if len(list_pmeths)==1: pmeth_bas=list_pmeths[0]
    else                  : pmeth_bas=mmlio.askselect('Select a histogram plotting method.',list_pmeths)
    pmeth=pmeth_bas
    # Get multiple
#    pmeth_typ=mmlio.askselect('What should be included in each plot.',['single typ/gal']+LIST_PLOTMULT)
#    if pmeth_typ in LIST_PLOTMULT: pmeth='{}_{}'.format(pmeth_bas,pmeth_typ)
#    else                         : pmeth=pmeth_bas
    # Return plot method
    return pmeth

####################################################################################################################################
# METHOD TO GET HISTDICT
def gethist(simstr,method,fpar,fname=None,overwrite=None,loadkeys=None,snapdict=None,
            fulloutput=None,skipsave=None,histmodule=None,**extra_kw):
    """
    Creates a histdict
    """
    # Pars input
    method=mmlpars.mml_pars(method,list=LIST_METHODS_SNP)
    ndim=int(float(method[0]))
    dims=[str(idim) for idim in range(1,ndim+1)]
    fpar=simfile.parsfpar('hist',method,fpar=fpar)
    fname=mmlpars.mml_pars(fname,default='simhist_{}'.format(method),type=str)
    overwrite=mmlpars.mml_pars(overwrite,default=False,type=bool)
    fulloutput=mmlpars.mml_pars(fulloutput,default=False,type=bool)
    skipsave=mmlpars.mml_pars(skipsave,default=False,type=bool)
    histmodule=mmlpars.mml_pars(histmodule,list=['numpy','pNbody'],default='numpy')
    if not mmlfiles.prep_write(fname,overwrite): return None
    # Load neccessary dictionary
    loadkeys=mmlpars.mml_pars(loadkeys,default={},type=dict)
    loadkeys=simstr.get_loadkeys(snapdict=snapdict,noidgal2=True,**loadkeys)
    if snapdict==None: snapdict=simstr.get_snapdict(**loadkeys)
    # Center accordingly
    comout=snapdict.get_center(centergal=fpar['cgal'],centertyp=fpar['ctyp'],method=fpar['cmeth'],
                                 move=True,vcomflag=True)
    # Get indices
    idxtot=snapdict.get_index(pgal=fpar['pgal'],ptyp=fpar['ptyp'])
    if np.any(idxtot):
        # Get limits and bins
        if fpar['y'] not in LIST_MODES_CLC: fpar['ylim']=snapdict.get_limits(fpar['y'],idxtot=idxtot)
        for idim in dims:
            # Limits
            if fpar['xlim'+idim][0]==fpar['xlim'+idim][1]:
                fpar['xlim'+idim]=snapdict.get_limits(fpar['x'+idim],scale=fpar['xscl'+idim],idxtot=idxtot)
            # Bins
            if   fpar['xscl'+idim]=='lin': fpar['xdat'+idim]=np.linspace(fpar['xlim'+idim][0],fpar['xlim'+idim][1],fpar['xbin'+idim])
            elif fpar['xscl'+idim]=='log': fpar['xdat'+idim]=mmlmath.logspace(fpar['xlim'+idim][0],fpar['xlim'+idim][1],fpar['xbin'+idim],addzero=True)
            else: raise Exception('Invalid scale method: {}'.format(fpar['xscl'+idim]))
        dparam=dict(overwrite=overwrite,fulloutput=True,askuser=False,**loadkeys)
        # Handle dependent variables
        if   fpar['y']=='toomreQ_calc':
            # Oort's constant
            Bfpar=copy.deepcopy(fpar) ; Bfpar['y']='Boort_calc'
            Bdict,snapdict=get_histsnap(simstr,method,histdict=Bfpar,snapdict=snapdict,**dparam)
            B=Bdict['ydat']
            # Radial velocity dispersion
            sfpar=copy.deepcopy(fpar) ; sfpar['y']='sigma_rxy_vR'
            sdict,snapdict=get_histsnap(simstr,method,histdict=sfpar,snapdict=snapdict,**dparam)
            sig_r=sdict['ydat']
            # Surface density
            dfpar=copy.deepcopy(fpar) ; dfpar['y']='surfdens'
            ddict,snapdict=get_histsnap(simstr,method,histdict=dfpar,snapdict=snapdict,**dparam)
            surfden=np.ma.masked_equal(ddict['ydat'],0.)
            # Toomre Q
            G=Bdict['gconst']
            fpar['ydat']=(sig_r*np.abs(B)/(G*surfden)).filled(0.)
            fpar['ynum']=Bdict['ynum']
            fpar['xdat1']=Bdict['xdat1']
        elif fpar['y']=='epicycle_calc':
            Vfpar=copy.deepcopy(fpar) ; Vfpar['y']='wcir'
            Vdict,snapdict=get_histsnap(simstr,method,histdict=Vfpar,snapdict=snapdict,**dparam)
            Rmid=Vdict['xdat1'][:-1]+np.diff(Vdict['xdat1'])/2.
            Nmid=Vdict['ynum'][:-1]+np.diff(Vdict['ynum'])/2.
            Vomega=Vdict['ydat']/Vdict['ynum']
            fpar['ydat']=mmlmath.epicycle(Rmid,Vomega)
            fpar['ynum']=Nmid #Nmid[:-1]+np.diff(Nmid)/2.
            fpar['xdat1']=Rmid
        elif fpar['y']=='dvcirdr_calc':
            Vfpar=copy.deepcopy(fpar) ; Vfpar['y']='vcir'
            Vdict,snapdict=get_histsnap(simstr,method,histdict=Vfpar,snapdict=snapdict,**dparam)
            Rmid=Vdict['xdat1'][:-1]+np.diff(Vdict['xdat1'])/2.
            Vcir=Vdict['ydat']/Vdict['ynum'].astype(float)
            fpar['ydat']=mmlmath.derivative(Vcir,Rmid)
            fpar['ynum']=Vdict['ynum']
            fpar['xdat1']=Vdict['xdat1']
#            print 'gethist: Rlim=({},{})'.format(fpar['xdat1'].min(),fpar['xdat1'].max())
        elif fpar['y']=='Boort_calc':
            Vfpar=copy.deepcopy(fpar) ; Vfpar['y']='vcir'
            Vdict,snapdict=get_histsnap(simstr,method,histdict=Vfpar,snapdict=snapdict,**dparam)
            Rmid=Vdict['xdat1'][:-1]+np.diff(Vdict['xdat1'])/2.
            Vcir=Vdict['ydat']/Vdict['ynum']
            fpar['ydat']=mmlmath.oorts_const(Rmid,Vcir,method='B')
            fpar['ynum']=Vdict['ynum']
            fpar['xdat1']=Vdict['xdat1']
        else:
            # Get data
            w=snapdict.get_var(fpar['y'],idxtot=idxtot)
#            print fpar['ylim']
#            print (w.min(),w.max())
#            mmlio.yorn('?')
            ndim=int(float(method[0])) ; xdat=[] ; bins=[]
            for idim in range(1,ndim+1): 
                xdat.append(snapdict.get_var(fpar['x{}'.format(idim)],idxtot=idxtot))
                bins.append(fpar['xdat{}'.format(idim)])
            # Create histogram data
            if len(xdat)==1: data=xdat[0]
            else           : data=np.vstack(tuple(xdat)).T
            # Calculate necessary dictionary
            if method=='2d' and fpar['histmodule']=='pNbody':
                import pNbody.libutil as libutil
                # Tranlate kws
                parkw={'mode':fpar['y'],'view':fpar['x1']+fpar['x2'],
                       'size':(fpar['xlim1'][1]-fpar['xlim1'][0],fpar['xlim2'][1]-fpar['xlim2'][0]),
                       'shape':(fpar['xbin1'],fpar['xbin2'])}
                # Decipher space
                spaceDict={'pos':[ 'x', 'y', 'z'],
                           'vel':['vx','vy','vz']}
                for ispace in spaceDict.keys():
                    if fpar['x1'] in spaceDict[ispace] and fpar['x2'] in spaceDict[ispace]: parkw['space']=ispace
                if 'space' not in parkw: raise Exception('Cannot identify space for modes {x1} and {x2}'.format(**fpar))
                # Create histogram
                pNparam=libutil.extract_parameters([],parkw,snapdict.defaultparameters)
                histw=snapdict.CombiMap(pNparam)
                pNparam['mode']='m'
                if fpar['y'] in ['n','m']: histn=copy.deepcopy(histw)
                else                     : histn=snapdict.CombiMap(pNparam)
            else:
                histw,binsw=mmlmath.myhist(data,bins=bins,weights=w)
                if fpar['y']=='n': histn=copy.deepcopy(histw)
                else             : histn,binsn=mmlmath.myhist(data,bins=bins)
            # Add data to parameter dictionary
            fpar['ydat']=histw
            fpar['ynum']=histn
        # Check ranges
        if fpar['ydat'].min()==fpar['ydat'].max(): mmlio.yorn('histogram wei null ({},{})'.format(fpar['ydat'].min(),fpar['ydat'].max()))
        if fpar['ynum'].min()==fpar['ynum'].max(): mmlio.yorn('histogram num null ({},{})'.format(fpar['ynum'].min(),fpar['ynum'].max()))
        # Add extra data
        fpar['fname' ]=snapdict.p_name[0]
        fpar['time'  ]=snapdict.atime
        fpar['gconst']=snapdict.get_gconst()
        for igal in snapdict.get_gallist():
            icom,ivcom=snapdict.get_center(centergal=igal,centertyp='all',method='mass',move=False,vcomflag=True)
            fpar[ 'com{}'.format(igal)]=icom
            fpar['vcom{}'.format(igal)]=ivcom
    # Save dictionary
    if not skipsave: savehist(fname,fpar)
    # Return dictionaries
    if fulloutput: return fpar,snapdict
    else         : return fpar
                                     
####################################################################################################################################
# METHOD TO SAVE A SIMULATION HIST FILE
def savehist(fname,histdict):
    """
    Save a histdict to a file
    """
    pickle.dump(histdict,open(fname,'wb'),pickle.HIGHEST_PROTOCOL)
    return

####################################################################################################################################
# METHOD TO LOAD A SIMULATION HIST FILE
def loadhist(fname):
    """
    Loads a histdict from a file
    """
    try:
        histdict=pickle.load(open(fname,'rb'))
        return histdict
    except:
        return None

####################################################################################################################################
# METHOD TO ANIMATE HISTOGRAMS
def animhist(simstr,histmeth,plotmeth=None,histdict=None,fparfile=None,
             verbose=None,overwrite=None,askuser=None,**anim_kw):
    """
    Animates histograms
    """
    # Pars input
    histmeth=mmlpars.mml_pars(histmeth,list=LIST_METHODS_SNP)
    askuser=mmlpars.mml_pars(askuser,default=False,type=bool)
    # Get plot method
    if not isinstance(plotmeth,str) and askuser: plotmeth=askplot(histmeth)
    plotmeth=mmlpars.mml_pars(plotmeth,type=str)
    # Get parameters
    histdict=simfile.parsfpar('hist',histmeth,fpar=histdict,askuser=askuser,tagstr=fparfile,plot=True)
    # Get file names
    fdict_plot=snapfiles(simstr,histmeth,histdict,plot=True,plotmeth=plotmeth)
    plotbase=fdict_plot['file']
    animfile=fdict_plot['anim']
    if not mmlfiles.prep_write(animfile,overwrite=overwrite,askuser=askuser): return None
    if len(glob.glob(plotbase))==0: 
        mmlio.verbose('No frames matching template:')
        print '    '+plotbase
        return None
    # Create animation
    movestr=mmlplot.ffmpeg(plotbase,animfile,overwrite=True,verbose=verbose,
                           rmtemp=True,**anim_kw)
    # Return output
    return movestr
    
####################################################################################################################################
# METHOD TO PLOT HISTOGRAM
def plothist(simstr,histmeth,plotmeth,histdict=None,snapdict=None,loadkeys=None,
             owplot=None,owhist=None,askuser=None,test=None,limdict=None,**input_kw):
    """
    Plots a histogram
    """
    # Pars input
    histmeth=mmlpars.mml_pars(histmeth,list=LIST_METHODS_SNP)
    histdict=mmlpars.mml_pars(histdict,default={},type=dict)
    limdict=mmlpars.mml_pars(limdict,default={},type=dict)
    owplot=mmlpars.mml_pars(owplot,default=False,type=bool)
    owhist=mmlpars.mml_pars(owhist,default=False,type=bool)
    askuser=mmlpars.mml_pars(askuser,default=False,type=bool)
    test=mmlpars.mml_pars(test,default=False,type=bool)
    loadkeys=mmlpars.mml_pars(loadkeys,default={},type=dict)
    loadkeys=simstr.get_loadkeys(snapdict=snapdict,noidgal2=True,askuser=askuser,**loadkeys)
    plot_kw,input_kw=mmlplot.parsopt(input_kw,outextra=True)
    # Get plot method
    if not isinstance(plotmeth,str) and askuser: plotmeth=askplot(histmeth)
    plotmeth=mmlpars.mml_pars(plotmeth,type=str)
    # Set parameters
    if test: owplot=True
    dims=[str(idim+1) for idim in range(int(histmeth[0]))]
    varlist=['y']+['x'+idim for idim in dims]
    limlist=['ylim']+['xlim'+idim for idim in dims]
    datlist=['ydat']+['xdat'+idim for idim in dims]
    # Get parameters
    histdict=simfile.parsfpar('hist',histmeth,fpar=histdict,askuser=askuser,plot=True)
    typList=[histdict['ptyp']]
    galList=[histdict['pgal']]
    if   histdict['clrvar']=='otr': pass
    elif histdict['clrvar']=='typ': typList=simlist.LIST_PTYPS
    elif histdict['clrvar']=='gal': galList=simlist.LIST_PGALS
    elif histdict['clrvar']=='galtyp':
        typList=simlist.LIST_PTYPS
        galList=simlist.LIST_PGALS
    else: raise Exception('Invalid color variable: {}'.format(histdict['clrvar']))
    pmeth=plotmeth
    plot_kw['clrvar']=histdict['clrvar']
    plot_kw['hmeth' ]=histdict['hmeth' ]
    # Get file names
    if not plot_kw['outplot']:
        fdict_plot=snapfiles(simstr,histmeth,histdict,fext=loadkeys['fext'],plot=True,plotmeth=plotmeth)
        plotfile=fdict_plot['file']
        if not mmlfiles.prep_write(plotfile,overwrite=owplot,askuser=askuser): return histdict,limdict,snapdict
        plot_kw['plotfile']=plotfile
    # Loop over parameters
    histlist=[]
    labllist=[]
    for ityp in typList:
        for igal in galList:
            # ith histdict
            ihistdict=copy.deepcopy(histdict)
            ihistdict['ptyp']=ityp
            ihistdict['pgal']=igal
            # Get histdict
            ihistdict,snapdict=get_histsnap(simstr,histmeth,histdict=ihistdict,askuser=askuser,
                                            fulloutput=True,snapdict=snapdict,overwrite=owhist,**loadkeys)
            for ikey in histdict.keys(): 
                if ikey not in ihistdict: ihistdict[ikey]=histdict[ikey]
            if histdict['ylim'][0]!=histdict['ylim'][1]: ihistdict['ylim']=histdict['ylim']
            # Set limits
            if 'ydat' in ihistdict:
                for idim in range(len(varlist)):
                    ivar=varlist[idim] ; ilim=limlist[idim] ; idat=datlist[idim]
#                    print 'plothist: {}={} in ({},{})'.format(ivar,ihistdict[ivar],ihistdict[idat].min(),ihistdict[idat].max())
                    if ilim not in ihistdict or ihistdict[ilim][0]==ihistdict[ilim][1]:
                        if ilim in limdict: ilimdef=limdict[ilim]
                        else: 
                            if ivar=='y': ix=ihistdict['ydat'][ihistdict['ynum']!=0]
                            else        : ix=ihistdict[idat]
                            ilimdef=(ix.min(),ix.max())
                        ihistdict[ilim]=mmlio.askquest('{} ({})?'.format(ilim,ihistdict[ivar]),dtype='tuple',default=ilimdef)
                        if ilim not in limdict: limdict[ilim]=ihistdict[ilim]
                # Append dictionary
                histlist.append(ihistdict)
    # Plot
    if   histmeth=='1d':
        if    pmeth=='plot1d' : out=plothist_plot1d(histlist,overwrite=owplot,**plot_kw)
        else: raise Exception('Invalid 1D plot method: {}'.format(pmeth))
    elif histmeth=='2d':
        if   pmeth=='ellipse': raise Exception('ellipse plot method not currently supported')
        elif pmeth=='imagpil': out=plothist_imagpil(histlist,overwrite=owplot,**plot_kw)
        elif pmeth=='imagshw': out=plothist_imagshw(histlist,overwrite=owplot,**plot_kw)
        elif pmeth=='contour': out=plothist_contour(histlist,overwrite=owplot,**plot_kw)
        else: raise Exception('Invalid 2D plot method: {}'.format(pmeth))
    else: raise Exception('Invalid histogram method: {}'.format(histmeth))
    # Return output/print info
    if plot_kw['outplot']: return out
    else                 : return histlist[0],limdict,snapdict
    
###################################################################################################################################
# METHOD TO PLOT 1D HISTOGRAMS
def plothist_plot1d(histlist,plotfile=None,overwrite=None,clrvar=None,hmeth=None,**plot_kw):
    """
    Plots 1D histograms
    """
    # Set constants
    fnameTAG='plothist_plot1d'
    fparLIST=list_fpar('1d',plot=True) ; fparDEF=list_fparDEF('1d',plot=True)
    clrtyp=mmlplot.get_typcolor(outstr=True)
    clrgal=mmlplot.get_galcolor()
    scldict={'lin':'linear','log':'log'}
    # Pars input
    options=mmlplot.parsopt(plot_kw)
    if isinstance(histlist,dict): histlist=[histlist]
    histlist=mmlpars.mml_pars(histlist,type=list)
    plotfile=mmlpars.mml_pars(plotfile,default=os.path.join(options['plotdir'],fnameTAG+options['plotext']),type=str)
    overwrite=mmlpars.mml_pars(overwrite,default=False,type=bool)
    clrvar=mmlpars.mml_pars(clrvar,default=fparDEF['clrvar'],list=fparLIST['clrvar'])
    hmeth=mmlpars.mml_pars(hmeth,default=fparDEF['hmeth'],list=fparLIST['hmeth'])
    if not mmlfiles.prep_write(plotfile,overwrite): return None
    # Set parameters
    nhist=len(histlist)
    if nhist>1: ncol=2
    else      : ncol=1
    nrow=1
    # Title string
    fstr=os.path.basename(histlist[0]['fname']).split('.')[0]
    tstr='t = {:5.2f}'.format(histlist[0]['time'])
    titstr='{}: {}'.format(fstr,tstr)
    # Initialize figure & axes
    options['pad']=2.0
    options['xlabpos']=(0.5,-0.06)
    options['ylabpos']=(-0.1,0.5)
    options['figsize']=(5.*ncol,5.*nrow)
    fig=plt.figure()
    fig.suptitle(titstr)
    axs=plt.subplot(nrow,ncol,1)
    if nhist>1: leg=plt.subplot(nrow,ncol,2)
    # Loop over add histograms
    for ihist in histlist:
        # Label and color
        ilab='{}{}'.format(ihist['ptyp'],ihist['pgal'])
        if   clrvar=='typ': iclr=mmlplot.get_mmlcolor(clrtyp[ihist['ptyp']])
        elif clrvar=='gal': iclr=mmlplot.get_mmlcolor(clrgal[ihist['pgal']])
        else              : iclr=options['fgclr']
        # Data
        if   hmeth=='dir': iy=np.ma.masked_where(ihist['ynum']==0,ihist['ydat'])
        elif hmeth=='avg': iy=ihist['ydat']/np.ma.masked_where(ihist['ynum']==0,ihist['ynum'])
        elif hmeth=='den': 
            if   ihist['x1']=='rxy' : ivol=np.pi*np.diff(ihist['xdat1']**2)
            elif ihist['x1']=='rxyz': ivol=(4./3.)*np.pi*np.diff(ihist['xdat1']**3)
            else                    : ivol=np.diff(ihist['xdat1'])
            iy=np.ma.masked_where(ihist['ynum']==0,ihist['ydat'])/ivol
        else: raise Exception('Invalid histogram method: {}'.format(hmeth))
        if len(ihist['xdat1'])==len(iy): ix=ihist['xdat1']
        else                           : ix=ihist['xdat1'][:-1]+np.diff(ihist['xdat1'])/2.
        # Plot
        ip=axs.plot(ix,iy,label=ilab,color=iclr)
        if nhist>1: il=leg.plot(ix[0],iy[0],label=ilab,color=iclr)
    # Set scales
    axs.set_xscale(scldict[ihist['xscl1']])
    axs.set_yscale(scldict[ihist['yscl' ]])
    # Set limits
    axs.set_xlim(ihist['xlim1'])
    axs.set_ylim(ihist['ylim' ])
    # Label axes
    axs.set_xlabel(ihist['x1'])
    axs.set_ylabel(ihist['y' ])
    # Create legend
    if nhist>1: 
        leg.legend()
        leg.set_xlim((-1,1))
        leg.set_ylim((-1,1))
    # Set axes properties
    mmlplot.set_figprop(fig,options)
    fig.tight_layout(pad=options['pad'])
    # Display plot
    if options['showplot']: fig.show()
    # Save/output plot
    if options['outplot']: return fig
    else:
        fig.savefig(plotfile,facecolor=options['bgclr'],edgecolor=options['bgclr'])
        return None

###################################################################################################################################
# METHOD TO PLOT 2D HISTOGRAMS AS IMAGPIL OBJECTS
def plothist_imagpil(histlist,plotfile=None,overwrite=None,annotate=None,outdata=None,
                     histlist2=None,residflag=None,clrvar=None,hmeth=None,**plot_kw):
    """
    Plots 2D histograms as imagpil objects
    """
    # Set constants
    fnameTAG='plothist_plot2d_{}'.format('imagpil')
    fparLIST=list_fpar('1d',plot=True) ; fparDEF=list_fparDEF('1d',plot=True)
    scldict={'lin':'linear','log':'log'}
    fntclr=255
    locfct=20
    sizDEF=512
    # Pars input
    options=mmlplot.parsopt(plot_kw)
    if isinstance(histlist,dict): histlist=[histlist]
    if isinstance(histlist2,dict): histlist2=[histlist2]
    histlist=mmlpars.mml_pars(histlist,type=list)
    if histlist2==None: residflag=False
    plotfile=mmlpars.mml_pars(plotfile,default=os.path.join(options['plotdir'],fnameTAG+options['plotext']),type=str)
    overwrite=mmlpars.mml_pars(overwrite,default=False,type=bool)
    annotate=mmlpars.mml_pars(annotate,default=True,type=bool)
    outdata=mmlpars.mml_pars(outdata,default=False,type=bool)
    residflag=mmlpars.mml_pars(residflag,default=False,type=bool)
    clrvar=mmlpars.mml_pars(clrvar,default=fparDEF['clrvar'],list=fparLIST['clrvar'])
    hmeth=mmlpars.mml_pars(hmeth,default=fparDEF['hmeth'],list=fparLIST['hmeth'])
    if not options['outplot'] and not outdata:
        if not mmlfiles.prep_write(plotfile,overwrite): return None
    # Set parameters
    nhist=len(histlist)
    if residflag:
        if len(histlist2)!=nhist: raise Exception('Second set of histograms must have same length as first for residual.')
    # Get sizes
    totdim=int(np.ceil(np.sqrt(nhist)))
    nypix=sizDEF
    nxpix=int(histlist[0]['xbin1']*float(sizDEF)/float(histlist[0]['xbin2']))
#    nxpix=histlist[0]['xbin1']
#    nypix=histlist[0]['xbin2']
    nxpix_tot=totdim*nxpix
    nypix_tot=totdim*nypix
    # Create total image
    piltot=Image.new("RGB",(nxpix_tot,nypix_tot))
    # Loop over histograms
    for idxhist in range(nhist):
        ixpix=int((idxhist%nhist)*nxpix)
        iypix=int(np.floor(idxhist/nhist)*nypix)
        if residflag:
            ihist=resid_hist2d(histlist[idxhist],histlist2[idxhist])
        else:
            ihist=histlist[idxhist]
        # Get mean for log scaled color
        if ihist['yscl']=='log': ihist['ymean']=10.**((np.log10(ihist['ylim'][1])+np.log10(ihist['ylim'][0]))/2.-1.)
        else                   : ihist['ymean']=0.
        # Get view window
        view_window=(ihist['xlim1'][1]-ihist['xlim1'][0],ihist['xlim2'][1]-ihist['xlim2'][0])
        # Get label
        ilab='{ptyp}{pgal}'.format(**ihist)
        # Get colormap
        if   clrvar=='otr': icmap=mmlplot.get_cmap('pnbody_light')
        elif clrvar=='gal': icmap=mmlplot.get_cmap(str(ihist['pgal']),cmaptype='mmlpgal')
        elif clrvar=='typ': icmap=mmlplot.get_cmap(str(ihist['ptyp']),cmaptype='mmlptyp')
        else: raise Exception('Invalid color variable: {}'.format(clrvar))
        # Data
        if   hmeth=='dir': iy=np.ma.masked_where(ihist['ynum']==0,ihist['ydat'])
        elif hmeth=='avg': iy=ihist['ydat']/np.ma.masked_where(ihist['ynum']==0,ihist['ynum'])
        elif hmeth=='den':
            iXX1,iXX2=np.meshgrid(np.diff(ihist['xdat1']),np.diff(ihist['xdat2']),indexing='ij')
            ivol=iXX1*iXX2
            iy=np.ma.masked_where(ihist['ynum']==0,ihist['ydat'])/ivol
        else: raise Exception('Invalid histogram method: {}'.format(ihist['hmeth']))
        isclkw=dict(scale=ihist['yscl'],clrmin=0,clip=True,mn=ihist['ylim'][0],mx=ihist['ylim'][1],cd=ihist['ymean'])
        iyint=mmlplot.map_array(iy,**isclkw).filled(0)
        if outdata: return iy.filled(0),iyint,view_window
        # Image
        iyint_T=np.transpose(iyint)
        ipil=Image.fromstring("P",iyint.shape,iyint_T.tostring())
        # Add colormap
        ipil.putpalette(mmlplot.cmap2palette(icmap))
        # Resize
        ipil=ipil.resize((nxpix,nypix))
        # Annotate
        if annotate:
            isiz=ipil.size
            idraw=ImageDraw.Draw(ipil)
            ifont=ImageFont.load_default()
            xlocpad=isiz[0]/locfct
            ylocpad=isiz[1]/locfct
            # Label specifying time
            tloc=(xlocpad,ylocpad)
            if residflag: tstr="t1 = {:5.2f}, t2 = {:5.2f}".format(ihist['time'],ihist2['time'])
            else        : tstr='t = {:5.2f}'.format(ihist['time'])
            idraw.text(tloc,tstr,fill=fntclr,font=ifont)
            # Label specifying dimensions
            dstr="{:3.3g} x {:3.3g}".format(view_window[0],view_window[1])
            dwid,dhgt=ifont.getsize(dstr)
            dloc=(isiz[0]-xlocpad-dwid,ylocpad)
            idraw.text(dloc,dstr,fill=fntclr,font=ifont)
        # Add to total image
        ipil.convert("RGB")
        piltot.paste(ipil,(ixpix,iypix))
    # Display image
    if options['showplot']: piltot.show()
    # Return output
    if options['outplot']: return piltot
    else                 : piltot.save(plotfile)
    return None
    
###################################################################################################################################
# METHOD TO PLOT 2D HISTOGRAMS
def plothist_plot2d(histlist,method,plotfile=None,overwrite=None,annotate=None,
                    histlist2=None,residflag=None,**plot_kw):
    """
    Plots 2D histograms
    """
    # Set constants
    fnameTAG='plothist_plot2d_{}'.format(method)
    scldict={'lin':'linear','log':'log'}
    # Pars input
    options=mmlplot.parsopt(plot_kw)
    if isinstance(histlist,dict): histlist=[histlist]
    if isinstance(histlist2,dict): histlist2=[histlist2]
    histlist=mmlpars.mml_pars(histlist,type=list)
    if histlist2==None: residflag=False
    plotfile=mmlpars.mml_pars(plotfile,default=os.path.join(options['plotdir'],fnameTAG+options['plotext']),type=str)
    overwrite=mmlpars.mml_pars(overwrite,default=False,type=bool)
    annotate=mmlpars.mml_pars(annotate,default=True,type=bool)
    residflag=mmlpars.mml_pars(residflag,default=False,type=bool)
    if not options['outplot']:
        if not mmlfiles.prep_write(plotfile,overwrite): return None
    # Set parameters
    nhist=len(histlist)
    # Initialize figure
    fig=plt.figure()
    fig.suptitle(os.path.basename(histlist[0]['fname']).split('.')[0])
    axs=[]
    leg=plt.subplot(1,nhist+1,nhist+1)
    # Loop over add histograms
    for idxhist in range(nhist):
        ihist=histlist[idxhist]
        # Axes
        iax=plt.subplot(1,nhist+1,idxhist+1)
        # Label and color
        ilab='{}{}'.format(ihist['ptyp'],ihist['pgal'])
        iax.set_title(ilab)
        # Data
        if   ihist['hmeth']=='dir': iy=np.ma.masked_where(ihist['ynum']==0,ihist['ydat'])
        elif ihist['hmeth']=='avg': iy=ihist['ydat']/np.ma.masked_where(ihist['ynum']==0,ihist['ynum'])
        else: raise Exception('Invalid histogram method: {}'.format(ihist['hmeth']))
        ix1=ihist['xdat1'][:-1]+np.diff(ihist['xdat1'])/2.
        ix2=ihist['xdat2'][:-1]+np.diff(ihist['xdat2'])/2.
        # Plot
        if   method=='imagshw': pass
        elif method=='contour': pass
        else: raise Exception('Invalid method: {}'.format(method))
        # Set scales
        iax.set_xscale(scldict[ihist['xscl1']])
        iax.set_yscale(scldict[ihist['xscl2']])
        # Set limits
        iax.set_xlim(ihist['xlim1'])
        iax.set_ylim(ihist['xlim2'])
        # Label axes
        iax.set_xlabel(ihist['x1'])
        iax.set_ylabel(ihist['x2'])
        # Add axes
        axs.append(iax)
    # Create legend
    if nhist>1: leg.legend()
    # Set axes properties
    mmlplot.set_figprop(fig,options)
    fig.tight_layout()
    # Display plot
    if options['showplot']: fig.show()
    # Save/output plot
    if options['outplot']: return fig
    else:
        fig.savefig(plotfile,facecolor=options['bgclr'],edgecolor=options['bgclr'])
        return None

###################################################################################################################################
###################################################################################################################################
# PROVIDE COMMAND LINE ACCESS
if __name__ == '__main__': main()


## Trash
def trash(x=False):
    if not x:
        if fpar['y'].endswith('_enc'):
            Ofpar=copy.deepcopy(fpar) ; Ofpar['y']=fpar['y'].split('_enc')[0]
            Odict,snapdict=get_histsnap(simstr,method,histdict=Ofpar,snapdict=snapdict,**dparam)
            if Odict['xdat1'][0]==0: Odict['xdat1']=Odict['xdat1'][1:]
            else:
                idxbin0=(np.digitize(snapdict.get_var(fpar['x1'])[idxtot],np.array([0.,Odict['xdat1'][0]]))==1)
                O0=snapdict.get_var(Ofpar['y'])[idxtot][idxbin0].sum()
                N0=np.sum(idxbin0)
                Odict['ydat']=np.array([O0]+list(Odict['ydat']))
                Odict['ynum']=np.array([N0]+list(Odict['ynum']))
            fpar['ydat']=np.cumsum(Odict['ydat'])
            fpar['ynum']=np.cumsum(Odict['ynum'])
            fpar['xdat1']=Odict['xdat1']
        elif fpar['y']=='spin':
            Mfpar=copy.deepcopy(fpar) ; Mfpar['y']='m_enc'
            Jfpar=copy.deepcopy(fpar) ; Jfpar['y']='Ltot_enc'
            Mdict,snapdict=get_histsnap(simstr,method,histdict=Mfpar,snapdict=snapdict,**dparam)
            Jdict,snapdict=get_histsnap(simstr,method,histdict=Jfpar,snapdict=snapdict,**dparam)
            Menc=Mdict['ydat']
            Jenc=Jdict['ydat']
            Renc=Mdict['xdat1']
            G=Mdict['gconst']
    #        G=snapdict.get_gconst()
            Venc=np.sqrt(G*Menc/Renc)
            fpar['ydat']=mmlmath.spin(Jenc,Menc,Venc,Renc)
            fpar['ynum']=Mdict['ynum']
            fpar['xdat1']=Mdict['xdat1']
