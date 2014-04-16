#!/usr/bin/python    
###################################################################################################################################
#
# MEAGAN LANG'S SIMCALC METHODS
#
####################################################################################################################################
import sys,os,shutil,glob,copy,pprint,scipy,math,itertools
import numpy as np
import pNbody
import cPickle
LIST_METHODS_RUN=[]
LIST_METHODS_SNP=['calcsnap','ellipfit','profile']
LIST_METHODS=LIST_METHODS_RUN+LIST_METHODS_SNP
from mmlutils import *
import main as mmlnbody
import simlist,simlyze,simplot,simfile

def main():
    """
    Provides command line access to CALC methods
    """
    out = mmlnbody.walk(mtype='calc')
    return out

####################################################################################################################################
####################################################################################################################################
# METHODS FOR GETTING LISTS
def list_methods2(ftype='all',fsize='all'):
    """
    Returns a list of supported CALC methods for simcalc.run
    """
    runLIST=[]
    snapLIST_sing=['calcsnap']
    snapLIST_mult=['ellipfit','profile']
    if   fsize=='all' : snapLIST=snapLIST_sing+snapLIST_mult
    elif fsize=='sing': snapLIST=snapLIST_sing
    elif fsize=='mult': snapLIST=snapLIST_mult
    else: raise Exception('Invalid file size: {}'.format(fsize))
    if   ftype=='all' : methLIST=runLIST+snapLIST
    elif ftype=='run' : methLIST=runLIST
    elif ftype=='snap': methLIST=snapLIST
    else: raise Exception('Invalid file type: {}'.format(ftype))
#    methLIST=['all','stat','calc','plot']
    return ['all']+methLIST
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
           'input'  :['shortin','longin'],
           'output' :['shortout','longout']}
    return fdict
def list_fpar(method):
    """
    Returns a list of supported parameters
    """
    method=mmlpars.mml_pars(method,list=LIST_METHODS)
    if   method=='calcsnap': pardict={}
    elif method=='profile' : pardict=simprof.list_fpar('profile')
    elif method=='ellipfit': pardict={'type'  :simlist.LIST_PTYPS,
                                      'mode'  :['pos_m','pos_vr','pos_gal','pos_scfpot','pos_scfpotm2',
                                                'EkLz_Lz','EpLz_Lz','EtotLz_Lz','EtotLtot_Ltot'],
                                      'galaxy':simlist.LIST_PGALS,
                                      'view'  :['xy','xz','yz']  }
#                                      'wind'  :tuple             ,
#                                      'size'  :int               }
    else: raise Exception('Invalid method: {}'.format(method))
    return pardict
def list_fparDEF(method):
    """
    Returns a dictionary of supported parameters defaults
    """
    method=mmlpars.mml_pars(method,list=LIST_METHODS)
    if   method=='calcsnap': pardict={}
    elif method=='profile' : pardict=simprof.list_fparDEF('profile')
    elif method=='ellipfit': pardict={'type'  :'visi' ,
                                      'mode'  :'pos_m',
                                      'galaxy':1      ,
                                      'view'  :'xy'   }
#                                      'wind'  :(0.,0.),
#                                      'size'  :512    }
    else: raise Exception('Invalid method: {}'.format(method))
    return pardict
def fpar2tag(method,inpar):
    """
    Returns a file tag given the set of parameters 
    """
    method=mmlpars.mml_pars(method,list=LIST_METHODS)
    outpar=simfile.parsfpar('calc',method,fpar=inpar)
    if   method=='calcsnap': tagstr=''
    elif method=='profile' : tagstr=simprof.fpar2tag('prof',inpar)
    elif method=='ellipfit': 
        if 'wind' in outpar: outpar['windstr']='{}x{}'.format(mmlstring.dec2str(outpar['wind'][0]),mmlstring.dec2str(outpar['wind'][1]))
        tagstr='{type}_{mode}_{galaxy}_{view}'.format(**outpar)
    else: raise Exception('Invalid method: {}'.format(method))
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
    elif method in LIST_METHODS_SNP: simlyze.analyze(simstr,method='calc_'+method,**method_kw)
    else: raise Exception('Invalid method: {}'.format(method))
    # Return output
    return out

####################################################################################################################################
####################################################################################################################################
# FILE CLASSES AND METHODS
                                                                    
####################################################################################################################################
# METHOD TO RETURN DICTIONARY OF SIMCALC FILES FOR CALCULATING SNAPSHOT STATISTICS
def snapfiles(simstr,method,fpar={},**exkw):
    """
    Returns a dictionary containing a snapshot calc directory and file base
    """
    snfdict=simfile.fpar2file(simstr,'calc',method,fpar,**exkw)
    return snfdict
        
####################################################################################################################################
# METHOD TO RETURN CALC FILES
def files(fdict):
    """
    Returns a dictionary of CALC files & directories
    """
    # Pars input
    fdict=mmlpars.mml_pars(fdict,type=dict)
    # Initialize dictionary with directory & default/static files
    fdict['calc']['calcstat']=os.path.join(fdict['rundir'],fdict['pfix_shrt']+'static_calc')
    fdict['calc']['subdirs']=[]
    # Return file dictionary
    return fdict

####################################################################################################################################
# METHOD TO RETURN LIST OF FILES
def get_filelist(fdict,keylist):
    """
    Returns a list of relavent CALC files
    """
    # Set constants
    ftypDICT=simfile.dict_ftype(fmeth='calc')
    fgrpDICT=simfile.dict_fgroup(fmeth='calc')
    ftypLIST=ftypDICT.keys()
    fgrpLIST=fgrpDICT.keys()
    # Pars input
    fdict=mmlpars.mml_pars(fdict,type=dict)
    keylist=mmlpars.mml_pars(keylist,type=list)
    # Add files
    filelist={}
    if 'static'   in keylist:
        filelist['static'  ]=[fdict['calcstat']]
    for ilong in ftypLIST['longout']:
        if ilong  in keylist:
            filelist[ilong]=[]
            fpar=list_fpar(ilong)
            if len(fpar.keys())==0: fpar={'none':[None]}
            iterpar=simfile.iter_fpar('calc',ilong,fpar=fpar)
            for ifpar in iterpar:
                filelist[ilong].append(snapfiles(fdict,ilong,ifpar)['file'])
    # Return output
    return filelist

####################################################################################################################################
# METHOD TO RETURN A CALCSNAP DICTIONARY FOR A SNAPSHOT
def get_calcsnap(simstr,fext=None,fname=None,idgal2=None,ftype=None,ptype=None,
                 snapdict=None,calcdict0=None,overwrite=None,owstat=None,fulloutput=None,**calc_kw):
    """
    Returns a calcsnap dictionary for a snapshot
    """
    # Pars input
    if not isinstance(fext,str) and not isinstance(fname,str):
        raise Exception('[simcalc.get_calcsnap] Must provide either fext or fname')
    if not isinstance(fext ,str): fext=simstr.get_fext(fname=fname,ftype=ftype,ptype=ptype)
    if not isinstance(fname,str): fname=simstr.get_snapname(fext,ftype=ftype,ptype=ptype)[0]
    loadkeys=dict(fname=fname,fext=fext,ftype=ftype,ptype=ptype,idgal2=idgal2)
    overwrite=mmlpars.mml_pars(overwrite,default=False,type=bool)
    owstat=mmlpars.mml_pars(owstat,default=False,type=bool)
    fulloutput=mmlpars.mml_pars(fulloutput,default=False,type=bool)
    # Set scaling for default snapshot
    if   simstr['runtyp'] == 'galsim': adjylim0=False
    elif simstr['runtyp'] == 'intsim': adjylim0=True
    else: raise Exception('Invalid run type: {}'.format(simstr['runtyp']))
    # Create file dictionary & get file extension
    fdict=simstr.fdict
    # Select files
    snapfile=loadkeys['fname']
    calcfile=fdict['calc']['csnapbase']+loadkeys['fext']
    snapfile0=simstr.get_snapfile0(**loadkeys)
    calcfile0=fdict['calc']['calcstat']
    loadkeys0=copy.deepcopy(loadkeys) ; loadkeys0['fname']=snapfile0 ; loadkeys0['idgal2']=None
    if not os.path.isfile(snapfile): raise Exception('[simcalc.get_calcsnap] Invalide file extension: {}'.format(loadkeys['fext']))
    # Create directories
    mmlfiles.mkdirs(os.path.dirname(calcfile0))
    mmlfiles.mkdirs(os.path.dirname(calcfile))
    # Intialize variables
    calcdict=None
    # Load calcsnap if file exists
    if os.path.isfile(calcfile) and not overwrite:
        calcdict=loadcalc(calcfile)
    # Otherwise create it
    if calcdict==None:
        # Load static calcfile if it exists
        if calcdict0==None and os.path.isfile(calcfile0) and not owstat:
            calcdict0=loadcalc(calcfile0)
        # Otherwise create it
        if calcdict0==None:
            snapdict0=simstr.get_pnbody(**loadkeys0)
            calcdict0=snapdict0.calc(calcfile0,adjustylim=adjylim0,overwrite=owstat,**calc_kw)
            del snapdict0
        # Load pNbody object if not provided
        if snapdict==None: snapdict=simstr.get_pnbody(**loadkeys)
        # Create calc dictionary
        calcdict=snapdict.calc(calcfile,overwrite=overwrite,calcdict=copy.deepcopy(calcdict0),**calc_kw)
    # Return dictionaries
    if fulloutput:
        return calcdict,calcdict0,snapdict
    else:
        return calcdict
            
####################################################################################################################################
# METHOD TO RETURN ELLIPSE FITS
def get_ellipfit(simstr,snapdict=None,calcdict=None,calcdict0=None,loadkeys=None,
                 mode=None,ptyp=None,pgal=None,view=None,efitfile=None,
                 askuser=None,verbose=None,overwrite=None,test=None,
                 semimajorlist=None,nsemimajor=None,params0=None,
                 interpmeth=None,nphi=None,errtol=None,dertol=None,miniter=None,maxiter=None,**input_kw):
    """
    Loads/returns ellipse fit to data
    """
    # Pars input
    askuser=mmlpars.mml_pars(askuser,default=False,type=bool)
    verbose=mmlpars.mml_pars(verbose,default=False,type=bool)
    overwrite=mmlpars.mml_pars(overwrite,default=False,type=bool)
    test=mmlpars.mml_pars(test,default=False,type=bool)
    loadkeys=mmlpars.mml_pars(loadkeys,default={},type=dict)
    loadkeys=simstr.get_loadkeys(snapdict=snapdict,noidgal2=True,outextra=False,askuser=askuser,**loadkeys)
    if test: overwrite=True
    # Create lists of things to calc/load
    input_kw['h1DList']=[]
    input_kw['h2DList']=[mode]
    input_kw['typList']=[ptyp]
    input_kw['galList']=[pgal]
    input_kw['rayList']=[view]
    # Get info on data
    input_kw['singflag']=True
    snapdict,calcdict,calckw=simcalc.askcalc(simstr,askuser=askuser,snapdict=snapdict,calcdict=calcdict,**input_kw)
    # Recover keys
    if mode==None and len(calcdict['h2DList'])>0: mode=calcdict['h2DList'][0]
    if ptyp==None and len(calcdict['typList'])>0: ptyp=calcdict['typList'][0]
    if pgal==None and len(calcdict['galList'])>0: pgal=calcdict['galList'][0]
    if view==None and len(calcdict['rayList'])>0: view=calcdict['rayList'][0]
    # Check file
    fpar=dict(type=ptyp,mode=mode,galaxy=pgal,view=view)
    snfdict=snapfiles(simstr,'ellipfit',fpar,fext=loadkeys['fext'],test=test)
    efitfile=mmlpars.mml_pars(efitfile,default=snfdict['file'],type=str)
    if not mmlfiles.prep_write(efitfile,overwrite=overwrite,askuser=askuser): return snapdict,calcdict,calcdict0
    # Get data
    snapdict,calcdict,calcdict0,calckw=getcalc(simstr,askuser=askuser,
                                               snapdict=snapdict,calcdict=calcdict,calcdict0=calcdict0,
                                               loadkeys=loadkeys,histpackage='numpy',**calckw)
    histopt=calcdict['hist2D'][mode][ptyp][pgal][view]
    mat=np.ma.masked_array(histopt['z'],mask=histopt['zmask'])
    imgdim=[list(histopt['xlim']),list(histopt['ylim'])]
    imgdat=mat.filled(0.)
    # Fit ellipses
    elliplist=mmlellipse.fitmultellip(semimajorlist=semimajorlist,nsemimajor=nsemimajor,params0=params0,
                                      imgdat=imgdat,imgdim=imgdim,interpmeth=interpmeth,
                                      nphi=nphi,errtol=errtol,dertol=dertol,miniter=miniter,maxiter=maxiter,verbose=verbose)
    # Save data
    mmlellipse.saveellipse(efitfile,elliplist,overwrite=overwrite,verbose=verbose)
    # Return data dictionaries
    return snapdict,calcdict,calcdict0
            
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
# METHOD TO INTIALIZE CALC DICTIONARY
def initcalc(calcdict=None,snapdict=None,h1DList=None,h2DList=None,typList=None,rayList=None,galList=None):
    """
    Initializes a calc dictionary
    """
    # Set constants
    hist1D_kw=simplot.parsopt_hist1D()
    hist2D_kw=simplot.parsopt_hist2D() ; hist2D_kw['zlimpad']=0.1
    typListDEF=simlist.LIST_PTYPS
    # Pars input
    h1DList=mmlpars.mml_pars(h1DList,default=[],type=list)
    h2DList=mmlpars.mml_pars(h2DList,default=[],type=list)
    typList=mmlpars.mml_pars(typList,default=[],type=list)
    rayList=mmlpars.mml_pars(rayList,default=[],type=list)
    galList=mmlpars.mml_pars(galList,default=[],type=list)
    # Remove empty types from typList
    if snapdict!=None:
        typDict={typListDEF[idx]:idx for idx in range(len(typListDEF))}
        for ityp in typListDEF:
            if ityp not in typList: continue
            if not snapdict.npart[typDict[ityp]]>0: typList.remove(ityp)
    # Initalize dictionary
    if calcdict==None: calcdict={}
    # Add missing keywords
    if 'statDict' not in calcdict: calcdict['statDict']=None
    if 'unitDict' not in calcdict: calcdict['unitDict']=None
    if 'h1DList' not in calcdict: calcdict['h1DList']=h1DList
    if 'h2DList' not in calcdict: calcdict['h2DList']=h2DList
    if 'typList' not in calcdict: calcdict['typList']=typList
    if 'galList' not in calcdict: calcdict['galList']=galList
    if 'rayList' not in calcdict: calcdict['rayList']=rayList
    if 'hist1D' not in calcdict: calcdict['hist1D']=None
    if 'hist2D' not in calcdict: calcdict['hist2D']=None
    # Fill in missing dictionaries
    if calcdict['hist1D']==None: calcdict['hist1D']=mmlclass.mydict_nd(hist1D_kw,h1DList,typList,galList)
    if calcdict['hist2D']==None: calcdict['hist2D']=mmlclass.mydict_nd(hist2D_kw,h2DList,typList,galList,rayList)
    # Fill in missing entries
    for ityp in typList:
        if ityp not in calcdict['typList']: calcdict['typList'].append(ityp)
        for igal in galList:
            if igal not in calcdict['galList']: calcdict['galList'].append(igal)
            for ih1D in h1DList:
                if ih1D not in calcdict['h1DList']: calcdict['h1DList'].append(ih1D)
                # Hist 1D mode
                if ih1D not in calcdict['hist1D']:
                    calcdict['hist1D'][ih1D]={ityp:{igal:copy.deepcopy(hist1D_kw)}}
                # Hist 1D type
                if ityp not in calcdict['hist1D'][ih1D]:
                    calcdict['hist1D'][ih1D][ityp]={igal:copy.deepcopy(hist1D_kw)}
                # Hist 1D gal
                if igal not in calcdict['hist1D'][ih1D][ityp]:
                    calcdict['hist1D'][ih1D][ityp][igal]=copy.deepcopy(hist1D_kw)
            for iray in rayList:
                if iray not in calcdict['rayList']: calcdict['rayList'].append(iray)
                for ih2D in h2DList:
                    if ih2D not in calcdict['h2DList']: calcdict['h2DList'].append(ih2D)
                    # Hist 2D mode
                    if ih2D not in calcdict['hist2D']:
                        calcdict['hist2D'][ih2D]={ityp:{igal:{iray:copy.deepcopy(hist2D_kw)}}}
                    # Hist 2D type
                    if ityp not in calcdict['hist2D'][ih2D]:
                        calcdict['hist2D'][ih2D][ityp]={igal:{iray:copy.deepcopy(hist2D_kw)}}
                    # Hist 2D gal
                    if igal not in calcdict['hist2D'][ih2D][ityp]:
                        calcdict['hist2D'][ih2D][ityp][igal]={iray:copy.deepcopy(hist2D_kw)}
                    # Hist 2D ray
                    if iray not in calcdict['hist2D'][ih2D][ityp][igal]:
                        calcdict['hist2D'][ih2D][ityp][igal][iray]=copy.deepcopy(hist2D_kw)
    # Return intialized calcdict
    return calcdict
    
####################################################################################################################################
# METHOD TO ASK FOR CALCDICT PARAMETERS
def askcalc(simstr,askuser=None,snapdict=None,calcdict=None,loadcalc=None,
            h1DList=None,h2DList=None,typList=None,galList=None,rayList=None,
            shape=None,center=None,centermeth=None,centertyp=None,centergal=None,
            limdict=None,singflag=None,**extrakw):
    """
    Asks user for info required to create a calcdict
    """
    # Set constants
    if   simstr['runtyp']=='intsim': galstrLIST=['both','primary','secondary']
    elif simstr['runtyp']=='galsim': galstrLIST=['both','primary']
    else: raise Exception('Invalid run type: {}'.format(simstr['runtyp']))
    listDictKeys=['h1DList','h2DList','typList','galList','rayList']
    listDictDEF={'h1DList':['rxyz_rho','z_n','vz_n','phi_n'],
                 'h2DList':list_fpar('ellipfit')['mode'],
                 'typList':list_fpar('ellipfit')['type'],
                 'galList':range(len(galstrLIST)),
                 'rayList':list_fpar('ellipfit')['view']}
    limdict=simstr.get_rundefaults()['limdict']
    centermethLIST=simlist.LIST_COMMETHS
    centertypLIST=simlist.LIST_COMPTYPS
    centermethDEF=centermethLIST[0]
    centertypDEF='disk'
    centergalDEF=1
    shapeDEF=512
    centerDEF=True
    # Pars input
    askuser=mmlpars.mml_pars(askuser,default=False,type=bool)
    loadcalc=mmlpars.mml_pars(loadcalc,default=False,type=bool)
    singflag=mmlpars.mml_pars(singflag,default=False,type=bool)
    # Ask for lists
    listDict={'h1DList':h1DList,'h2DList':h2DList,'typList':typList,'galList':galList,'rayList':rayList}
    for ikey in listDictKeys:
        if isinstance(listDict[ikey],str): listDict[ikey]=[listDict[ikey]]
        if ikey=='rayList' and len(listDict['h2DList'])==0:
            listDict['rayList']=[]  ; continue
        if ikey=='galList' and simstr['runtyp']!='intsim':
            listDict['galList']=[1] ; continue
        if not isinstance(listDict[ikey],list):
            if askuser:
                listDict[ikey]=[]
                stopflag=False
                while not stopflag:
                    if ikey=='galList': ians=mmlio.askselect('Select {} list member.'.format(ikey),galstrLIST+['done'])
                    else              : ians=mmlio.askselect('Select {} list member.'.format(ikey),listDictDEF[ikey]+['done'])
                    if ians=='done': stopflag=True
                    else           :
                        if ikey=='galList': listDict[ikey].append(galstrLIST.index(ians))
                        else              : listDict[ikey].append(ians)
                    if singflag: stopflag=True
            else:
                if singflag:
                    listDict[ikey]=[listDictDEF[ikey][0]]
                else:
                    listDict[ikey]=listDictDEF[ikey]
        else:
            rmele={}
            for iele in listDict[ikey]:
                if iele not in listDictDEF[ikey]:
                    if askuser:
                        rmele[iele]=mmlio.askselect('Invalid {} list member {}. Select another or remove:'.format(ikey,iele),listDictDEF[ikey]+['remove'])
                    else:
                        rmele[iele]='remove'
            for iele in rmele.keys():
                listDict[ikey].remove(iele)
                if rmele[iele]!='remove': listDict[ikey].append(rmele[iele])
    # Properties for getting histogram
    if loadcalc:
        propdict={}
    else:
        # Initalize calcdict
        calcdict=initcalc(calcdict,snapdict=snapdict,**listDict)
        # Loop over 2D histogram stuff
        for ih2D in listDict['h2DList']:
            ismode,iwmode=ih2D.split('_')
            for ityp in listDict['typList']:
                for igal in listDict['galList']:
                    for iray in listDict['rayList']:
                        iopt=calcdict['hist2D'][ih2D][ityp][igal][iray]
                        # Default limits
#                        print limdict
#                        mmlio.yorn(ismode)
#                        mmlio.yorn(ityp)
#                        mmlio.yorn(str(igal))
                        if ismode not in limdict: limdict[ismode]={}
                        if ityp not in limdict[ismode]: limdict[ismode][ityp]={}
                        if igal not in limdict[ismode][ityp]: limdict[ismode][ityp][igal]=(0.,0.)
                        if iwmode not in limdict: limdict[iwmode]={}
                        if ityp not in limdict[iwmode]: limdict[iwmode][ityp]={}
                        if igal not in limdict[iwmode][ityp]: limdict[iwmode][ityp][igal]=(0.,0.)
                        # Histogram size
                        if iopt['xlim'][0]==iopt['xlim'][1]:
                            if askuser: iopt['xlim']=mmlio.askquest('{} {} {} {} x variable limits?'.format(ismode,ityp,igal,iray),dtype='tuple',default=limdict[ismode][ityp][igal])
                            else      : iopt['xlim']=limdict[ismode][ityp][igal]
                        if iopt['ylim'][0]==iopt['ylim'][1]:
                            if askuser: iopt['ylim']=mmlio.askquest('{} {} {} {} y variable limits?'.format(ismode,ityp,igal,iray),dtype='tuple',default=limdict[ismode][ityp][igal])
                            else      : iopt['ylim']=limdict[ismode][ityp][igal]
                        # Histogram limits
                        if iopt['zlim'][0]==iopt['zlim'][1]:
                            if askuser: iopt['zlim']=mmlio.askquest('{} {} {} {} z variable limits?'.format(iwmode,ityp,igal,iray),dtype='tuple',default=limdict[iwmode][ityp][igal])
                            else      : iopt['zlim']=limdict[iwmode][ityp][igal]
                        # Histogram dimensions
                        if not isinstance(shape,int) or shape<=0:
                            if askuser: shape=mmlio.askquest('How many bins should the histogram include?',dtype='int',default=shapeDEF)
                            else      : shape=shapeDEF
                        iopt['xnbins']=shape
                        iopt['ynbins']=shape
                        # Reassign options
                        calcdict['hist2D'][ih2D][ityp][igal][iray]=iopt
        # Loop over 1D histogram stuff
        for ih1D in listDict['h1DList']:
            ismode,iwmode=ih1D.split('_')
            for ityp in listDict['typList']:
                for igal in listDict['galList']:
                    iopt=calcdict['hist1D'][ih1D][ityp][igal]
                    # Default limits
                    if ismode not in limdict: limdict[ismode]={}
                    if ityp not in limdict[ismode]: limdict[ismode][ityp]={}
                    if igal not in limdict[ismode][ityp]: limdict[ismode][ityp][igal]=(0.,0.)
                    if iwmode not in limdict: limdict[iwmode]={}
                    if ityp not in limdict[iwmode]: limdict[iwmode][ityp]={}
                    if igal not in limdict[iwmode][ityp]: limdict[iwmode][ityp][igal]=(0.,0.)
                    # Histogram size
                    if iopt['xlim'][0]==iopt['xlim'][1]:
                        if askuser: iopt['xlim']=mmlio.askquest('{} {} {} x variable limits?'.format(ismode,ityp,igal),dtype='tuple',default=limdict[ismode][ityp][igal])
                        else      : iopt['xlim']=limdict[ismode][ityp][igal]
                    # Histogram limits
                    if iopt['ylim'][0]==iopt['ylim'][1]:
                        if askuser: iopt['ylim']=mmlio.askquest('{} {} {} y variable limits?'.format(iwmode,ityp,igal),dtype='tuple',default=limdict[iwmode][ityp][igal])
                        else      : iopt['ylim']=limdict[iwmode][ityp][igal]
                    # Histogram dimensions
                    if not isinstance(shape,int) or shape<=0:
                        if askuser: shape=mmlio.askquest('How many bins should the histogram include?',dtype='int',default=shapeDEF)
                        else      : shape=shapeDEF
                    iopt['xnbins']=shape
                    # Reassign options
                    calcdict['hist1D'][ih1D][ityp][igal]=iopt
        # Ask for center info
        if not isinstance(center,bool):
            if askuser: center=mmlio.yorn('Should the particles be recentered?')
            else      : center=centerDEF
        if center:
            if centermeth not in centermethLIST:
                if askuser: centermeth=mmlio.askselect('How should the particles be centered?',centermethLIST)
                else      : centermeth=centermethDEF
            if centermeth not in ['massbytype','densitybytype']:
                if centertyp not in centertypLIST:
                    if askuser: centertyp=mmlio.askselect('What type of particles should the image be centered on?',centertypLIST)
                    else      : centertyp=centertypDEF
            if centergal not in range(len(galstrLIST)):
                if simstr['runtyp']=='intsim':
                    if askuser:
                        centergalstr=mmlio.askselect('Which galaxy should the image be centered on?',galstrLIST)
                        centergal=galstrLIST.index(centergalstr)
                    else:
                        centergal=centergalDEF
                else:
                    centergal=1
        # Set property dict
        propdict={'center':center,'centermeth':centermeth,'centergal':centergal,'centertyp':centertyp,'shape':shape}
    # Return dictionaries
    kwdict=dict(listDict,**propdict)
    return snapdict,calcdict,dict(kwdict,**extrakw)

####################################################################################################################################
# METHOD TO GET CALCDICT
def getcalc(simstr,askuser=None,snapdict=None,calcdict=None,calcdict0=None,loadcalc=None,loadkeys=None,
            histpackage=None,limdict=None,**extra_kw):
    """
    Creates a calcdict
    """
    # Pars input
    loadcalc=mmlpars.mml_pars(loadcalc,default=False,type=bool)
    loadkeys=mmlpars.mml_pars(loadkeys,default={},type=dict)
    loadkeys=simstr.get_loadkeys(snapdict=snapdict,noidgal2=True,askuser=askuser,**loadkeys)
    askuser=mmlpars.mml_pars(askuser,default=False,type=bool)
    # Load neccessary dictionary
    if loadcalc:
        if calcdict==None:
            calcdict,calcdict0,snapdict=simstr.get_calcdict(snapdict=snapdict,calcdict0=calcdict0,fulloutput=True,
                                                            histpackage=histpackage,**loadkeys)
    else:
        if snapdict==None:
            snapdict=simstr.get_snapdict(**loadkeys)
    # Get keywords
    snapdict,calcdict,kwdict=askcalc(simstr,askuser=askuser,snapdict=snapdict,calcdict=calcdict,loadcalc=loadcalc,limdict=limdict,**extra_kw)
    # Calculate necessary dictionary
    if not loadcalc:
        calcdict=snapdict.calc(calcdict=calcdict,skipsave=True,histpackage=histpackage,**kwdict)
    # Return dictionaries
    return snapdict,calcdict,calcdict0,kwdict
                                     
####################################################################################################################################
# METHOD TO LOAD A SIMULATION CALC FILE
def loadcalc(fname):
    """
    Loads a calcdict from a file
    """
#    import mmlutils.mmlclass as mmlclass
    try:
        fid=open(fname,'r')
        calcdict=cPickle.load(fid)
        return calcdict
    except:
        return None

###################################################################################################################################
###################################################################################################################################
# PROVIDE COMMAND LINE ACCESS
if __name__ == '__main__': main()
