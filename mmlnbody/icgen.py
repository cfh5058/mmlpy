#!/usr/bin/python
#ext='_%(snapnum1)03d' % {'snapnum1':snapnum1}
####################################################################################################################################
#
# MEAGAN LANG'S ICGEN METHODS
#
####################################################################################################################################
import sys,os,shutil,glob,copy,pprint,scipy,math
from time import gmtime, strftime
import numpy as np
import matplotlib as mplib
import pNbody
LIST_METHTYP_MKIC=['buildgal','gadget']
LIST_METHODS=['clean']+LIST_METHTYP_MKIC
LIST_RUNTYPS=['galsim','intsim']
LIST_PARMETH=['galmod','mmlsim']+LIST_RUNTYPS
DICT_METHODS_MKIC={'buildgal':['mmlmodel','existing'],
                   'gadget'  :['buildgal','snapshot','restart','interact']}
from mmlutils import *
from mmlastro import mmlprofiles
import main as mmlnbody
import simlist,simcalc,simplot,simfile,mmlbuildgal,mmlgadget,mmlscf,mmlgalfit

def main():
    """
    Provides command line access to ICGEN methods
    """
    out = mmlnbody.walk(mtype='icgen')
    return out

####################################################################################################################################
####################################################################################################################################
# METHODS FOR GETTING LISTS
def dict_ftype():
    """
    Returns a dictionary listing files by types
    """
    fdict={'static':['static']}
#    fdict={'static':['gdic','gdpm','gdmklog','bgic','bgpm','bgmklog']}
    return fdict
def dict_fgroup():
    """
    Returns a dictionary listing file types by group
    """
    fdict={'clean':['static']}
    return fdict
def list_mklogkeys(methtype,method,nohead=None):
    """
    Returns a list of valid mklog keys for a given method (header included)
    """
    # Pars input
    methtype=mmlpars.mml_pars(methtype,list=LIST_METHTYP_MKIC)
    nohead=mmlpars.mml_pars(nohead,default=False,type=bool)
    # Handle header
    if nohead: logkeys=[]
    else     : logkeys=list_mklogkeys_head()
    # Add method specific keys
    logkeys+=list_mklogkeys_meth(methtype,method)
    # Return keys
    return logkeys
def list_mklogkeys_head():
    """
    Returns a list of valid mklog header keys
    """
    logkeys=['runtyp','runtag','addtag','timestamp','makecomp','method','baserun','ntot']
    return logkeys
def list_mklogkeys_meth(mtype,method):
    """
    Returns a list of valid mklog keys for a given method
    """
    mtype=mmlpars.mml_pars(mtype,list=LIST_METHTYP_MKIC)
    method=mmlpars.mml_pars(method,list=DICT_METHODS_MKIC[mtype])
    if   mtype=='buildgal':
        if   method=='mmlmodel': logkeys=list_fpar('galsim')['keylist']
        else: raise Exception('Invalid buildgal method: {}'.format(method))
    elif mtype=='gadget'  :
        if   method=='buildgal': logkeys=list_fpar('galsim')['keylist']+['nproc','rsoft','typsoft','softs']
        elif method=='snapshot': logkeys=['snaprun','snapext','snaptime','resettime',
                                          'change','heatfact','changesoft','softfact','softtype','softs']
        elif method=='restart' : logkeys=['restrun','nproc']
        elif method=='interact': logkeys=list_fpar('intsim')['keylist']+['snapext1','snapext2','snaptime1','snaptime2','idgal2','softs1','softs2']
        else: raise Exception('Invalid buildgal method: {}'.format(method))
    else: raise Exception('Invalid method type: {}'.format(mtype))
    return logkeys

####################################################################################################################################
####################################################################################################################################
# FPAR METHODS
def list_fpar(method):
    """
    Returns a dictionary defining valid parameters
    """
    method=mmlpars.mml_pars(method,list=LIST_PARMETH)
    if   method=='galmod': 
        fpar={'keylist'  :['mvir','concen','halo_prof'],
              'mvir'     :float,
              'concen'   :float,
              'halo_prof':mmlbuildgal.LIST_HALOPROFS}
        for ityp in mmlbuildgal.LIST_PTYPS:
            parlist=['mtot','rscl','rmax','zscl','zmax']
            for ipar in parlist:
                fpar['keylist'].append('{}_{}'.format(ityp,ipar))
                fpar['{}_{}'.format(ityp,ipar)]=float
    elif method=='mmlsim': fpar={'keylist'    :['runtag','runtyp','subtag','memcomp','runcomp','nprocbg','nprocgd','nprocscf',
                                                'inclbgtag','inclgdtag','inclscftag','newstyletag'],
                                 'runtag'     :str                       ,
                                 'runtyp'     :LIST_RUNTYPS              ,
                                 'subtag'     :str                       ,
                                 'memcomp'    :mmlinfo.LIST_COMPUTERS    ,
                                 'runcomp'    :mmlinfo.LIST_COMPUTERS    ,
                                 'nprocbg'    :int                       ,
                                 'nprocgd'    :int                       ,
                                 'nprocscf'   :int                       ,
                                 'inclbgtag'  :bool                      ,
                                 'inclgdtag'  :bool                      ,
                                 'inclscftag' :bool                      ,
                                 'newstyletag':bool                      }
    elif method=='galsim': fpar={'keylist'    :['addtag','model','ntot','mvir','haloprof','concen','fhalo','fdisk','fbulge','fgas','bulgespin','toomreQ'],
                                 'model'      :simfile.fpar_taglist('icgen','galmod'),
                                 'ntot'       :long                      ,
                                 'mvir'       :float                     ,
                                 'haloprof'   :mmlbuildgal.LIST_HALOPROFS,
                                 'concen'     :float                     ,
                                 'fhalo'      :float                     ,
                                 'fdisk'      :float                     ,
                                 'fbulge'     :float                     ,
                                 'fgas'       :float                     ,
                                 'bulgespin'  :float                     ,
                                 'toomreQ'    :float                     ,
                                 'addtag'     :str                       }
    elif method=='intsim': fpar={'keylist'    :['addtag','gal1','gal2','rperi','ecc','incl'],
                                 'gal1'       :simfile.fpar_taglist('icgen','galsim'),
                                 'gal2'       :simfile.fpar_taglist('icgen','galsim'),
                                 'rperi'      :float                     ,
                                 'ecc'        :float                     ,
                                 'incl'       :float                     ,
                                 'addtag'     :str                       }
    else: raise Exception('Invalid parameter method: {}'.format(method))
    return fpar
def list_fparDEF(method):
    """
    Returns a dictionary defining valid parameter defaults
    """
    method=mmlpars.mml_pars(method,list=LIST_PARMETH)
    if   method=='galmod':
        fpar={'mvir':0.0,'concen':0.0,'halo_prof':'LH'}
        for ityp in mmlbuildgal.LIST_PTYPS:
            parlist=['mtot','rscl','rmax','zscl','zmax']
            for ipar in parlist: fpar['{}_{}'.format(ityp,ipar)]=0.0
    elif method=='mmlsim': fpar={'runtag'     :'galsim'        ,
                                 'runtyp'     :'galsim'        ,
                                 'subtag'     :''              ,
                                 'memcomp'    :'bender'        ,
                                 'runcomp'    :'accre'         ,
                                 'nprocbg'    :100             ,
                                 'nprocgd'    :256             ,
                                 'nprocscf'   :50              ,
                                 'inclbgtag'  :False           ,
                                 'inclgdtag'  :False           ,
                                 'inclscftag' :False           ,
                                 'newstyletag':True            }
    elif method=='galsim': fpar={'model'      :'MW'            ,
                                 'fbulge'     :1.0             ,
                                 'fdisk'      :1.0             ,
                                 'fhalo'      :1.0             ,
                                 'haloprof'   :'LH'            ,
                                 'concen'     :11.0            ,
                                 'ntot'       :long(1e7)       ,
                                 'mvir'       :1.0e12          ,
                                 'fgas'       :0.0             ,
                                 'bulgespin'  :0.0             ,
                                 'toomreQ'    :0.0             ,
                                 'addtag'     :''              }
    elif method=='intsim': fpar={'gal1'       :'MW_7M12LH11FG0',
                                 'gal2'       :'MW_7M12LH11FG0',
                                 'rperi'      :0.1             ,
                                 'ecc'        :1.0             ,
                                 'incl'       :0.0             ,
                                 'addtag'     :''              }
    else: raise Exception('Invalid parameter method: {}'.format(method))
    return fpar
def fpar2tag(method,inpar):
    """
    Creates a string tag from parameters
    """
    method=mmlpars.mml_pars(method,list=LIST_PARMETH)
    outpar=simfile.parsfpar('icgen',method,fpar=inpar)
    if   method=='galmod':
        validfile=False
        while not validfile:
            tagstr=mmlio.askquest('Entre a new name for this model:',default='None',dtype='str')
            fname=simfile.fpar2fname('icgen',method,tagstr=tagstr)
            if os.path.isfile(fname): mmlio.verbose('That file already exists...')
            else                    : validfile=True
    elif method=='mmlsim':
        tagstr=outpar['runtag']
    elif method=='galsim':
        # Create short strings
        logm=mmlstring.decimal(np.log10(outpar['mvir']),formstr='.4g')
        gmod=outpar['model'].upper()
        if outpar['fhalo' ]!=1: gmod+='{}H'.format(mmlstring.decimal(outpar['fhalo' ]))
        if outpar['fdisk' ]!=1: gmod+='{}D'.format(mmlstring.decimal(outpar['fdisk' ]))
        if outpar['fbulge']!=1: gmod+='{}B'.format(mmlstring.decimal(outpar['fbulge']))
        logn=mmlstring.decimal(np.log10(outpar['ntot']))
        gpro=outpar['haloprof'].upper()
        conc=mmlstring.decimal(outpar['concen'])
        fgas=mmlstring.decimal(100*outpar['fgas'])
        # Combine strings
        tagstr='{}_{}M{}{}{}FG{}'.format(gmod,logn,logm,gpro,conc,fgas)
        if outpar['toomreQ'  ]!=0: tagstr+='Q_{}'.format(mmlstring.decimal(outpar['toomreQ']))
        if outpar['bulgespin']!=0: tagstr+='BS_{}'.format(mmlstring.decimal(outpar['bulgespin']))
        # Add addtag
        if len(outpar['addtag'])!=0: tagstr+=outpar['addtag']
    elif method=='intsim':
        galsim1=simfile.loadfpar('icgen','galsim',outpar['gal1'])
        galsim2=simfile.loadfpar('icgen','galsim',outpar['gal2'])
        # Create short strings for model
        gmod1 = galsim1['model'].upper()
        if galsim1['fhalo' ]!=1: gmod1+='{}H'.format(mmlstring.decimal(galsim1['fhalo' ]))
        if galsim1['fdisk' ]!=1: gmod1+='{}D'.format(mmlstring.decimal(galsim1['fdisk' ]))
        if galsim1['fbulge']!=1: gmod1+='{}B'.format(mmlstring.decimal(galsim1['fbulge']))
        gmod2 = galsim2['model'].upper()
        if galsim2['fhalo' ]!=1: gmod2+='{}H'.format(mmlstring.decimal(galsim2['fhalo' ]))
        if galsim2['fdisk' ]!=1: gmod2+='{}D'.format(mmlstring.decimal(galsim2['fdisk' ]))
        if galsim2['fbulge']!=1: gmod2+='{}B'.format(mmlstring.decimal(galsim2['fbulge']))
        # Create short strings
        ntot1 = mmlstring.decimal(np.log10(galsim1['ntot']))
        ntot2 = mmlstring.decimal(np.log10(galsim2['ntot']))
        qmas  = mmlstring.decimal(galsim1['mvir']/galsim2['mvir'])
        rper  = mmlstring.decimal(outpar['rperi'])
        ecc   = mmlstring.decimal(outpar['ecc'])
        incl  = mmlstring.decimal(outpar['incl']*180./np.pi) # Angle in degrees  
        # Combine strings
        tagstr='{}_{}_Q{}RP{}E{}I{}'.format(gmod1+ntot1,gmod2+ntot2,qmas,rper,ecc,incl)
        # Add addtag
        if len(outpar['addtag'])!=0: tagstr+=outpar['addtag']
    else: raise Exception('Invalid parameter method: {}'.format(method))
    # Return output
    return tagstr

####################################################################################################################################
####################################################################################################################################
# HIGH LEVEL METHODS REQUIRING MMLSIM

####################################################################################################################################
# METHOD FOR RUNNING DIFFERENT ICGEN METHODS
def run(methtype,**method_kw):
    """
    Provides interface for running different ICGEN operations
    """
    # Pars input
    method=mmlpars.mml_pars(methtype,list=LIST_METHODS)
    # Initialize default output
    out=None
    # Proceed based on method
    if   methtype=='clean':
        raise Exception('Clean disabled for now')
##         if mmlio.yorn('Are you sure you want to clean up static files?'):
##             method_kw['ftype']='clean'
##             method_kw['memcomp']=os.environ['HOSTNAME']
##             simfile.rmfiles('icgen',simstr=simstr,**method_kw)
##         else:
##             out='fail'
    elif methtype in LIST_METHTYP_MKIC: mkout=mkic(methtype=methtype,**method_kw)
    else: raise Exception('Option for method {} needs to be added to icgen.run.'.format(methtype))
    # Return output
    return out

####################################################################################################################################
####################################################################################################################################
# FILE CLASSES AND METHODS

####################################################################################################################################
# METHOD TO CREATE DICTIONARY OF ICGEN FILES
def files(fdict):
    """
    Returns a list of ICGEN files & directories
    """
    # Pars input
    fdict=mmlpars.mml_pars(fdict,type=dict)
    pfix_shrt=fdict['pfix_shrt']
    icdir=fdict['icgen']['dir']
    # Buildgal static files
    fdict['icgen']['bgic']=os.path.join(icdir,pfix_shrt+'static_bgic')
    fdict['icgen']['bgpm']=os.path.join(icdir,pfix_shrt+'static_bgpm')
    fdict['icgen']['bgmklog']=os.path.join(icdir,pfix_shrt+'static_bgmklog')
    # Gadget static files
    fdict['icgen']['gdic']=os.path.join(icdir,pfix_shrt+'static_gdic')
    fdict['icgen']['gdpm']=os.path.join(icdir,pfix_shrt+'static_gdpm')
    fdict['icgen']['gdmklog']=os.path.join(icdir,pfix_shrt+'static_gdmklog')
    # Add subdirs key
    fdict['icgen']['subdirs']=[]
    # Return file dictionary
    return fdict

####################################################################################################################################
# METHOD TO RETURN LIST OF RELAVENT FILES
def get_filelist(fdict,keylist):
    """
    Returns a list of relevant ICGEN files
    """
    # Get groups and lists
    ftypDICT=simfile.dict_ftype(fmeth='icgen')
    fgrpDICT=simfile.dict_fgroup(fmeth='icgen')
    ftypLIST=ftypDICT.keys()
    fgrpLIST=fgrpDICT.keys()
    # Pars input
    fdict=mmlpars.mml_pars(fdict,type=dict)
    keylist=mmlpars.mml_pars(keylist,type=list)
    # Add files
    filelist={}
    if 'static'   in keylist:
        fkeys=['gdic','gdpm','gdmklog','bgic','bgpm','bgmklog']
        filelist['static'  ]=[fdict[ikey] for ikey in fkeys]
    # Return filelist
    return filelist

####################################################################################################################################
####################################################################################################################################
# METHODS TO MAKE IC FILES
def mkic(methtype=None,method=None,runtag=None,**method_kw):
    """
    Creates IC files
    """
    # Set constants
    runLIST=simfile.fpar_taglist('icgen','mmlsim')
    # Pars input
    methtype=mmlpars.mml_pars(methtype,list=LIST_METHTYP_MKIC)
    methLIST=DICT_METHODS_MKIC[methtype]
    if isinstance(method,str):
        method=mmlpars.mml_pars(method,list=methLIST)
    else:
        method=mmlio.askselect('What method should be used to created the GADGET IC files?',methLIST)
    lkeyLIST=list_mklogkeys(methtype,method)
    # Get key for mklog fdict entry
    if   methtype=='buildgal': mklogkey='bgmklog'
    elif methtype=='gadget'  : mklogkey='gdmklog'
    else: raise Exception('Invalid method type: {}'.format(methtype))
    # Initialize mklog
    mklog={}
    if runtag in runLIST:
        simobj=mmlnbody.runtag2simobj(runtag)
        files=simobj.fdict
        mklogfile=files['icgen'][mklogkey]
        if os.path.isfile(mklogfile):
            mklog=rw_makelog('R',mklogfile,methtype,method=method,rmhead=False)
        else:
            mklog['runtag']=simobj['runtag']
            mklog['runtyp']=simobj['runtyp']
            mklog['addtag']=simobj.infodict['addtag']
    for ilkey in lkeyLIST:
        if ilkey not in mklog: mklog[ilkey]=None
    # Update mklog based on methtyp
    if   methtype=='buildgal': mklog=mkbgic(method,'mklog',mklog,**method_kw)
    elif methtype=='gadget'  : mklog=mkgdic(method,'mklog',mklog,**method_kw)
    else: raise Exception('Invalid method type: {}'.format(methtype))
    if isinstance(mklog,str):
        print mklog
        return
    # Generate mmlsim file
    mklog=mkic_mmlsimfile(mklog)
    if isinstance(mklog,str):
        print mklog
        return
    # Get file info & make directories
    files=mmlnbody.runtag2simobj(mklog['runtag']).fdict
    mmlfiles.mkdirs(files['icgen']['dir'])
    fmake=files['icgen'][mklogkey]
    # Create files and add additional info to mklog
    if   methtype=='buildgal': mklog=mkbgic(method,'mkfiles',mklog,**method_kw)
    elif methtype=='gadget'  : mklog=mkgdic(method,'mkfiles',mklog,**method_kw)
    else: raise Exception('Invalid method type: {}'.format(methtype))
    if isinstance(mklog,str):
        print mklog
        return
    # Create make log
    if os.path.isfile(fmake): makelogflag=mmlio.yorn('The specified make log already exists. Overwrite?')
    else                    : makelogflag=True
    if makelogflag:
        rw_makelog('W',fmake,methtype,method=method,meth_kw=mklog,runtag=mklog['runtag'],overwrite=True)
    # Return control
    return
    
####################################################################################################################################
# METHOD TO CREATE A NEW MMLSIM FILE
def mkic_mmlsimfile(mklog):
    """
    Creates a mmlsim file for a new IC & parameter file
    """
    # Set constants
    runtagLIST=simfile.fpar_taglist('icgen',mklog['runtyp'])
    # Alter affected parameters
    isimdct={}
    if isinstance(mklog['baserun'],str):
        irunfile=simfile.fpar2fname('icgen','mmlsim',tagstr=mklog['baserun'])
        if os.path.isfile(irunfile): isimdct=simfile.loadfpar('icgen','mmlsim',irunfile)
    isimdct['runtyp']=mklog['runtyp']
    isimdct['runtag']=mklog['runtag']
    # Ask if file should be overwritten
    frunfile=simfile.fpar2fname('icgen','mmlsim',tagstr=mklog['runtag'])
    if os.path.isfile(frunfile): makerapflag=mmlio.yorn('Specified mmlsim file already exists. Overwrite?')
    else: makerapflag=True
    # Create mmlsim file
    if makerapflag:
        # Save mmlsim file
        fsimdct=simfile.parsfpar('icgen','mmlsim',fpar=isimdct,tagstr=isimdct['runtag'],askuser=True,overwrite=makerapflag)
        # Ask user to check it
        if not mmlio.yorn('Check the mmlsim file {} now! Does it look ok?'.format(frunfile)):
            os.remove(frunfile)
            return 'fail'
        else: 
            mklog['runtag']=fsimdct['runtag']
    # Return mklog
    return mklog

####################################################################################################################################
####################################################################################################################################
# METHODS TO MAKE BUILDGAL IC FILES
def mkbgic(method,runmeth,mklog,**method_kw):
    """
    Creates BUILDGAL IC files
    """
    # Pars input
    method=mmlpars.mml_pars(method,list=DICT_METHODS_MKIC['buildgal'])
    # Update mklog based on method
    if   method=='mmlmodel': mklog=mkbgic_mmlmodel(runmeth,mklog,**method_kw)
    elif method=='existing': mklog=mkbgic_existing(runmeth,mklog,**method_kw)
    else: raise Exception('Option for method {} needs to be added to icgen.mkbgic.'.format(method))
    # Return updated mklog
    return mklog

####################################################################################################################################
# METHOD TO CREATE BUILDGAL ICFILES USING MMLMODEL PREMISE
def mkbgic_mmlmodel(runmeth,mklog,**extra_kw):
    """
    Creates BUILDGAL IC files using the mmlmodel premise
    """
    # Pars input
    runmeth=mmlpars.mml_pars(runmeth,list=['mklog','mkfiles'])
    # Fill in initial mklog params
    if runmeth=='mklog':
        mklog['runtyp']='galsim'
        # Get standard galaxy parameters
        mklog['mkinfo']=simfile.parsfpar('icgen',mklog['runtyp'],fpar=mklog,askuser=True)
        for ikey in mklog['mkinfo'].keys(): mklog[ikey]=mklog['mkinfo'][ikey]
        # Get runtag
        mklog['runtag']=simfile.fpar2tag('icgen',mklog['runtyp'],mklog['mkinfo'])
    # Create files
    elif runmeth=='mkfiles':
        # Get mmlsim object & file names
        simstr=mmlnbody.runtag2simobj(mklog['runtag'])
        fparm=simstr.fdict['icgen']['bgpm']
        fmkbg=simstr.fdict['icgen']['bgmklog']
        # Create direcotires
        mmlfiles.mkdirs(os.path.dirname(fparm))
        # Check file overwrite
        if os.path.isfile(fparm): makeparmflag=mmlio.yorn('Specified parameter file already exists. Overwrite?')
        else                    : makeparmflag=mmlio.yorn('Specified parameter file does not exist. Create?')
        # Create parameter file
        if makeparmflag:
            print fparm
            fparams=mmlbuildgal.mkparam(simstr,overwrite=True,parfile=fparm)
    # Handle error
    else: raise Exception('Invalid run method: {}'.format(runmeth))
    # Return make log
    return mklog
    
####################################################################################################################################
####################################################################################################################################
# METHODS TO MAKE GADGET IC FILES
def mkgdic(method,runmeth,mklog,**method_kw):
    """
    Creates GADGET IC files
    """
    # Set constants
    methLIST=DICT_METHODS_MKIC['gadget']
    lkeyLIST=list_mklogkeys('gadget',method)
    # Pars input
    method=mmlpars.mml_pars(method,list=methLIST)
    # Update mklog based on method
    if   method=='buildgal': mklog=mkgdic_buildgal(runmeth,mklog,**method_kw)
    elif method=='snapshot': mklog=mkgdic_snapshot(runmeth,mklog,**method_kw)
    elif method=='restart' : mklog=mkgdic_restart( runmeth,mklog,**method_kw)
    elif method=='interact': mklog=mkgdic_interact(runmeth,mklog,**method_kw)
    else: raise Exception('Option for method {} needs to be added to icgen.mkgdic.'.format(method))
    # Return updated mklog
    return mklog

####################################################################################################################################
# METHOD TO CREATE A NEW GALAXY IC FILE
def mkgdic_buildgal(runmeth,mklog,**extra_kw):
    """
    Creates a GADGET IC file from BUILDGAL output
    """
    # Set constants
    runlist=simfile.fpar_taglist('icgen','galsim')
    # Pars input
    runmeth=mmlpars.mml_pars(runmeth,list=['mklog','mkfiles'])
    # Fill in initial mklog params
    if runmeth=='mklog':
        # Get runtag if not provided
        if mklog['runtag'] not in runlist: mklog['runtag']=mmlio.askselect('For what run should a GADGET IC file be generated from BUILDGAL output?',runlist)
        # Load mklog for selected run
        simstr=mmlnbody.runtag2simobj(mklog['runtag'])
        mklog0=rw_makelog('R',simstr.fdict['icgen']['bgmklog'],'buildgal')
        mklog.update(mklog0)
        mklog['baserun']=mklog['runtag']
        mklog['mkinfo' ]=simfile.parsfpar('icgen','galsim',fpar=simstr.infodict,tagstr=mklog['runtag'],askuser=True)
        mklog['nproc'  ]=simstr['nprocbg']
    # Create files
    elif runmeth=='mkfiles':
        # Get initial file info
        isimstr=mmlnbody.runtag2simobj(mklog['runtag'])
        isnap=isimstr.fdict['icgen']['bgic']
        iparm=isimstr.fdict['icgen']['bgpm']
        imkbg=isimstr.fdict['icgen']['bgmklog']
        iunit=isimstr.fdict['buildgal']['input']['units']
        # Get file info
        fsimstr=isimstr
        fsnap=fsimstr.fdict['icgen']['gdic']
        fparm=fsimstr.fdict['icgen']['gdpm']
        fmkbg=fsimstr.fdict['icgen']['gdmklog']
        funit=fsimstr.fdict['gadget']['pmdeft']
        # Load file
        nb_bg=isimstr.get_pnbody(fname=isnap,ftype='buildgal')
        # Convert units
        iudct=mmlio.rwdict('R',iunit)
        fudct=mmlio.rwdict('R',funit)
        units={'M':iudct['UnitMass_in_g'           ]/fudct['UnitMass_in_g'           ],
               'L':iudct['UnitLength_in_cm'        ]/fudct['UnitLength_in_cm'        ],
               'V':iudct['UnitVelocity_in_cm_per_s']/fudct['UnitVelocity_in_cm_per_s']}
        # Create GADGET IC file
        nb_gd=pNbody.Nbody(p_name=fsnap,ftype='gadget',status='new',unitsfile=funit,
                           mass=nb_bg.mass*units['M'],
                           pos =nb_bg.pos *units['L'],
                           vel =nb_bg.vel *units['V'],
                           num =nb_bg.num            ,
                           tpe =nb_bg.tpe            )
        nb_gd.atime=0.0
        # Get softenings
        if mklog['rsoft']==None:
            mklog['rsoft']=mmlio.askquest('At what radius (in kpc) should the interparticle separation be estimated?',dtype='float',default=mdict['disk_rscl'])
        if mklog['typsoft']==None or mklog['softs']==None:
            if mmlio.yorn('Print softening info?'):
                mdict=mmlbuildgal.model2param(isimstr.infodict)
                sdict={}
                print 'Determining softenings:'
                for ityp in mmlbuildgal.LIST_PTYPS:
                    cylindric=(ityp=='disk')
                    sdict[ityp]=nb_gd.get_intersep(mklog['rsoft'],ptyp=ityp,cylindric=cylindric)
                    print '    {} = {} kpc (prof)'.format(ityp,mdict[ityp+'_eps']/1000.)
                    print '    {} = {} kpc (part)'.format(ityp,sdict[ityp])
            else: sdict={ityp:0. for ityp in mmlbuildgal.LIST_PTYPS}
        if mklog['typsoft']==None:
            mklog['typsoft']=mmlio.askselect('What particle type would you like to set the softening from?',mmlbuildgal.LIST_PTYPS)
        if mklog['softs']==None:
            mklog['softs']={}
            mklog['softs'][mklog['typsoft']]=mmlio.askquest('What should the softening be? (kpc)',dtype='float',default=sdict[mklog['typsoft']])
            for ityp in simlist.LIST_PTYPBASE:
                if ityp in mmlbuildgal.LIST_PTYPS:
                    mklog['softs'][ityp]=mklog['softs'][mklog['typsoft']]*((mdict[ityp+'_mp']/mdict[mklog['typsoft']+'_mp'])**(1./3.))
                else:
                    mklog['softs'][ityp]=0.
        # Create direcotires
        mmlfiles.mkdirs(os.path.dirname(fsnap))
        # Check file overwrite
        if os.path.isfile(fsnap): makesnapflag=mmlio.yorn('Specified IC file already exists. Overwrite?')
        else                    : makesnapflag=mmlio.yorn('Specified IC file does not exist. Create?')
        if os.path.isfile(fparm): makeparmflag=mmlio.yorn('Specified parameter file already exists. Overwrite?')
        else                    : makeparmflag=mmlio.yorn('Specified parameter file does not exist. Create?')
        # Create IC file
        if makesnapflag:
            print fsnap
            nb_gd.write()
        # Create parameter file
        if makeparmflag:
            print fparm
            fparams=mmlgadget.mkparam(fsimstr,overwrite=True,parstat=funit,parfile=fparm,TimeBegin=0.0,softdict=mklog['softs'])
    # Handle error
    else: raise Exception('Invalid run method: {}'.format(runmeth))
    # Return log info
    return mklog

####################################################################################################################################
# METHOD TO CREATE A NEW IC FILE FROM A PREVIOUS GADGET SNAPSHOT
def mkgdic_snapshot(runmeth,mklog,**extra_kw):
    """
    Creates a GADGET IC file from a GADGET snapshot
    """
    # Set contents
    partypLIST=simlist.LIST_PTYPBASE
    galtypLIST=simlist.LIST_RUNTYPS
    changeLIST=['nothing','disktemp','softenings','nproc']
    changesoftLIST=['soften particles by a factor','change all softenings']
    # Pars input
    runmeth=mmlpars.mml_pars(runmeth,list=['mklog','mkfiles'])
    # Fill in initial mklog params
    if runmeth=='mklog':
        # Ask for simulation type
        if mklog['runtyp']==None: mklog['runtyp']=mmlio.askselect('What type of simulation should the snapshot be taken from?',galtypLIST)
        runtagLIST=simfile.fpar_taglist('icgen',mklog['runtyp'])
        # Ask for simulation name
        if mklog['snaprun']==None: mklog['snaprun']=mmlio.askselect('Which {} simulation should the snapshot be taken from?'.format(mklog['runtyp']),runtagLIST)
        # Get info on selected run
        isimstr=mmlnbody.runtag2simobj(mklog['snaprun'])
        iparm=isimstr.fdict['gadget']['input']['param']
        isofts=mmlgadget.read_softenings(iparm)
        mklog['mkinfo']=isimstr.infodict
        # Ask for snapshot extension
        if mklog['snapext']==None:
            validsnap=False
            while not validsnap:
                mklog['snapext']=mmlio.askquest('What is the extension for the snapshot that should be used?',default='_000ic',dtype='str')
                isnap,isnap_icflag=isimstr.get_snapname(mklog['snapext'],ftype='gadget')
                if os.path.isfile(isnap):
                    validsnap=True
                else:
                    if not mmlio.yorn('That file does not exist. Continue?'): return 'fail'
        # Ask if the time should be reset
        if mklog['resettime']==None: mklog['resettime']=mmlio.yorn('Should the time be reset to 0 in the new files?')
        # Ask if anything should be changed & get corresponding info
        if mklog['change']==None: mklog['change']=mmlio.askselect('Do you want to change anything?',changeLIST)
        if   mklog['change']=='nothing'   : addtagDEF=''
        # Heat the disk
        elif mklog['change']=='disktemp'  :
            if mklog['heatfact']==None:
                mklog['heatfact']=mmlio.askquest('What factor do you want to heat the disk by?',default=1.0,dtype='float')
            addtagDEF='_DH{}'.format(mmlstring.decimal(mklog['heatfact']))
        # Soften particles
        elif mklog['change']=='softenings':
            if mklog['changesoft']==None: mklog['changesoft']=mmlio.askselect('How do you want to change the particle softenings?',changesoftLIST)
            if mklog['changesoft']=='soften particles by a factor':
                if mklog['softtype']==None: mklog['softtype']=mmlio.askselect('What particles do you want to soften?',['all']+partypLIST)
                if mklog['softfact']==None: mklog['softfact']=mmlio.askquest('What factor should {} particles be softened by?'.format(mklog['softtype']),default=1.0,dtype='float')
                if mklog['softtype']=='all': addtagDEF='_SFT{}'.format(mmlstring.decimal(mklog['softfact']))
                else                       : addtagDEF='_SFT{}{}'.format(mklog['softtype'].upper(),mmlstring.decimal(mklog['softfact']))
                if mklog['softtype']=='all': softfacttypes=partypLIST
                else                       : softfacttypes=[mklog['softtype']]
                mklog['softs']=isofts
                for ityp in softfacttypes:
                    mklog['softs'][ityp]*=mklog['softfact']
            elif mklog['changesoft']=='change all softenings'     :
                if mklog['softs']==None:
                    print 'What should be the new softening length (in kpc) for:'
                    mklog['softs']={}
                    for ityp in partypLIST:
                        mklog['softs'][ityp]=mmlio.askquest('    {} particles?'.format(ityp),default=isofts[ityp],dtype='float')
                addtagDEF='_SFTNEW'
            else: raise Exception('Invalid method for changing softenings: {}'.format(mklog['changesoft']))
        # Change the # of processors
        elif change=='nproc'     : raise Exception('Work in progress')
        else: raise Exception('Invalid method for changing a snapshot: {}'.format(mklog['change']))
        # Check addtag
        addtag=mmlio.askquest('What additional tag should be added to the runtag?',dtype='str',default=addtagDEF)
        # Add some things for creating run
        mklog['addtag'  ]=isimstr.infodict['addtag']+addtag
        mklog['runtag'  ]=mklog['snaprun']+addtag
        mklog['baserun' ]=mklog['snaprun']
    # Create files
    elif runmeth=='mkfiles':
        # Get initial file info
        isimstr=mmlnbody.runtag2simobj(mklog['snaprun'])
        isnap,isnap_icflag=isimstr.get_snapname(mklog['snapext'],ftype='gadget')
        iparm=isimstr.fdict['gadget']['input']['param']
        imkbg=isimstr.fdict['icgen']['bgmklog']
        isofts=mmlgadget.read_softenings(iparm)
        if mklog['softs']==None: mklog['softs']=isofts
        # Get file info
        fsimstr=mmlnbody.runtag2simobj(mklog['runtag'])
        fsnap=fsimstr.fdict['icgen']['gdic']
        fparm=fsimstr.fdict['icgen']['gdpm']
        fmkbg=fsimstr.fdict['icgen']['bgmklog']
        # Create new directories
        mmlfiles.mkdirs(os.path.dirname(fsnap))
        # Check if files should be overwritten
        if os.path.isfile(fsnap): makesnapflag=mmlio.yorn('The specified IC file already exists. Overwrite?')
        else                    : makesnapflag=mmlio.yorn('The specified IC file does not exist. Create?')
        if os.path.isfile(fparm): makeparmflag=mmlio.yorn('The specified parameter file already exists. Overwrite?')
        else                    : makeparmflag=mmlio.yorn('The specified parameter file does not exist. Create?')
        # Generate new IC file
        if makesnapflag:
            nb=isimstr.get_pnbody(fname=isnap,ftype='gadget')
            mklog['snaptime']=nb.atime
            if mklog['resettime']: nb.atime=0.0
            if mklog['change']=='heatdisk'  : nb.heatdisk(mklog['heatfact'])
            nb.rename(fsnap)
            nb.update_flags(flag_ic=True)
            nb.write()
            print fsnap
        # Generate new parameter file
        if makeparmflag:
            fparams=mmlgadget.mkparam(fsimstr,overwrite=True,parstat=iparm,parfile=fparm,
                                      TimeBegin=nb.atime,softfact=mklog['softfact'],softdict=mklog['softs'])
            print fparm
        # Copy BUILDGAL mklog
        if makesnapflag or makeparmflag:
            shutil.copy2(imkbg,fmkbg)
    # Handle error
    else: raise Exception('Invalid run method: {}'.format(runmeth))
    # Return log info
    return mklog
    
####################################################################################################################################
# METHOD TO CREATE A NEW INTERACTION IC FILE
def mkgdic_interact(runmeth,mklog,**extra_kw):
    """
    Creates an interaction IC file from two GADGET snapshots for evolved isolated galaxies
    """
    # Set contents
    partypLIST=simlist.LIST_PTYPBASE
    runlist=simfile.fpar_taglist('icgen','intsim')
    # Pars input
    runmeth=mmlpars.mml_pars(runmeth,list=['mklog','mkfiles'])
    # Fill in initial mklog params
    if runmeth=='mklog':
        mklog['runtyp']='intsim'
        # Get standard interaction parameters
        mklog['mkinfo']=simfile.parsfpar('icgen','intsim',fpar=mklog,askuser=True)
        for ikey in mklog['mkinfo'].keys(): mklog[ikey]=mklog['mkinfo'][ikey]
        mklog['runtag']=simfile.fpar2tag('icgen','intsim',mklog['mkinfo'])
        # Get info on selected run
        isimstr1=mmlnbody.runtag2simobj(mklog['gal1'])
        isimstr2=mmlnbody.runtag2simobj(mklog['gal2'])
        iparm1=isimstr1.fdict['gadget']['input']['param']
        iparm2=isimstr2.fdict['gadget']['input']['param']
        # Get snapshot extensions
        if mklog['snapext1']==None:
            validsnapext1=False
            while not validsnapext1:
                mklog['snapext1']=mmlio.askquest('What snapshot should galaxy 1 be loaded from?',default='_000ic',dtype='str')
                isnap1,isnap_icflag1=isimstr1.get_snapname(mklog['snapext1'],ftype='gadget')
                if os.path.isfile(isnap1): validsnapext1=True
                else: 
                    if not mmlio.yorn('That file does not exist. Continue?'): return 'fail'
        if mklog['snapext2']==None:
            validsnapext2=False
            while not validsnapext2:
                mklog['snapext2']=mmlio.askquest('What snapshot should galaxy 2 be loaded from?',default='_000ic',dtype='str')
                isnap2,isnap_icflag2=isimstr2.get_snapname(mklog['snapext2'],ftype='gadget')
                if os.path.isfile(isnap2): validsnapext2=True
                else: 
                    if not mmlio.yorn('That file does not exist. Continue?'): return 'fail'
        mklog['baserun']=mklog['runtag']
    # Create files
    elif runmeth=='mkfiles':
        # Get initial file info
        isimstr1=mmlnbody.runtag2simobj(mklog['gal1'])
        isimstr2=mmlnbody.runtag2simobj(mklog['gal2'])
        isnap1,isnap_icflag1=isimstr1.get_snapname(mklog['snapext1'],ftype='gadget')
        isnap2,isnap_icflag2=isimstr2.get_snapname(mklog['snapext2'],ftype='gadget')
        iparm1=isimstr1.fdict['gadget']['input']['param']
        iparm2=isimstr2.fdict['gadget']['input']['param']
        isofts1=mmlgadget.read_softenings(iparm1)
        isofts2=mmlgadget.read_softenings(iparm2)
        # Get final file info
        fsimstr=mmlnbody.runtag2simobj(mklog['runtag'])
        fsnap=fsimstr.fdict['icgen']['gdic']
        fparm=fsimstr.fdict['icgen']['gdpm']
        fsofts={}
        for ikey in isofts1.keys():
            if isofts1[ikey]==isofts2[ikey]: fsofts[ikey]=isofts1[ikey]
            else: raise Exception('The softening lengths are not equal. Decide how to treat this.')
        mklog['softs1']=isofts1
        mklog['softs2']=isofts2
        # Create new directories
        mmlfiles.mkdirs(os.path.dirname(fsnap))
        # Check if files should be overwritten
        print isnap1
        print isnap2
        if os.path.isfile(fsnap): makesnapflag=mmlio.yorn('The specified IC file already exists. Overwrite?')
        else                    : makesnapflag=mmlio.yorn('The specified IC file does not exist. Create?')
        if os.path.isfile(fparm): makeparmflag=mmlio.yorn('The specified parameter file already exists. Overwrite?')
        else                    : makeparmflag=mmlio.yorn('The specified parameter file does not exist. Create?')
        # Create snapshot
        if makesnapflag:
            # Load galaxies from snapshots
            nb1=isimstr1.get_pnbody(fname=isnap1,ftype='gadget')
            nb2=isimstr2.get_pnbody(fname=isnap2,ftype='gadget')
            # Get make log params
            mklog['idgal2']=nb1.get_idgal2()
            mklog['snaptime1']=nb1.atime
            mklog['snaptime2']=nb2.atime
            # Center galaxies by center of mass
            nb1.cmcenter()
            nb2.cmcenter()
            # Get mass in solar masses
            M1=nb1.mass_tot*nb1.localsystem_of_units.convertionFactorTo(pNbody.units.Unit_Msol)
            M2=nb2.mass_tot*nb2.localsystem_of_units.convertionFactorTo(pNbody.units.Unit_Msol)
            # Get maximum radii in pc
        #    R1=nb1.rxyz().max()*nb1.localsystem_of_units.convertionFactorTo(pNbody.units.Unit_pc)
        #    R2=nb2.rxyz().max()*nb2.localsystem_of_units.convertionFactorTo(pNbody.units.Unit_pc)
            R1=None ; R2=None
            # Put secondary galaxy on specified orbit
            nb2.add_orbit(M1,M2,R1_pc=R1,R2_pc=R2,**mklog)
            # Combine galaxies
            nb=nb1+nb2
            nb.atime=0.0
            nb.rename(fsnap)
            nb.update_flags(flag_ic=True)
            nb.write()
            print fsnap
        # Create new parameter file
        if makeparmflag:
            fparams=mmlgadget.mkparam(fsimstr,overwrite=True,parstat=iparm1,parfile=fparm,
                                      TimeBegin=nb.atime,softdict=fsofts)
            print fparm
    # Handle error
    else: raise Exception('Invalid run method: {}'.format(runmeth))
    # Return log info
    return mklog

####################################################################################################################################
####################################################################################################################################
# METHODS TO READ/WRITE ICGEN FILES

####################################################################################################################################
# METHOD TO READ/WRITE ICGEN MAKELOG FILE
def rw_makelog(rwid,fname,methtype,method=None,meth_kw=None,runtag=None,overwrite=None,rmhead=None):
    """
    Reads/writes information on how an IC file was generated
    """
    # Set constants
    headkeys=list_mklogkeys_head()
    # Pars input
    rwid=mmlpars.mml_pars(rwid,list=['R','W'])
    fname=mmlpars.mml_pars(fname,type=str)
    methtype=mmlpars.mml_pars(methtype,list=LIST_METHTYP_MKIC)
    rmhead=mmlpars.mml_pars(rmhead,default=False,type=bool)
    # Handle read
    if   rwid=='R':
        if not os.path.isfile(fname): raise Exception('Provided makelog file does not exist: {}'.format(fname))
        # Read
        mklog=mmlio.rwdict('R',fname)
        # Fill in missing header values
        for ikey in headkeys:
            if ikey not in mklog: mklog[ikey]=None
        # Get method from header
        method=mmlpars.mml_pars(mklog['method'],list=DICT_METHODS_MKIC[methtype])
        methkeylist=list_mklogkeys(methtype,method)
        # Fill in missing body values
        for ikey in methkeylist:
            if ikey not in mklog: mklog[ikey]=None
        # Remove header keys if necessary
        if rmhead:
            for ikey in headkeys: del mklog[ikey]
        # Return output
        return mklog
    # Handle write
    elif rwid=='W':
        method=mmlpars.mml_pars(method,list=DICT_METHODS_MKIC[methtype])
        methkeylist=list_mklogkeys(methtype,method)
        meth_kw=mmlpars.mml_pars(meth_kw,type=dict)
        simobj=mmlnbody.runtag2simobj(runtag)
        # Create dictionary
        mklog_head=dict(runtag    = simobj['runtag'],
                        runtyp    = simobj['runtyp'],
                        timestamp = strftime("%Y%b%d-%H:%M:%S", gmtime()),
                        makecomp  = mmlinfo.hostname2compid(),
                        method    = method,
                        addtag    = meth_kw['addtag'],
                        baserun   = meth_kw['baserun'],
                        ntot      = simobj.get_ntot()
                        )
        mklog={}
        for ikey in methkeylist:
            if ikey in headkeys:
                if ikey in mklog_head: mklog[ikey]=mklog_head[ikey]
                else: mklog[ikey]=None
            else:
                if ikey in meth_kw: mklog[ikey]=meth_kw[ikey]
                else: mklog[ikey]=None
            if mklog[ikey]==None: mmlio.verbose('({}) Key {} is None.'.format(runtag,ikey))
        # Add keylist
        mklog['keylist']=methkeylist
        # Write dictionary to file
        mmlio.rwdict('W',fname,mklog,overwrite=overwrite)
        # Return
        return
    return

####################################################################################################################################
####################################################################################################################################
# PROVIDE COMMAND LINE ACCESS
if __name__ == '__main__': main()
