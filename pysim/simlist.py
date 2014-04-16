#!/usr/bin/python
####################################################################################################################################
#
# MEAGAN LANG'S SIMFILE METHODS
#
####################################################################################################################################

import numpy as np
from mmlutils import *
from pysim import files

####################################################################################################################################
####################################################################################################################################
# METHODS FOR GETTING LISTS
LIST_PROGMAKE=['buildgal']
LIST_PROGEVOL=['gadget']
LIST_PROGPOST=['scf']
LIST_PROG=LIST_PROGMAKE+LIST_PROGEVOL+LIST_PROGPOST

LIST_RUNTYPS=['galsim','intsim','cossim']
LIST_YORN=['Y','N']
LIST_PTYPS=['gas','halo','disk','bulge','star']
LIST_PGALS=['galaxy{}'.format(i) for i in range(1,3)]
LIST_COMMETHS=['mass','dens','pot','hyb','com']#'dens_avg','pot_avg']
LIST_HISTMETHS=['dir','avg','den']
LIST_DISPVARS=['otr','gal','typ','galtyp']
LIST_SCALES=['lin','log','symlog']
DICT_PGALS={0:'All',1:'Primary',2:'Secondary'}
LIST_HALOPROFS=['LH','NF','IS']

LIST_METHTYPS=['general']+LIST_PROG
DICT_METHODS={}

####################################################################################################################################
####################################################################################################################################
# PARAMETER PARSING OPTIONS
def listpar(partyp,default=None,**exkw):
    plot=exkw.get('plot',False)
    partyp=mmlpars.mml_pars(partyp,type=str)
    ################################################################################################################################
    # RUN INFO PASSING
    if   partyp=='galmod':
        par={'keylist'  :['mvir','concen','halo_prof'],
             'mvir'     :{'def':0.0 ,'form':float         },
             'concen'   :{'def':0.0 ,'form':float         },
             'halo_prof':{'def':'LH','form':LIST_HALOPROFS}}
        bgpar=['mtot','rscl','rmax','zscl','zmax']
        for ityp in files.buildgal.LIST_PTYPS:
            for ipar in bgpar:
                par['keylist'].append('{}_{}'.format(ityp,ipar))
                par['{}_{}'.format(ityp,ipar)]={'form':float,'def':0.}
    elif partyp=='mmlsim': par={'keylist'    :['runtag','runtyp','subtag','memcomp','mkcomp','runcomp','scfcomp',
                                               'nprocmk','nprocev','nprocscf','inclmktag','inclevtag','inclscftag'],
                                'runtag'     :{'def':'galsim' ,'form':str                       },
                                'runtyp'     :{'def':'galsim' ,'form':LIST_RUNTYPS              },
                                'subtag'     :{'def':''       ,'form':str                       },
                                'memcomp'    :{'def':'bender' ,'form':mmlinfo.LIST_COMPUTERS    },
                                'mkcomp'     :{'def':'accre'  ,'form':mmlinfo.LIST_COMPUTERS    },
                                'runcomp'    :{'def':'accre'  ,'form':mmlinfo.LIST_COMPUTERS    },
                                'scfcomp'    :{'def':'accre'  ,'form':mmlinfo.LIST_COMPUTERS    },
                                'nprocev'    :{'def':256      ,'form':int                       },
                                'nprocmk'    :{'def':100      ,'form':int                       },
                                'nprocscf'   :{'def':50       ,'form':int                       },
                                'inclevtag'  :{'def':False    ,'form':bool                      },
                                'inclmktag'  :{'def':False    ,'form':bool                      },
                                'inclscftag' :{'def':False    ,'form':bool                      }}
#                                'runcomp'    :{'def':'accre'  ,'form':mmlinfo.LIST_COMPUTERS    },
    elif partyp=='galsim': par={'keylist'    :['addtag','model','ntot','mvir','haloprof','concen','fhalo','fdisk','fbulge','fgas','bulgespin','toomreQ'],
                                'addtag'     :{'def':''       ,'form':str                       },
                                'model'      :{'def':'MW'     ,'form':mmlparam.par_taglist('pysim.simlist','galmod')},
                                'ntot'       :{'def':long(1e7),'form':long                      },
                                'mvir'       :{'def':1.0e12   ,'form':float                     },
                                'haloprof'   :{'def':'LH'     ,'form':LIST_HALOPROFS            },
                                'concen'     :{'def':11.0     ,'form':float                     },
                                'fhalo'      :{'def':1.0      ,'form':float                     },
                                'fdisk'      :{'def':1.0      ,'form':float                     },
                                'fbulge'     :{'def':1.0      ,'form':float                     },
                                'fgas'       :{'def':0.0      ,'form':float                     },
                                'bulgespin'  :{'def':0.0      ,'form':float                     },
                                'toomreQ'    :{'def':0.0      ,'form':float                     }}
    elif partyp=='intsim': par={'keylist'    :['addtag','gal1','gal2','rperi','ecc','incl'],
                                'addtag'     :{'def':''       ,'form':str                       },
                                'gal1'       :{'def':'galsim' ,'form':mmlparam.par_taglist('pysim.simlist','galsim')},
                                'gal2'       :{'def':'galsim' ,'form':mmlparam.par_taglist('pysim.simlist','galsim')},
                                'rperi'      :{'def':0.1      ,'form':float                     },
                                'ecc'        :{'def':1.0      ,'form':float                     },
                                'incl'       :{'def':0.0      ,'form':float                     }}
    ################################################################################################################################
    # PLOTTING PARAMETERS
    elif partyp=='hist1d':
        par={'keylist':['family','galaxy','x','xscl','xlim','w','wscl','wlim','a','cmeth'],
             'family':{'def':''     ,'form':['']+exkw.get('famlist',LIST_PTYPS)},
             'galaxy':{'def':''     ,'form':['']+exkw.get('gallist',LIST_PGALS)},
             'x'     :{'def':'r'    ,'form':str},
             'xscl'  :{'def':'lin'  ,'form':LIST_SCALES},
             'xlim'  :{'def':(0.,0.),'form':tuple},
             'w'     :{'def':''     ,'form':str},
             'wscl'  :{'def':'log'  ,'form':LIST_SCALES},
             'wlim'  :{'def':(0.,0.),'form':tuple},
             'a'     :{'def':''     ,'form':str},
             'cmeth' :{'def':'hyb'  ,'form':['']+exkw.get('comlist',LIST_COMMETHS)}}
    elif partyp in ['pilimg','hist2d']:
        par={'keylist':['family','galaxy','x','xscl','xlim','y','yscl','ylim','w','wscl','wlim','a','cmeth'],
             'family':{'def':''     ,'form':['']+exkw.get('famlist',LIST_PTYPS)},
             'galaxy':{'def':''     ,'form':['']+exkw.get('gallist',LIST_PGALS)},
             'x'     :{'def':'x'    ,'form':str},
             'xscl'  :{'def':'lin'  ,'form':LIST_SCALES},
             'xlim'  :{'def':(0.,0.),'form':tuple},
             'y'     :{'def':'y'    ,'form':str},
             'yscl'  :{'def':'lin'  ,'form':LIST_SCALES},
             'ylim'  :{'def':(0.,0.),'form':tuple},
             'w'     :{'def':''     ,'form':str},
             'wscl'  :{'def':'log'  ,'form':LIST_SCALES},
             'wlim'  :{'def':(0.,0.),'form':tuple},
             'a'     :{'def':''     ,'form':str},
             'cmeth' :{'def':'pot'  ,'form':['']+exkw.get('comlist',LIST_COMMETHS)}}
    elif partyp in ['map2d','ellipse']:
        par=listpar('hist2d',**exkw)
        if plot:
            par['plotmeth']={'def':'contourf','form':['contourf','pilimg']}
            par['keylist']=['plotmeth']+par['keylist']
        par['histmeth']={'def':'numpy','form':['numpy','pNbody']}
        par['keylist']=['histmeth']+par['keylist']
    elif partyp in ['multipanel']:
        par={'keylist' :['name','plotlist','ncol'],
             'name'    :{'def':'misc','form':str },
             'plotlist':{'def':[]    ,'form':list},
             'ncol'    :{'def':3     ,'form':int }}

    ################################################################################################################################
    # ORBIT MEASUREMENT
    elif partyp=='orbit':
        par={'keylist':['ftype','cmeth'],
             'ftype':{'def':'gadget','form':LIST_PROG    },
             'cmeth':{'def':'pot'   ,'form':LIST_COMMETHS}}

    ################################################################################################################################
    # BAR MEASUREMENT
    elif partyp=='barprop':
        par={'keylist':['galaxy','family','rlim','rscl','cmeth','nbins'],
             'galaxy' :{'def':''         ,'form':['']+LIST_PGALS   },
             'family' :{'def':''         ,'form':['']+LIST_PTYPS   },
             'rlim'   :{'def':(0.01,30.0),'form':tuple             },
             'rscl'   :{'def':'log'      ,'form':LIST_SCALES       },
             'cmeth'  :{'def':'pot'      ,'form':['']+LIST_COMMETHS},
             'nbins'  :{'def':100        ,'form':int               }}
    elif partyp=='barprop_Am2':
        par=listpar('barprop',**exkw)
        par['maxmeth']={'def':'amp','form':['amp','phi']}
        par['phierr' ]={'def':10.  ,'form':float        }
        par['keylist']+=['maxmeth','phierr']
    elif partyp=='barprop_scfAm2':
        par=listpar('barprop',**exkw)
        par['phierr' ]={'def':10.  ,'form':float        }
        par['ampmin' ]={'def':0.04 ,'form':float        }
        par['keylist']+=['phierr','ampmin']
    elif partyp=='barprop_derAm2':
        par=listpar('barprop',**exkw)
        par['phierr' ]={'def':10.  ,'form':float        }
        par['keylist']+=['phierr']
    elif partyp=='barprop_rot':
        par=listpar('barprop',**exkw)
    elif partyp=='barprop_ellipse':
        par=listpar('barprop',**exkw)
        par['phierr']={'def':10. ,'form':float}
        par['ellerr']={'def':0.1 ,'form':float}
        par['nellip']={'def':100 ,'form':int  }
        par['ellmin']={'def':0.25,'form':float}
        par['keylist']+=['phierr','ellerr','nellip','ellmin']
    ################################################################################################################################
    # ERROR HANDLING AND RETURN
    else: raise Exception('Invalid parameter type: {}'.format(partyp))
    return par

####################################################################################################################################
# TAG CREATION
def par2tag(partyp,inpar,**exkw):
    """
    Creates tags
    """
    from mmlutils.mmlstring import lim2str
    partyp=mmlpars.mml_pars(partyp,type=str)
    plot=exkw.get('plot',False)
    outpar=mmlparam.parspar('pysim.simlist',partyp,inpar,**exkw)
    if   partyp=='galmod':
        validfile=False
        while not validfile:
            tagstr=mmlio.askquest('Entre a new name for this model:',default='None',dtype='str')
            fname=mmlparam.par2fname('pysim.simlist',partyp,tagstr=tagstr)
            if os.path.isfile(fname): mmlio.verbose('That file already exists...')
            else                    : validfile=True
    elif partyp=='mmlsim': tagstr=outpar['runtag']
    elif partyp=='galsim':
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
    elif partyp=='intsim':
        galsim1=mmlparam.loadpar('pysim.simlist','galsim',outpar['gal1'])
        galsim2=mmlparam.loadpar('pysim.simlist','galsim',outpar['gal2'])
        # Create short strings for model & components
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
        incl  = mmlstring.decimal(outpar['incl'])#*180./np.pi) # Angle in degrees
        # Combine strings
        tagstr='{}_{}_Q{}RP{}E{}I{}'.format(gmod1+ntot1,gmod2+ntot2,qmas,rper,ecc,incl)
        # Add addtag
        if len(outpar['addtag'])!=0: tagstr+=outpar['addtag']
    # 1D Histogram
    elif partyp=='hist1d':
        tagstr=''
        if outpar['family']: tagstr+='{family}_'
        if outpar['galaxy']: tagstr+='{galaxy}'
        for var in ['x','w','a']:
            varstr='' if not outpar[var] else outpar[var]
            if var=='a':
                limstr=''
                sclstr=''
            else:
                limstr='' if outpar[var+'lim'][0]==outpar[var+'lim'][1] else lim2str(outpar[var+'lim'])
                sclstr='' if outpar[var+'scl']=='lin' else outpar[var+'scl']
            if len(varstr+limstr+sclstr)!=0:
                if len(tagstr)!=0 and not tagstr.endswith('_'): tagstr+='_'
                tagstr+='{'+var+'str}'
                outpar[var+'str']=sclstr+varstr+limstr
        tagstr+='_{cmeth}Center'
        tagstr=tagstr.format(**outpar)
    # 2D Histogram
    elif partyp in ['hist2d','pilimg']:
        tagstr=''
        if outpar['family']: tagstr+='{family}_'
        if outpar['galaxy']: tagstr+='{galaxy}'
        for var in ['x','y','w','a']:
            varstr='' if not outpar[var] else outpar[var]
            if var=='a':
                limstr=''
                sclstr=''
            else:
                limstr='' if outpar[var+'lim'][0]==outpar[var+'lim'][1] else lim2str(outpar[var+'lim'])
                sclstr='' if outpar[var+'scl']=='lin' else outpar[var+'scl']
            if len(varstr+limstr+sclstr)!=0:
                if len(tagstr)!=0 and not tagstr.endswith('_'): tagstr+='_'
                tagstr+='{'+var+'str}'
                outpar[var+'str']=sclstr+varstr+limstr
        tagstr+='_{cmeth}Center'
        tagstr=tagstr.format(**outpar)
    elif partyp in ['map2d','ellipse']:
        tagstr=outpar['histmeth']+'_'
        if plot: tagstr+=outpar['plotmeth']+'_'
        tagstr+=par2tag('hist2d',outpar,**exkw)
    elif partyp in ['multipanel']: tagstr=outpar['name']
    # Orbit
    elif partyp=='orbit':
        tagstr='{ftype}_{cmeth}Center'.format(**outpar)
    # Bar properties
    elif partyp=='barprop':
        tagstr=''
        if outpar['family']: tagstr+='{family}_'
        if outpar['galaxy']: tagstr+='{galaxy}'
        varstr='r'
        limstr='' if outpar['rlim'][0]==outpar['rlim'][1] else lim2str(outpar['rlim'])
        sclstr='' if outpar['rscl']=='lin' else outpar['rscl']
        if len(varstr+limstr+sclstr)!=0:
            if len(tagstr)!=0 and not tagstr.endswith('_'): tagstr+='_'
            tagstr+='{rstr}'
            outpar['rstr']=sclstr+varstr+limstr
        tagstr+='_{cmeth}Center'
        tagstr=tagstr.format(**outpar)
    elif partyp=='barprop_Am2':
        tagstr=par2tag('barprop',outpar,**exkw)
        if outpar['maxmeth']!='amp': tagstr+='_'+outpar['maxmeth']
        tagstr+='_phi'+mmlstring.decimal(outpar['phierr'])
        tagstr+='_'+mmlstring.decimal(outpar['nbins'])+'bins'
    elif partyp=='barprop_derAm2':
        tagstr=par2tag('barprop',outpar,**exkw)
        #if outpar['maxmeth']!='amp': tagstr+='_'+outpar['maxmeth']
        tagstr+='_phi'+mmlstring.decimal(outpar['phierr'])
        tagstr+='_'+mmlstring.decimal(outpar['nbins'])+'bins'
    elif partyp=='barprop_scfAm2':
        tagstr=par2tag('barprop',outpar,**exkw)
        #if outpar['maxmeth']!='amp': tagstr+='_'+outpar['maxmeth']
        tagstr+='_phi'+mmlstring.decimal(outpar['phierr'])
        tagstr+='_'+mmlstring.decimal(outpar['nbins' ])+'bins'
        tagstr+='_'+mmlstring.decimal(outpar['ampmin'])+'Amin'
    elif partyp=='barprop_ellipse':
        tagstr=par2tag('barprop',outpar,**exkw)
        tagstr+='_ell'+mmlstring.decimal(outpar['ellerr'])
        tagstr+='_phi'+mmlstring.decimal(outpar['phierr'])
        tagstr+='_'+mmlstring.decimal(outpar['nbins' ])+'bins'
        tagstr+='_'+mmlstring.decimal(outpar['nellip'])+'ellp'
    # Unknown type
    else: raise Exception('Invalid parameter type: {}'.format(partyp))
    return tagstr
