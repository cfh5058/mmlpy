#!/usr/bin/python
####################################################################################################################################
#
# MEAGAN LANG'S SIMFILE METHODS
#
####################################################################################################################################

import numpy as np
from mmlutils import *

####################################################################################################################################
####################################################################################################################################
# METHODS FOR GETTING LISTS
LIST_PROGMAKE=['buildgal']
LIST_PROGEVOL=['gadget']
LIST_PROGPOST=['scf']
LIST_PROG=LIST_PROGMAKE+LIST_PROGEVOL+LIST_PROGPOST

LIST_RUNTYPS=['galsim','intsim']
LIST_YORN=['Y','N']
LIST_PTYPGRPS=['all','visi']
LIST_PTYPBASE=['gas','halo','disk','bulge','stars','bndry']
LIST_PTYPS=LIST_PTYPGRPS+LIST_PTYPBASE
LIST_PTYPVISI=['disk','bulge','stars']
LIST_PTYPBG=['halo','disk','bulge','gas']
LIST_PGALS=[0,1,2]
LIST_COMMETHS=['mass','dens','pot','dens_avg','pot_avg']
LIST_COMPTYPS=LIST_PTYPS+['v']
LIST_COMPGALS=LIST_PGALS+['v']
LIST_HISTMETHS=['dir','avg','den']
LIST_DISPVARS=['otr','gal','typ','galtyp']
LIST_SCALES=['lin','log']
DICT_PGALS={0:'All',1:'Primary',2:'Secondary'}
LIST_HALOPROFS=['LH','NF','IS']

LIST_METHTYPS=['general']+LIST_PROG
DICT_METHODS={}

def listpar(partyp,default=None,**exkw):
    partyp=mmlpars.mml_pars(partyp,type=str)
    default=mmlpars.mml_pars(default,type=bool,default=False)
    if   partyp=='galmod':
        par={'keylist'  :['mvir','concen','halo_prof'],
             'mvir'     :{'def':0.0 ,'form':float         },
             'concen'   :{'def':0.0 ,'form':float         },
             'halo_prof':{'def':'LH','form':LIST_HALOPROFS}}
        bgpar=['mtot','rscl','rmax','zscl','zmax']
        for ityp in LIST_PTYPBG:
            for ipar in bgpar:
                par['keylist'].append('{}_{}'.format(ityp,ipar))
                par['{}_{}'.format(ityp,ipar)]={'form':float,'def':0.}
    elif partyp=='mmlsim': par={'keylist'    :['runtag','runtyp','subtag','memcomp','nprocmk','nprocev','nprocscf','inclmktag','inclevtag','inclscftag'],
                                'runtag'     :{'def':'galsim' ,'form':str                       },
                                'runtyp'     :{'def':'galsim' ,'form':LIST_RUNTYPS              },
                                'subtag'     :{'def':''       ,'form':str                       },
                                'memcomp'    :{'def':'bender' ,'form':mmlinfo.LIST_COMPUTERS    },
                                'runcomp'    :{'def':'bender' ,'form':mmlinfo.LIST_COMPUTERS    },
                                'nprocev'    :{'def':256      ,'form':int                       },
                                'nprocmk'    :{'def':100      ,'form':int                       },
                                'nprocscf'   :{'def':50       ,'form':int                       },
                                'inclevtag'  :{'def':False    ,'form':bool                      },
                                'inclmktag'  :{'def':False    ,'form':bool                      },
                                'inclscftag' :{'def':False    ,'form':bool                      }}
#                                'runcomp'    :{'def':'accre'  ,'form':mmlinfo.LIST_COMPUTERS    },
    elif partyp=='galsim': par={'keylist'    :['addtag','model','ntot','mvir','haloprof','concen','fhalo','fdisk','fbulge','fgas','bulgespin','toomreQ'],
                                'addtag'     :{'def':''       ,'form':str                       },
                                'model'      :{'def':'MW'     ,'form':mmlparam.par_taglist('mmlyt.simlist','galmod')},
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
                                'gal1'       :{'def':'galsim' ,'form':mmlparam.par_taglist('mmlyt.simlist','galsim')},
                                'gal2'       :{'def':'galsim' ,'form':mmlparam.par_taglist('mmlyt.simlist','galsim')},
                                'rperi'      :{'def':0.1      ,'form':float                     },
                                'ecc'        :{'def':1.0      ,'form':float                     },
                                'incl'       :{'def':0.0      ,'form':float                     }}
    else: raise Exception('Invalid parameter type: {}'.format(partyp))
    return par

def par2tag(partyp,inpar,**exkw):
    partyp=mmlpars.mml_pars(partyp,type=str)
    outpar=mmlparam.parspar('mmlyt.simlist',partyp,inpar)
    if   partyp=='galmod':
        validfile=False
        while not validfile:
            tagstr=mmlio.askquest('Entre a new name for this model:',default='None',dtype='str')
            fname=mmlparam.par2fname('mmlyt.simlist',partyp,tagstr=tagstr)
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
        galsim1=mmlparam.loadpar('mmlyt.simlist','galsim',outpar['gal1'])
        galsim2=mmlparam.loadpar('mmlyt.simlist','galsim',outpar['gal2'])
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
        incl  = mmlstring.decimal(outpar['incl']*180./np.pi) # Angle in degrees
        # Combine strings
        tagstr='{}_{}_Q{}RP{}E{}I{}'.format(gmod1+ntot1,gmod2+ntot2,qmas,rper,ecc,incl)
        # Add addtag
        if len(outpar['addtag'])!=0: tagstr+=outpar['addtag']
    else: raise Exception('Invalid parameter type: {}'.format(partyp))
    return tagstr
