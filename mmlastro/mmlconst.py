#!/usr/bin/python
from mmlutils import mmlpars
import numpy as np

####################################################################################################################################
####################################################################################################################################
# METHODS TO RETURN LISTS
LIST_CONSTMETH=['astro','cosmo','phys','units']
LIST_COSMOMETH=['gadget2','LCDM','WMAP5','PLANK1']
LIST_UNITSMETH=['cgs','mks','galaxy','gadget2','buildgal','SI']
LIST_ASTROMETH=LIST_UNITSMETH
LIST_PHYSMETH=LIST_UNITSMETH

####################################################################################################################################
# METHOD TO HANDLE ALL TYPES OF CONSTANTS
def main(idtag,subtag=None):
    """
    Returns a dictionary of constants
    """
    # Pars input
    idtag=mmlpars.mml_pars(idtag,list=LIST_CONSTMETH)
    # Proceed based on idtag
    if   idtag=='units': const=const_units(subtag)
    elif idtag=='astro': const=const_astro(subtag)
    elif idtag=='cosmo': const=const_cosmo(subtag)
    elif idtag=='phys' : const=const_phys(subtag)
    else: raise Exception('Invalid constant ID string: {}'.format(idtag))
    # Return constants
    return const
    
####################################################################################################################################
# METHOD TO RETURN DICTIONARY OF UNITS
def const_units(subtag=None):
    """
    Returns a dictionary of units (in cgs)
    """
    # Pars input
    subtag=mmlpars.mml_pars(subtag,list=LIST_UNITSMETH,default=LIST_UNITSMETH[0])
    # Set constants
    if   subtag=='cgs':
        const={'M':1.0        , # g
               'L':1.0        , # cm
               'V':1.0        , # cm/s
               'Q':1.0        , # statC (erg*cm)**0.5
               'B':1.0        } # G (statC/(cm**2))
    elif subtag in ['mks','SI']:
        const={'M':1.0e3      , # kg
               'L':1.0e2      , # m
               'V':1.0e2      , # m/s
               'Q':2.9979246e9, # C
               'B':1.0e4      } # T
    elif subtag=='gadget2':
        const={'M':1.989e43   , # 10**10 Msol
               'L':3.085678e21, # kpc
               'V':1.0e5      } # km/s
    elif subtag=='galaxy':
        const={'M':1.989e33   , # Msol
               'L':3.085678e18, # pc
               'V':1.0e5      } # km/s
    elif subtag=='buildgal': const=const_units('galaxy')
    else: raise Exception('Invalid units method: {}'.format(subtag))
    const['T']=const['L']/const['V']
    const['K']=1.0
    const['E']=const['M']*(const['V']**2)
    if 'Q' not in const: const['Q']=np.sqrt(const['E']*const['L'])
    if 'B' not in const: const['B']=const['Q']/(const['L']**2)
    # Return constants
    return const

####################################################################################################################################
# METHOD TO RETURN DICTIONARY OF PHYS CONSTANTS
def const_phys(subtag=None):
    """
    Returns a dicitonary of physics constants
    """
    # Pars input
    subtag=mmlpars.mml_pars(subtag,list=LIST_PHYSMETH,default=LIST_PHYSMETH[0])
    # Set constant base and conversion
    const_cgs={'c'    :2.99792458e10 , 'unit_c'    :'{V}'                        , # cm/s
               'h'    :6.6260755e-27 , 'unit_h'    :'{E}*{T}'                    , # erg*s
               'hbar' :1.05457266e-27, 'unit_hbar' :'{E}*{T}'                    , # erg*s
#               'G'    :6.67259e-8    , 'unit_G'    :'{L}*{V}*{V}/{M}'            , # cm*((cm/s)**2)/g
               'G'    :6.672e-8      , 'unit_G'    :'{L}*{V}*{V}/{M}'            , # cm*((cm/s)**2)/g
               'e'    :4.8032068e-10 , 'unit_e'    :'{Q}'                        , # statC
               'me'   :9.1093897e-28 , 'unit_me'   :'{M}'                        , # g
               'mp'   :1.6726231e-24 , 'unit_mp'   :'{M}'                        , # g
               'mn'   :1.6749286e-24 , 'unit_mn'   :'{M}'                        , # g
               'mh'   :1.6733e-24    , 'unit_mh'   :'{M}'                        , # g
               'amu'  :1.6605402e-24 , 'unit_amu'  :'{M}'                        , # g
               'Na'   :6.0221367     , 'unit_Na'   :'1'                          , # number
               'k'    :1.380658e-16  , 'unit_k'    :'{E}/{K}'                    , # erg/K
               'eV'   :1.6021772e-12 , 'unit_eV'   :'{E}'                        , # erg
               'a'    :7.5646e-15    , 'unit_a'    :'{E}/(({L}**3)*({K}**4))'    , # erg/((cm**3)*(K**4))
               'sigma':5.67051e-5    , 'unit_sigma':'{E}/({T}*({L}**2)*({K}**4))', # erg/(s*(cm**2)*(K**4))
               'alpha':7.29735308e-3 , 'unit_alpha':'1'                          , # unitless
               'Rh'   :2.1798741e-11 , 'unit_Rh'   :'{E}'                        } # erg
    # Get units
    units_cgs=const_units(subtag)
    # Add scaled constants
    const={}
    for iconst in const_cgs:
        if iconst.startswith('unit_'): continue
        const[iconst]=const_cgs[iconst]/eval(const_cgs['unit_'+iconst].format(**units_cgs))
    # Return constants
    return const

####################################################################################################################################
# METHOD TO RETURN DICTIONARY OF ASTRO CONSTANTS
def const_astro(subtag=None):
    """
    Returns a dicitonary of astronomy constants
    """
    # Pars input
    subtag=mmlpars.mml_pars(subtag,list=LIST_ASTROMETH,default=LIST_ASTROMETH[0])
    # Set constant base and conversion
    const_cgs={'AU'  :1.496e13, 'unit_AU'  :'{L}'    , # cm
               'pc'  :3.086e18, 'unit_pc'  :'{L}'    , # cm
               'ly'  :9.463e17, 'unit_ly'  :'{L}'    , # cm
               'Msol':1.99e33 , 'unit_Msol':'{M}'    , # g
               'Rsol':6.96e10 , 'unit_Rsol':'{L}'    , # g
               'Lsol':3.9e33  , 'unit_Lsol':'{E}/{T}', # erg/s
               'Tsol':5.780e3 , 'unit_Tsol':'{K}'    , # K
               'yr'  :3.1557e7, 'unit_yr'  :'{T}'    } # s
    # Get units
    units_cgs=const_units(subtag)
    # Add scaled constants
    const={}
    for iconst in const_cgs:
        if iconst.startswith('unit_'): continue
        const[iconst]=const_cgs[iconst]/eval(const_cgs['unit_'+iconst].format(**units_cgs))
    # Return constants
    return const

####################################################################################################################################
# METHOD TO RETURN COSMOLOGICAL CONSTANTS
def const_cosmo(subtag=None):
    """
    Returns a dictionary of cosmological constants
    """
    # Pars input
    subtag=mmlpars.mml_pars(subtag,list=LIST_COSMOMETH,default=LIST_COSMOMETH[0])
    # Set constants
    if   subtag=='gadget2': const={'h'     :0.7 ,
                                   'omegaR':0.0 ,
                                   'omegaM':0.3 ,
                                   'omegaB':0.04,
                                   'omegaK':0.0 ,
                                   'omegaL':0.7 }
    elif subtag=='LCDM'   : const={'h'     :0.7 ,
                                   'omegaR':0.0 ,
                                   'omegaM':0.3 ,
                                   'omegaB':0.04,
                                   'omegaK':0.0 ,
                                   'omegaL':0.7 ,
                                   'sigma8':0.8 ,
                                   'ncosmo':1.0 ,
                                   'wdarkE':None}
    elif subtag=='WMAP5'  : const={}
    # Planck 2013 results. XVI. Cosmological parameters
    elif subtag=='PLANK1' : const={'h'     :0.673  , # 0.012
                                   'omegaR':None   ,
                                   'omegaM':0.315  , # 0.017
                                   'omegaK':-0.0005, # 0.00655
                                   'omegaL':0.685  , # 0.017
                                   'sigma8':0.829  , # 0.012
                                   'ncosmo':0.9603 , # 0.0073
                                   'wdarkE':-1.109 } # 0.24
    else: raise Exception('Invalid constant ID: {}'.format(subtag))
    # Add dependent constants
    if subtag=='PLANK1':
        const['omegaB']=0.02205/(const['h']**2) # 0.00028
        const['omegaC']=0.11990/(const['h']**2) # 0.00270
    # Return constants
    return const

####################################################################################################################################
# METHOD TO RETURN PHYSICAL CONSTANTS

if __name__ == '__main__':
    main()
