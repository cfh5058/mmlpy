#!/usr/bin/python
import sys,os,shutil,glob,copy,pprint,scipy,math
import numpy as np
from mmlutils import *
#import main as mmlmtree

LIST_METHTYP=['general','simstat','simcalc','simplot']
LIST_METHTYP_NOF=[]
DICT_INTTYPS={'ALL'    :[0,1,2,3,4,5],
              'MERGER' :[0,1,2,3    ],
              'FLYBY'  :[      3,  5],
              'FLYSLOW':[      3    ],
              'FLYFAST':[          5]}
LIST_INTTYPS=sorted(DICT_INTTYPS.keys())
LIST_ICMETHS=['ZA','2LPT']

import simstat

DICT_METHODS={'general':['get_simobj','get_snapshot','get_halosnap','get_groups','get_intlist','plot_intlist'],
              'simstat':simstat.LIST_METHODS}

def listpar(partyp,default=None,**exkw):
    partyp=mmlpars.mml_pars(partyp,type=str)
    default=mmlpars.mml_pars(default,type=bool,default=False)
    if default:
        if   partyp=='cossim': par={'user'    :'sinham'       ,
                                    'boxsize' :50.            ,
                                    'npart3'  :512            ,
                                    'zstart'  :249.           ,
                                    'icmeth'  :LIST_ICMETHS[0],
                                    'memcomp' :'bender'       ,
                                    'grpdir'  :'/data1/sinham/isolatedhalos/50mpc/512/z249/groups'}
        elif partyp=='intpar': par={'tag'     :3              ,
                                    'iz'      :0.             ,
                                    'fz'      :0.             ,
                                    'isnap'   :0              ,
                                    'fsnap'   :0              ,
                                    'haloid1' :0L             ,
                                    'haloid2' :0L             ,
#                                    'm1_z0'   :0.             ,
                                    'm1'      :0.             ,
                                    'rvir1'   :0.             ,
                                    'rvir1_phys':0.             ,
                                    'q'       :0.             ,
                                    'minrsep' :0.             }
        else: raise Exception('Invalid parameter type: {}'.format(partyp))
    else:
        if   partyp=='cossim': par={'user'    :str            ,
                                    'boxsize' :float          ,
                                    'npart3'  :int            ,
                                    'zstart'  :float          ,
                                    'icmeth'  :LIST_ICMETHS   ,
                                    'memcomp' :mmlinfo.LIST_COMPUTERS,
                                    'grpdir'  :str            }
        elif partyp=='intpar': par={'tag'     :DICT_INTTYPS['ALL'],
                                    'iz'      :float          ,
                                    'fz'      :float          ,
                                    'isnap'   :int            ,
                                    'fsnap'   :int            ,
                                    'haloid1' :long           ,
                                    'haloid2' :long           ,
#                                    'm1_z0'   :float          ,
                                    'm1'      :float          ,
                                    'rvir1'   :float          ,
                                    'rvir1_phys':float          ,
                                    'q'       :float          ,
                                    'minrsep' :float          }
        else: raise Exception('Invalid parameter type: {}'.format(partyp))
    return par

LIST_HALOREAD=['mtot','Vtot','htot','spin_tot','mgas','N','Lx','Ly','Lz','G']
#LIST_HALOCALC=['Nvir','mvir','Vvir','hvir','spin_vir']
LIST_HALOCALC=['mvir','Vvir','hvir','spin_vir']
LIST_HALOPROP=LIST_HALOREAD+LIST_HALOCALC
LIST_INTREAD=mmlparam.listpar('mmlmtree.simlist','intpar')['keylist']
LIST_INTCALC=['tint','Lphidiff']+[ih+'diff' for ih in LIST_HALOPROP]
LIST_INTHALO=['i'+ih for ih in LIST_HALOPROP]+['f'+ih for ih in LIST_HALOPROP]
LIST_INTPROP=LIST_INTREAD+LIST_INTCALC+LIST_INTHALO
#LIST_HALOPROP=['cm','cmv','mtot','mgas','spin','tiso']#,'iz','fz','nint']
LIST_HALOVINT=['halo_'+ihkey for ihkey in LIST_HALOPROP]+['int_'+iikey for iikey in LIST_INTREAD]

LIST_TABKEYS=['Snapshot','GroupNum','Groupid','Mtot','LastMergerRed','Formationz','Nmergers','TotNmergers',
              'ParentLev','Nsub','FofHaloid','ContainerNum','Rvir','Rvir_anyl','Rhalf','Nflybys','TotNflybys']
LIST_TABVARS=['N','m','v','h','spin']
LIST_TABRADS=['Rvir','Rvir_anyl','Rhalf']
LIST_SNAPTABKEYS=LIST_TABKEYS
LIST_INTVARS=['haloid1','haloid2','tag','minrsep','q']
LIST_LASTINTKEYS=LIST_SNAPTABKEYS+LIST_INTVARS

for v in LIST_TABVARS:
    for r in LIST_TABRADS:
        LIST_SNAPTABKEYS.append(v+'_'+r)

def par2tag(partyp,inpar,**exkw):
    partyp=mmlpars.mml_pars(partyp,type=str)
    outpar=mmlparam.parspar('mmlmtree.simlist',partyp,inpar)
    if partyp=='cossim':
        outpar['boxstr']=mmlstring.dec2str(outpar['boxsize'])
        outpar['zstr'  ]=mmlstring.dec2str(outpar['zstart' ])
        tagstr='{user}_{boxstr}mpc_{npart3}_z{zstr}_{icmeth}'.format(**outpar)
    elif partyp=='intpar': 
        if outpar['tag'] in [3,5]: outpar['type']='FLYBY'
        else                     : outpar['type']='MERGER'
        tagstr='{type}_{haloid1:06d}_{haloid2:06d}_{isnap:03d}'.format(**outpar)
    else: raise Exception('Invalid parameter type: {}'.format(partyp))
    return tagstr





        
