#!/usr/bin/python
####################################################################################################################################
#
# MEAGAN LANG'S SIMFILE METHODS
#
####################################################################################################################################

####################################################################################################################################
####################################################################################################################################
# METHODS FOR GETTING LISTS
LIST_YORN=['Y','N']
LIST_RUNTYPS=['galsim','intsim']
LIST_PTYPGRPS=['all','visi']
LIST_PTYPBASE=['gas','halo','disk','bulge','stars','bndry']
LIST_PTYPS=LIST_PTYPGRPS+LIST_PTYPBASE
LIST_PTYPVISI=['disk','bulge','stars']
LIST_PGALS=[0,1,2]
LIST_COMMETHS=['mass','dens','pot','dens_avg','pot_avg','com']
LIST_COMPTYPS=LIST_PTYPS+['v']
LIST_COMPGALS=LIST_PGALS+['v']
LIST_HISTMETHS=['dir','avg','den']
LIST_DISPVARS=['otr','gal','typ','galtyp']
LIST_SCALES=['lin','log']
DICT_PGALS={0:'All',1:'Primary',2:'Secondary'}

LIST_METHTYP=['general','icgen','buildgal','gadget','scf','galfit','stat','calc','hist','plots','prof','compare','file']
LIST_METHTYP_NOF=['general','file']

import icgen,mmlbuildgal,mmlgadget,mmlscf,mmlgalfit,simstat,simcalc,simhist,simplot,simprof,compare,simfile

DICT_METHODS={'general' :['getsimobj','getpnbody','printinfo','checkrepeat','fitprofiles','displayplot','scatterplot'],
              'icgen'   :icgen.LIST_METHODS      ,
              'buildgal':mmlbuildgal.LIST_METHODS,
              'gadget'  :mmlgadget.LIST_METHODS  ,
              'scf'     :mmlscf.LIST_METHODS     ,
              'galfit'  :mmlgalfit.LIST_METHODS  ,
              'stat'    :simstat.LIST_METHODS    ,
              'calc'    :simcalc.LIST_METHODS    ,
              'hist'    :simhist.LIST_METHODS    ,
              'plots'   :simplot.LIST_METHODS    ,
              'prof'    :simprof.LIST_METHODS    ,
              'compare' :compare.LIST_METHODS    ,
              'file'    :simfile.LIST_METHODS    }
