#!/usr/bin/python
import pysim
import os
from mmlutils import mmlplot,mmlio,mmlmath,mmlplot
from matplotlib import pyplot as plt

# Define run and plot
runkey='q1'
meth='ellipse'
#meth='scf'

# Load run info
runs={'q1' :'MW0B7_MW0B7_Q1RP0p1E1I0',
      'q10':'MW0B7_MW0B6_Q10RP0p1E1I0'}
if   runkey=='q1' : fext=205
elif runkey=='q10': fext=245
x=pysim.loadsim(runs[runkey])

# Get plot keyword
if meth=='scf':
    plotkw={'annotate':False,'textsize':18,'ticksize':15,
            'cmap':plt.get_cmap('jet'),'fgcolor':[0,0,0],
            'xlabel':'X [kpc]','ylabel':'Y [kpc]','clabel':'Relative m=2 Amplitude',
            'cbtickloc':[-1.0,-0.1,-0.01,-0.001,0,0.001,0.01,0.1,1.0],
            'cbtickfmt':"% .3f"}
    inkw={'askuser':True,'overwrite':True,
          'plotkw':plotkw,'loadkw':{'fext':fext,'ftype':'scf'},
          'plotfile':os.path.expanduser('~/'+runkey+'_'+meth)+'{:03d}'.format(fext)+'.eps',
          'pmeth':'map2d','histmeth':'numpy','plotmeth':'contourf',
          'tagstr':'numpy_disk_galaxy1_x_y_symlogAm20to1_N_potCenter'}
          #'tagstr':'numpy_disk_galaxy2_x_y_symlogAm20to1_N_potCenter'}
          #'tagstr':'numpy_disk_galaxy1_x_y_symlogAm20p01to100_hybCenter'
elif meth=='ellipse':
    plotkw={'annotate_text':False,'textsize':18,'ticksize':15,
            'arwcolor':[0,1,0],'arwsize':400,'arwtail':0.1,'arwwidth':3.0,
            'xlabel':'X [kpc]','ylabel':'Y [kpc]','clabel':'log(Mass) [M$_{\odot}$]',
            'cbtickloc':[5.0,5.5,6.0,6.5,7.0,7.5,8.0,8.5],
            'cbtickfmt':'% .1f'}
    inkw={'askuser':True,'overwrite':True,
          'plotkw':plotkw,'loadkw':{'fext':fext,'ftype':'gadget'},
          'plotfile':os.path.expanduser('~/'+runkey+'_'+meth)+'{:03d}'.format(fext)+'.eps',
          'pmeth':'ellipse','histmeth':'numpy','plotmeth':'contourf',
          'tagstr':'numpy_disk_galaxy1_x_y_logmass_hybCenter'}
          #'tagstr':'numpy_disk_galaxy2_x_y_logmass_hybCenter'}

# Plot
out=x.plotsnap(**inkw)
#out.save(inkw['plotfile'])
print '    '+inkw['plotfile']
