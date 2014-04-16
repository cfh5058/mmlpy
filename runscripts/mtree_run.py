#!/usr/bin/python 
import mmlmtree
from mmlutils import mmlfiles

#inttype='flyby'
inttype='merger'

mtype='simstat'
#method='plotcut'
method='snaptabs'

runtag='sinham_50mpc_512_z249_ZA'
cuttag='ALL_mass10to500_q1to100_rsep0to100_z0to0'

# Cut parameters for plotting
cutpar={'x'   :'fmvir',
        'xlim':(0.,0.),
        'xscl':'log'  ,
        'y'   :'otr'  ,
        'ylim':(0.,0.),
        'yscl':'log'  }
for key in ['clr','mrk','siz']:
    cutpar[key+'var']='otr'
    cutpar[key+'lim']=(0.,0.)
    cutpar[key+'scl']='lin'

# Mask options
mask={'all'   :{'tag':[]},
      'flyby' :{'tag':[-2,-1,-999,-666,99,0,1,2,4]},
      'merger':{'tag':[-2,-1,-999,-666,99,3,5]}}

# Variables to make histograms of
varbase=['mtot','htot','spin_tot',
         'mvir','hvir','spin_vir']
varlist=[]
for var in varbase: varlist+=['i'+var,'f'+var,var+'diff']


if method=='plotcut':
    print ' '
    for var in varlist:
        cutpar['x']=var
        xold=mmlmtree.main.walk(mtype=mtype,method=method,
                                runtag=runtag,cuttag=cuttag,cutpar=cutpar,mask=mask[inttype],
                                verbose=False,overwrite=True,askuser=True)
        xnew=xold.replace('ALL',inttype.upper())
        mmlfiles.cpfiles(xold,xnew,move=True)
        print 'get '+xnew
    print ' '
elif method=='snaptabs':
    z=mmlmtree.main.walk(mtype=mtype,method=method,runtag=runtag,overwrite=True)
elif method=='lastints':
    par={'cuttag':cuttag,'nback':2}
    tag='lastint_n2_ALL_mass10to500_q1to100_rsep0to100_z0to0'
    intlist=simstr.get_intlist(cuttag=par['cuttag'])
    xold=mmlmtree.main.walk(mtype=mtype,method=method,runtag=runtag,
                            par=par,tag=tag,intlist=intlist,
                            verbose=False,overwrite=True,askuser=True)
