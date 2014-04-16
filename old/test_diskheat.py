#!/usr/bin/python
import copy
import numpy as np
import matplotlib.pyplot as plt

####################################################################################################################################
# METHOD
####################################################################################################################################
def test(mean=0.,sig0=1.,npart=1000,nbins=100):
    '''
    NAME:
        
    PURPOSE:
        
    CALLING:
        
    ARGUMENTS:
        
    KEYWORDS:
        
    OUTPUT:
        
    '''
    # Create velocity distributions
    distrib={}
    vel_sig1p0=np.random.normal(loc=mean,scale=sig0,size=(npart,1))
    vel_sig1p1=np.random.normal(loc=mean,scale=1.1*sig0,size=(npart,1))
    vmax=abs(vel_sig1p0).max() ; print 'vmax={}'.format(vmax)
    range=(-vmax,vmax)
    distrib['sig1p0']={'label':'sigma=1.0'  ,'color':'-k' ,'vel':vel_sig1p0                   }
    distrib['sig1p1']={'label':'sigma=1.1'  ,'color':'--k','vel':vel_sig1p1                   }
    distrib['ln_flt']={'label':'b*v (b=1.1)','color':'-b' ,'vel':heat_linefact(vel_sig1p0,1.1)}
#    distrib['ln_sig']={'label':'(av+b)*v (a=0.1/sig(v),b=1)','color':'-r' ,'vel':heat_linefact(vel_sig1p0,1.,0.1/sig0)}
#    distrib['ln_max']={'label':'(av+b)*v (b=0.1/max(v),b=1)','color':'-g' ,'vel':heat_linefact(vel_sig1p0,1.,0.1/vmax)}
    keyord=['ln_flt','ln_sig','ln_max','sig1p0','sig1p1']

    # Create histograms
    for k in distrib.keys():
        distrib[k]['hist'],distrib[k]['bins']=np.histogram(distrib[k]['vel'],nbins,range=range,normed=True)

    # Plot
    fig=plt.figure()
    axObj=fig.add_subplot(1,1,1)
    for k in keyord:
        if distrib.has_key(k):
            midbins=distrib[k]['bins'][:-1]+(distrib[k]['bins'][1:]-distrib[k]['bins'][:-1])/2.
            axObj.plot(midbins,distrib[k]['hist'],distrib[k]['color'],label=distrib[k]['label'])
    axObj.legend(bbox_to_anchor=(0., 1.02, 1., .102), loc=3,
                 ncol=3, mode="expand", borderaxespad=0.)
#    axObj.legend(bbox_to_anchor=(1.05, 1), loc=2, borderaxespad=0.)
    axObj.set_xlabel('v')
    axObj.set_ylabel('N(v)/Ntot')
#    fig.tight_layout()
    fig.savefig('test_diskheat.png')
    plt.close(fig)

def heat_flatfact(velin,factor):
    vel=copy.deepcopy(velin)
    vel*=factor
    return vel
def heat_linefact(velin,intercept=1.,slope=0.):
    vel=copy.deepcopy(velin)
    vel*=(slope*abs(vel)+intercept)
    return vel
def heat_polyfact(velin,pfactors):
    vel=copy.deepcopy(velin)
    for ip in range(len(pfactors)):
        factor=pfactors[ip]*(abs(vel)**ip)
    vel*=factor
    return vel

if __name__ == '__main__':
    main()
