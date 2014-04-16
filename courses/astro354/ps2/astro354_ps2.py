#!/usr/bin/python
import numpy as np
import matplotlib.pyplot as plt
import scipy.integrate as integrate
from mmlutils import mmlio

def main():
    print "Astro 354: PS2"
    # Set constants
    datafile='/home/langmm/utils/python/mine/courses/astro354/ps2/SDSS_DR7.dat'
    keystr=['ra','dec','redshift','Mg','Mr']
    typstr=len(keystr)*['float']
    grcut=0.75
    # Load data
    data=mmlio.rwtable('R',datafile,keystr=keystr,typstr=typstr)
    for ikey in data['keylist']:
        data[ikey]=np.array(data[ikey])
    data['g-r']=data['Mg']-data['Mr']
    Ngal=len(data['redshift'])
    # Determine volume
    zlim=0.1
    omega=2.295 # steradians
    volume=z2vol(zlim,omega)

    # PART A
    print "Part a)"
    # Plot
    plotname1='/home/langmm/utils/python/mine/courses/astro354/ps2/astro354_ps2_parta.png'
    fig1=plt.figure()
    ax11=plt.subplot(2,1,1)
    ax11.plot(data['ra'],data['dec'],'.')
    ax11.set_xlabel('Right Ascension (deg)')
    ax11.set_ylabel('Declination (deg)')
    ax12=plt.subplot(2,1,2)
    ax12.plot(data['redshift'],data['Mr'],'.')
    ax12.set_xlabel('Redshift')
    ax12.set_ylabel('Mr (mag)')
    ax12.invert_yaxis()
    fig1.tight_layout()
    fig1.savefig(plotname1)
    plt.close(fig1)
    print '    {}'.format(plotname1)

    # PART B
    print "Part b)"
    grhist,grbins=np.histogram(data['g-r'],bins=100)
    bluefrac=get_bluefrac(data['g-r'],grcut)
    print '    Ngal (tot): {}'.format(Ngal)
    print '    blue frac:  {}'.format(bluefrac)
    print '    red  frac:  {}'.format(1.-bluefrac)
    # Plot
    plotname2='/home/langmm/utils/python/mine/courses/astro354/ps2/astro354_ps2_partb.png'
    fig2=plt.figure()
    ax21=plt.subplot(1,1,1)
    ax21.plot(grbins[:-1]+0.5*(grbins[1:]-grbins[:-1]),grhist)
    ax21.set_xlabel('g-r Color')
    ax21.set_ylabel('N(g-r)')
    fig2.tight_layout()
    fig2.savefig(plotname2)
    plt.close(fig2)
    print '    {}'.format(plotname2)

    # PART C
    print "Part c)"
    print "    {} Mpc^3".format(volume)
    # Compute histogram
    bins_Mr,n_Mr=get_lumfunc(data['Mr'],volume)
    # Plot
    plotname3='/home/langmm/utils/python/mine/courses/astro354/ps2/astro354_ps2_partc.png'
    fig3=plt.figure()
    ax31=plt.subplot(1,1,1)
    ax31.plot(bins_Mr,np.ma.log10(n_Mr))
    ax31.set_xlabel('Mr (mag)')
    ax31.set_ylabel('log(dn/dMr) (h^3/Mpc^3)')
    fig3.tight_layout()
    fig3.savefig(plotname3)
    plt.close(fig3)
    print '    {}'.format(plotname3)

    # PROBLEM 3
    print 'Part d)'
    # Get g-r and M-M0 at specified times 
    Mrlimlist=[-20.,-19.,-18.]
    Mrlimdict={}
    labellist=[]
    for iMrlim in Mrlimlist:
        idxcut_mag=(data['Mr']<=iMrlim)
        iMrlimdict={}
        # Calculate quantities for part d
        iMrlimdict['Mrlim'   ]=iMrlim
        iMrlimdict['dmax'    ]=Mr2dmax(iMrlim)
        iMrlimdict['zmax'    ]=d2z(iMrlimdict['dmax'])
        iMrlimdict['zmin'    ]=0.
        idxcut_red=(data['redshift']<=iMrlimdict['zmax'])
        idxcut=np.logical_and(idxcut_mag,idxcut_red)
        iMrlimdict['vol'     ]=d2vol(iMrlimdict['dmax'],omega)
        iMrlimdict['Ngal'    ]=len(data['Mr'][idxcut])
        iMrlimdict['bluefrac']=get_bluefrac(data['g-r'][idxcut],grcut)
        # Calculate luminosity function for part e
        iMrlimdict['lfunc_bins'],iMrlimdict['lfunc_dn']=get_lumfunc(data['Mr'][idxcut],iMrlimdict['vol'])
        # Add dictionary for this magnitude cut
        iMrlimstr=repr(iMrlim)
        labellist.append(iMrlimstr)
        Mrlimdict[iMrlimstr]=iMrlimdict
        # Print info to screen
        print '    Mrlim = {}'.format(iMrlimdict['Mrlim'])
        print '      zlim     = [{},{}]'.format(iMrlimdict['zmin'],iMrlimdict['zmax'])
        print '      volume   = {} Mpc^3/h^3'.format(iMrlimdict['vol'])
        print '      Ngal     = {}'.format(iMrlimdict['Ngal'])
        print '      bluefrac = {}'.format(iMrlimdict['bluefrac'])
#        print '      dmax     = {} Mpc/h'.format(iMrlimdict['dmax'])

    # PART E
    print 'Part e)'
    # Plot
    plotname4='/home/langmm/utils/python/mine/courses/astro354/ps2/astro354_ps2_parte.png'
    colorlist=['b','g','r']
    fig4=plt.figure()
    ax41=plt.subplot(1,1,1)
    for ilab,iclr in zip(labellist,colorlist):
        ax41.plot(Mrlimdict[ilab]['lfunc_bins'],np.ma.log10(Mrlimdict[ilab]['lfunc_dn']),iclr,label=ilab)
    ax41.set_xlabel('Mr (mag)')
    ax41.set_ylabel('log(dn/dMr) (h^3/Mpc^3)')
    plt.legend()
    fig4.tight_layout()
    fig4.savefig(plotname4)

    plt.close(fig4)
    print '    {}'.format(plotname4)

    # PART F
    print 'Part f)'
    # Determine weigths
    dmax=Mr2dmax(data['Mr'])
    vmax=d2vol(dmax,omega)
    bins_vmax,lfunc_vmax=get_lumfunc(data['Mr'],vmax)
    # Plot
    plotname5='/home/langmm/utils/python/mine/courses/astro354/ps2/astro354_ps2_partf.png'
    fig5=plt.figure()
    ax51=plt.subplot(1,1,1)
    ax51.plot(bins_vmax,lfunc_vmax)
    ax51.set_xlabel('Mr (mag)')
    ax51.set_ylabel('log(dn/dMr) (h^3/Mpc^3)')
    fig5.tight_layout()
    fig5.savefig(plotname5)
    plt.close(fig5)
    print '    {}'.format(plotname5)
    
    return

# Convert redshift to distance
def z2d(z):
    c=3.0e5 # km/s
    H0=100.0 # h*km/s/Mpc
    d=c*z/H0 # Mpc/h 
    return d
def d2z(d):
    c=3.0e5 # km/s
    H0=100.0 # h*km/s/Mpc
    z=H0*d/c
    return z
# Convert redshift to volume
def z2vol(zlim,omega):
    dlim=z2d(zlim)
    volume=d2vol(dlim,omega)
    return volume
# Convert distance to volume
def d2vol(dlim,omega):
    volume=omega*(dlim**3)/3
    return volume
# Determine the fraction of 'blue' galaxies based on cut
def get_bluefrac(grcolor,grcut):
    Ngal=len(grcolor)
    Ngal_blue=len(grcolor[grcolor<grcut])
    return float(Ngal_blue)/float(Ngal)
# Determine luminosity function
def get_lumfunc(Mr,volume):
    min_Mr=Mr.min()
    max_Mr=Mr.max()
    binsize=0.1
    nbins_Mr=int(np.ceil((max_Mr-min_Mr)/binsize))
    bins_Mr=min_Mr+binsize*np.array(range(nbins_Mr+1))
    if isinstance(volume,float):
        N_Mr,bins_Mr=np.histogram(Mr,bins=bins_Mr)
        n_Mr=N_Mr/volume
    else:
        n_Mr,bins_Mr=np.histogram(Mr,bins=bins_Mr,weights=1./volume)
    bins_mid=bins_Mr[:-1]+0.5*(bins_Mr[1:]-bins_Mr[:-1])
    return bins_mid,n_Mr
# Determine maximum distance at which magnitude can be observed
def Mr2dmax(Mr):
    rlim=17.77
    dmax=(10.0**((rlim-Mr+5.)/5.))/(1.0e6)
    return dmax

if __name__ == '__main__':
    main()
