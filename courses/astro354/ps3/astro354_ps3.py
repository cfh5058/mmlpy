#!/usr/bin/python
import os
import numpy as np
import matplotlib.pyplot as plt
import scipy.integrate as integrate
from mmlutils import mmlio

def main():
    print "Astro 354: PS3"
    # Set constants
    fname_SDSS_Mr21_rspace='/home/langmm/utils/python/mine/courses/astro354/ps3/data/SDSS_Mr21_rspace.dat'
    fname_SDSS_Mr20_rspace='/home/langmm/utils/python/mine/courses/astro354/ps3/data/SDSS_Mr20_rspace.dat'
    fname_SDSS_Mr20_zspace='/home/langmm/utils/python/mine/courses/astro354/ps3/data/SDSS_Mr20_zspace.dat'
    fname_SDSS_rand='/home/langmm/utils/python/mine/courses/astro354/ps3/data/SDSS_random.dat'
    fname_DM_data='/home/langmm/utils/python/mine/courses/astro354/ps3/data/DM.dat'
    fname_DM_rand='/home/langmm/utils/python/mine/courses/astro354/ps3/data/DM_random.dat'
    fcorr_SDSS_Mr21_rspace='/home/langmm/utils/python/mine/courses/astro354/ps3/corr/corr_SDSS_Mr21_rspace.dat'
    fcorr_SDSS_Mr20_rspace='/home/langmm/utils/python/mine/courses/astro354/ps3/corr/corr_SDSS_Mr20_rspace.dat'
    fcorr_SDSS_Mr20_zspace='/home/langmm/utils/python/mine/courses/astro354/ps3/corr/corr_SDSS_Mr20_zspace.dat'
    fcorr_DM='/home/langmm/utils/python/mine/courses/astro354/ps3/corr/corr_DM.dat'
    rbins=np.logspace(-1.0,1.301,15)
    rbins_mid=rbins[:-1]+(rbins[1:]-rbins[:-1])/2.
    # Load data
    dict_SDSS_Mr21_rspace=read_SDSS(fname_SDSS_Mr21_rspace)
    dict_SDSS_Mr20_rspace=read_SDSS(fname_SDSS_Mr20_rspace)
    dict_SDSS_Mr20_zspace=read_SDSS(fname_SDSS_Mr20_zspace)
    dict_SDSS_rand=read_SDSS(fname_SDSS_rand)
    dict_DM_data=read_DM(fname_DM_data)
    dict_DM_rand=read_DM(fname_DM_rand)

    # PART 1
    print "Part 1"
    # Compute correlation functions
    xi_SDSS_Mr21_rspace,xierr_SDSS_Mr21_rspace=corrfunc_2pnt(dict_SDSS_Mr21_rspace,dict_SDSS_rand,rbins,fcorr_SDSS_Mr21_rspace)
    xi_SDSS_Mr20_rspace,xierr_SDSS_Mr20_rspace=corrfunc_2pnt(dict_SDSS_Mr20_rspace,dict_SDSS_rand,rbins,fcorr_SDSS_Mr20_rspace)
    xi_SDSS_Mr20_zspace,xierr_SDSS_Mr20_zspace=corrfunc_2pnt(dict_SDSS_Mr20_zspace,dict_SDSS_rand,rbins,fcorr_SDSS_Mr20_zspace)
    print "    Finished correlation functions"
    # Plot
    plotname1='/home/langmm/utils/python/mine/courses/astro354/ps3/astro354_ps3_part1.png'
    fig1=plt.figure()
    # Different magnitude limits
    ax11=plt.subplot(2,1,1)
    ax11.set_xscale('log')
    ax11.set_yscale('log')
    ax11.plot(rbins_mid,xi_SDSS_Mr21_rspace,'b',label='Mr < -21 (real space)')
    ax11.plot(rbins_mid,xi_SDSS_Mr20_rspace,'r',label='Mr < -20 (real space)')
    ax11.set_xlabel('r (Mpc/h)')
    ax11.set_ylabel('xi(r)')
    ax11.set_title('Magnitude Limit Comparison')
    ax11.legend()
    # Different spaces
    ax12=plt.subplot(2,1,2)
    ax12.set_xscale('log')
    ax12.set_yscale('log')
    ax12.plot(rbins_mid,xi_SDSS_Mr20_rspace,'b',label='Mr < -20 (real space)')
    ax12.plot(rbins_mid,xi_SDSS_Mr20_zspace,'r',label='Mr < -20 (redshift space)')
    ax12.set_xlabel('r (Mpc/h)')
    ax12.set_ylabel('xi(r)')
    ax12.set_title('Real vs. Redshift Space Comparison')
    ax12.legend()
    # Save figure
    fig1.tight_layout()
    fig1.savefig(plotname1)
    plt.close(fig1)
    print '    {}'.format(plotname1)

    # PART 2
    print 'Part 2'
    # Compute correlation function
    xi_DM,xierr_DM=corrfunc_2pnt(dict_DM_data,dict_DM_rand,rbins,fcorr_DM)
    print "    Finished correlation function"
    # Compute bias functions
    bias_SDSS_Mr21_rspace=biasfunc(xi_SDSS_Mr21_rspace,xi_DM)
    bias_SDSS_Mr20_rspace=biasfunc(xi_SDSS_Mr20_rspace,xi_DM)
    # Prepare to plot
    plotname2='/home/langmm/utils/python/mine/courses/astro354/ps3/astro354_ps3_part2.png'
    fig2=plt.figure()
    # Plot dark matter correlation function
    ax21=plt.subplot(2,1,1)
    ax21.set_xscale('log')
    ax21.set_yscale('log')
    ax21.plot(rbins_mid,xi_DM)
    ax21.set_xlabel('r (Mpc/h)')
    ax21.set_ylabel('xi(r)')
    ax21.set_title('Dark Matter Correlation Function')
    # Plot comparison of bias functions
    ax22=plt.subplot(2,1,2)
    ax22.set_xscale('log')
    ax22.plot(rbins_mid,bias_SDSS_Mr21_rspace,'b',label='Mr < -21 (real space)')
    ax22.plot(rbins_mid,bias_SDSS_Mr20_rspace,'r',label='Mr < -20 (real space)')
    ax22.set_xlabel('r (Mpc/h)')
    ax22.set_ylabel('b(r)')
    ax22.set_title('Bias Comparison')
    ax22.legend()
    # Save figure
    fig2.tight_layout()
    fig2.savefig(plotname2)
    plt.close(fig2)
    print '    {}'.format(plotname2)

    # Part 3
    print 'Part 3'
    # Plot
    plotname3='/home/langmm/utils/python/mine/courses/astro354/ps3/astro354_ps3_part3.png'
    fig3=plt.figure()
    # Different magnitude limits
    ax31=plt.subplot(2,1,1)
    ax31.set_xscale('log')
    ax31.set_yscale('log')
    ax31.errorbar(rbins_mid,xi_SDSS_Mr21_rspace,yerr=xierr_SDSS_Mr21_rspace,label='Mr < -21 (real space)')
    ax31.errorbar(rbins_mid,xi_SDSS_Mr20_rspace,yerr=xierr_SDSS_Mr20_rspace,label='Mr < -20 (real space)')
    ax31.set_xlabel('r (Mpc/h)')
    ax31.set_ylabel('xi(r)')
    ax31.set_title('Magnitude Limit Comparison (with errors)')
    ax31.legend()
    # Different spaces
    ax32=plt.subplot(2,1,2)
    ax32.set_xscale('log')
    ax32.set_yscale('log')
    ax32.errorbar(rbins_mid,xi_SDSS_Mr20_rspace,yerr=xierr_SDSS_Mr20_rspace,label='Mr < -20 (real space)')
    ax32.errorbar(rbins_mid,xi_SDSS_Mr20_zspace,yerr=xierr_SDSS_Mr20_zspace,label='Mr < -20 (redshift space)')
    ax32.set_xlabel('r (Mpc/h)')
    ax32.set_ylabel('xi(r)')
    ax32.set_title('Real vs. Redshift Space Comparison (with errors)')
    ax32.legend()
    # Save figure
    fig3.tight_layout()
    fig3.savefig(plotname3)
    plt.close(fig3)
    print '    {}'.format(plotname3)
    
    
    return

# Bias function
def biasfunc(xi_gal,xi_DM):
    bias=np.sqrt(xi_gal/xi_DM)
    return bias
# 2 Point correlation function
def corrfunc_2pnt(D,R,rbins,outfile):
    keystr=['rmid','xi','err']
    typstr=len(keystr)*['float']
    # Load file if it exists
    if os.path.isfile(outfile):
        # Load data from file
        outdict=mmlio.rwtable('R',outfile,keystr=keystr,typstr=typstr)
    # Otherwise create it
    else:
        # Get pair count densities
        DD,DDerr=paircount(D,D,rbins)
        RR,RRerr=paircount(R,R,rbins)
        Nr=float(len(R['x']))
        Nd=float(len(D['x']))
        # Compute xi
        xi=((Nr/Nd)**2)*(DD/RR)-1.
        # Compute xi error
        xierr_DD=((Nr/Nd)**2)*(DDerr/RR)
        # xierr_RR=-((Nr/Nd)**2)*(RRerr/(RR**2))
        xierr_RR=0.
        xierr=np.sqrt((xierr_DD**2)+(xierr_RR**2))
        xierr=xi*np.sqrt((1./DD)+(1./RR))
        # Write data to file
        rmid=rbins[:-1]+(rbins[1:]-rbins[:-1])/2
        outdict=dict(xi=list(xi),err=list(xierr),rmid=list(rmid),keylist=keystr)
        mmlio.rwtable('W',outfile,outdict,keystr=keystr,typstr=typstr)
    # Return xi & error
    return np.array(outdict['xi']),np.array(outdict['err'])
# Count pairs
def paircount(data1,data2,rbins):
    # Get bin volumes
    binvols=(4.*np.pi/3.)*(rbins[1:]**3 - rbins[:-1]**3)
    # Initialize histogram
    tothist=np.zeros(len(rbins)-1,dtype=float)
    # Loop over positions in 1st dataset
    for ipnt in range(len(data1['x'])):
        # Get location of ith position
        ix0=data1['x'][ipnt]
        iy0=data1['y'][ipnt]
        iz0=data1['z'][ipnt]
        # Get separations of 2nd data set from ith position
        irsep=np.sqrt((data2['x']-ix0)**2 + (data2['y']-iy0)**2 + (data2['z']-iz0)**2)
        # Calc histogram
        ihist,ibins=np.histogram(irsep,bins=rbins)
        if ihist.max()==0.: print 'histogram empty'
        # Sum histogram
        tothist+=ihist.astype(float)
    if tothist.max()==0.: print 'total histogram empty'
    # Get error
    toterr=np.sqrt(tothist)
    # Return histogram & error as number density
    return tothist/binvols,toterr/binvols
# Read SDSS files
def read_SDSS(fname):
    # Define file contents
    keystr=['ra','dec','z']
    typstr=len(keystr)*['float']
    # Load data from file
    data=mmlio.rwtable('R',fname,keystr=keystr,typstr=typstr)
    for ikey in data['keylist']:
        data[ikey]=np.array(data[ikey])
    # Convert to spherical coordinates
    data['r']=z2d(data['z']) # Mpc/h
    data['phi']=(np.pi/180.)*data['ra'] # [0,2pi] rad
    data['theta']=(np.pi/180.)*(90.-data['dec']) # [0,pi] rad
    # Convert to cartesian coordinates
    data['x']=data['r']*np.sin(data['theta'])*np.cos(data['phi'])
    data['y']=data['r']*np.sin(data['theta'])*np.sin(data['phi'])
    data['z']=data['r']*np.cos(data['theta'])
    # Return data
    return data
# Read DM files
def read_DM(fname):
    # Define file contents
    keystr=['x','y','z']
    typstr=len(keystr)*['float']
    # Load data from file
    data=mmlio.rwtable('R',fname,keystr=keystr,typstr=typstr)
    for ikey in data['keylist']:
        data[ikey]=np.array(data[ikey])
    # Return data
    return data
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

if __name__ == '__main__':
    main()
