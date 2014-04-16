#!/usr/bin/python
import os
import numpy as np
import matplotlib.pyplot as plt
import scipy.integrate as integrate
from mmlutils import mmlio

def main():
    print "Astro 354: PS4"
    # Set values
    zrange=[0.,3.]
    nz=100
    z=np.linspace(zrange[0],zrange[1],nz)
    h=[1.,1.,1.,1.]
    omega_m=[1.00,0.25,0.25,0.25]
    omega_l=[0.00,0.75,0.75,0.75]
    omega_k=[0.00,0.00,0.00,0.00]
    w=[0.,-1.,-0.8,-1.2]
    labels=['No dark energy','"cosmological constant"','high w','low w']
    # Initialize plots
    plotname='/home/langmm/utils/python/mine/courses/astro354/ps4/astro354_ps4_part1.png'
    axdim=(3,2)
    fig=plt.figure()
    ylabel=['Dc','Dl','Da','Vc','tL']
    yunits=['Gpc/h','Gpc/h','Gpc/h','Gpc^3/h^3','Gyr']
    axdict={} ; datadict={}
    for idx in range(len(ylabel)):
        ilab=ylabel[idx] ; iuni=yunits[idx]
        axdict[ilab]=plt.subplot(axdim[0],axdim[1],idx+1)
        axdict[ilab].set_xlabel('z')
        axdict[ilab].set_ylabel('{} ({})'.format(ilab,iuni))
        datadict[ilab]=[]
    axdict['leg']=plt.subplot(axdim[0],axdim[1],axdim[0]*axdim[1])
    # Calculate variables
    for ih,iomega_m,iomega_l,iomega_k,iw,ilab in zip(h,omega_m,omega_l,omega_k,w,labels):
        # Comoving distance
        iDc=z2Dc(z,ih,iomega_m,iomega_k,iomega_l,iw)
        axdict['Dc'].plot(z,iDc,label=ilab)
        datadict['Dc'].append(iDc)
        # Luminosity distance
        iDl=z2Dl(z,ih,iomega_m,iomega_k,iomega_l,iw)
        axdict['Dl'].plot(z,iDl,label=ilab)
        datadict['Dl'].append(iDl)
        # Angular diameter distance
        iDa=z2Da(z,ih,iomega_m,iomega_k,iomega_l,iw)
        axdict['Da'].plot(z,iDa,label=ilab)
        datadict['Da'].append(iDa)
        # Comoving volume
        iVc=z2Vc(z,ih,iomega_m,iomega_k,iomega_l,iw)
        axdict['Vc'].plot(z,iVc,label=ilab)
        datadict['Vc'].append(iVc)
        # Lookback time
        itL=z2tL(z,ih,iomega_m,iomega_k,iomega_l,iw)
        axdict['tL'].plot(z,itL,label=ilab)
        datadict['tL'].append(itL)
        # Legend
        axdict['leg'].plot(z,z,label=ilab)

    axdict['leg'].legend()
    fig.tight_layout()
#    axdict['leg'].set_visible(False)
    axdict['leg'].set_axis_off()
    fig.savefig(plotname)
    plt.close(fig)
    print '    {}'.format(plotname)

    
    return

def E_z(z,h,omega_m0,omega_k0,omega_l0,w):
    E=np.sqrt(omega_m0*((1.+z)**3)+omega_k0*((1.+z)**2)+omega_l0*((1.+z)**(3.*(1.+w))))
    return E
# Hubble distance
def H02Dh(h):
    H0_kmsMpc=100.*h
    c=3.0e5 # km/s
    Dh=c/H0_kmsMpc # Mpc
    return Dh/1000. # Gpc
# Hubble time
def H02tH(h):
    H0_kmsMpc=100.*h
    Mpc2km=3.09e19
    Gyr2s=3.16e16
    H0_s=H0_kmsMpc/Mpc2km
    H0_Gyr=H0_s*Gyr2s
    tH=1./H0_Gyr
    return tH
# Proper distance Dp
def z2Dp(z,h,omega_m0,omega_k0,omega_l0,w):
    from scipy import integrate
    Dh=H02Dh(h)
    dDp=lambda zp: Dh/((1.+zp)*E_z(zp,h,omega_m0,omega_k0,omega_l0,w))
    Dp=np.zeros(z.shape) ; Dp_err=np.zeros(z.shape)
    for idx in range(len(z)): Dp[idx],Dp_err[idx]=intergrate.quad(dDp,0.,z[idx])
    return Dp
# Comoving distance Dc
def z2Dc(z,h,omega_m0,omega_k0,omega_l0,w):
    from scipy import integrate
    Dh=H02Dh(h)
    dDc=lambda zp: Dh/E_z(zp,h,omega_m0,omega_k0,omega_l0,w)
    Dc=np.zeros(z.shape) ; Dc_err=np.zeros(z.shape)
    for idx in range(len(z)): Dc[idx],Dc_err[idx]=integrate.quad(dDc,0.,z[idx])
    return Dc
# Transverse comoving distance Dm
def z2Dm(z,h,omega_m0,omega_k0,omega_l0,w):
    Dh=H02Dh(h)
    Dc=z2Dc(z,h,omega_m0,omega_k0,omega_l0,w)
    if   omega_k0 >  0.: Dm=Dh*np.sinh(np.sqrt(omega_k0)*Dc/Dh)/np.sqrt(omega_k0)
    elif omega_k0 == 0.: Dm=Dc
    elif omega_k0 <  0.: Dm=Dh*np.sin(np.sqrt(omega_k0)*Dc/Dh)/np.sqrt(np.abs(omega_k0))
    return Dm
# Luminosity distance Dl
def z2Dl(z,h,omega_m0,omega_k0,omega_l0,w):
    Dm=z2Dm(z,h,omega_m0,omega_k0,omega_l0,w)
    Dl=(1.+z)*Dm
    return Dl
# Angular diameter distance Da
def z2Da(z,h,omega_m0,omega_k0,omega_l0,w):
    Dm=z2Dm(z,h,omega_m0,omega_k0,omega_l0,w)
    Da=Dm/(1.+z)
    return Da
# Comoving volue Vc
def z2Vc(z,h,omega_m0,omega_k0,omega_l0,w,solidang=4.*np.pi):
    if omega_k0==0.:
        Dm=z2Dm(z,h,omega_m0,omega_k0,omega_l0,w)
        Vc=(4.*np.pi/2.)*(Dm**3)
    else:
        from scipy import integrate
        Dh=H02Dh(h)
        Da=z2Da(z,h,omega_m0,omega_k0,omega_l0,w)
        dVc=lambda zp: solidang*Dh*((1.+zp)**2)*(Da**2)/E_z(zp,h,omega_m0,omega_k0,omega_l0,w)
        Vc=np.zeros(z.shape) ; Vc_err=np.zeros(z.shape)
        for idx in range(len(z)): Vc[idx],Vc_err[idx]=integrate.quad(dVc,0.,z[idx])
    return Vc
# Lookback time tL
def z2tL(z,h,omega_m0,omega_k0,omega_l0,w):
    from scipy import integrate
    tH=H02tH(h)
    dtL=lambda zp: tH/((1.+zp)*E_z(zp,h,omega_m0,omega_k0,omega_l0,w))
    tL=np.zeros(z.shape) ; tL_err=np.zeros(z.shape)
    for idx in range(len(z)): tL[idx],tL_err[idx]=integrate.quad(dtL,0.,z[idx])
    return tL

if __name__ == '__main__':
    main()
