#!/usr/bin/python
import numpy as np
import matplotlib.pyplot as plt
import scipy.integrate as integrate
from mmlutils import mmlio
import scipy.special as ss

def main():
    print "Astro 310: PS7"
    # Set constants
    yr2s=3.16e7 # s
    Mpc2cm=3.0857e24
    pc2cm=3.09e18
    jansky2cgs=1.0e-23
    c_cgs=3.0e10  # cm/s
    h_cgs=6.6e-27 # erg/s
    k_cgs=1.4e-16 # erg/K
    sigma_cgs=5.67e-5 # cgs
    sigmaT_cgs=0.665e-25 # cm^2
    mu0_cgs=4.*np.pi/c_cgs
#    qe_cgs=-10. # statC
    statT_T=2.9979e6 # Tesla
    tesla_gauss=1.0e4 # Gauss
    me_cgs=9.109e-28 # g
    mp_cgs=1.67e-24 # g
    qe_cgs=4.8e-10 # statC

    # Problem 1
    print 'Problem 1'
    B=1.0e-4 # Gauss
    sinalpha=2.0/3.0
    R=2.5e4*pc2cm # cm
    Ltot=1.0e45 # erg
    nu_typ=1.0e9 # Hz
    Tgas=1.0e7 # K
    print 'R={} cm'.format(R)
    # Part a
    print 'Part a)'
    uB=(B**2)/(8.*np.pi)
    print '    uB={} erg/cm^3'.format(uB)
    up=(B**2)/(6.*np.pi)
    print '    up={} erg/cm^3'.format(up)
    # Part b
    print 'Part b)'
    nu_ctyp=nu_typ/0.29
    print '    nu_ctyp={} Hz'.format(nu_ctyp)
    gamma_typ=np.sqrt(2.*np.pi*me_cgs*c_cgs*nu_ctyp/(qe_cgs*B))
    print '    gamma_typ={}'.format(gamma_typ)
    # Part c
    print 'Part c)'
    E_typ=gamma_typ*me_cgs*(c_cgs**2)
    print '    E_typ={} erg'.format(E_typ)
    n_p=up/E_typ
    print '    n_p={} cm^-3'.format(n_p)
    # Part d
    print 'Part d)'
    tcool=3.*me_cgs*c_cgs/(4.*sigmaT_cgs*gamma_typ*uB)
    print '    tcool={} s'.format(tcool)
    print '    tcool={} yrs'.format(tcool/yr2s)
    # Part e
    print 'Part e)'
    c_sound=np.sqrt(k_cgs*Tgas/me_cgs)
    print '    c_sound={} cm/s'.format(c_sound)
    t_cross=2.*R/c_sound
    print '    t_cross={} s'.format(t_cross)
    print '    t_cross={} yrs'.format(t_cross/yr2s)

    # Problem 2
    print 'PROBLEM 2'
    Sx=1.0e12 # erg/cm^2/s/arcmin
    T_eV=7.0e3 # eV
    k_eV=8.6e-5 # eV/K
    T_cgs=T_eV/k_eV
    nu_cmb=9.0e7 # Hz
    dI_I=5.0e-5
    thetaR=3.0 # arcmin
    thetaR_rad=thetaR*np.pi/(180.*60.)
    print '    T={} K'.format(T_cgs)
    print '    thetaR {} rad'.format(thetaR_rad)
    urad=Sx*np.pi*(thetaR**2)/c_cgs
    print '    urad={} erg/cm^3'.format(urad)
    part1a=me_cgs*mp_cgs*(c_cgs**2)/(4.*k_cgs*(me_cgs+mp_cgs)*sigmaT_cgs)
    part1b=np.log(dI_I+1.)
    part1=(part1a*part1b)**2
    part2=(1.4e-27)*1.2/(3.*np.pi*thetaR_rad*Sx*(T_cgs**(3./2.)))
    d_cm=part1*part2
    print '    d={} cm'.format(d_cm)
    print '    d={} pc'.format(d_cm/pc2cm)
    R_cm=d_cm*thetaR_rad
    print '    R={} pc'.format(R_cm/pc2cm)




if __name__ == '__main__':
    main()
