#!/usr/bin/python
import numpy as np
import matplotlib.pyplot as plt
import scipy.integrate as integrate
from mmlutils import mmlio
import scipy.special as ss

def main():
    print "Astro 310: PS6"
    # Set constants
    yr2s=3.16e7 # s
    Mpc2cm=3.0857e24
    pc2cm=3.09e18
    jansky2cgs=1.0e-23
    c_cgs=3.0e10  # cm/s
    h_cgs=6.6e-27 # erg/s
    k_cgs=1.4e-16 # erg/K
    sigma_cgs=5.67e-5 # cgs
    mu0_cgs=4.*np.pi/c_cgs
#    qe_cgs=-10. # statC
    statT_T=2.9979e6 # Tesla
    tesla_gauss=1.0e4 # Gauss

    # Problem 1
    print 'Problem 1'
    me_cgs=9.109e-28 # g
    qe_cgs=(1.602e-19/10.) # abcoulombs
    B=2.0e-6 # Gauss
    sinalpha=2.0/3.0
    Fv=1.0*jansky2cgs
    R=2.0e5*pc2cm
    d=9.0e7*pc2cm
    print 'R={} cm'.format(R)
    print 'd={} cm'.format(d)
    nu=1.0e9 # Hz
    # Part a
    print 'Part a)'
    pidx=3.6
    print '    p index={}'.format(pidx)
    # Part b
    print 'Part b)'
    Lv=4.0*np.pi*(d**2)*Fv
    print '    Lv={} erg/s/Hz'.format(Lv)
    # Part c
    print 'Part c)'
    C0_part1=3.*Lv/(4.*np.pi*(R**3))
    C0_part2=2.*np.pi*me_cgs*(c_cgs**2)*(pidx+1.)/(np.sqrt(3)*(qe_cgs**3)*B*sinalpha)
    C0_part3=1./(ss.gamma(pidx/4.+19./12.)*ss.gamma(pidx/4.-1./12.))
    C0_part4=(2.*np.pi*me_cgs*c_cgs*nu/(2.*qe_cgs*B*sinalpha))**(-(pidx-1.)/2.)
    C0=C0_part1*C0_part2*C0_part3*C0_part4
    print '    C0={}'.format(C0)
    nu_c1=nu/0.29
    gamma_c1=np.sqrt(4.*np.pi*me_cgs*c_cgs*nu_c1/(3.*qe_cgs*B*sinalpha))
    print '    gamma_c={}'.format(gamma_c1)
    # Part d
    print 'Part d)'
    E_e=-me_cgs*(c_cgs**2)*C0*(gamma_c1**(2.-pidx))/(2.-pidx)
    print '    E_e={}'.format(E_e)
    u_e=3.*E_e/(4.*np.pi*(R**3))
    print '    u_e={}'.format(u_e)
    # Part e
    print 'Part e)'
    nu_c2=(1.0e7)/0.29
    gamma_c2=np.sqrt(4.*np.pi*me_cgs*c_cgs*nu_c2/(3.*qe_cgs*B*sinalpha))
    print '    gamma_c2={}'.format(gamma_c2)
    E_e2=me_cgs*(c_cgs**2)*C0*(gamma_c1**(2.-pidx)-gamma_c2**(2.-pidx))/(2.-pidx)
    print '    E_e2={}'.format(E_e2)
    u_e2=2.*E_e2/(4.*np.pi*(R**3))
    print '    u_e2={}'.format(u_e2)
    # Part f
    print 'Part f)'
    u_B=(B**2)/(8.*np.pi)
    print '    u_B={}'.format(u_B)
    # Part g
    print 'Part g)'
    t_1=gamma_c1*me_cgs*(c_cgs**2)/Lv
    print '    t_1={}'.format(t_1)

    # Problem 3
    print 'Problem 3'
    d_casA=(3.4e3)*pc2cm # cm
    theta_casA=(4./60.)*(np.pi/180) # radians
    print '    d = {} cm'.format(d_casA)
    print '    theta = {}'.format(theta_casA)
    R_casA=d_casA*np.tan(theta_casA/2.)
    print '    R = {} cm'.format(R_casA)
    nu1_casA=0.01 # GHz
    nu2_casA=100. # GHz
    Fv_1=2.7e-14 # erg/s/cm^2/GHz
    Fv_m=-0.77
    # Part a
    print 'Part a)'
    F_casA=Fv_1*((nu2_casA**(Fv_m+1.))-(nu1_casA**(Fv_m+1.)))/(Fv_m+1.)
    print '    Ftot (casA)={} erg/s/cm^2'.format(F_casA)
    L_casA=4.*np.pi*(d_casA**2)*F_casA
    print '    Ltot (casA)={} erg/s'.format(L_casA)
    # Part b
    print 'Part b)'
    k=40.
    pidx_casA=-2.*Fv_m+1.
    print '    pidx (casA)={}'.format(pidx_casA)
    Bmin_ratio=(nu2_casA**(2.-pidx_casA)-nu1_casA**(2.-pidx_casA))/(nu2_casA**(3.-pidx_casA)-nu1_casA**(3.-pidx_casA))
    Bmin=((5.75e11)*L_casA*Bmin_ratio*9.*(1.+k)*(3.-pidx_casA)/(2.*(R**3)*(2.-pidx_casA)))**(2./7.)
    print '    Bmin (casA)={}'.format(Bmin)
    # Part c
    print 'Print c)'
    Wtotmin=(Bmin**2)*(R**3)/6. + (5.75e11)*L_casA*Bmin_ratio*(1.+k)*(3.-pidx_casA)*(Bmin**(-3./2.))/(2.-pidx_casA)
    print '    Wtot (casA)={}'.format(Wtotmin)
    K2=1.91e4
    Bminalt=(9.*K2/(2.*(R**3)))**(2./7.)
    print '    Bmin (alt)={}'.format(Bminalt)
    Wtotalt=K2*(Bminalt**(-3./2.))+(Bminalt**2)*(R**3)/6.
    print '    Wtot (alt)={}'.format(Wtotalt)
    # Part d
    print 'Print d)'
    We=(5.75e11)*L_casA*Bmin_ratio*(3.-pidx_casA)*(Bmin**(-3./2.))/(2.-pidx_casA)
    print '    We (casA)={} erg'.format(We)
    Wealt=K2*(Bminalt**(-3./2.))/(1.+k)
    print '    We (alt)={} erg'.format(Wealt)
    # Part e
    print 'Part e)'
    tage=We/L_casA
    tagealt=Wealt/L_casA
    print '    tage (casA)={} s'.format(tage)
    print '    tage (casA)={} yrs'.format(tage/yr2s)
    print '    tage (alt)={} s'.format(tagealt)
    print '    tage (alt)={} yrs'.format(tagealt/yr2s)
    rsn=1000./(tage/yr2s)
    rsnalt=1000./(tagealt/yr2s)
    print '    rate (casA)={} per year'.format(rsn)
    print '    rate (alt)={} per year'.format(rsnalt)


if __name__ == '__main__':
    main()
