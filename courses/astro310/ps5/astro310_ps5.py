#!/usr/bin/python
import numpy as np
import matplotlib.pyplot as plt
import scipy.integrate as integrate
from mmlutils import mmlio

def main():
    print "Astro 310: PS5"
    # Set constants
    c_cgs=3.0e10  # cm/s
    h_cgs=6.6e-27 # erg/s
    k_cgs=1.4e-16 # erg/K
    sigma_cgs=5.67e-5 # cgs
    mu0_cgs=4.*np.pi/c_cgs
    me_cgs=9.11e-28 # g
    qe_cgs=-10. # statC
    statT_T=2.9979e6 # Tesla
    tesla_gauss=1.0e4 # Gauss

    # Problem 1
    print 'Problem 1'
    deltaPhi=-1.883 # radians
    deltaT=1.5576 # s
    Bmean_esu=me_cgs*(c_cgs**2)*deltaPhi/(qe_cgs*deltaT) # statT
    Bmean_si=statT_T*Bmean_esu
    Bmean_cgs=tesla_gauss*Bmean_si
    print '    <Bpara>={} statT'.format(Bmean_esu)
    print '    <Bpara>={} T'.format(Bmean_si)
    print '    <Bpara>={} Gauss'.format(Bmean_cgs)

    # Problem 2
    print 'Problem 2'
    plotname2='/home/langmm/utils/python/mine/courses/astro310/ps5/astro310_ps5_problem2.png'
    gff=1.2
    pc_cm=3.09e18 #cm
    nlam=100
    lam_i=1.0e-2 # cm
    lam_f=1.0e+2 # cm
    lamlist=np.logspace(np.log10(lam_i),np.log10(lam_f),nlam)
    cgs2jansky=1.0e-23 
    nu_i=c_cgs/lam_i
    nu_f=c_cgs/lam_f
    ne_cgs=2.0e3 # cm^-3
    Te=8.0e3 # K
    Tp=Te
    R=1.0*pc_cm #cm
    d=400.0*pc_cm #cm
    rat_i=h_cgs*nu_i/(k_cgs*Te)
    rat_f=h_cgs*nu_f/(k_cgs*Te)
    print '    hv/kT (lambda 1)={}'.format(rat_i)
    print '    hv/kT (lambda 2)={}'.format(rat_f)
    lam_knee=h_cgs*c_cgs/(k_cgs*Te)
    print '    lambda @ knee = {}'.format(lam_knee)
    epsi_ff=epsilon_ff(Te,1.,ne_cgs,ne_cgs,lamlist,gff=gff)
    alph_ff=alpha_ff(Te,1.,ne_cgs,ne_cgs,lamlist,gff=gff)
    flux_highE=cgs2jansky*(R**3)*epsi_ff/(3*(d**2))
    flux_lowE=cgs2jansky*((R/d)**2)*epsi_ff/(4*alph_ff)
    fig=plt.figure()
    axs=plt.subplot(1,1,1)
    axs.loglog(lamlist,flux_highE,label='high E')
    axs.loglog(lamlist,flux_lowE ,label='low  E')
    plt.legend()
    fig.tight_layout()
    fig.savefig(plotname2)
    plt.close(fig)
    print '    '+plotname2

    # Problem 3
    print 'Problem 3'
    G_cgs=6.67e-8 # cgs
    mh_cgs=1.67e-24 # g
    yr2sec=3.16e7 # s
    Msol2g=1.99e33 # g
    gB=2.0
    Lx=5.0e44 # erg/s
    R=0.5*(1.0e6)*pc_cm # cm
    T=1.0e8 # K
    A=(1.4e-27)*(3.6e14)*4.*np.pi/3.
    print '    A={}'.format(A)
    rho=np.sqrt(Lx/(A*(R**3)*(T**(-0.5))*gB))
    print '    rho={} g/cm^3'.format(rho)
    Mtot_cgs=(4.*np.pi*(R**3)*rho/3.)
    Mtot_sol=Mtot_cgs/Msol2g
    print '    Mtot={} Msol'.format(Mtot_sol)
    Ntot=Mtot_cgs/mh_cgs
    Etherm=3.*Ntot*k_cgs*T/2.
    print '    Etherm={} erg'.format(Etherm)
    tcool_s=Etherm/Lx
    tcool_yrs=tcool_s/yr2sec
    print '    tcool={} yrs'.format(tcool_yrs)
    Mclust_cgs=3.*k_cgs*T*R/(mh_cgs*G_cgs)
    Mclust_sol=Mclust_cgs/Msol2g
    print '    Mclust={} Msol'.format(Mclust_sol)

def alpha_ff(T,Z,ne,ni,lamb,gff=1.2):
    c_cgs=3.0e10  # cm/s
    h_cgs=6.6e-27 # erg/s
    k_cgs=1.4e-16 # erg/K
    nu=c_cgs/lamb
    exppart=(1.0-np.exp(-h_cgs*c_cgs/(lamb*k_cgs*T)))
    a_ff=(3.7e8)*(T**(-0.5))*(Z**2)*ne*ni*(nu**(-3.0))*exppart*gff
    return a_ff

def epsilon_ff(T,Z,ne,ni,lamb,gff=1.2):
    c_cgs=3.0e10  # cm/s
    h_cgs=6.6e-27 # erg/s
    k_cgs=1.4e-16 # erg/K
    nu=c_cgs/lamb
    exppart=np.exp(-h_cgs*c_cgs/(lamb*k_cgs*T))
    e_ff=(6.8e-38)*(T**(-0.5))*(Z**2)*ne*ni*exppart*gff
    return e_ff
    


if __name__ == '__main__':
    main()
