#!/usr/bin/python
import numpy as np
import matplotlib.pyplot as plt
import scipy.integrate as integrate
from mmlutils import mmlio

def main():
    print "Astro 310: PS4"
    # Set constants
    c_cgs=3.0e10  # cm/s
    h_cgs=6.6e-27 # erg/s
    k_cgs=1.4e-16 # erg/K
    sigma_cgs=5.67e-5 # cgs
    mu0_cgs=4.*np.pi/c_cgs

    # Problem 1
    print 'Problem 1'
    lamb=np.array([3302.98,5895.94])
    widt=np.array([0.067  ,0.560  ])
    fact=np.array([0.0049 ,0.3250 ])
    logW=np.log10(widt/lamb)
    print '    log(W/lambda)={}'.format(logW)
    lgNf=np.array([12.5   ,14.5   ])
    Nfla=5000.0*(10.**lgNf)
    N_ang=Nfla/(fact*lamb)
    N_cm=N_ang
#    print '    N={} (ang^-2)'.format(N_ang)
#    ang2cm=1.0e-8
#    N_cm=N_ang/(ang2cm**2)
    print '    N={} (cm^-2)'.format(N_cm)
    print '    mean(N)={} (cm^-2)'.format(np.mean(N_cm))

    # Problem 2
    print 'Problem 2'
    day2sec=24.*60.*60.
    msol2g=1.99e33
    km2cm=1.0e5
    pc2cm=3.09e18
    P_s=0.033 # s
    dP_dt_sday=91.0e-9 #s/day
    dP_dt_nounit=dP_dt_sday/day2sec
    d_pc=2000. # pc
    d_cm=d_pc*pc2cm
    M_msol=1.3 # Msol
    M_gram=M_msol*msol2g
    R_km=10. # km
    R_cm=R_km*km2cm
    print '    Part a)'
    tau_day=P_s/(2*dP_dt_sday)
    tau_yrs=tau_day/365.
    print '        dP_dt={}'.format(dP_dt_nounit)
    print '        tau={} days'.format(tau_day)
    print '        tau={} years'.format(tau_yrs)
    bday=2013.-tau_yrs
    print '        bday={}'.format(bday)
    print '    Part b)'
    I_cgs=(2./5.)*M_gram*(R_cm**2)
    w_cgs=2.*np.pi/P_s
    dw_dt_cgs=-2.*np.pi*dP_dt_nounit/(P_s**2)
    d_cgs=np.sqrt(-3.*I_cgs*dw_dt_cgs*(c_cgs**3)/(2.*(w_cgs**3)))
    print '        d={} (cgs)'.format(d_cgs)
    # sqrt(mrrrrrsss/sssss)
    # sqrt(mrrrrr/ss)
    # (rr/s)sqrt(mr)=rsqrt(r)E
    # sqrt(Errr)
    # r/
    print '    Part c)'
    B_cgs=mu0_cgs*d_cgs/(2.*np.pi*(R_cm**3))
    print '        B={} (cgs)'.format(B_cgs)
    print '    Part d)'
    u_cgs=(w_cgs**4)*(d_cgs**2)/(4.*np.pi*(R_cm**2)*(c_cgs**5))
    print'         u={} (erg/cm^3)'.format(u_cgs)
    print '    Part e)'
    flux=(w_cgs**4)*(d_cgs**2)/(4.*np.pi*(d_cm**2)*(c_cgs**3))
    print '        f={} (erg/s/cm^2)'.format(flux)
    print '    Part f)'
    R_earth=6.38e8
    r_earth=pc2cm
    area=np.pi*(R_earth**2)
    flux_earth=(w_cgs**4)*(d_cgs**2)/(4.*np.pi*(r_earth**2)*(c_cgs**3))
    lum_earth=area*flux_earth
    print '        L={} (erg/s)'.format(lum_earth)
    T=(flux/(4.*sigma_cgs))**(1./4.)
    print '        T={} (K)'.format(T)


if __name__ == '__main__':
    main()
