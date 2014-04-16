#!/usr/bin/python
import numpy as np
import matplotlib.pyplot as plt
import scipy.integrate as integrate
from mmlutils import mmlio

def main():
    print "Astro 314: PS2"
    # Set constants
    c=3.0e10  # cm/s
    h=6.6e-27 # erg/s
    k=1.4e-16 # erg/K
    sigma=5.67e-5 # cgs

    # PROBLEM 3
    print "Problem 3"
    # PART A
    print "Part a)"
    T=314. # K
    nlambda=100
    lambda_min=1.0e-06 # cm
    lambda_max=8.0e-03 # cm
##     # Compute black body
##     lambda_list=np.linspace(lambda_min,lambda_max,nlambda)
##     B_lambda=plank(T,lambda_list,space='lambda')
##     # Plot
##     plotname1='/home/langmm/utils/python/mine/courses/astro310/ps2/astro310_ps2_problem1a.png'
##     fig1=plt.figure()
##     ax11=plt.subplot(1,1,1)
##     ax11.plot(lambda_list,B_lambda,'.')
##     ax11.set_xlabel('Wavelength (cm)')
##     ax11.set_ylabel('Specific Intensity (erg/(s*cm^2*str*Hz))')
##     fig1.tight_layout()
##     fig1.savefig(plotname1)
##     plt.close(fig1)
##     print '    {}'.format(plotname1)
    # PART C
    print "Part c)"
    import scipy.optimize as optimize
    zero=optimize.newton(lambda x: (5.*np.exp(-x)+x-5.),2.82)
#    zero=optimize.newton(lambda x: ((3.-x)*np.exp(x)-3.),2.82)
#    zero=optimize.newton(lambda x: ((x-5.)*np.exp(x)+5.),2.82)
    zero_scaled=h*c/(zero*k*T)
    print '    zero,unscaled = {}'.format(zero)
    print '    zero,scaled   = {} cm'.format(zero_scaled)

    print "Problem 4"
    print "Part a)"
    T_eye=310.
    r_eye=1.
    lambda_eye = h*c/(zero*k*T_eye)
    L_eye = 4*np.pi*(r_eye**2)*sigma*(T_eye**4)
    E_eyephot = h*c/lambda_eye
    photrate_eye = L_eye/E_eyephot
    print "    lambda_eye  = {} cm".format(lambda_eye)
    print "    L_eye       = {} erg/s".format(L_eye)
    print "    E_photon    = {} erg".format(E_eyephot)
    print "    photon rate = {} s^-1".format(photrate_eye)
    print "Part b)"
    L_bulb=1.0e9
    lambda_bulb=5.0e-8
    r_pupil=100.
    R_pupil=0.5
    E_bulbphot=h*c/lambda_bulb
    photrate_bulb = L_bulb/E_bulbphot
    photrate_pupil = photrate_bulb*(np.pi*(R_pupil**2))/(4*np.pi*(r_pupil**2))
    print "    bulb photon rate (tot)   = {} z^-1".format(photrate_bulb)
    print "    bulb photon rate (pupil) = {} z^-1".format(photrate_pupil)

# Blackbody
def plank(T,x,space='lambda'):
    c=3.0e10  # cm/s
    h=6.6e-27 # erg/s
    k=1.4e-16 # erg/K
    if space=='lambda':
        B=(2.*h*(c**2)/(x**5))/(np.exp(h*c/(x*k*T))-1.)
    elif space=='nu':
        B=(2.*h*(x**3)/(c**2))/(np.exp(h*x/(k*T))-1)
    else: raise Exception('[plank] Invalid space: {}'.format(space))
    return B

if __name__ == '__main__':
    main()
