#!/usr/bin/python
import numpy as np
from mmlutils import mmlpars
from mmlastro import mmlconst
from scipy import integrate

####################################################################################################################################
# METHOD TO RETURN DELTA VIR FOR A GIVEN REDSHIFT
def deltavir(z=None,cosmo=None,**exkw):
    """
    Returns delta vir for a given redshift (Bryan & Norman,1998)
    """
    # Pars input
    z=mmlpars.mml_pars(z,type=[float,np.ndarray],default=0.)
    dcosmo=mmlconst.main('cosmo',cosmo)
    # Get deltavir
    if dcosmo['omegaK']!=0:
        if dcosmo['omegaL']==0:
            om=omegaz('M',z,cosmo=cosmo)
            x=om-1.
            return (18.*(np.pi**2)+60.*x-32.*(x**2))
        else:
            raise NotImplementedError("can't compute deltavir with cosmological constant and omega!=1.0")
    else:
        om=omegaz('M',z,cosmo=cosmo)
        x=om-1.
        return (18.*(np.pi**2)+82.*x-39.*(x**2))

####################################################################################################################################
# METHOD TO CALCULATE OMEGAS FOR A GIVEN REDSHIFT
def omegaz(var,z=None,cosmo=None,**exkw):
    # Pars input
    var=mmlpars.mml_pars(var,list=['R','M','K','L'])
    z=mmlpars.mml_pars(z,type=[float,np.ndarray],default=0.)
    dcosmo=mmlconst.main('cosmo',cosmo)
    # Get supporting params
    a=1./(1.+z)
    rho0=rhocrit(0.,cosmo=cosmo)
    rhoz=rhocrit(z ,cosmo=cosmo)
    dpow={'R':4,'M':3,'K':2,'L':0}
    # Get omega
    omega=(rho0/rhoz)*dcosmo['omega'+var]*(a**-dpow[var])
    # Return output
    return omega

####################################################################################################################################
# METHOD TO RETURN CRITICAL DENSITY OF THE UNIVERSE
def rhocrit(z=None,units=None,cosmo=None,**exkw):
    """
    Returns the critical density of the universe as a function of redshift
    """
    # Pars input
    z=mmlpars.mml_pars(z,type=[float,np.ndarray],default=0.)
    dcosmo=mmlconst.main('cosmo',cosmo)
    dphys =mmlconst.main('phys' ,units)
    # Get hubble constant
    H=hubble(z,units=units,cosmo=cosmo)
    # Get rho crit
    rhoc=2.*(H**2.)/(8.*np.pi*dphys['G'])
    # Return output
    return rhoc

####################################################################################################################################
# METHOD TO RETURN THE TIME BETWEEN TWO REDSHIFTS
def t_lookback(zi=None,zf=None,cosmo=None):
    """
    Computes the lookback time (in Gyr) between two redshifts
    """
    singflag_i=False ; singflag_f=False ; arrayflag=False
    # Most distant redshift (higher number)
    if   isinstance(zi,np.ndarray): zi=list(zi) ; arrayflag=True
    elif isinstance(zi,list      ): pass
    else                          : zi=[mmlpars.mml_pars(zi,type=float,min=0.)] ; singflag_i=True # past
    # Most recent redshift (lower number)
    if   isinstance(zf,np.ndarray): zf=list(zf) ; arrayflag=True
    elif isinstance(zf,list      ): pass
    else                          : zf=[mmlpars.mml_pars(zf,type=float,default=0.,min=0.,max=zi)] ; singflag_f=True # today
    if singflag_i and sinflag_f  : singflag=True
    elif singflag_i: zi*=len(zf) ; singflag=False
    elif singflag_f: zf*=len(zi) ; singflag=False
    else                         : singflag=False
    if len(zi)!=len(zf): raise Exception('zi and zf must be the same length.')
    # Get cosmology and units
    dcosmo=mmlconst.main('cosmo',cosmo)
    dunits0={'L':mmlconst.main('astro','cgs')['pc']*(1.0e6),'V':1.0e5}
    H0_uni='{V}/{L}'
    H0_cgs=100.*dcosmo['h']*eval(H0_uni.format(**dunits0))
    # Get hubble time
    tH=(1./H0_cgs)/(mmlconst.main('astro','cgs')['yr']*(1.0e9))
    # Integrate
    E_z=lambda z: np.sqrt(dcosmo['omegaM']*((1.+z)**3)+dcosmo['omegaK']*((1.+z)**2)+dcosmo['omegaL'])
    tL=[]
    for idx in range(len(zi)): 
        if zf[idx]==0: zf[idx]=1.0e-6
        if zi[idx]==0: zi[idx]=1.0e-6
        if zf[idx]<0 or zi[idx]<0: tL.append(-1)
        else                     : tL.append(tH*integrate.quad(lambda z: 1./((1.+z)*E_z(z)),zf[idx],zi[idx])[0])
    # Return output
    if   singflag : return tL[0]
    elif arrayflag: return np.array(tL)
    else          : return tL

####################################################################################################################################
# METHOD TO RETURN HUBBLE PARAMETER AS A FUNCTION OF TIME
def hubble(z=None,units=None,cosmo=None,**exkw):
    """
    Returns the Hubble parameter as a function of redshift
    """
    # Pars input
    z=mmlpars.mml_pars(z,type=[float,np.ndarray],default=0.)
    dunits=mmlconst.main('units',units)
    dcosmo=mmlconst.main('cosmo',cosmo)
    # Get constants
    dunits0={'L':mmlconst.main('astro','cgs')['pc']*(1.0e6),'V':1.0e5}
    H0_uni='{V}/{L}'
    H0_cgs=100.*dcosmo['h']*eval(H0_uni.format(**dunits0))
    H0=H0_cgs/eval(H0_uni.format(**dunits))
    omegaR=dcosmo['omegaR']
    omegaM=dcosmo['omegaM']
    omegaK=dcosmo['omegaK']
    omegaL=dcosmo['omegaL']
    # Calculate Hubble constant
    a=1./(1.+z)
    H=H0*np.sqrt(omegaR/(a**4.)+omegaM/(a**3.)+omegaK/(a**2.)+omegaL)
    # Return output
    return H

if __name__ == '__main__':
    main()
