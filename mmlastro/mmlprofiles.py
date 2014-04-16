#!/usr/bin/python
from mmlutils import mmlpars,mmlio,mmlmath
import mmlcosmo
import math
import numpy as np
from scipy import integrate,optimize,special

####################################################################################################################################
####################################################################################################################################
# METHODS TO RETURN LISTS
LIST_PROFILES_SPH=['NFW','HERNQUIST','ISOTHERMAL','SPH2POW','JAFFE','EINASTO']
LIST_PROFILES_CYL=['HIYASHI','FLATDISK']
LIST_PROFILES=LIST_PROFILES_SPH+LIST_PROFILES_CYL
DICT_PROFABRV={'LH':'HERNQUIST','NF':'NFW','LH':'HERNQUIST'}
LIST_PROFMETH=['rho','mass','pot']
DICT_PROFPAR={'ISOTHERMAL':['ms','rs','rt'],
              'HIYASHI'   :['ms','rs','zs'],
              'EINASTO'   :['ms','rs','a'],
              'SPH2POW'   :['ms','rs','a','b'],
              'OTHER'     :['ms','rs']}

####################################################################################################################################
####################################################################################################################################
# METHODS TO RETURN PROFILE INFO

####################################################################################################################################
# METHOD TO RETURN GENERAL PROFILE
def main(profid,r=None,z=None,a=None,b=None,**inkw):
    """
    Provides command line access to profiles
    """
    # Pars input
    if profid.upper() in DICT_PROFABRV: profid=DICT_PROFABRV[profid.upper()]
    profid=mmlpars.mml_pars(profid.upper(),list=LIST_PROFILES)
    # Proceed based on profile
    if   profid=='SPH2POW'   : outvar=profile_sph2pow(r,a,b,**inkw)
    elif profid=='JAFFE'     : outvar=profile_jaffe(r,**inkw)
    elif profid=='NFW'       : outvar=profile_nfw(r,**inkw)
    elif profid=='HERNQUIST' : outvar=profile_hernquist(r,**inkw)
    elif profid=='ISOTHERMAL': outvar=profile_isothermal(r,**inkw)
    elif profid=='HIYASHI'   : outvar=profile_hiyashi(r,z,**inkw)
    elif profid=='FLATDISK'  : outvar=profile_flatdisk(r,**inkw)
    elif profid=='EINASTO'   : outvar=profile_einasto(r,a,**inkw)
    else: raise Exception('Invalid profile ID: {}'.format(profid))
    # Return output
    return outvar

####################################################################################################################################
# METHOD TO FIT PROFILE
def fit(profid,y,r=None,z=None,method=None,G=1.0,plotfile=None,**inkw):
    """
    Returns fit parameters for a profile
    """
    # Pars input
    y=mmlpars.mml_pars(y,type=np.ndarray)
    if profid.upper() in DICT_PROFABRV: profid=DICT_PROFABRV[profid.upper()]
    profid=mmlpars.mml_pars(profid.upper(),list=LIST_PROFILES)
    method=mmlpars.mml_pars(method,list=LIST_PROFMETH,default='mass')
    if profid in DICT_PROFPAR.keys(): parlist=DICT_PROFPAR[profid]
    else                            : parlist=DICT_PROFPAR['OTHER']
    # Get initial guess
    p0=[] ; bounds=[]
    for ipar in parlist: 
        if ipar in inkw: p0.append(inkw[ipar])
        else           : p0.append(1.)
        if ipar in ['ms','rs','rt','zs']: bounds.append((0.  ,None))
        else                            : bounds.append((None,None))
    # Handle function bassed on profile
    if  profid=='HIYASHI'   : 
        if r==None and z==None: raise Exception('Either r or z must be provided.')
        elif r==None: defval=z.max()
        elif z==None: defval=r.max()
        if method=='rho': vardef=np.zeros(y.shape)
        else            : vardef=defval*np.ones(y.shape)
        r=mmlpars.mml_pars(r,default=vardef,type=np.ndarray,shape=y.shape)
        z=mmlpars.mml_pars(z,default=vardef,type=np.ndarray,shape=y.shape)
        inargs=(r,z,y)
        fitfunc=lambda p,ri,zi: main(profid,method=method,r=ri,z=zi,G=G,**{parlist[idx]:p[idx] for idx in range(len(parlist))})
        errfunc_lsq=lambda p,ri,zi,yi: fitfunc(p,ri,zi)-yi
        errfunc_min=lambda p,ri,zi,yi: sum(errfunc_lsq(p,ri,zi,yi)**2)
    else: 
        r=mmlpars.mml_pars(r,type=np.ndarray,shape=y.shape)
        inargs=(r,y)
        fitfunc=lambda p,ri: main(profid,method=method,r=ri,G=G,**{parlist[idx]:p[idx] for idx in range(len(parlist))})
        errfunc_lsq=lambda p,ri,yi: fitfunc(p,ri)-yi
        errfunc_min=lambda p,ri,yi: sum(errfunc_lsq(p,ri,yi)**2)
    # Fit
    if hasattr(optimize,'minimize'):
        res=optimize.minimize(errfunc_min,p0[:],args=inargs,bounds=bounds)
        p1=res.x ; success=res.success
    else:
        p1,success=mmlmath.leastsq(errfunc_lsq,p0[:],args=inargs,bounds=bounds)
    # Plot
    if plotfile:
        import matplotlib.pyplot as plt
        plt.close('all')
        plt.plot(r,y,'r.')
        plt.plot(r,fitfunc(p1,r),'b-')
        plt.savefig(plotfile)
    # Create dictionary of parameters
    outpar={parlist[idx]:p1[idx] for idx in range(len(parlist))}
    return outpar

####################################################################################################################################
# METHOD TO RETURN EINASTO PROFILE
def profile_einasto(r,a,method=None,ms=None,rs=None,G=1.0):
    """
    Returns enclosed mass or density at a given radius for a two-power profile
       [By default scale mass and radius are set to unity]
    """
    # Pars input
    r=mmlpars.mml_pars(r,min=0.,type=[float,np.ndarray])
    a=mmlpars.mml_pars(a,type=float)
    method=mmlpars.mml_pars(method,list=LIST_PROFMETH,default=LIST_PROFMETH[0])
    ms=mmlpars.mml_pars(ms,min=0.,type=float,default=1.0)
    rs=mmlpars.mml_pars(rs,min=0.,type=float,default=1.0)
    # Set scale density & dimensionless parameter
    rho0=ms/(4.*np.pi*(rs**3.))
    s=r/rs
    # Calculate profile
    if   method=='rho' : out=rho0*np.exp(-(2./a)*(s**a-1.))
    elif method=='mass': out=integrate.quad(lambda x: 4.*np.pi*(x**2.)*profile_einasto(x,a,method='rho',ms=ms,rs=rs,G=G),0.,r)[0]
    elif method=='pot' : out=integrate.quad(lambda x: -(G/x**2.)*profile_einasto(x,a,method='mass',ms=ms,rs=rs,G=G),r,float('inf'))[0]
    else: raise Exception('Invalid profile method: {}'.format(method))
    # Return output
    return out

####################################################################################################################################
# METHOD TO RETURN SPHERICAL TWO-POWER PROFILE
def profile_sph2pow(r,a,b,method=None,ms=None,rs=None,G=1.0,**exkw):
    """
    Returns enclosed mass or density at a given radius for a two-power profile
       [By default scale mass and radius are set to unity]
    """
    # Pars input
    r=mmlpars.mml_pars(r,min=0.,type=[float,np.ndarray])
    a=mmlpars.mml_pars(a,type=float)
    b=mmlpars.mml_pars(b,type=float)
    method=mmlpars.mml_pars(method,list=LIST_PROFMETH,default=LIST_PROFMETH[0])
    ms=mmlpars.mml_pars(ms,min=0.,type=float,default=1.0)
    rs=mmlpars.mml_pars(rs,min=0.,type=float,default=1.0)
    # Set scale density & dimensionless parameter
    rho0=ms/(4.*np.pi*(rs**3.))
    s=r/rs
    # Calculate profile
    if   method=='rho' : out=rho0/((s**a)*((1.+s)**(b-a)))
    elif method=='mass': out=integrate.quad(lambda x: 4.*np.pi*(x**2.)*profile_sph2pow(x,a,b,method='rho',ms=ms,rs=rs,G=G),0.,r)[0]
    elif method=='pot' : out=integrate.quad(lambda x: -(G/x**2.)*profile_sph2pow(x,a,b,method='mass',ms=ms,rs=rs,G=G),r,float('inf'))[0]
    else: raise Exception('Invalid profile method: {}'.format(method))
    # Return output
    return out

####################################################################################################################################
# METHOD TO RETURN NFW PROFILE
def profile_nfw(r,method=None,ms=None,rs=None,G=1.0,**exkw):
    """
    Returns enclosed mass or density at a given radius for an NFW profile
       [By default scale mass and radius are set to unity]
    """
    # Pars input
    r=mmlpars.mml_pars(r,min=0.,type=[float,np.ndarray])
    method=mmlpars.mml_pars(method,list=LIST_PROFMETH,default=LIST_PROFMETH[0])
    ms=mmlpars.mml_pars(ms,min=0.,type=float,default=1.0)
    rs=mmlpars.mml_pars(rs,min=0.,type=float,default=1.0)
    # Set dimensionless parameter
    rho0=ms/(4.*np.pi*(rs**3.))
    s=r/rs
    # Calculate profile
    if   method=='rho' : out=rho0/(s*((1.0+s)**2.0))
    elif method=='mass': out=ms*(np.log(1.+s)-s/(1.+s))
    elif method=='pot' : out=-(G*ms/rs)*(np.log(1.+s)/s)
    else: raise Exception('Invalid profile method: {}'.format(method))
    # Return output
    return out

####################################################################################################################################
# METHOD TO RETURN HERNQUIST PROFILE
def profile_hernquist(r,method=None,ms=None,rs=None,G=1.0,**exkw):
    """
    Returns enclosed mass or density at a given radius for a Hernquist profile
       [By default scale mass and radius are set to unity]
    """
    # Pars input
    r=mmlpars.mml_pars(r,min=0.,type=[float,np.ndarray])
    method=mmlpars.mml_pars(method,list=LIST_PROFMETH,default=LIST_PROFMETH[0])
    ms=mmlpars.mml_pars(ms,min=0.,type=float,default=1.0)
    rs=mmlpars.mml_pars(rs,min=0.,type=float,default=1.0)
    # Set dimensionless parameter
    rho0=ms/(4.*np.pi*(rs**3.))
    s=r/rs
    # Calculate profile
    if   method=='rho' : out=rho0/(s*((1.+s)**3.))
    elif method=='mass': out=ms*(s**2.)/(2.*((1.+s)**2.))
    elif method=='pot' : out=-(G*ms/rs)/(2.*(1+s))
    else: raise Exception('Invalid profile method: {}'.format(method))
    # Return output
    return out

####################################################################################################################################
# METHOD TO RETURN JAFFE PROFILE
def profile_jaffe(r,method=None,ms=None,rs=None,G=1.0,**exkw):
    """
    Returns enclosed mass or density at a given radius for a Jaffe profile
       [By default scale mass and radius are set to unity]
    """
    # Pars input
    r=mmlpars.mml_pars(r,min=0.,type=[float,np.ndarray])
    method=mmlpars.mml_pars(method,list=LIST_PROFMETH,default=LIST_PROFMETH[0])
    ms=mmlpars.mml_pars(ms,min=0.,type=float,default=1.0)
    rs=mmlpars.mml_pars(rs,min=0.,type=float,default=1.0)
    # Set dimensionless parameter
    rho0=ms/(4.*np.pi*(rs**3.))
    s=r/rs
    # Calculate profile
    if   method=='rho' : out=rho0/((s**2.)*((1.+s)**2.))
    elif method=='mass': out=ms*s/(1.+s)
    elif method=='pot' : out=-(G*ms/rs)*np.log(1.+1./s)
    else: raise Exception('Invalid profile method: {}'.format(method))
    # Return output
    return out

####################################################################################################################################
# METHOD TO RETURN ISOTHERMAL PROFILE
def profile_isothermal(r,method=None,ms=None,rs=None,rt=None,G=1.0,**exkw):
    """
    Returns enclosed mass or density at a given radius for an isothermal profile
       [By default scale mass and radius are set to unity]
    """
    # Pars input
    r=mmlpars.mml_pars(r,min=0.,type=[float,np.ndarray])
    method=mmlpars.mml_pars(method,list=LIST_PROFMETH,default=LIST_PROFMETH[0])
    ms=mmlpars.mml_pars(ms,min=0.,type=float,default=1.0)
    rs=mmlpars.mml_pars(rs,min=0.,type=float,default=1.0)
    rt=mmlpars.mml_pars(rt,min=0.,type=float,default=1.0)
    # Set dimensionless parameter
    q=rs/rt ; s=r/rs ; t=r/rt
    alpha=1./(1.-np.sqrt(np.pi)*q*np.exp(q**2)*(1.-special.erf(q)))
    rho0=ms*alpha/(2.*(np.pi**1.5)*(rs**2)*rt)
    # Calculate profile
    if   method=='rho' : out=rho0*np.exp(-(t**2.))/(1.+(s**2.))
    elif method=='mass': out=integrate.quad(lambda x: 4.*np.pi*(x**2.)*profile_isothermal(x,method='rho',ms=ms,rs=rs,rt=rt,G=G),0.,r)[0]
    elif method=='pot' : out=integrate.quad(lambda x: -(G/x**2.)*profile_isothermal(x,method='mass',ms=ms,rs=rs,rt=rt,G=G),r,float('inf'))[0]
    else: raise Exception('Invalid profile method: {}'.format(method))
    # Return output
    return out

####################################################################################################################################
# METHOD TO RETURN HIYASHI DISK PROFILE
def profile_hiyashi(r,z,method=None,ms=None,rs=None,zs=None,G=1.0,**exkw):
    """
    Returns enclosed mass or density at a given radius for a Hiyashi disk
    """
    # Pars input
    r=mmlpars.mml_pars(r,min=0.,type=[float,np.ndarray])
    z=mmlpars.mml_pars(z,type=[float,np.ndarray])
    method=mmlpars.mml_pars(method,list=LIST_PROFMETH,default=LIST_PROFMETH[0])
    ms=mmlpars.mml_pars(ms,min=0.,type=float,default=1.0)
    rs=mmlpars.mml_pars(rs,min=0.,type=float,default=1.0)
    zs=mmlpars.mml_pars(zs,min=0.,type=float,default=1.0)
    # Set dimensionless parameter and scale density
    s=r/rs ; t=z/zs
    rho0=ms/(4.*np.pi*zs*(rs**2))
    # Calculate profile
    if   method=='rho' : out=rho0*np.exp(-s)*4./((np.exp(-t)+np.exp(t))**2)
    # elif method=='mass':
    #     nele=1
    #     if isinstance(r,np.ndarray): nele=max(nele,len(r))
    #     if isinstance(z,np.ndarray): nele=max(nele,len(z))
    #     if isinstance(r,float): r=np.array(nele*[r])
    #     if isinstance(z,float): z=np.array(nele*[z])
    #     out=np.zeros(nele,dtype=float)
    #     for idx in range(nele):
    #         out[idx]=integrate.dblquad(lambda x1,x2: 2.*np.pi*x1*profile_hiyashi(x1,x2,method='rho',ms=ms,rs=rs,zs=zs,G=G),
    #                                    -z[idx],z[idx],lambda x: 0.,lambda x:r[idx])[0]
    #     if len(out)==1: out=out[0]
    elif method=='mass': out=ms*(1.-(s+1.)*np.exp(-s))
    elif method=='pot' : out=integrate.quad(lambda x: -(G/x**2.)*profile_hiyashi(x,z,method='mass',ms=ms,rs=rs,zs=zs,G=G),r,float('inf'))[0]
    else: raise Exception('Invalid profile method: {}'.format(method))
    # Return output
    return out

####################################################################################################################################
# METHOD TO RETURN FLAT DISK PROFILE
def profile_flatdisk(r,method=None,ms=None,rs=None,G=1.0,**exkw):
    """
    Returns enclosed mass or density at a given radius for a flat disk
    """
    # Pars input
    r=mmlpars.mml_pars(r,min=0.,type=[float,np.ndarray])
    method=mmlpars.mml_pars(method,list=LIST_PROFMETH,default=LIST_PROFMETH[0])
    ms=mmlpars.mml_pars(ms,min=0.,type=float,default=1.0)
    rs=mmlpars.mml_pars(rs,min=0.,type=float,default=1.0)
    # Set dimensionless parameter and scale density
    s=r/rs
    rho0=ms/(np.pi*(rs**2))
    # Calculate profile
    if   method=='rho' : out=rho0*np.exp(-s)
    elif method=='mass': out=ms*(1.-(s+1.)*np.exp(-s))
    elif method=='pot' : out=integrate.quad(lambda x: -(G/x**2.)*profile_flatdisk(x,z,method='mass',ms=ms,rs=rs,zs=zs,G=G),r,float('inf'))[0]
    else: raise Exception('Invalid profile method: {}'.format(method))
    # Return output
    return out

####################################################################################################################################
####################################################################################################################################
# METHODS TO CALCULATE PROFILE STUFF

####################################################################################################################################
# RETURN RADIUS OF PEAK VELOCITY
def rvpeak(profile,**inkw):
    """
    Returns the radius at which the velocity peaks
    """
    vel=lambda x: -np.sqrt((1./x)*main(profile,r=x,method='mass',**inkw))
    r=optimize.fmin(vel,2.,disp=False)[0]
    return r

####################################################################################################################################
####################################################################################################################################
# METHODS CONVERT BETWEEN PARAMETERS
def convert_rs(prof1,prof2,rs1,profkw1={},profkw2={}):
    """
    Converts between scale radii by placing vpeak at the same radius
    """
    rv1=rvpeak(prof1,**profkw1) # in rs1
    rv2=rvpeak(prof2,**profkw2) # in rs2
    rs2=rs1*rv1/rv2
    rv1c=rvpeak(prof1,rs=rs1,**profkw1)
    rv2c=rvpeak(prof2,rs=rs2,**profkw2)
    if not np.allclose(rv1c,rv2c): raise Exception('Error in converting scale radius')
    return rs2
def convert_ms(prof1,prof2,ms1,profkw1={},profkw2={}):
    """
    Converts between scale masses by equating average densities inside the vpeak radius
    """
    method='mass'
    rs1=1.0
    rs2=convert_rs(prof1,prof2,rs1,profkw1=profkw1,profkw2=profkw2)
    rv1=rvpeak(prof1,rs=rs1,**profkw1) # in rs1
    rv2=rvpeak(prof2,rs=rs2,**profkw2) # in rs2
    menc1=main(prof1,method=method,r=rv1,rs=rs1,**profkw1) # in ms1
    menc2=main(prof2,method=method,r=rv2,rs=rs2,**profkw2) # in ms2
    ms2=ms1*menc1/menc2
    menc1c=main(prof1,method=method,r=rv1,rs=rs1,ms=ms1,**profkw1)
    menc2c=main(prof2,method=method,r=rv2,rs=rs2,ms=ms2,**profkw2)
    if not np.allclose(menc1c,menc2c): raise Exception('Error in converting scale mass')
    return ms2
def menc2ms(profile,menc,**profkw):
    """
    Converts enclosed mass to scale mass
    """
    menc=mmlpars.mml_pars(menc,type=float,min=0.)
    profkw['ms']=1.0
    menc0=main(profile,method='mass',**profkw)
    ms=menc/menc0
    profkw['ms']=ms
    mtot=main(profile,method='mass',**profkw)
    if not np.allclose(menc,mtot): raise Exception('Error in converting enclosed mass to scale mass')
    return ms
def mvir2ms(profile,mvir,c=None,rs=None,**profkw):
    """
    Converts virial mass to scale mass for a given profile
    """
    mvir=mmlpars.mml_pars(mvir,type=float,min=0.)
    if isinstance(c,float): c=mmlpars.mml_pars(c,type=float,min=0.)
    else:
        rs=mmlpars.mml_pars(rs,type=float,min=0.)
        rs_nfw=convert_rs(profile,'NFW',rs,profkw1=profkw)
        rvir=mvir2rvir(mvir,**profkw)
        c=rvir/rs_nfw
    ms_nfw=mvir/(np.log(1.+c)-c/(1.+c))
    if profile in ['NFW','NF']: ms=ms_nfw
    else                      : ms=convert_ms('NFW',profile,ms_nfw,profkw2=profkw)
    return ms
def ms2mvir(profile,ms,c=None,rs=None,**profkw):
    """
    Converts scale mass for a given profile to virial mass
    """
    ms=mmlpars.mml_pars(ms,type=float,min=0.)
    if isinstance(c,float): c=mmlpars.mml_pars(c,type=float,min=0.)
    else:
        rs=mmlpars.mml_pars(rs,type=float,min=0.)
        rs_nfw=convert_rs(profile,'NFW',rs,profkw1=profkw)
        rvir=mvir2rvir(mvir,**profkw)
        c=rvir/rs_nfw
    ms_nfw=convert_ms(profile,'NFW',ms,profkw1=profkw)
    mvir=ms_nfw*(np.log(1.+c)-c/(1.+c))
    return mvir
def rvir2mvir(rvir,rhoc=None,deltavir=None,units='galaxy',**inkw):
    """
    Converts virial radius to virial mass
    """
    # Pars input
    rvir=mmlpars.mml_pars(rvir,type=float,min=0.)
    rhoc=mmlpars.mml_pars(rhoc,type=float,min=0.,default=mmlcosmo.rhocrit(units=units,**inkw))
    deltavir=mmlpars.mml_pars(deltavir,type=float,min=0.,default=mmlcosmo.deltavir(**inkw))
    # Calculate Mvir
    mvir=deltavir*rhoc*(4./3.)*np.pi*(rvir**3)
    return mvir
def mvir2rvir(mvir,rhoc=None,deltavir=None,units='galaxy',**inkw):
    """
    Converts virial mass to virial radius
    """
    # Pars input
    mvir=mmlpars.mml_pars(mvir,type=float)
    rhoc=mmlpars.mml_pars(rhoc,type=float,min=0.,default=mmlcosmo.rhocrit(units=units,**inkw))
    deltavir=mmlpars.mml_pars(deltavir,type=float,min=0.,default=mmlcosmo.deltavir(**inkw))
    # Calculate Rvir
    rvir=(mvir/(deltavir*rhoc*(4./3.)*np.pi))**(1./3.)
    return rvir

if __name__ == '__main__':
    main()
