"""

derived
=======

Holds procedures for creating new arrays from existing ones, e.g. for
getting the radial position.
<http://code.google.com/p/pynbody/wiki/AutomagicCalculation>

"""

from . import snapshot
from snapshot import SimSnap
from . import array
from . import analysis
from . import sph
from . import config
from . import units
import numpy as np
import sys
from mmlutils import mmlmath

@SimSnap.derived_quantity
def N(self):
    """Just ones for histograms"""
    return array.SimArray(np.ones(len(self),dtype=float),1)

@SimSnap.derived_quantity
def G(self):
    """Gravitational constant in sim units"""
    G=SimArray(6.672e-8,'cm cm**2 s**-2 g**-1')
    return G.in_unitsys([self['mass'].units,self['vel'].units,self['pos'].units])

@SimSnap.derived_quantity
def Am2(self):
    return self['pm2']/self['phi']

@SimSnap.derived_quantity
def r(self):
    return ((self['pos']**2).sum(axis = 1))**(1,2)

@SimSnap.derived_quantity
def rxy(self):
    return ((self['pos'][:,0:2]**2).sum(axis = 1))**(1,2)

@SimSnap.derived_quantity
def R(self): return self['rxy']

@SimSnap.derived_quantity
def vr(self):
    """Radial velocity"""
    return (self['pos']*self['vel']).sum(axis=1)/self['r']


@SimSnap.derived_quantity
def v2(self) :
    """Squared velocity"""
    return (self['vel']**2).sum(axis=1)

@SimSnap.derived_quantity
def vt(self) :
    """Tangential velocity"""
    return np.sqrt(self['v2']-self['vr']**2)


@SimSnap.derived_quantity
def ke(self) :
    """Specific kinetic energy"""
    return 0.5*(self['vel']**2).sum(axis=1)

@SimSnap.derived_quantity
def te(self) :
    """Specific total energy"""
    return self['ke']+self['phi']

@SimSnap.derived_quantity
def Ek(self) :
    """Kinetic energy"""
    return self['mass']*self['ke']

@SimSnap.derived_quantity
def Ep(self) :
    """Potential energy"""
    return self['mass']*self['phi']

@SimSnap.derived_quantity
def Et(self) :
    """Total energy"""
    return self['Ek']+self['Ep']

@SimSnap.derived_quantity
def L(self) :
    """Angular momentum"""
    return self['mass']*self['j']

@SimSnap.derived_quantity
def Ltot(self) :
    """Total angular momentum"""
    return np.sqrt(np.sum(np.sum(self['L'],axis=0)**2))

@SimSnap.derived_quantity
def j(self) :
    """Specific angular momentum"""
    angmom = np.cross(self['pos'], self['vel']).view(array.SimArray)
    angmom.units = self['pos'].units*self['vel'].units
    return angmom

@SimSnap.derived_quantity
def j2(self) :
    """Square of the specific angular momentum"""

    return (self['j']**2).sum(axis=1)

@SimSnap.derived_quantity
def jz(self):
    """z-component of the angular momentum"""
    return self['j'][:,2]

@SimSnap.derived_quantity
def spin(self):
    """dimensionless spin parameter"""
    L=np.array(self['L'])
    M=np.array(self['mass'])
    R=np.array(self['rxyz'])
    G=float(self['G'])
    var=mmlmath.spinvec(L,M,R,G)

@SimSnap.derived_quantity
def spintot(self):
    """Total dimensionless spin parameter"""
    L=float(self['Ltot'])
    M=float(self['mass'].sum())
    R=float(self['rxyz'].max())
    G=float(self['G'])
    V=np.sqrt(G*M/R)
    return mmlmath.spinvec(L,M,V,R)

@SimSnap.derived_quantity
def vrxy(self):
    return (self['pos'][:,0:2]*self['vel'][:,0:2]).sum(axis=1)/self['rxy']

@SimSnap.derived_quantity
def vcxy(self) :
    f = (self['x']*self['vy']-self['y']*self['vx'])/self['rxy']
    f[np.where(f!=f)]=0
    return f

@SimSnap.derived_quantity
def vomega(self):
    f = (self['x']*self['vy']-self['y']*self['vx'])/(self['rxy']**2)
    f[np.where(f!=f)]=0
    return f

@SimSnap.derived_quantity
def v_mean(self):
    import sph
    
    sph.build_tree(self)
    
    nsmooth = config['sph']['smooth-particles']
    
    if config['verbose']: print 'Calculating mean velocity with %d nearest neighbours' % nsmooth

    sm = array.SimArray(np.empty((len(self['pos']),3)), self['vel'])
    self.kdtree.populate(sm, 'v_mean', nn=nsmooth, smooth=self['smooth'], rho=self['rho'])
    if config['verbose']: print 'Mean velocity done.'

    return sm 

@SimSnap.derived_quantity
def v_disp(self):
    import sph

    sph.build_tree(self)
    nsmooth = config['sph']['smooth-particles']
    self['rho']
    
    if config['verbose']: print 'Calculating velocity dispersion with %d nearest neighbours' % nsmooth

    sm = array.SimArray(np.empty(len(self['pos'])), self['vel'].units)
    self.kdtree.populate(sm, 'v_disp', nn=nsmooth, smooth=self['smooth'], rho=self['rho']) 
    if config['verbose']: print 'Velocity dispersion done.'

    return sm 

@SimSnap.derived_quantity
def age(self) :
    return self.properties['time'].in_units(self['tform'].units, **self.conversion_context()) - self['tform']

bands_available = ['u','b','v','r','i','j','h','k','U','B','V','R','I',
                   'J','H','K']

for band in bands_available :
    X = lambda s, b=str(band): analysis.luminosity.calc_mags(s,band=b)
    X.__name__ = band+"_mag"
    SimSnap.derived_quantity(X)

    X = lambda s, b=str(band): (10**(-0.4*s[b+"_mag"]))*s['rho']/s['mass']
    X.__name__ = band+"_lum_den"
    SimSnap.derived_quantity(X)

@SimSnap.derived_quantity
def theta(self) :
	"""Angle from the z axis, from [0:2pi]"""
	return np.arccos(self['z']/self['r'])

@SimSnap.derived_quantity
def alt(self) :
	"""Angle from the horizon, from [-pi/2:pi/2]"""
	return np.pi/2 - self['theta']

@SimSnap.derived_quantity
def az(self) :
	"""Angle in the xy plane from the x axis, from [-pi:pi]"""
	return np.arctan2(self['y'],self['x'])

@SimSnap.derived_quantity
def modaz(self):
    """Angle in the xy plane from the x axes folded onto [0:pi]"""
    az0=self['az']
    modaz=np.zeros(len(az0),dtype=float)
    modaz[az0>0]=np.array(az0[az0>0])
    modaz[az0<0]=np.array(az0[az0<0])+np.pi
    return array.SimArray(modaz,az0.units)

@SimSnap.derived_quantity
def cs(self):
    mu = np.zeros(len(self))
    mu[np.where(self['temp']>=1e4)[0]] = 0.59
    mu[np.where(self['temp']<1e4)[0]] = 1.3
    return np.sqrt(5.0*units.k*self['temp'] / mu/units.m_p)



@SimSnap.derived_quantity
def zeldovich_offset(self) :
    """The position offset in the current snapshot according to
    the Zel'dovich approximation applied to the current velocities.
    (Only useful in the generation or analysis of initial conditions.)"""
    from . import analysis
    bdot_by_b = analysis.cosmology.rate_linear_growth(self, unit='km Mpc^-1 s^-1')/analysis.cosmology.linear_growth_factor(self)

    a = self.properties['a']
    
    offset = self['vel']/(a*bdot_by_b)
    offset.units=self['vel'].units/units.Unit('km Mpc^-1 s^-1 a^-1')
    return offset
