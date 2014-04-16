"""

halo
====

Functions for dealing with and manipulating halos in simulations.


"""

from .. import filt, util, config, snapshot, array
from . import cosmology
import numpy as np
import math

def list_centers_by_var(sim,**kwargs):
    """
    Return a list of halo centers
    """
    gals=kwargs.pop('gals',sim.galaxies())
    # Determine class of snapshot
    pNbody = False if issubclass(sim.__class__,snapshot.SimSnap) else True
    # Loop over galaxies
    c=[]
    for g in gals:
        if pNbody: isim=sim.select(g)
        else     : isim=sim[g]
        c.append(center_by_var(isim,**kwargs))
    # Return centers
    return c

def center_by_var(sim, pvars=['x','y','z'], wvars=['mass']) :
    """
    Return the center of the sim in certain vars
    """
    # Determine class of snapshot
    pNbody = False if issubclass(sim.__class__,snapshot.SimSnap) else True
    # Remove empty variables
    pvars=[x for x in pvars if x!='']
    wvars=[x for x in wvars if x!='']
    # Get weights and positions
    ws=[] ; ps=[]
    if pNbody:
        for wvar in wvars: ws.append(array.SimArray(sim.get_var(wvar),sim.get_units(wvar)))
        for pvar in pvars: ps.append(array.SimArray(sim.get_var(pvar),sim.get_units(pvar)))
        Nw=sim.nbody
    else:
        for wvar in wvars: ws.append(sim[wvar])
        for pvar in pvars: ps.append(sim[pvar])
        Nw=len(sim)
    # Get weights
    w=np.ones(Nw,dtype=float)
    for iw in ws: w*=iw
    wtot=w.sum()
    # Add 'positions'
    p=[]
    for ip in ps:
        iptot=array.SimArray(np.sum(w*ip)/wtot,ip.units)
        p.append(iptot)
    # Return
    return tuple(p)


def center_of_mass(sim, var="pos") : 
    """

    Return the centre of mass of the SimSnap

    """
    mtot = sim["mass"].sum()
    p = np.sum(sim["mass"]*sim[var].transpose(), axis=1)/mtot

    p.units = sim[var].units # otherwise behaviour is numpy version dependent

    return p # only return position to be consistent with other functions in halo.py

def center_of_mass_velocity(sim) :
    """

    Return the center of mass velocity of the SimSnap

    """
    return center_of_mass(sim, var="vel")
    # mtot = sim["mass"].sum()
    # v = np.sum(sim["mass"]*sim["vel"].transpose(), axis=1)/mtot
    # v.units = sim["vel"].units # otherwise behaviour is numpy version dependent

    # return v

def shrink_sphere_center(sim, r=None, shrink_factor = 0.7, min_particles = 100, verbose=False) :
    """
    
    Return the center according to the shrinking-sphere method of
    Power et al (2003)
    
    """
    x = sim

    if r is None :
        # use rough estimate for a maximum radius
        # results will be insensitive to the exact value chosen
        r = (sim["x"].max()-sim["x"].min())/2
    com=None
    while len(x)>min_particles or com is None :
        com = center_of_mass(x)#, cov = center_of_mass(x)
        r*=shrink_factor
        x = sim[filt.Sphere(r, com)]
        if verbose:
            print com,r,len(x)
    return com

def virial_radius(sim, cen=None, overden=178, r_max=None) :
    """
    
    Calculate the virial radius of the halo centerd on the given
    coordinates.

    This is here defined by the sphere centerd on cen which contains a
    mean density of overden * rho_c_0 * (1+z)^3.

    """

    if r_max is None :
        r_max = (sim["x"].max()-sim["x"].min())
    else :
        if cen is not None :
            sim = sim[filt.Sphere(r_max,cen)]
        else :
            sim = sim[filt.Sphere(r_max)]

    r_min = 0.0

    if cen is not None :
        sim['pos']-=cen
        
    # sim["r"] = ((sim["pos"]-cen)**2).sum(axis=1)**(1,2)

    rho = lambda r : sim["mass"][np.where(sim["r"]<r)].sum()/(4.*math.pi*(r**3)/3)
    target_rho = overden * sim.properties["omegaM0"] * cosmology.rho_crit(sim, z=0) * (1.0+sim.properties["z"])**3

    result = util.bisect(r_min, r_max, lambda r : target_rho-rho(r), epsilon=0, eta=1.e-3*target_rho, verbose=True)
    if cen is not None :
        sim['pos']+=cen

    return result

def potential_minimum(sim, var="pos", r='10 kpc') :
    cen_a = center_of_mass(sim)
    idxsph = filt.Sphere(r,cen_a)
    i = sim[idxsph]["phi"].argmin()
    return sim[idxsph][var][i].copy()

def hybrid_center(sim, r='3 kpc', **kwargs) :
    """

    Determine the center of the halo by finding the shrink-sphere
    -center inside the specified distance of the potential minimum

    """

    try:
        cen_a = potential_minimum(sim)
    except KeyError:
        cen_a = center_of_mass(sim)
    return shrink_sphere_center(sim[filt.Sphere(r, cen_a)], **kwargs)

def index_center(sim, ind = None, var = "pos", **kwargs) :
    """

    Determine the center of mass based on specific particles.

    Supply a list of indices using the ``ind`` keyword.

    """

    if 'ind' is not None :
        return center_of_mass(sim[ind], var=var)
    else :  
        raise RuntimeError("Need to supply indices for centering")
    

def vel_center(sim, cen_size = "1 kpc", cen = None, retcen=False, **kwargs) :
    """

    Use stars from a spehre to calculate center of velocity. The size
    of the sphere is given by the ``cen_size`` keyword and defaults to
    1 kpc.


    """

    if cen==None: cen=array.SimArray([0.,0.,0.],sim['vel'].units)
    if config['verbose'] :
        print "Finding halo velocity center..."
    simcen = sim.star[filt.Sphere(cen_size,cen)]
    if len(simcen)<5 :
        # fall-back to DM
        simcen = sim.dm[filt.Sphere(cen_size,cen)]
    if len(simcen)<5 :
        # fall-back to gas
        simcen = sim.gas[filt.Sphere(cen_size,cen)]
    if len(simcen)<5 :
        # fall-back to all partilces
        simcen=sim[filt.Sphere(cen_size,cen)]
    if len(simcen)<5 :
        # very weird snapshot, or mis-centering!
        print "Insufficient particles ({}) around center to get velocity".format(len(simcen))
        simcen=sim
        #raise ValueError, "Insufficient particles ({}) around center to get velocity".format(len(cen))

    vcen = (simcen['vel'].transpose()*simcen['mass']).sum(axis=1)/simcen['mass'].sum()
    vcen.units = simcen['vel'].units
    if config['verbose'] :
        print "vcen=",vcen

    if retcen:  return vcen
    else:  sim.ancestor["vel"]-=vcen

def center(sim, mode=None, retcen=False, vel=None, velmeth='mass', **kwargs) :
    """

    Determine the center of mass of the given particles using the
    specified mode, then recenter the particles (of the entire
    ancestor snapshot) accordingly

    Accepted values for *mode* are

      *pot*: potential minimum

      *com*: center of mass

      *ssc*: shrink sphere center

      *ind*: center on specific particles; supply the list of particles using the ``ind`` keyword.

      *hyb*: for sane halos, returns the same as ssc, but works faster by
             starting iteration near potential minimum

    or a function returning the COM.

    **Other keywords:**

    *retcen*: if True only return the center without centering the
     snapshot (default = False)

    *ind*: only used when *mode=ind* -- specifies the indices of
     particles to be used for centering

    *vel*: if True, translate velocities so that the velocity of the
    central 1kpc is zeroed
    """

    global config
    if mode is None:
        mode=config['centering-scheme']

    if retcen: velDEF=False
    else     : velDEF=True
    if vel==None: vel=velDEF

    # Identify function
    try:
        fn = {'pot': potential_minimum,
              'com': center_of_mass,
              'mass': center_of_mass,
              'ssc': shrink_sphere_center,
              'hyb': hybrid_center,
              'ind': index_center}[mode]
    except KeyError :
        fn = mode

    # Center of positions
    cen_pos = fn(sim, **kwargs)

    # Center of velocities
    if vel:
        if velmeth=='sphere': cen_vel=vel_center(sim,cen_size="1 kpc",cen=cen_pos,retcen=True)
        elif velmeth=='mass': cen_vel=center_of_mass_velocity(sim)
        
    # Return or shift positions/velocities
    if retcen:
        if vel: return cen_pos,cen_vel
        else  : return cen_pos
    else:
        sim.ancestor["pos"]-=cen_pos
        if vel: sim.ancestor["vel"]-=cen_vel
