#!/usr/bin/python
import numpy as np

####################################################################################################################################
def main(npart=1000,nbin=100,nmax=10,rmax=None,disp=4):
    """
    Method to create my universe from command line
    """
    dphi=np.random.randint(-disp,disp,npart).astype('float')
    phi0=np.linspace(0.,2.*np.pi,npart,endpoint=False).astype('float')
    phi=phi0+(phi0[1]-phi0[0])*dphi/float(disp**2)
    circ=float(npart)
    if not isinstance(rmax,float): rmax=circ/(2.*np.pi)
    phibins=np.linspace(0.,2.*np.pi,nbin,endpoint=False).astype('float')
    vel0=

def calc_r(phi,bins,rmax=None):
    """
    Method to calculate 'radius'
    """
    mass,obins=np.histogram(phi,bins=bins)
    if not isinstance(rmax,float): rmax=float(mass.max())
    r=rmax-mass.astype('float')
    return r

def calc_rad(r):
    R=r.max()
    return R

def calc_pos(phi,R):
    pos=R*phi
    return pos

def calc_grav(phi,R):
    grav=np.zeros(phi.shape)
    for iphi,iidx in enumerate(phi):
        phisep=np.ma.masked_equal(phi-iphi,0.)
        idxrapneg=np.where(phisep > +np.pi)
        idxrappos=np.where(phisep < -np.pi)
        phisep[idxrapneg]+=(-2.*np.pi)
        phisep[idxrappos]+=(+2.*np.pi)
        xsep=calc_pos(phisep,R)
        xacc=1./(xsep**2)
        grav[idx]=xacc.sum()
    return grav

def calc_surf(phi,r):
    r_half=r[:-1]+np.diff(r)/2.
    dphi=np.diff(phi)
    surf=r_half*dphi
    return surf.sum()

def calc_area(phi,r):
    r_half=r[:-1]+np.diff(r)/2.
    dphi=np.diff(phi)
    area=(r_half**2)*dphi/2.
    return area.sum()
    
    

if __name__ == '__main__':
    main()
