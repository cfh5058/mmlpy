#!/usr/bin/python
####################################################################################################################################
#
# MEAGAN LANG'S MATHEMATICAL CLASSES & METHODS
#
####################################################################################################################################
import numpy as np
import mmlpars,mmlio
import copy,math

RAN2_IDUM2=123456789
RAN2_IV=32*[0]
RAN2_IY=0
HELLO=123456789


def logN(x,base=10.):
    """
    Compute logarithm of arbitrary base
    """
    return np.log(x)/np.log(base)

def smooth1d(arr,window=5,method='mean',wrap=False):
    """
    Smooth array
    """
    if wrap and (len(arr)%window)!=0:
        arr=np.hstack(arr,arr[:len(arr)%window])
    new=np.array([np.mean(arr[i*window:(i+1)*window]) for i in len(arr)/window])
    return new

####################################################################################################################################
# INTERSECTION OF SORTED ARRAYS
def intersect1d(A,B,sort=True):
    """
    Returns intersection of two arrays
    """
    if sort:
        A=np.sort(A)
        B=np.sort(B)
    x=[] ; i=0 ; j=0
    while i<len(A) and j<len(B):
        if   A[i]>B[j]: j+=1
        elif B[j]>A[i]: i+=1
        else:
            x.append(A[i])
            i+=1
            j+=1
    return np.array(x)

####################################################################################################################################
# NO REPEATS
def index_norepeats(arr,cut=0.,axis=0):
    """
    Return indices of unrepeated elements in first dimension
    """
    # Check input
    shp=arr.shape
    dim=len(shp)
    num=arr.shape[axis]
    if dim>2: raise Exception('Unsure how this behaves for more than 2 dimensions')
    # Get indices for sorting along nth dimension
    idxa=[slice(0,arr.shape[d]) for d in range(dim)]
    idxT=range(dim)
    idxT.remove(axis)
    idxT.append(axis)
    a=np.lexsort(np.transpose(arr,idxT))
    idxa[axis]=a
    # Get square of difference along sorted axis
    diff=np.diff(arr[a],axis=axis)**2
    # Sum differences along other dimensions
    for d in range(dim)[::-1]:
        if d==axis: continue
        diff=diff.sum(axis=d)
    diff=np.sqrt(diff)
    if len(diff)!=(num-1):
        raise Exception('arr.shape={} but diff.shape={}'.format(shp,diff.shape))
    diff=np.insert(diff,0,1.+cut)
    # Return sorted indices that are not repeats
    idx=np.sort(a[diff>cut])
    return idx

####################################################################################################################################
# RUNNING AVERAGES
def runavg1d(arr,window=5,weights=None):
    """
    Compute average of 1d histogram
    """
    if weights is None: weights=np.ones(arr.shape)
    nx=len(arr)
    avg=np.zeros(arr.shape)
    for x in range(nx):
        xslc=slice(max(x-window,0),min(x+window,nx))
        avg[x]=np.sum(arr[xslc]*weights[xslc])/np.sum(weights[xslc])
    return avg
    
def runavg2d(arr,window=5,weights=None):
    """
    Compute average of 2d histogram
    """
    if weights is None: weights=np.ones(arr.shape)
    nx,ny=arr.shape
    avg=np.zeros((nx,ny))
    for x in range(nx):
        xslc=slice(max(x-window,0),min(x+window,nx))
        for y in range(ny):
            yslc=slice(max(y-window,0),min(y+window,ny))
            avg[x,y]=np.sum(arr[xslc,yslc]*weights[xslc,yslc])/np.sum(weights[xslc,yslc])
    return avg

def wrapder(x,y,wrpwin,derwin=1,avgwin=5,retder2=False):
    """
    Calculate derivative allowing wrapping
    """
    xder1=(x[derwin:]+x[:-derwin])/2.
    yder1=np.diff(y,n=derwin)
    idxwrp=(abs(yder1)>wrpwin/2.)
    yder1[idxwrp]-=np.sign(yder1[idxwrp])*wrpwin
    yder1=runavg1d(yder1,window=avgwin)/np.diff(x,n=derwin)
    if retder2:
        xder2=(xder1[derwin:]+xder1[:-derwin])/2.
        yder2=np.diff(yder1,n=derwin)/np.diff(xder1,n=derwin)
        return xder1,yder1,xder2,yder2
    else:
        return xder1,yder1

def findruns(arr0,err,val='abs',wrap=None,nMin=2,noOverlap=False,idxvalid=None,
             retval=False):
    """
    Find longest run of values within error
    """
    # Select valid arrays
    n0=len(arr0)
    if isinstance(idxvalid,np.ndarray):
        flg=np.arange(n0)[idxvalid]
        arr=arr0[idxvalid]
    else:
        flg=None
        arr=arr0
    n=len(arr)
    # Pars input
    if val is None: val='abs'
    if   nMin<2 : raise Exception('Why look for runs shorter than 2? (nMin={})'.format(nMin))
    elif nMin>n0: return []
    #elif nMin>n0: raise Exception('Cannot have a run longer than the array. (nMin={},N={})'.format(nMin,n0))
    elif nMin>n : return []
    # Loop over indices looking for runs
    pairs=[]
    idx=0
    while idx<n:
        ibeg=idx ; iend=idx
        # Get values to find run in
        if   val == 'abs': ifit=arr-arr[idx]
        elif val == 'rel': ifit=np.hstack([np.diff(arr),err])
        else             : ifit=arr-val
        if wrap:
            idxwrp=(abs(ifit)>wrap/2.)
            ifit[idxwrp]-=np.sign(ifit[idxwrp])*wrap
        # Continue if started on a bad value
        if abs(ifit[idx])>err:
            idx+=1
            continue
        # Retreat beginning until error exceeded or no more array
        while True:
            ibeg-=1
            if ibeg<0: break
            if val=='rel': ibegfit=arr[min(ibeg+1,n-1)]-arr[ibeg]
            else         : ibegfit=ifit[ibeg]
            if abs(ibegfit)>err:
                ibeg+=1
                break
        if noOverlap and len(pairs)>0:
            if ibeg<=pairs[-1][1]: 
                ibeg+=1
        ibeg=max(ibeg,0)
        # Advance end until error exceeded or no more array
        while True:
            iend+=1
            if iend>=n: break
            if val=='rel': iendfit=arr[iend]-arr[max(iend-1,0)]
            else         : iendfit=ifit[iend]
            if abs(iendfit)>err:
                iend-=1
                break
        iend=min(iend,n-1)
        # Add pair to list if longer than minimum length
        if (iend-ibeg)>=(nMin-1): 
            pairs.append([ibeg,iend])
            if noOverlap: idx=iend
        # Advance index
        idx+=1
    # Return correct indices
    if isinstance(idxvalid,np.ndarray): pairs=[[flg[p[0]],flg[p[1]]] for p in pairs]
    # Include array values rather than indices
    if retval: pairs=[[arr0[p[0]]-err/2.,arr0[p[1]]+err/2.] for p in pairs]
    # Return pairs
    return pairs

def findrun(arr,err,**kws):
    """Find longest run of a single value/difference"""
    runs=findruns(arr,err,**kws)
    if len(runs)==0: beg,end=[0,0]
    else:
        imax=np.argmax(np.array([p[1]-p[0] for p in runs]))
        beg,end=runs[imax]
    avg=np.mean(arr[beg:(end+1)])
    return avg,beg,end

####################################################################################################################################
# POLAR SQUARE
def polar_square(phi0,side=1.0):
    phi=phi0 % (2.*np.pi)
    if phi<0: phi+=2.*np.pi
    if phi<=(np.pi/4.) or phi>=(7.*np.pi/4.):
        x=side/2.
        y=x*np.tan(phi)
    elif phi>(np.pi/4.) and phi<=(3.*np.pi/4.):
        y=side/2.
        x=y/np.tan(phi)
    elif phi>(3.*np.pi/4.) and phi<=(5.*np.pi/4.):
        x=-side/2.
        y=x*np.tan(phi)
    elif phi>(5.*np.pi/4.) and phi<=(7.*np.pi/4.):
        y=-side/2.
        x=y/np.tan(phi)
    return (x,y)

####################################################################################################################################
# METHOD TO SIMULATE BUILDGAL'S RAN2 FUNCTION
def ran2(idum):
    """
    Re-creates BUILDGAL's ran2 function
    """
    IM1=2147483563 ; IM2=2147483399 ; AM=1./IM1 ; IMM1=IM1-1
    IA1=40014 ; IA2=40692 ; IQ1=53668 ; IQ2=52774 ; IR1=12211 ; IR2=3791
    NTAB=32 ; NDIV=1+IMM1/NTAB ; EPS=1.2e-7 ; RNMX=1.-EPS

    global RAN2_IDUM2,RAN2_IV,RAN2_IY
    idum2=RAN2_IDUM2
    iv=RAN2_IV
    iy=RAN2_IY
    # defdict=dict(IDUM2=123456789,RAN2_IV=NTAB*[0],RAN2_IY=0)
    # for ikey in defdict.keys():
    #     if ikey not in globals():
    #         exec('global {}'.format(ikey))
    #         exec('{}={}'.format(ikey,str(defdict[ikey])))
#    global IDUM2,RAN2_IV,RAN2_IY


    if idum <= 0:
        idum=max(-idum,1)
        idum2=idum
        for j in reversed(range(1,NTAB+8+1)):
            k=idum/IQ1
            idum=IA1*(idum-k*IQ1)-k*IR1
            if (idum<0): idum=idum+IM1
            if (j<=NTAB): iv[j-1]=idum
        iy=iv[0]
    k=idum/IQ1
    idum=IA1*(idum-k*IQ1)-k*IR1
    if (idum<0): idum=idum+IM1
    k=idum2/IQ2
    idum2=IA2*(idum2-k*IQ2)-k*IR2
    if (idum2<0): idum2=idum2+IM2
    j=1+iy/NDIV
    iy=iv[j-1]-idum2
    iv[j-1]=idum
    if (iy<1): iy=iy+IMM1
    ran2=float(min(AM*iy,RNMX))

    RAN2_IDUM2=idum2
    RAN2_IV=iv
    RAN2_IY=iy

    return ran2

####################################################################################################################################
# METHOD TO HACK BOUNDS ON LEASTSQ
def leastsq(errfunc0,x0,bounds=None,**exkw):
    """
    Hack to allow leastsq to take bounds
    """
    from scipy import optimize
    lres=1.0e6
    def errfunc(x,*args):
        if bounds!=None:
            for idx in range(len(x)):
                if bounds[idx][0]!=None and x[idx]<bounds[idx][0]: return lres
                if bounds[idx][1]!=None and x[idx]>bounds[idx][1]: return lres
        return errfunc0(x,*args)
    return optimize.leastsq(errfunc,x0,**exkw)

####################################################################################################################################
# METHOD TO CONVERT CARTESIAN COORDINATES TO SPHERICAL COORDINATES
def xyz2sph(xyz):
    """
    Converts positions from cartesian to spherical coordinates
    """
    r=xyz2rxyz(xyz)
    phi=xyz2phi(xyz)
    theta=xyz2theta(xyz)
    np.vstack((Lx,Ly,Lz)).T
    sph=np.vstack((r,phi,theta)).T
    return sph
def sph2xyz(sph):
    """
    Converts positions from spherical to cartesian coordinates
    """
    x=sph[:,0]*np.cos(sph[:,1])*np.sin(sph[:,2])
    y=sph[:,0]*np.sin(sph[:,1])*np.sin(sph[:,2])
    z=sph[:,0]*np.cos(sph[:,2])
    xyz=np.vstack((x,y,z)).T
    return xyz
def sph2cyl(sph):
    """
    Converts positions from spherical to cylindrical coordinates
    """
    xyz=sph2xyz(sph)
    cyl=xyz2cyl(xyz)
    return cyl
####################################################################################################################################
# METHOD TO CONVERT CARTESIAN COORDINATES TO CYLINDRICAL COORDINATES
def xyz2cyl(xyz):
    """
    Converts positions from cartesian to cylindrical coordinates
    """
    R=xyz2rxy(xyz)
    phi=xyz2phi(xyz)
    z=xyz[:,2]
    cyl=np.vstack((R,phi,z)).T
    return cyl
def cyl2xyz(cyl):
    """
    Converts positions from cylindrical to cartesian coordinates
    """
    x=cyl[:,0]*np.cos(cyl[:,1])
    y=cyl[:,0]*np.sin(cyl[:,1])
    z=cyl[:,2]
    xyz=np.vstack((x,y,z)).T
    return xyz
def cyl2sph(cyl):
    """
    Converts positions from cylindrical to spherical coordinates
    """
    xyz=cyl2xyz(cyl)
    sph=xyz2sph(xyz)
    return sph

####################################################################################################################################
# METHOD TO CONVERT INDIVIDUAL COORDINATES
def xyz2rxy  (xyz): return np.sqrt(xyz[:,0]**2 + xyz[:,1]**2)
def xyz2rxyz (xyz): return np.sqrt(xyz[:,0]**2 + xyz[:,1]**2 + xyz[:,2]**2)
def xyz2phi  (xyz): return np.arctan(xyz[:,1]/xyz[:,0])
def xyz2theta(xyz): return np.arccos(xyz[:,2]/xyz2rxyz(xyz))
####################################################################################################################################
# METHODS TO GET TIME DERIVATIVE OF COORDINATES
def drxyz_dt (pos,vel): return (pos[:,0]*vel[:,0]+pos[:,1]*vel[:,1]+pos[:,2]*vel[:,2])/xyz2rxyz(pos)
def drxy_dt  (pos,vel): return (pos[:,0]*vel[:,0]+pos[:,1]*vel[:,1])/xyz2rxy(pos)
def dphi_dt  (pos,vel): return (pos[:,0]*vel[:,1]-pos[:,1]*vel[:,0])/(xyz2rxy(pos)**2.)
def dtheta_dt(pos,vel): 
    r=xyz2rxyz(pos)
    rdot=drxyz_dt(pos,vel)
    return (pos[:,2]*rdot-vel[:,2]*r)/(r*np.sqrt(r**2. - pos[:,2]**2.))
####################################################################################################################################
# METHOD TO GET VRAD FROM XYZ POSITIONS AND VELOCITIES
def vrad(pos,vel,form='deriv',coord='sph'):
    """
    Converts xyz positions and velocities to v in the radial direction
    """
    if   form=='angle': 
        if   coord=='cyl': vrad=xyz2rxy(vel)*np.cos(xyz2phi(vel)-xyz2phi(pos))
        elif coord=='sph': raise Exception('Work in progress')
    elif form=='deriv':
        if   coord=='cyl': vrad=drxy_dt (pos,vel)
        elif coord=='sph': vrad=drxyz_dt(pos,vel)
        else: raise Exception('Invalid coordinate system: {}'.format(coord))
    else: raise Exception('Invalid form: {}'.format(form))
    return vrad
####################################################################################################################################
# METHOD TO GET VPHI FROM XYZ POSITIONS AND VELOCITIES
def vphi(pos,vel,form='deriv',coord='sph'):
    """
    Converts xyz positions and velocities to v in phi direction
    """
    if   form=='angle': vphi=xyz2rxy(vel)*np.sin(xyz2phi(vel)-xyz2phi(pos))
    elif form=='deriv': 
        if   coord=='cyl': vphi=xyz2rxy (pos)*dphi_dt(pos,vel)
        elif coord=='sph': vphi=xyz2rxyz(pos)*dphi_dt(pos,vel)
        else: raise Exception('Invalid coordinate system: {}'.format(coord))
    else: raise Exception('Invalid form: {}'.format(form))
    return vphi
####################################################################################################################################
# METHOD TO GET VTHETA FROM XYZ POSITIONS AND VELOCITIES
def vtheta(pos,vel,**exkw):
    """
    Converts xyz positions and velocities to v in the theta direction
    """
    vtheta=xyz2rxyz(pos)*np.sin(xyz2theta(pos))*dtheta_dt(pos,vel)
    return vtheta
    
####################################################################################################################################
# METHOD TO GET OORT'S CONSTANTS
def oorts_const(R=None,v=None,method=None):
    """
    Calculates Oort's constants
        R: cylindrical radius
        v: rotational velocity
    """
    # Pars input
    npart=len(R) ; sshp=(npart,) ; vshp=(npart,3)
    v=mmlpars.mml_pars(v,type=np.ndarray,shape=sshp)
    R=mmlpars.mml_pars(R,type=np.ndarray,shape=sshp)
    method=mmlpars.mml_pars(method,default='both',list=['A','B','both'])
    # Get derivative
    dv_dR=derivative(v,R)
    # Get constants
    A=0.5*(v/R - dv_dR)
    B=-0.5*(v/R + dv_dR)
    # Return output
    if   method=='both': return A,B
    elif method=='A'   : return A
    elif method=='B'   : return B
    else: raise Exception('Invalid method: {}'.format(method))

####################################################################################################################################
# METHOD TO DETERMINE EPICYCLIC FREQUENCY
def epicycle(r,v,approx=None):
    """
    Returns epicyclic frequency as function of radius
    """
    # Pars input
    npart=len(r) ; sshp=(npart,) ; vshp=(npart,3) 
    r=mmlpars.mml_pars(r,type=np.ndarray,shape=sshp)
    v=mmlpars.mml_pars(v,type=np.ndarray,shape=sshp)
    approx=mmlpars.mml_pars(approx,default=True,type=bool)
    # Epicyclic approximation (nearly circular orbits)
    if approx:
        # Get Oort's constant
        B=oorts_const(r,v,method='B')
        # Mask values
        r=np.ma.masked_less_equal(r,0.)
        # Calculate kappa
        kappa=np.sqrt(-4.*B*v/r)
    # Direct calculation
    else:
        w=v/r
        # Get derivative
        dwr2_dr=derivative(r*r*w,r)
        # Mask values
        r=np.ma.masked_less_equal(r,0.)
        # Calculated epicyclic frequency
        kappa=np.sqrt(2.*w*dwr2_dr/r)
    # Return output
    return kappa.filled(0.)

####################################################################################################################################
# METHOD TO ESTIMATE DERIVATIVE
def derivative(f,x):
    """
    Uses finite differences to estimate derivative
    """
    # Pars input
    npart=len(x) ; sshp=(npart,) ; vshp=(npart,3) 
    f=mmlpars.mml_pars(f,type=np.ndarray,shape=sshp)
    x=mmlpars.mml_pars(x,type=np.ndarray,shape=sshp)
    # Get indices
    idxsort=np.argsort(x)
    revsort=np.argsort(idxsort)
    xsrt=x[idxsort]
    fsrt=f[idxsort]
    if np.array_equal(x,xsrt): print 'x variable sorted'
    # Get derivative
    df0=np.ma.masked_equal(np.diff(fsrt),0)
    dx0=np.ma.masked_equal(np.diff(xsrt),0)
    print 'Nf masked = {}'.format(np.sum(df0.mask))
    print 'Nx masked = {}'.format(np.sum(dx0.mask))
    df_dx0=(df0/dx0).filled(0.)
    # Interpolate
    der_beg=[df_dx0[0]]
    der_mid=list((df_dx0[1:]+df_dx0[:-1])/2.)
    der_end=[df_dx0[-1]]
    df_dx=np.array(der_beg+der_mid+der_end)[revsort]
    # Return output
    return df_dx

####################################################################################################################################
# METHOD TO RETURN CENTER OF MASS
def center_mass(m,pos,vel=None,vflag=None):
    """
    Returns center of mass
    """
    # Pars
    npart=len(m) ; sshp=(npart, ) ; vshp=(npart,3)
    m=mmlpars.mml_pars(m,type=np.ndarray,shape=sshp)
    pos=mmlpars.mml_pars(pos,type=np.ndarray,shape=vshp)
    vflag=mmlpars.mml_pars(vflag,default=False,type=bool)
    if vflag: vel=mmlpars.mml_pars(vel,type=np.ndarray,shape=vshp)
    # Compute center
    mtot=m.sum()
    mx=m*pos[:,0]
    my=m*pos[:,1]
    mz=m*pos[:,2]
    com=np.array([mx.sum(),my.sum(),mz.sum()])/mtot
    # Compute velocity center
    if vflag:
        mvx=m*vel[:,0]
        mvy=m*vel[:,1]
        mvz=m*vel[:,1]
        vcom=np.array([mvx.sum(),mvy.sum(),mvz.sum()])/mtot
    # Return output
    if vflag: return com,vcom
    else    : return com

####################################################################################################################################
# METHOD TO RETURN CENTER OF DENSITY
def center_density(m,pos,vel=None,vflag=None,navg=None,nbin=None,window=None):
    """
    Return center based on density
    """
    # Pars
    npart=len(m) ; sshp=(npart, ) ; vshp=(npart,3)
    m=mmlpars.mml_pars(m,type=np.ndarray,shape=sshp)
    pos=mmlpars.mml_pars(pos,type=np.ndarray,shape=vshp)
    vflag=mmlpars.mml_pars(vflag,default=False,type=bool)
    if vflag: vel=mmlpars.mml_pars(vel,type=np.ndarray,shape=vshp)
    navg=mmlpars.mml_pars(navg,default=0,min=0,type=int)
    nbin=mmlpars.mml_pars(nbin,default=100,min=2,type=int)
    windowDEF=(pos.max()-pos.min())/2.
    window=mmlpars.mml_pars(window,default=windowDEF,type=float)
    # Get center of mass to search on
    com0=center_mass(m,pos)
    # Create histogram of density
    bins=[]
    for icom0 in list(com0): bins.append(np.linspace(icom0-window,icom0+window,nbin))
    hist,bins=myhist(pos,bins=bins,weights=m)
    # Compute maximum of density
    idx0_flat=np.argmax(hist)
    idx0_mult=np.unravel_index(idx0_flat,hist.shape)
    bin0=[]
    bincen=[]
    for ibin,iidx in zip(bins,idx0_mult):
        bin0.append(ibin[iidx:iidx+2])
        bincen.append(np.mean(bin0[-1]))
    bincen=np.array(bincen)
    # Return the location of maximum density if navg=0
    if navg==0:
        com=bincen
        if vflag:
            vhist,vbins=myhist(pos,bins=bin0,weights=m*vel)
            vcom=vhist[0,0,0,:]/m.sum()
            out=com,vcom
        else: out=com
    # Get particles close to maximum density
    else:
        pos0=pos-bincen
        r0=xyz2rxyz(pos0)
        idx0=np.argsort(r0)[:navg]
        # Get center
        mb=m[idx0] ; posb=pos[idx0]
        if vflag: velb=vel[idx0]
        else    : velb=None
        out=center_mass(mb,posb,vel=velb,vflag=vflag)
    # Print
    mmlio.verbose('navg={}'.format(navg))
    mmlio.verbose('window={}'.format(window))
    mmlio.verbose('com (mass) = {}'.format(com0))
    if vflag: com,vcom=out
    else    : com=out
    mmlio.verbose('com (dens) = {}'.format(com))
    mmlio.yorn('Look OK?')
    # Return output
    return out
   
####################################################################################################################################
# METHOD TO RETURN CENTER OF POTENTIAL
def center_potential(pot,m,pos,vel=None,vflag=None,navg=None):
    """
    Returns center based on potential
    """
    # Pars
    npart=len(pot) ; sshp=(npart, ) ; vshp=(npart,3)
    if npart==0: raise Exception('No particles provided.')
    pot=mmlpars.mml_pars(pot,type=np.ndarray,shape=sshp)
    m=mmlpars.mml_pars(m,type=np.ndarray,shape=sshp)
    pos=mmlpars.mml_pars(pos,type=np.ndarray,shape=vshp)
    vflag=mmlpars.mml_pars(vflag,default=False,type=bool)
    if vflag: vel=mmlpars.mml_pars(vel,type=np.ndarray,shape=vshp)
    navg=mmlpars.mml_pars(navg,default=0,min=0,type=int)
    if navg==0: navg=1
    # Compute minimum of potential
    idxmax=np.argmin(pot)
    com0=pos[idxmax,:]
    # Find navg particles at center of potential
    pos0=pos-com0
    r0=xyz2rxyz(pos0)
    idx0=np.argsort(r0)[:navg]
    # Get center
    mb=m[idx0] ; posb=pos[idx0]
    if vflag: velb=vel[idx0]
    else    : velb=None
    out=center_mass(mb,posb,vel=velb,vflag=vflag)
    # Return output
    return out

####################################################################################################################################
# METHOD TO COMPUTE HISTOGRAM WITH MULTIDIMENSIONAL WEIGHTS
def myhist(x,weights=None,**histkw):
    """
    Multidimensional histogram with option for multidimensional weights
    """
    # Pars input
    x=mmlpars.mml_pars(x,type=np.ndarray)
    xshp=x.shape ; N=xshp[0]
    # Histogram w/o weights
    if weights==None: hist,bins=np.histogramdd(x,**histkw)
    # Histogram w/ weights
    else:
        wshp=weights.shape
        # Histogram w/ scaler weights
        if wshp==(N,): hist,bins=np.histogramdd(x,weights=weights,**histkw)
        # Histogram w/ vector weights
        else:
            # Handle errors
            if wshp[0]!=N: raise Exception('Weights must have same size first dimension as data. (data.shape={},weights.shape={})'.format(xshp,wshp))
            if len(wshp)>2: raise Exception('Weights with more than 2 dimensions not supported. (weights.shape={})'.format(wshp))
            # Get histograms
            owshp=tuple(list(histlist[0].shape)+[wshp[1]])
            hist=np.zeros(owshp)
            for iw in range(wshp[1]):
                ihist,bins=np.histogramdd(x,weights=weights[:,iw],**histkw)
                if iw==0: histkw['bins']=bins
                if   len(xshp)==1: hist[:,iw]=ihist
                elif len(xshp)==2: hist[:,:,iw]=ihist
                elif len(xshp)==3: hist[:,:,:,iw]=ihist
                elif len(xshp)==4: hist[:,:,:,:,iw]=ihist
                else: raise Exception('Multidimensional weights only supported for data with <4 dimensions. (data.shape={})'.format(xshp))
    # Return output
    return hist,bins

####################################################################################################################################
# METHOD TO ENSURE LIMITS ARE LOG FRIENDLY
def loglim(inlim,minord=3):
    '''
    NAME:
        mmlmath.loglim
    PURPOSE:
        To ensure limits are log friendly (ignoring sign)
    CALLING
        outlim=loglim(inlim[,minord=3])
    ARGUMENTS:
        inlim:  Tuple with two elements limits to make log friendly
    KEYWORDS:
        minord: Int # of orders of magnitude to scale min to from max if a zero present
    OUTPUT:
        outlim: Tuple with two elements log friendly limits
    '''
    inlim=mmlpars.mml_pars(inlim,type=tuple,nelements=2)
    minord=mmlpars.mml_pars(minord,default=3,type=int,min=1)
    if inlim[0]==0.:
        if inlim[1]==0.:
            raise Exception('Both limits can"t be zero.')
        else:
            lim1=inlim[1]/(10.**minord)
    else:
        lim1=inlim[0]
    if inlim[1]==0.:
        lim2=inlim[0]/(10.**minord)
    else:
        lim2=inlim[1]
    return (lim1,lim2)

####################################################################################################################################
# METHOD TO CREATE BROKEN LOG SPACE 
####################################################################################################################################
def logspace(imin,imax,n,minord=None,addzero=False,mirror=False):
    '''
    '''
    # Mirror limits
    if mirror:
        iabsmax=max(abs(imin),abs(imax))
        imin=-iabsmax
        imax=+iabsmax
    # Remove zeros
    (logmin,logmax)=loglim((imin,imax),minord)
    # Find negative/positive limits
    if logmin<0 and logmax<0:
        n1=n ; n2=0
        lim1=(logmin,logmax)
        lim2=None
        midpnt=None
        if addzero:
            n1=n-1 ; midpnt=0.
    elif logmin<0 and logmax>0:
        n1=n/2 ; n2=n/2
        lim1=loglim((logmin,0.),minord)
        lim2=loglim((0.,logmax),minord)
        if n%2 == 0: midpnt=None
        else:        midpnt=0.
    elif logmin>0 and logmax<0:
        n1=n/2 ; n2=n/2
        lim1=loglim((logmin,0.),minord)
        lim2=loglim((0.,logmax),minord)
        if n%2 == 0: midpnt=None
        else:        midpnt=0.
    elif logmin>0 and logmax>0:
        n1=0 ; n2=n
        lim1=None
        lim2=(logmin,logmax)
        midpnt=None
        if addzero:
            n2=n-1 ; midpnt=0.
    else:
        raise Exception('Not sure what\'s happening... ({},{})'.format(logmin,logmax))
    # Create vectors
    if lim1 != None:
        space1=np.sign(lim1[0])*np.logspace(np.log10(abs(lim1[0])),np.log10(abs(lim1[1])),n1)
    else:
        space1=np.array([])
    if lim2 != None:
        space2=np.sign(lim2[0])*np.logspace(np.log10(abs(lim2[0])),np.log10(abs(lim2[1])),n2)
    else:
        space2=np.array([])
    if midpnt !=None:
        space_mid=np.array([midpnt])
    else:
        space_mid=np.array([])
    # Concatenate vectors
    space=np.concatenate((space1,space_mid,space2),axis=0)
    return space


####################################################################################################################################
# METHOD TO ORDER OF MAGNITUDE ROUND
####################################################################################################################################
def oom(xin,nsig=None,method=None):
    '''
    NAME:
        mmlmat.oom
    PURPOSE:
        To round a number to the nearest order of magnitude.
    CALLING:
        xout=oom(xin[,nsig=,method=])
    ARGUMENTS:
        xin:    Float to round to nearest order of magnitude.
    KEYWORDS:
        nsig:   Number of significant figures to keep [DEFAULT=1]
        method: String identifying method of rounding
          "ROUND": [DEFAULT] Round to closest order of magnitude
          "FLOOR": Round down to closest order of magnitude
          "CEIL":  Round up to closest order of magnitude
    OUTPUT:
        xout:   xin rounded to the nearest order of magnitude
    '''
    # Pars input
    methLIST=['ROUND','FLOOR','CEIL']
    xin=mmlpars.mml_pars(xin,type=float)
    nsig=mmlpars.mml_pars(nsig,type=int,min=0,default=1)
    if isinstance(method,str): method=method.upper()
    method=mmlpars.mml_pars(method,list=methLIST,default=methLIST[0],type=str)
    # Find zeros
    if xin != 0.0:
        # Get power
        pow=math.floor(math.log10(abs(xin)))-float(nsig-1)
        # Round
        if method == 'ROUND':   fact=round(xin/(10.0**pow))
        elif method == 'FLOOR': fact=math.floor(xin/(10.0**pow))
        elif method == 'CEIL':  fact=math.ceil(xin/(10.0**pow))
        else: raise Exception('{} is invalid for the keyword method.'.format(method))
        # Handle round to zero
        if fact == 0.0: fact=0.1
        else: pass
        # Return output
        xout=fact*(10.0**pow)
    else:
        xout=xin
    # Return answer
    return xout

####################################################################################################################################
# METHOD TO TAKE AMPLITUDE OF VECTORS
def absvect(x):
    """
    Returns amplitude of vectors
    """
    npart=len(x) ; vshp=(npart,3)
    x=mmlpars.mml_pars(x,type=np.ndarray,shape=vshp)
    absx=np.sqrt(np.sum(x*x,axis=1))
    return absx

####################################################################################################################################
# METHOD TO CALCULATE ANGULAR MOMENTUM
def angmom(m,r,v):
    """
    Returns angular momentum vector
    """
    npart=len(m) ; sshp=(npart,) ; vshp=(npart,3)
    if npart==0: raise Exception('No particles provided.')
    m=mmlpars.mml_pars(m,type=np.ndarray,shape=sshp)
    r=mmlpars.mml_pars(r,type=np.ndarray,shape=vshp)
    v=mmlpars.mml_pars(v,type=np.ndarray,shape=vshp)
    rxv=np.cross(r,v)
    Lx=m*rxv[:,0]
    Ly=m*rxv[:,1]
    Lz=m*rxv[:,2]
    angmom=np.vstack((Lx,Ly,Lz)).T
    return angmom

####################################################################################################################################
# METHOD TO RETURN SPIN
def spin(*args):
    """
    Returns spin
    """
    if   len(args)==3: h,v,r=args
    elif len(args)==4: 
        j,m,v,r=args
        h=j/m
    else: raise Exception('Only 3 or 4 arguments accepted.')
    spin=h/(np.sqrt(2.)*v*r)
    return spin

####################################################################################################################################
# METHOD TO RETURN SPIN FOR PARTICLE POSITIONS
def spinvec(j,m,r,G):
    """
    Returns spin for particle positions
    """
    # Pars input
    npart=len(m) ; sshp=(npart,) ; vshp=(npart,3)
    j=mmlpars.mml_pars(j,type=np.ndarray,shape=vshp)
    m=mmlpars.mml_pars(m,type=np.ndarray,shape=sshp)
    r=mmlpars.mml_pars(r,type=np.ndarray,shape=sshp)
    G=mmlpars.mml_pars(G,type=float,min=0.)
    # Get indices
    idxsort=np.argsort(r)
    revsort=np.argsort(idxsort)
    # Get mass enclosed and circular velocity
    menc=np.cumsum(m[idxsort])[revsort]
    venc=np.sqrt(G*menc/r)
    # Get angular momentum enclosed
    jencxyz=np.cumsum(j[idxsort,:],axis=0)[revsort,:]
    jenc=absvect(jencxyz)
    # Get spin
    senc=spin(jenc,menc,venc,r)
    # Return output
    return senc

####################################################################################################################################
# METHOD TO COMPUTE SURFACE DENSITY
def surfdens(m,r,nwin=1000):
    """
    Estimates surface mass density at each particle's location
    """
    # Pars input
    npart=len(m) ; sshp=(npart,) ; vshp=(npart,3)
    m=mmlpars.mml_pars(m,type=np.ndarray,shape=sshp)
    r=mmlpars.mml_pars(r,type=np.ndarray,shape=sshp)
    # Get radius & sorting
    idxsort=np.argsort(r)
    revsort=np.argsort(idxsort)
    rsrt=r[idxsort]
    msrt=m[idxsort]
    sigma=np.zeros(r.shape)
    # Loop over windows
    imin=0 ; imax=0 ; irmin=0. ; irmax=0.
    while imin<npart:
        imax=min(npart,imin+nwin)
        iidx=range(imin,imax)
        if rsrt[iidx].min()==rsrt[iidx].max():
            irmin=irmax
        else:
            irmin=rsrt[iidx].min()
        irmax=rsrt[iidx].max()
        iarea=np.pi*(irmax**2 - irmin**2)
        imass=msrt[iidx].sum()
        sigma[iidx]=imass/iarea
        imin=imax
    # Return output
    return sigma[revsort]

####################################################################################################################################
# METHOD TO COMPUTE PARTICLE AVERAGED DISPERSION 
def dispersion(x,y,nwin=1000):
    """
    Estimates dispersion of y at each x
    """
    # Pars input
    npart=len(x) ; sshp=(npart,) ; vshp=(npart,3)
    x=mmlpars.mml_pars(x,type=np.ndarray,shape=sshp)
    y=mmlpars.mml_pars(y,type=np.ndarray,shape=sshp)
    # Get radius & sorting
    idxsort=np.argsort(x)
    revsort=np.argsort(idxsort)
    xsrt=x[idxsort]
    ysrt=y[idxsort]
    sigma=np.zeros(x.shape)
    # Loop over windows
    imin=0 ; imax=0
    while imin<npart:
        imax=min(npart,imin+nwin)
        iidx=range(imin,imax)
        sigma[iidx]=np.sqrt(np.mean((ysrt[iidx]-np.mean(ysrt[iidx]))**2))
        imin=imax
    # Return output
    return sigma[revsort]

####################################################################################################################################
# METHOD TO COMPUTE TOOMRE Q PARAMETER
def toomreq(vdisp,kappa,sdens,G):
    """
    Returns the toomre Q parameter for particle info
    """
    # Pars input
#    mmlio.yorn('vdisp: {}'.format(type(vdisp)))
#    mmlio.yorn('kappa: {}'.format(type(kappa)))
#    mmlio.yorn('sdens: {}'.format(type(sdens)))
    if isinstance(vdisp,np.ndarray):
        npart=len(vdisp) ; sshp=(npart,) ; vshp=(npart,3)
        vdisp=mmlpars.mml_pars(vdisp,type=np.ndarray,shape=sshp)
        kappa=mmlpars.mml_pars(kappa,type=np.ndarray,shape=sshp)
        sdens=mmlpars.mml_pars(sdens,type=np.ndarray,shape=sshp)
    else:
        vdisp=mmlpars.mml_pars(vdisp,type=float)
        kappa=mmlpars.mml_pars(kappa,type=float)
        sdens=mmlpars.mml_pars(sdens,type=float)
    G=mmlpars.mml_pars(G,type=float,min=0.)
    # Mask values
    sdens=np.ma.masked_equal(sdens,0)
    # Compute Toomre Q
    Q=vdisp*kappa/(3.36*G*sdens)
    # Return output
    return Q.filled(0.)

####################################################################################################################################
# METHOD TO RETURN ENCLOSED VARIABLE FOR PARTICLE POSITIONS
def xenc(r,x):
    """
    Returns variable enclosed for particle positions
    """
    # Pars input
    npart=len(r) ; sshp=(npart,) ; vshp=(npart,3) ; xshp=x.shape
    if len(xshp)==1:
        vecflag=False
        xshp=sshp
    elif len(xshp)==2:
        vecflag=True
        xshp=vshp
    else: raise Exception('Only 2D variables currently supported (x.shape={})'.format(xshp))
    x=mmlpars.mml_pars(x,type=np.ndarray,shape=xshp)
    r=mmlpars.mml_pars(r,type=np.ndarray,shape=sshp)
    # Get indices
    idxsort=np.argsort(r)
    revsort=np.argsort(idxsort)
    # Get variable enclosed
    if vecflag: xenc=np.cumsum(x[idxsort,:])[revsort,:]
    else      : xenc=np.cumsum(x[idxsort  ])[revsort  ]
    # Return output
    return xenc

####################################################################################################################################
# METHOD TO DETERMINE INITIAL CARTESIAN POSITIONS/VELOCITIES FOR AN ORBIT
def orbit(M1,M2,rinit,rperi,ecc,gconst,verbose=False):
    # Pars input
    M1=mmlpars.mml_pars(M1,type=float,min=0.)
    M2=mmlpars.mml_pars(M2,type=float,min=0.)
    rinit=mmlpars.mml_pars(rinit,type=float,min=0.)
    rperi=mmlpars.mml_pars(rperi,type=float,min=0.,max=rinit)
    ecc=mmlpars.mml_pars(ecc,type=float,min=0.)
    gconst=mmlpars.mml_pars(gconst,type=float,min=0.)
    if rinit<=rperi: raise Exception('Pericentric distance is larger than or equal to initial separation.')

    # Set useful values
    M=M1+M2
    mu=(M1*M2)/M
    vfrac=np.sqrt(1.+ecc)
    vperi_c=np.sqrt(gconst*M/rperi)
    
    # Get velocity at pericenter from orbit equation
    vperi=vperi_c*vfrac
  
    # Get tangential initial velocity from conservation of angular momentum
    vinit_phi=rperi*vperi/rinit
  
    # Get total initial velocity from conservation of energy
    vinit=np.sqrt(gconst*M*((ecc/rperi)+(1.0/rinit)))
  
    # Get radial initial velocity from definition of total velocity
    vinit_r=-np.sqrt(vinit*vinit - vinit_phi*vinit_phi)
  
    # Get angle of rotation required to put peri on negative x-axis 
    hinit=rinit*vinit_phi
    phiinit=np.arccos(-(1./ecc)*(hinit*hinit/(gconst*M*rinit) - 1.))
  
    # Set output 
    vorb=np.array([vinit_r,vinit_phi,0.])

    # Print some info
    if verbose:
        if ecc == 0:              print 'CIRCULAR orbit'
        elif ecc > 0 and ecc < 1: print 'ELLIPTICAL orbit'
        elif ecc == 1:            print 'PARABOLIC orbit'
        elif ecc > 1:             print 'HYPERBOLIC orbit'
        print '    ecc   = {}'.format(ecc)
        print '    vfrac = {}'.format(vfrac)
        print '    vcirc = {}'.format(np.sqrt(gconst*M1/rinit))
        print '    vel   = {}'.format(vorb)
        print '    phi   = {}'.format(phiinit)
        print '    M1    = {}'.format(M1)
        print '    rinit = {}'.format(rinit)
        print '    G     = {}'.format(gconst)
  
    return vorb,phiinit
