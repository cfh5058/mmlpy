#!/usr/bin/python
from mmlutils import *
import mmlutils
import numpy as np
import os,pprint,copy,time,timeit

FTEST='/home/langmm/code/python/mine/mmlastro/ellipse_test.pickle'

PARAM0={'I0'   :0.0, # 1.  Mean intensity
        'B3'   :0.0, # 2.  3rd even Fourier coefficient
        'der'  :0.0, # 3.  Derivative of intensity
        'a'    :1.0, # 4.  Semi-major axis
        'rms'  :0.0, # 5.  RMS residual of intensity
        'B4'   :0.0, # 6.  4th even Fourier coefficient
        'niter':0  , # 7.  Number of iterations
        'nphi' :0  , # 8.  Number of points used for fit
        'B1'   :0.0, # 9.  1st even Fourier coefficient
        'A1'   :0.0, # 10. 1st odd Fourier coefficient
        'B2'   :0.0, # 11. 2nd even Fourier coefficient
        'A2'   :0.0, # 12. 2nd odd Fourier coefficient
        'ellip':0.0, # 13. Ellipticity
        'pa'   :0.0, # 14. Position angle
        'x0'   :0.0, # 15. X center
        'y0'   :0.0, # 16. Y center
        'A3'   :0.0, # 17. 3rd odd Fourier coefficient
        'A4'   :0.0, # 18. 4th odd Fourier coefficient
        'err'  :1.0, # Error
        'hmax' :'' , # Maximum harmonic
        'b'    :1.0, # Semi-minor axis length
        'c'    :2.0, # Disky/Boxy exponent
        'fitMeth':''}

kwsDEF={
    # nphi
    'nphimin':5,
    # mk_funcs
    'avgaz':True,'avgref':False,'avgwin':0.,
    # 'meth_interp':'linear',
    # 'meth_interp':'nearest',
    'meth_interp':'RectBivariateSpline',
    'nsect':64,'radwid':0.05,
    #'meth_interp':'linear','radwid':0.05,
    'meth_derivt':'rings' ,'derwid':0.05,
    # getharmonics
    'meth_error':'rms', 
    # fitharm
    'meth_fharm':'matrix','center':True,'fitI0':True,
    # maxharm
    'constrain':False, 
    # checkfit
    'errtol':None,'miniter':8,'maxiter':25,
    'dertol':0.5,'rmstol':1.,
    # adjparams 
    #    sign method performs better for higher adjfact, but requires more iterations to converge
    #    full method performs better for lower adjfact at fewer iterations
    'meth_adjust':'full','adjfact':1.01,'adjharm':True,'adjiter':True,'flipellip':False
    }

kwsTEST={
    # INTERPOLATION
    # Interpolation methods
    'meth_interp':['linear','cubic','RectBivariateSpline'],
    # Number of sectors
    'nsect'      :[32,64,128],
    # Width of sector rings radiall
    'radwid'     :[0.01,0.05,0.1],
    # Minimum number of phi points sampled
    'nphimin'    :[3,5,10],
    # DERIVATIVE
    # Derivative method
    'meth_derivt':['rings','major','minor','max'],
    # Width of ring derivative
    'derwid'     :[0.01,0.05,0.1],
    # Tolerance below which derivative is unreliable
    'dertol'     :[0.1,0.5,1.0],
    # FITTING
    # Error method
    'meth_error' :['amp','rms'],
    # Error tolerance
    'errtol'     :[0.01,0.04,0.1], # For rms, 0.001 for amp
    # RMS tolerance for tiny derivative
    'rmstol'     :[0.1,1.0,10.0],
    # Numbers of iterations
    'miniter'    :[5,8,10],
    'maxiter'    :[25,50,100],
    # OPTIONS
    # Fit center of ellipses
    'center'     :[True,False],
    # Effect of averaging azimuthally
    'avgaz'      :[True,False],
    # Effect of averaging under reflection
    'avgref'     :[True,False],
    # Force concentric ellipses
    'concentric' :[False,True],
    # Adjust harmonics to fit constraints
    'meth_adjust':['full','sign'],
    'constrain'  :[True,False],
    'adjharm'    :[True,False],
    'adjfact'    :[1.0,1.01,1.05,1.1],
    'adjiter'    :[True,False],
    'flipellip'  :[True,False]
    }


####################################################################################################################################
# TEST
def test(method,vals=None,atest=3.,plotFlag=True,**kws):
    """Test ellipse fitting"""
    import time,pickle
    # Pars input
    kws.setdefault('verbose',True)
    if vals is None: vals=kwsTEST[method]
    # File names
    histfile=kws.pop('histfile',FTEST)
    plotfile=kws.pop('plotfile',None)
    if plotfile is None:
        plotfile=os.path.join(os.path.dirname(histfile),'testellipse_{}.png'.format(method))
    # Load data
    if not kws.get('imgdat',None):
        f=open(histfile,'r')
        kws['imgdat'],kws['imgdim']=parshist(pickle.load(f))
        f.close()
    # Loop over vals and run
    ell=[] ; times=[] ; lablist=[] ; plist=[]
    for v in vals:
        kws[method]=v
        tic=time.clock()
        ip,piter=fit(atest,retIter=True,**kws)
        plist.append(ip)
        ell.append(piter)
        toc=time.clock()
        times.append(toc-tic)
        lablist.append('{} = {}'.format(method,v))
    # Plot
    if plotFlag:
        import matplotlib.pyplot as plt
        clrs=['b','r','g','m']
        lw=4.
        f,axs=plt.subplots(2,2,squeeze=False)
        # Plot image and final fits
        plot(plist,imgdat=kws['imgdat'],imgdim=kws['imgdim'],axs=axs[0][0],
             label=lablist,color=clrs,linewidth=lw)
        for i,v in enumerate(vals):
            label=lablist[i]
            # Print
            print label
            print '    time    = {}'.format(times[i])
            print '    niter   = {niter}'.format(**plist[i])
            print '    err     = {err}'.format(**plist[i])
            print '    rms     = {rms}'.format(**plist[i])
            print '    fitMeth = {fitMeth}'.format(**plist[i])
            x=np.arange(len(ell[i]))+1
            #x=np.array([p['niter'] for p in ell[i]])
            # Plot err vs niter
            iy=np.array([p['err'] for p in ell[i]])
            axs[1][0].plot(x,iy,clrs[i],label=label,lw=lw)
            axs[1][0].set_xlabel('niter')
            axs[1][0].set_ylabel('err')
            # Plot rms vs niter
            iy=np.array([p['rms'] for p in ell[i]])
            axs[1][1].plot(x,iy,clrs[i],label=label,lw=lw)
            axs[1][1].set_xlabel('niter')
            axs[1][1].set_ylabel('rms')
        axs[1][1].legend()
        # Save figure
        plt.savefig(plotfile)
        print ' '+plotfile
    # Return
    return method,vals,ell,times


####################################################################################################################################
# FITTING METHODS
def fitmult(imgdat,imgdim=None,p0=None,filename=None,overwrite=False,
            alist=None,Na=10,amin=None,amax=None,ascl='loglin',abase=0.1,
            concentric=False,reverse=False,verbose=False,vv=False,**kws):
    """Fit multiple ellipses to an image at different radii"""
    # Load and return filename if it exists
    if isinstance(filename,str):
        if os.path.isfile(filename) and not overwrite:
            return read(filename,retlist=True)
    # Initialize parameters
    for k in kwsDEF: kws.setdefault(k,kwsDEF[k])
    Nx,Ny=imgdat.shape
    if imgdim is None: imgdim=((-float(Nx)/2.,float(Nx)/2.),
                               (-float(Ny)/2.,float(Ny)/2.))
    ximg=imgdim[0][1]-imgdim[0][0]
    yimg=imgdim[1][1]-imgdim[1][0]
    xbin=ximg/float(Nx)
    ybin=yimg/float(Ny)
    xybin=max(xbin,ybin)
    if p0 is None: p0=copy.deepcopy(PARAM0)
    # Create list of radii
    if isinstance(alist,list):
        amin=min(alist)
        amax=max(alist)
        Na=len(alist)
    else:
        if amin==None: amin=1.*xybin
        if amax==None: amax=min(ximg/2.,yimg/2.)
        amin=10.**np.log10(amin)
        amax=10.**np.log10(amax)
        if   ascl=='lin' : alist=list(np.linspace(amin,amax,Na))
        elif ascl=='log' : alist=list(np.logspace(np.log10(amin),np.log10(amax),Na))
        elif ascl=='log' : alist=list(np.logspace(mmlmath.logN(amin/amin,1.+abase),
                                                  mmlmath.logN(amax/amin,1.+abase),
                                                  num=Na,base=1.+abase))
        elif ascl=='loglin':
            alistlog=list(np.logspace(np.log10(amin),np.log10(amax),Na))#np.ceil(float(Na)*logfrac)))
            alistlin=list(np.linspace(alistlog[-1],amax,Na+1))#int(Na*logfrac))+1)
            alist=alistlog#+alistlin[1:]
        elif ascl=='fill': 
            Na=int(np.ceil(np.log(amax/amin)/np.log(1.+abase)))+1
            alist=list(amin*((1.+abase)**np.arange(Na)))
        else: raise Exception('Invalid ascl: {}'.format(ascl))
    # Reverse
    if reverse: alist=alist[::-1]
    # Functions
    finterp,fderivt=mk_funcs(imgdat,imgdim=imgdim,**kws)
    # Loop over radii
    plist=[] ; ip=p0
    for ia in alist:
        # Break if limit exceeded
        if   ia<amin: continue
        elif ia>amax: continue
        # Get number of phi bins
        nphi=int(2.*np.pi*ia/xybin)
        # Fit
        ip=fit(ia,p0=copy.deepcopy(ip),finterp=finterp,fderivt=fderivt,nphi=nphi,
               verbose=verbose,vv=False,**kws)
        # Break for overlap
        if len(plist)!=0 and concentric:
            if chkoverlap(ip,plist[-1],err=xybin): continue
        # Add parameters to list
        #pprint.pprint(ip)
        plist.append(ip)
    # Save
    if filename: save(filename,plist,overwrite=overwrite)
    # Return fit parameters
    return plist


def fit(a,imgdat=None,p0=None,finterp=None,fderivt=None,retIter=False,verbose=False,vv=False,**kws):
    """
    Fits single ellipse
    """
    for k in kwsDEF: kws.setdefault(k,kwsDEF[k])
    # Pars input
    if kws['errtol']==None:
        if   kws['meth_error']=='rms': kws['errtol']=0.04
        elif kws['meth_error']=='amp': kws['errtol']=0.001
        else: raise Exception('Invalid meth_error: {}'.format(kws['meth_error']))
    # Initialize parameters
    if p0 is None: p0=copy.deepcopy(PARAM0)
    p0['a']=float(a)
    p0['nphi']=calc_nphi(a,imgdat=imgdat,verbose=verbose,**kws)
    p=copy.deepcopy(p0)
    p['err']=1.0
    p['niter']=0
    pmin=None
    # Initialize functions
    if finterp==None or fderivt==None: finterp,fderivt=mk_funcs(imgdat,**kws)
    # Get list of sampling points
    phiLIST=np.linspace(0.,2.*np.pi,p0['nphi'],endpoint=False)
    # Get position of maxima
    if p['pa']==0:
        cp=copy.deepcopy(p)
        cp['ellip']=0.
        cimg=finterp(cp,phiLIST)
        p['pa']+=phiLIST[np.argmax(cimg)]
        p['pa']=p['pa']%np.pi
    # Iterate to solution
    fitFlag=False ; pIter=[]
    while not fitFlag:
        p['niter']+=1
        # Get image and derivative
        iimg=finterp(p,phiLIST)
        ider=fderivt(p,phiLIST)
        # Fit harmonics
        p=getharmonics(p,phiLIST,iimg,ider,verbose=vv,**kws)
        # Get minimum overall
        if pmin==None or abs(p[p['hmax']])<abs(pmin[pmin['hmax']]):
            pmin=copy.deepcopy(p)
        # Check if parameters meet requirements
        pIter.append(copy.deepcopy(p))
        fitFlag,p=checkfit(p,pmin,verbose=verbose,**kws)
        # Break if fit found
        if fitFlag:
            if verbose: mmlio.verbose('Breaking...')
        # Adjust parameters
        else: p=adjparams(p,p['hmax'],p[p['hmax']],verbose=vv,**kws)
    # Fit higher harmonics
    fharm3,frms3=fitharm(phiLIST,iimg,ord_beg=3,ord_end=3,**kws)
    fharm4,frms4=fitharm(phiLIST,iimg,ord_beg=4,ord_end=4,**kws)
    p['A3']=fharm3['A3'] ; p['B3']=fharm3['B3']
    p['A4']=fharm4['A4'] ; p['B4']=fharm4['B4']
    # Modulate position angle
    if p['pa']<0: p['pa']+=(2.*np.pi)
    p['pa']=p['pa']%np.pi
    # Pars parameters
    for k in PARAM0.keys(): p.setdefault(k,PARAM0[k])
    p['b']=p['a']*(1.-p['ellip'])
    # Return params
    if verbose:
        mmlio.verbose('Parameters:')
        pprint.pprint(p)
    if retIter: return p,pIter
    else      : return p

def checkfit(p,pmin,miniter=8,maxiter=20,errtol=None,dertol=0.5,rmstol=1.,verbose=False,meth_error='rms',**kws):
    """Determine if parameters meet fit requirements"""
    p['fitMeth']='None'
    fitFlag=False
    # Pars input
    if errtol==None:
        if   meth_error=='rms': errtol=0.04
        elif meth_error=='amp': errtol=0.001
        else: raise Exception('Invalid meth_error: {}'.format(meth_error))
    # Handle error requirement
    if not fitFlag and (abs(p['err']) < errtol and p['niter'] > miniter):
        if verbose: mmlio.verbose('Error tolerance reached [err={},errtol={}].'.format(p['err'],errtol))
        p['fitMeth']='Error tolerance'
        fitFlag=True
    # Handle tiny derivative
    if not fitFlag and (p['der']==0 or (abs(p['der']) < dertol*p['rms'] and p['rms'] < rmstol)):
        if verbose: mmlio.verbose('Slope tolerance reached [der={},dertol={},rms={}].'.format(p['der'],dertol,p['rms']))
        p=copy.deepcopy(pmin)
        p['fitMeth']='derivative'
        fitFlag=True
    # Handle max iterations
    if not fitFlag and p['niter'] >= maxiter:
        if verbose: mmlio.verbose('Maximum number of iterations reached [niter={}].'.format(p['niter']))
        p=copy.deepcopy(pmin)
        p['fitMeth']='maxiter'
        fitFlag=True
    # Return
    return fitFlag,p


def adjparams(p,hmax,coef,verbose=False,meth_adjust='full',adjharm=False,adjfact=1.5,adjiter=True,
              flipellip=False,**exkw):
    """Adjust params to negate harmonics"""
    if verbose: mmlio.verbose('Maximum harmonic: {}={}'.format(hmax,coef))
    #mmlio.yorn('deltaa={} or {} [hmax={}]'.format(-coef/p['der'],np.sign(coef)/(adjfact**(p['niter']-1)),coef))
    if   meth_adjust=='full': deltaa=-coef/p['der']
    elif meth_adjust=='sign': deltaa=np.sign(coef)*abs(coef)/10. #; adjiter=True
    else: raise Exception('Invalid meth_adjust: {}'.format(meth_adjust))
    if adjiter: deltaa/=(adjfact**(p['niter']-1))
    deltab=deltaa*(1.-p['ellip'])
    # Adjust harmonics to prevent invalid corrections
    if adjharm:
        if verbose: mmlio.verbose('Adjusting harmonic to meet criteria.')
        if   hmax=='A1': 
            while deltab>abs(p['b']): deltab/=adjfact
        elif hmax=='B1': 
            while deltab>abs(p['b']): deltab/=adjfact
        elif hmax=='A2': 
            # while abs(-2.*deltab/(p['a']*((1.-p['ellip'])**2 - 1.)))>=np.pi/2.: deltab/=adjfact
            pass
        elif hmax=='B2': 
            # while p['ellip']+2.*deltab/p['a']< 0: deltab/=adjfact
            # while p['ellip']+2.*deltab/p['a']>=1: deltab/=adjfact
            pass
        else: raise Exception('Invalid index for maximum harmonic: {}={}'.format(hmax,coef))
    # Adjust parameters based on hmax
    if verbose: mmlio.verbose('Adjusting corresponding parameter.')
    if   hmax=='A1': # A1: minor axis position
        p['x0']+=deltab*np.cos(p['pa']+np.pi/2.)
        p['y0']+=deltab*np.sin(p['pa']+np.pi/2.)
    elif hmax=='B1': # B1: major axis position
        p['x0']+=deltaa*np.cos(p['pa'])
        p['y0']+=deltaa*np.sin(p['pa'])
    elif hmax=='A2': # A2: position angle
        p['pa']+=-2.*deltab/(p['a']*((1.-p['ellip'])**2 - 1.))
    elif hmax=='B2': # B2: ellipticity
        p['ellip']+=2.*deltab/p['a']
        # Ellipticity above 1 (negative b)
        if p['ellip']>1 and adjharm: 
            if verbose: mmlio.verbose('Ellipticity above 1. Adjusting.')
            p['ellip']=1.-abs(1.-p['ellip'])
        # Ellipticity below 0 (b exceeds a)
        if p['ellip']<0 and adjharm: 
            if flipellip:
                newe=p['ellip']-2.*deltab/p['a']
                newe+=(-2.*deltab/p['a'])*(0.05)
                p['pa'   ]+=np.pi/2.
            else: 
                newe=1.-1./(1.-p['ellip'])
                p['pa'   ]+=np.pi/2.
            if verbose or flipellip:
                mmlio.verbose('Ellipticity below 0 ({}). Adjusting to {}.'.format(p['ellip'],newe))
            p['ellip']=newe
    else: raise Exception('Invalid index for maximum harmonic: {}={}'.format(hmax,coef))
    p['b']=p['a']*(1.-p['ellip'])
    # Ensure correct parameters
    # Return adjusted params
    return p

####################################################################################################################################
# METHOD FOR GETTING HARMONIC TERM
def mk_fharm(m):
    order=float(m[-1])
    if   m=='I0': fharm=lambda x: np.ones(len(x))
    elif m[0]=='A': fharm=lambda x: np.sin(order*x)
    elif m[0]=='B': fharm=lambda x: np.cos(order*x)
    else: raise Exception('Invalid order: {}'.format(m))
    return fharm
def ferr_harm(v,x,y,I0,ordlist,var=1.):
    v=list(v)
    np.seterr(all='raise')
    nord=len(ordlist)
    if   len(v)==2*nord+1: I0=v.pop(0)
    elif len(v)==2*nord  : pass
    else: raise Exception('Invalid coefficient length for {} orders.'.format(nord))
    I=I0
    for idx in range(nord):
        I+=v[2*idx  ]*np.sin(float(ordlist[idx])*x)
        I+=v[2*idx+1]*np.cos(float(ordlist[idx])*x)
    err=I-y
    return var*err
def frms_harm(*args,**kwargs):
    rms=np.sqrt(np.mean(np.power(ferr_harm(*args,**kwargs),2.)))
    return rms

####################################################################################################################################
# METHODS FOR FITTING HARMONICS
def getharmonics(p,phi,img,der,meth_error='rms',verbose=False,**exkw):
    """Determine error from harmonic fit"""
    # Assign parameters
    p['I0' ]=np.mean(img)
    p['der']=der
    # Fit harmonics
    harm,rms=fitharm(phi,img,**exkw)
    p.update(**harm)
    p['rms' ]=rms
    # Get maximum harmonic
    p['hmax']=maxharm(p,harm,verbose=verbose,**exkw)
    # Get error
    if   meth_error=='amp': p['err']=abs(harm[p['hmax']])
    elif meth_error=='rms':
        if p['rms']==0: p['err']=float('Inf')
        else          : p['err']=abs(harm[p['hmax']]/p['rms'])
    else: raise Exception('Invalid meth_error: {}'.format(meth_error))
    # Return params
    return p

def fitharm(phi,img,meth_fharm='matrix',ord_beg=None,ord_end=2,center=True,fitI0=True,warn=True,**exkw):
    """Fit harmonics to image along ellipse"""
    from numpy.linalg import lstsq
    # Pars input
    I0=np.mean(img)
    if np.isnan(I0): raise Exception('Image data contains a NaN')
    if ord_beg is None:
        if center: ord_beg=1
        else     : ord_beg=2
    ordlist=range(ord_beg,ord_end+1)
    # Create list of coefficients
    if fitI0: ordcoef=['I0'] ; nimg=0.
    else    : ordcoef=[]     ; nimg=I0
    for iord in ordlist: ordcoef.extend(['A{}'.format(iord),'B{}'.format(iord)])
    # Matrix solutions
    if meth_fharm=='matrix':
        # Create functions
        ordfunc=[mk_fharm(c) for c in ordcoef]
        nrow=len(ordfunc)
        # Allocate
        B=np.zeros((nrow))
        A=np.zeros((nrow,nrow))
        # Populate matrices
        for irow in range(nrow):
            B[irow]=np.sum((img-nimg)*ordfunc[irow](phi))
            for icol in range(nrow):
                A[irow,icol]=np.sum(ordfunc[irow](phi)*ordfunc[icol](phi))
        # Division
        out=lstsq(A,B,rcond=1.)
        v=out[0]
        success=1
    # Least squares fitting
    elif meth_fharm=='leastsq':
        from scipy.optimize import leastsq
        # Get initial guess
        if vdict0 is None: 
            vdict0={k:0. for k in ordcoef}
            vdict0['I0']=I0
        v0=[vdict0[k] for k in ordcoef]
        # Fit
        v,success=leastsq(ferr_harm,v0,args=(phi,img,I0,ordlist,1.))
    # Error
    else: raise Exception('Invalid meth_fharm: {}'.format(meth_fharm))
    # Check for NaNs
    vnanchk=np.isnan(v)
    if any(vnanchk): 
        raise Exception('Detected NaN coefficients: {}'.format(np.array(ordcoef)[vnanchk]))
    # Calculate error
    rms=frms_harm(v,phi,img,I0,ordlist,var=1.)
    # Distribute coefficients in dictionary
    vdict={ordcoef[idx]:v[idx] for idx in range(nrow)}
    if not fitI0: vdict['I0']=I0
    # Print things if there are errors
    if success > 4: mmlio.verbose('fit success: {}'.format(success))
    if rms==0 and img.min()!=img.max() and warn:
        mmlio.verbose('fit success: {}'.format(success))
        mmlio.verbose('rms={}'.format(rms))
        mmlio.verbose('v={}'.format(v))
        mmlio.verbose('img (min,max)=({},{})'.format(img.min(),img.max()))
        mmlio.verbose('I0={}'.format(I0))
        pprint.pprint(vdict)
        # mmlio.yorn('Continue?')
    # Return harmonics and rms error
    return vdict,rms

def maxharm(p,harm,center=True,constrain=False,verbose=False,**exkw):
    """Get maximum harmonic"""
    p['b']=p['a']*(1.-p['ellip'])
    # Sort harmonics in order of amplitude
    hord=sorted(harm.keys(),key=lambda k: abs(harm[k]),reverse=True)
    # Handle 0 derivative
    if p['der']==0.: return hord[0]
    # Skip in invalid values
    imax=-1 ; valid=False
    while not valid:
        imax+=1 ; ivalid=True
        hmax=hord[imax]
        # Continue if I0 somehow included
        if   hmax=='I0': ivalid=False
        # Dont use A1 or B1 if center not fitted
        elif hmax in ['A1','B1'] and not center: ivalid=False
        # Don't use A2 if ellipticity is zero
        elif hmax=='A2' and p['ellip']==0: ivalid=False
        # Just for completeness
        elif hmax=='B2': pass
        # Handle constraints
        if ivalid and constrain:
            # Get adjustment
            padj=adjparams(copy.deepcopy(p),hmax,harm[hmax],verbose=verbose,**exkw)
            #deltab=padj['b']-p['b']
            # Dont use A1 or B1 if center changes by more than semi-minor length
            if hmax in ['A1','B1']:
                dr=np.sqrt((padj['x0']-p['x0'])**2.+(padj['y0']-p['y0'])**2.)
                if dr > p['b']: ivalid=False
            # Don't use A2 if pa changes by more than pi/2
            elif hmax=='A2':
                pass
                #if abs(padj['pa']-p['pa'])>=np.pi/2.: ivalid=False
            # Don't use B2 if ellipticity becomes <0 or >1
            elif hmax=='B2':
                if   padj['ellip']<0: ivalid=False
                elif padj['ellip']>1: ivalid=False
        # Set validity
        valid=ivalid
        if not ivalid and verbose:
            mmlio.verbose('Skipping {}={}'.format(hmax,hord[hmax]))
    # Return valid max harmonic
    return hmax


def calcharmonics(phi,img,ord_beg=None,ord_end=2,center=True,fitI0=True,vdict0=None):
    """Calculate harmonics based on matrix algebra Athanassoula et al. (1990)"""
    from numpy.linalg import lstsq
    # Set constants
    minord=1
    maxord=4
    # Pars input
    if np.isnan(np.mean(img)): raise Exception('Image data contains a NaN')
    if ord_beg is None:
        if center: ord_beg=1
        else     : ord_beg=2
    ordlist=range(ord_beg,ord_end+1)
    # Create list of functions
    if fitI0: 
        ordstrg=['I0']
        ordfunc=[lambda x: np.ones(len(x))]
        I0=0.
    else    : 
        ordstrg=[]
        ordfunc=[]
        I0=np.mean(img)
    for iord in ordlist:
        astr='A{}'.format(iord)
        bstr='B{}'.format(iord)
        ordfunc.extend([mk_fharm(astr),mk_fharm(bstr)])
        ordstrg.extend([astr,bstr])
    nrow=len(ordfunc)
    # Allocate
    B=np.zeros((nrow))
    A=np.zeros((nrow,nrow))
    # Populate matrices
    for irow in range(nrow):
        B[irow]=np.sum((img-I0)*ordfunc[irow](phi))
        for icol in range(nrow):
            A[irow,icol]=np.sum(ordfunc[irow](phi)*ordfunc[icol](phi))
    # Division
    out=lstsq(A,B,rcond=1.)
    v=out[0]
    rms=frms_harm(v,phi,img,I0,ordlist,var=1.)
    # Distribute coefficients
    vdict={ordstrg[idx]:v[idx] for idx in range(nrow)}
    if not fitI0: vdict['I0']=I0
    return vdict,rms


#def fitharm(phi,img,meth_fharm='matrix',ord_beg=None,ord_end=2,center=True,fitI0=True,warn=True,**exkw):
# def fitharmonics(phi,img,vdict0=None,center=True,ord_beg=None,ord_end=2,
#                  fitI0=True,verbose=False,warn=True):
#     """Fit harmonics to image residual along ellipse"""
#     from scipy.optimize import leastsq
#     # Set constants
#     minord=1
#     maxord=4
#     astr=['A{}'.format(iord) for iord in range(minord,maxord+1)]
#     bstr=['B{}'.format(iord) for iord in range(minord,maxord+1)]
#     vstrall=['I0']
#     for iord in range(len(astr)): vstrall+=[astr[iord],bstr[iord]]
#     # Parse input
#     if np.isnan(np.mean(img)): raise Exception('Image data contains a NaN')
#     if vdict0==None: 
#         vdict0={ivstr:0. for ivstr in vstrall}
#         vdict0['I0']=np.mean(img)
#     if ord_beg is None:
#         if center: ord_beg=1
#         else     : ord_beg=2
#     # Initialize vector of coefficients
#     ordlist=range(ord_beg,ord_end+1)
#     aOrd=[] ; bOrd=[]
#     if fitI0: vOrd=['I0']
#     else    : vOrd=[]
#     for iord in ordlist:
#         iAstr='A{}'.format(iord)
#         iBstr='B{}'.format(iord)
#         aOrd.append(iAstr)
#         bOrd.append(iBstr)
#         vOrd.extend([iAstr,iBstr])
#     v0=[vdict0[ivstr] for ivstr in vOrd]
#     # Fit
#     fitargs=(phi,img,vdict0['I0'],ordlist,1.)
#     v,success=leastsq(ferr_harm,v0,args=fitargs)
#     rms=np.sqrt(np.power(ferr_harm(v,*fitargs),2).sum()/float(len(phi)))
#     rms2=frms_harm(v,*fitargs)
#     if success > 4:
#         mmlio.verbose('fit success: {}'.format(success))
#     if rms==0 and img.min()!=img.max() and warn:
#     # if np.allclose(rms,0.):
#         mmlio.verbose('fit success: {}'.format(success))
#         mmlio.verbose('rms={}'.format(rms))
#         mmlio.verbose('v={}'.format(v))
#         mmlio.verbose('img (min,max)=({},{})'.format(img.min(),img.max()))
#         mmlio.verbose('I0={}'.format(vdict0['I0']))
#         #mmlio.yorn('Continue?')
#     # Check for NaNs
#     vnanchk=np.isnan(v)
#     if any(vnanchk):
#         raise Exception('Detected NaN coefficients: {}'.format(np.array(vOrd)[vnanchk]))
#     # Set excluded coefficients to 0
#     vdict={vOrd[idx]:v[idx] for idx in range(len(vOrd))}
#     if not fitI0: vdict['I0']=vdict0['I0']
#     for ikey in vstrall:
#         if ikey not in vdict:
#             pass
# #            vdict[ikey]=0.
# #            if verbose: mmlio.verbose('{}={} (excluded)'.format(ikey,vdict[ikey]))
#         else:
#             if verbose: mmlio.verbose('{}={} (included)'.format(ikey,vdict[ikey]))
#     # Return fit parameters
#     return vdict,rms

def maxharmonic(p,harm,imgcen=(0.,0.),center=True,adjharm=False,verbose=False):
    """Get maximum harmonic"""
    p['b']=p['a']*(1.-p['ellip'])
    # Sort harmonics in order of amplitude
    hord=sorted(harm.keys(),key=lambda k: abs(harm[k]),reverse=True)
    # Handle 0 derivative
    if p['der']==0.: return hord[0]
    # Skip in invalid values
    imax=-1 ; valid=False
    while not valid:
        imax+=1 ; ivalid=True
        hmax=hord[imax]
        # Get changes to a & b
        deltaa=-harm[hmax]/p['der']
        deltab=deltaa*(1.-p['ellip'])
        # Continue if I0 somehow included
        if hmax=='I0': ivalid=False
        # Dont use A1 if center changes by more than semi-minor length
        if hmax=='A1':
            if not adjharm:
                ix0=p['x0']+deltab*np.cos(p['pa']+np.pi/2)
                iy0=p['y0']+deltab*np.sin(p['pa']+np.pi/2)
                idr=np.sqrt((ix0-imgcen[0])**2 + (iy0-imgcen[1])**2)
                if idr > p['b']: ivalid=False
            if not center: ivalid=False
        # Don't use B1 if center changes by more than semi-minor length
        if hmax=='B1':
            if not adjharm:
                ix0=p['x0']+deltaa*np.cos(p['pa'])
                iy0=p['y0']+deltaa*np.sin(p['pa'])
                idr=np.sqrt((ix0-imgcen[0])**2 + (iy0-imgcen[1])**2)
                if idr > p['b']: ivalid=False
            if not center: ivalid=False
        # Don't use A2 if ellipticity is zero
        if hmax=='A2':
            # if not adjharm:
            #     if abs(-2.*deltab/(p['a']*((1.-p['ellip'])**2 - 1.)))>=np.pi/2.: ivalid=False
            if p['ellip']==0: ivalid=False 
        # Don't use B2 if ellipticity becomes smaller than zero
        if hmax=='B2':
            if not adjharm: 
                if   p['ellip']+2.*deltab/p['a']< 0: ivalid=False
                # elif p['ellip']+2.*deltab/p['a']> 1: ivalid=False
        # Set validity
        valid=ivalid
        if not ivalid and verbose:
            mmlio.verbose('Skipping {}={} (deltab={})'.format(hmax,hord[hmax],deltab))
    # Return valid max harmonic
    return hmax

####################################################################################################################################
# METHOD FOR CREATING FUNCTIONS
def mk_funcs(imgdat,imgdim=None,meth_interp='linear',meth_derivt='rings',
             derwid=0.05,radwid=0.05,nsect=64,avgaz=True,avgref=False,
             avgwin=0.,**kws):
    """
    Create function for interpolating image and derivatives
        avgaz : average image azimuthally
        avgref: average image under reflection across max in phi
    """
    import scipy.interpolate as interpolate
    # Create xy coordinate system
    Nx,Ny=imgdat.shape
    if not imgdim: imgdim=((-float(Nx)/2.,float(Nx)/2.),
                           (-float(Ny)/2.,float(Ny)/2.))
    xbin=(imgdim[0][1]-imgdim[0][0])/float(Nx)
    ybin=(imgdim[1][1]-imgdim[1][0])/float(Ny)
    x=np.linspace(imgdim[0][0]+xbin/2.,imgdim[0][1]-xbin/2.,Nx)
    y=np.linspace(imgdim[1][0]+ybin/2.,imgdim[1][1]-ybin/2.,Ny)
    # Average
    avgbin=int(avgwin/max(xbin,ybin))
    if avgbin>0: imgdat=mmlmath.runavg2d(imgdat,window=avgbin)
    # Grid xy and flatten arrays
    yy,xx=np.meshgrid(y,x)
    fxx=xx.flatten()
    fyy=yy.flatten()
    fzz=imgdat.flatten()
    # Create interpolation subfunction
    if meth_interp=='RectBivariateSpline':
        #from scipy.interpolate import RectBivariateSpline
        from hacks.fitpack2 import RectBivariateSpline
        finterp_sub0=RectBivariateSpline(x,y,imgdat)
        finterp_sub=lambda pxy: finterp_sub0.ev(*pxy)
    else:
        from scipy.interpolate import griddata
        finterp_sub=lambda pxy: griddata((fxx,fyy),fzz,pxy,method=meth_interp,fill_value=0.)
    # Create interpolation function, averaging azimuthally
    def finterp(p,phi,nsect=nsect,radwid=radwid,maxpix=20,verbose=False):
        try             : nphi=len(phi)
        except TypeError: nphi=1
        phi=phi%(2.*np.pi)
        if avgaz: phi1=np.append(phi,(phi+np.pi)%(2.*np.pi))
        else    : phi1=phi
        if nphi<=64:
        #if (p['a']/max(xbin,ybin))<maxpix:
            ixy1=phi2xy(phi1,**p)
            iz1=finterp_sub(ixy1)
        else:
            if verbose: mmlio.verbose('nphi = {}: Using sectors'.format(len(phi)))
            swid=2.*np.pi/float(nsect)
            zsect=sectors(p,fxx,fyy,fzz,nsect=nsect,radwid=radwid)
            isect1=np.array(phi1/swid,dtype=int)%nsect
            iz1=np.array([zsect[s] for s in isect1],dtype=float)
        # Average azimuthally
        if avgaz: iz=(iz1[:nphi]+iz1[nphi:])/2.
        else    : iz=iz1
        # Average under reflection
        if avgref: iz=(iz+iz[::-1])/2.
            # phimax1=phi[np.argmax(iz)]
            # phimax2=(phimax1+np.pi)%(2.*np.pi)
        return iz
    # Create derivative function
    # Ring method of derivative
    if meth_derivt=='rings':
        def fderivt(p,phi):
            p_inn=copy.deepcopy(p) ; p_inn['a']*=1.-derwid
            p_out=copy.deepcopy(p) ; p_out['a']*=1.+derwid
            iz_inn=finterp(p_inn,phi,verbose=False)
            iz_out=finterp(p_out,phi,verbose=False)
            izder=np.mean(iz_out-iz_inn)/(p_out['a']-p_inn['a'])
            return izder
    # Linear along semi major axis
    elif meth_derivt=='major':
        def fderivt(p,phi):
            phi=np.array(p['pa'])
            p_inn=copy.deepcopy(p) ; p_inn['a']*=1.-derwid
            p_out=copy.deepcopy(p) ; p_out['a']*=1.+derwid
            iz_inn=finterp(p_inn,phi,verbose=False)
            iz_out=finterp(p_out,phi,verbose=False)
            izder=np.mean(iz_out-iz_inn)/(p_out['a']-p_inn['a'])
            return izder
    # Linear along semi minor axis
    elif meth_derivt=='minor':
        def fderivt(p,phi):
            phi=np.array(p['pa'])+np.pi/2.
            p_inn=copy.deepcopy(p) ; p_inn['a']*=1.-derwid
            p_out=copy.deepcopy(p) ; p_out['a']*=1.+derwid
            iz_inn=finterp(p_inn,phi,verbose=False)
            iz_out=finterp(p_out,phi,verbose=False)
            izder=np.mean(iz_out-iz_inn)/(p_out['a']-p_inn['a'])
            return izder
    # Linear along maximum
    elif meth_derivt=='max':
        def fderivt(p,phi):
            p_inn=copy.deepcopy(p) ; p_inn['a']*=1.-derwid
            p_out=copy.deepcopy(p) ; p_out['a']*=1.+derwid
            iz_inn=finterp(p_inn,phi,verbose=False)
            iz_out=finterp(p_out,phi,verbose=False)
            izder=(iz_out-iz_inn)/(p_out['a']-p_inn['a'])
            return izder[np.argmax(np.abs(izder))]
    # Error
    else: raise Exception('Invalid method: {}'.format(method))
    #mmlio.yorn(meth_derivt)
    # Return functions
    return finterp,fderivt

####################################################################################################################################
# Utility functions
def calc_nphi(a,imgdat=None,imgdim=None,nphi=None,nphimin=5,verbose=False,**exkw):
    """Determine how many times to bin in phi"""
    if nphi==None:
        if imgdim==None: pixsiz=1.
        else           : pixsiz=max((imgdim[0][1]-imgdim[0][0])/imgdat.shape[0],
                                    (imgdim[1][1]-imgdim[1][0])/imgdat.shape[1])
        nphi=int(2.*np.pi*float(a)/float(pixsiz))
    nphi=max(nphi,nphimin)
    if nphi==0: raise Exception('nphi cannot be 0')
    else      : 
        if verbose: mmlio.verbose('a = {:5.2f} & nphi = {:3d}'.format(a,nphi))
    return nphi

def sectors(p,x,y,z,nsect=64,radwid=0.1):
    """Iterpolate by sector in phi"""
    swid=2.*np.pi/float(nsect)
    rad0=np.sqrt((x-p['x0'])**2+(y-p['y0'])**2)
    phi0=np.arctan2((y-p['y0']),(x-p['x0'])) ; phi0[phi0<0]+=(2.*np.pi)
    # Get inner/outer radii for each phi position
    pmin=copy.deepcopy(p) ; pmin['a']*=(1.-radwid)
    pmax=copy.deepcopy(p) ; pmax['a']*=(1.+radwid)
    radmin=phi2rad(phi0,**pmin)
    radmax=phi2rad(phi0,**pmax)
    radidx=np.logical_and(np.greater_equal(rad0,radmin),np.less(rad0,radmax))
    # Average in each sector
    zs=np.zeros((nsect,),float)
    for idx in range(nsect):
        phimin=swid*float(idx) ; phimax=swid*float(idx+1)
        phiidx=np.logical_and(phi0>=phimin,phi0<phimax)
        totidx=np.logical_and(radidx,phiidx)
        if np.sum(totidx)==0:
            # mmlio.verbose('sphi = {}'.format((phimin+phimax)/2.))
            # mmlio.verbose('    nRadidx = {}'.format(np.sum(radidx)))
            # mmlio.verbose('    nPhiidx = {}'.format(np.sum(phiidx)))
            # mmlio.verbose('    nTotidx = {}'.format(np.sum(totidx)))
            zs[idx]=0.
        else: zs[idx]=np.mean(z[totidx])
    # Return
    return zs

def phi2xy(phi,**p):
    """Returns xy positions for given phi values"""
    phi=phi%(2.*np.pi)
    phi[phi<0]+=(2.*np.pi)
    # Postions relative to major/minor axes
    p['b']=p['a']*(1.-p['ellip'])
    # xp=p['a']*np.cos(phi)
    # yp=p['b']*np.sin(phi)
    # rp=np.sqrt((xp**2.)+(yp**2.))
    # Get positions
    x=p['x0']+p['a']*np.cos(p['pa'])*(np.cos(phi)**(2./p['c']))-p['b']*np.sin(p['pa'])*(np.sin(phi)**(2./p['c']))
    y=p['y0']+p['a']*np.sin(p['pa'])*(np.cos(phi)**(2./p['c']))+p['b']*np.cos(p['pa'])*(np.sin(phi)**(2./p['c']))
    # Return positions
    return x,y

def phi2rad(phi,**p):
    """Returns radii (relative to ellipse center) for given phi values"""
    x,y=phi2xy(phi,**p)
    r=np.sqrt((x-p['x0'])**2+(y-p['y0'])**2)
    return r

def chkoverlap(p1,p2,err=0.,nphi=10):
    """Checks for overlap between two ellipses"""
    phiLIST=np.linspace(0.,2*np.pi,nphi)
    x1,y1=phi2xy(phiLIST,**p1) ; r1=np.sqrt(x1**2+y1**2)
    x2,y2=phi2xy(phiLIST,**p2) ; r2=np.sqrt(x2**2+y2**2)
    return not (np.all(np.greater(r1+err,r2)) or np.all(np.greater(r2+err,r1)))

####################################################################################################################################
# METHODS FOR FILE HANDLING
def plist2dict(plist,sortkey=None,array=False):
    """Transform list of parameter dictionaries to dictionary of parameter lists"""
    # Parse input
    if   isinstance(plist,list): pass
    elif isinstance(plist,dict):
        if isinstance(plist['a'],list): return plist
        else                          : plist=[plist]
    else: raise Exception('Invalid plist type: {}'.format(type(plist)))
    # Get keys
    klist=sorted(PARAM0.keys())
    # Create dictionary
    pdict={k:[] for k in klist}
    for idx in range(len(plist)):
        for k in klist: pdict[k].append(plist[idx][k])
    # Sort 
    if sortkey in klist:
        idxsort=np.argsort(np.array(pdict[sortkey]))
        for k in klist: pdict[k]=list(np.array(pdict[k])[idxsort])
    # Convert to arrays
    if array:
        for k in klist: pdict[k]=np.array(pdict[k])
    # Return dictionary of lists
    return pdict

def pdict2list(pdict):
    """Transform dictionary of parameter lists to list of parameter dictionaries"""
    # Get keys
    klist=sorted(PARAM0.keys())
    # Parse input
    if   isinstance(pdict,list): return pdict
    elif isinstance(pdict,dict):
        if   isinstance(pdict['a'],list      ): pass
        elif isinstance(pdict['a'],np.ndarray): pdict={k:np.array(pdict[k]) for k in klist} 
        else                                  : return [pdict]
    else: raise Exception('Invalid pdict type: {}'.format(type(pdict)))
    # Create list
    plist=[]
    for idx in range(len(pdict[klist[0]])):
        ip={k:pdict[k][idx] for k in klist}
        plist.append(ip)
    # Return list of dictionaries
    return plist

def save(fname,plist,overwrite=False,verbose=False):
    """Save list of parameters to file"""
    # Covert list to dictionary
    pdict=plist2dict(plist)
    # Check overwrite
    if not mmlfiles.prep_write(fname,overwrite): return
    # Write to file
    mmlio.rwtable('W',fname,pdict,overwrite=overwrite)
    if verbose: mmlio.verbose('Wrote ellipse fits to file:')
    if verbose: print '    '+fname
    # Return control
    return

def read(fname,retlist=False):
    """Read dictionary of parameters from file"""
    pdict=mmlio.rwtable('R',fname)
    if retlist: return pdict2list(pdict)
    else      : return pdict

def plot(plist,imgdat=None,imgdim=None,axs=None,label=None,
         nphi=100,color='m',linewidth=4.0,cmap='pnbody_light',
         fname=None,overwrite=False,verbose=False,**exkw):
    """Plot ellipses"""
    import matplotlib.pyplot as plt
    # Pars input
    plist=pdict2list(plist)
    if isinstance(color,list): clist=color
    else                     : clist=len(plist)*[color]
    if label is None: label=['a = {}'.format(ip['a']) for ip in plist]
    # Initialize plot
    if axs is None:
        axs=plt.subplot(1,1,1)
    # Plot image
    if imgdat!=None:
        Nx,Ny=imgdat.shape
        if imgdim==None: imgdim=((-float(Nx)/2.,float(Nx)/2.),
                                 (-float(Ny)/2.,float(Ny)/2.))
        extent=(imgdim[0][0],imgdim[0][1],imgdim[1][0],imgdim[1][1])
        mmlplot.hist2imshow(imgdat,cmapname=cmap,extent=extent,axs=axs)
    axs.set_autoscale_on(False)
    # Plot ellipses
    phi=np.linspace(0.,2.*np.pi,nphi,endpoint=False)
    for i,ip in enumerate(plist):
        ix,iy=phi2xy(phi,**ip)
        axs.plot(ix,iy,color=clist[i],linewidth=linewidth,label=label[i])
    # Set limits
    plt.xlim=imgdim[0]
    plt.ylim=imgdim[1]
    # Save
    fig=plt.gcf()
    fig.tight_layout()
    if isinstance(fname,str):
        if mmlfiles.prep_write(fname,overwrite): 
            fig.savefig(fname)
            plt.close(fig)
            return None
    return fig

def parshist(hist):
    """Covert histogram input to ellipse parameters"""
    # Image dimensions
    imgdim=[hist['xlim'],hist['ylim']]
    # Scaled image data
    clrkw=dict(clip=True,mn=float(hist['wlim'][0]),mx=float(hist['wlim'][1]),scale=hist['wscl'])
    if clrkw['scale']=='log': clrkw['cd']=10.**((np.log10(clrkw['mx'])+np.log10(clrkw['mn']))/2.-1.)
    imgdat=mmlplot.map_array(hist['ws'],**clrkw)
    if clrkw['scale']=='log': imgdat[hist['ws']<=0]=0
    # Return results
    return imgdat,imgdim

