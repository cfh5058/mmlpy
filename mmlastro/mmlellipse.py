#!/usr/bin/python
from mmlutils import *
import mmlutils
import numpy as np
import scipy.interpolate as interpolate
from scipy.optimize import leastsq
import matplotlib.pyplot as plt
import os,pprint,copy

####################################################################################################################################
# METHOD TO PLOT MULTIPLE ELLIPSES
def plotellipse(paramlist,imgdat=None,imgdim=None,nphi=None,fname=None,
                overwrite=None,verbose=None,ellipcolor=None,palette=None):
    """
    Plots multiple ellipses over an image
    """
    # Set constants
    nphiDEF=100
    # Pars input
    nphi=mmlpars.mml_pars(nphi,type=int,default=nphiDEF,min=0)
    if isinstance(paramlist,dict): paramlist=[paramlist]
    paramlist=mmlpars.mml_pars(paramlist,type=list)
    overwrite=mmlpars.mml_pars(overwrite,default=False,type=bool)
    verbose=mmlpars.mml_pars(verbose,default=False,type=bool)
    ellipcolor=mmlpars.mml_pars(ellipcolor,default='b')
    ellipwidth=mmlpars.mml_pars(ellipwidth,default=4.0)
    palette=mmlpars.mml_pars(palette,default='pnbody_light',type=str)
    # Plot image
    fig=plt.figure()
    ax1=plt.subplot(1,1,1)
    if imgdat !=None:
        extent=(imgdim[0][0],imgdim[0][1],imgdim[1][0],imgdim[1][1])
        mmlplot.hist2imshow(imgdat,cmapname=palette,extent=extent)
    ax1.set_autoscale_on(False)
    # Loop over parameters adding ellipses
    phiLIST=np.linspace(0.,2*np.pi,nphi)
    for iparam in paramlist:
        ix,iy=phi2xy_ellipse(phiLIST,**iparam)
        ax1.plot(ix,iy,color=ellipcolor,linewidth=10.0,label='a = {}'.format(iparam['a']))#,lw=ellipwidth)
    plt.xlim=(imgdim[0])
    plt.ylim=(imgdim[1])
#    ax1.set_xlim=(imgdim[0])
#    ax1.set_ylim=(imgdim[1])
    # Save or display image
    fig.tight_layout()
    if fname != None:
        fname=mmlpars.mml_pars(fname,type=str)
        if overwrite or not os.path.isfile(fname):
            fig.savefig(fname)
        plt.close(fig)
        return None
    else:
        return fig
#        fig.show()

####################################################################################################################################
def chkoverlap(par0,par1,nphi=100):
    """
    Check for overlap between two ellipses
    """
    phiLIST=np.linspace(0.,2*np.pi,nphi)
    x0,y0=phi2xy_ellipse(phiLIST,**par0) ; r0=np.sqrt(x0**2+y0**2)
    x1,y1=phi2xy_ellipse(phiLIST,**par1) ; r1=np.sqrt(x1**2+y1**2)
    return not (np.all(np.greater(r0,r1)) or np.all(np.greater(r1,r0)))
#    return np.any(np.greater_equal(r0,r1))

####################################################################################################################################
# METHOD TO SAVE ELLIPSE DATA TO FILE
def saveellipse(fname,paramlist,overwrite=None,verbose=None):
    """
    Saves fitted ellipse info to file
    """
    # Pars input
    fname=mmlpars.mml_pars(fname,type=str)
    paramdict=paramlist2dict(paramlist)
    overwrite=mmlpars.mml_pars(overwrite,default=False,type=bool)
    verbose=mmlpars.mml_pars(verbose,default=False,type=bool)
    # Check if file exists
    if not mmlfiles.prep_write(fname,overwrite): return
    # Write output
    mmlio.rwtable('W',fname,paramdict,overwrite=overwrite)
    if verbose: mmlio.verbose('Wrote ellipse fits to file:')
    if verbose: print '    '+fname
    # Return
    return

####################################################################################################################################
# METHOD TO READ ELLIPSE DATA FROM FILE
def readellipse(fname):
    """
    Reads ellipse fit info from a file
    """
    # Pars input
    fname=mmlpars.mml_pars(fname,type=str,isfile=True)
    # Read file
    paramdict=mmlio.rwtable('R',fname)
    # Return
    return paramdict2list(paramdict)
    
####################################################################################################################################
# METHOD TO CONVERT LIST OF PARAMETERS DICTS TO DICT OF PARAMETER LISTS
def paramlist2dict(paramlist):
    """
    Converts list of parameter dictionaries into dictionary of parameter lists
    """
    # Set constants
    keylist=list_ellipseparam()
    # Pars input
    if isinstance(paramlist,dict): paramlist=[paramlist]
    paramlist=mmlpars.mml_pars(paramlist,type=list)
    # Create table
    paramdict={ikey:[] for ikey in keylist}
    for iparam in paramlist:
        iparam=pars_ellipseparam(iparam)
        for ikey in keylist:
            paramdict[ikey].append(iparam[ikey])
    paramdict['keylist']=keylist
    return paramdict

####################################################################################################################################
# METHOD TO CONVERT DICT OF PARAMETER LIST TO LIST OF PARAMETER DICTS
def paramdict2list(paramdict):
    """
    Converts dictionary of parameter lists to list of parameter dictionaries
    """
    # Set constants
    keylist=list_ellipseparam()
    # Pars input
    paramdict=mmlpars.mml_pars(paramdict,type=dict,keylist=keylist)
    lenparam=len(paramdict[keylist[0]])
    # Create list
    paramlist=[]
    for idx in range(lenparam):
        iparam={}
        for ikey in keylist:
            iparam[ikey]=paramdict[ikey][idx]
        paramlist.append(iparam)
    return paramlist

####################################################################################################################################
# METHOD TO FIT MULTIPLE ELLIPSES
def fitmultellip(semimajorlist=None,nsemimajor=None,isemimajor=None,fsemimajor=None,
                 semimajorscl='mult',semimajormult=0.5,
                 params0=None,concentric=False,
                 imgdat=None,imgdim=None,interpmeth=None,
                 nphi=None,errtol=None,dertol=None,derwid=None,miniter=None,maxiter=None,
                 verbose=False,reverse=False,nproc=10,**exkw):
    """
    Fits multiple ellipses to an image
    """
    # Set constants
    nsemimajorDEF=10
    # Pars input
    params0=mmlpars.mml_pars(params0,default={},type=dict)
    params0=pars_ellipseparam(params0)
    imgdat,imgdim,imgsiz,imgbin,imgcen=pars_imagedata(imgdat,imgdim=imgdim)
    params0['x0']=imgcen[0]
    params0['y0']=imgcen[1]
#    mmlio.yorn('center={}'.format(imgcen))
    xbinsize,ybinsize=imgbin
    # Get aliases
    chkalias=['alist','Na','semimajorscale']
    for ichk in chkalias:
        if ichk in exkw: raise Exception('{} outdated'.format(ichk))
    # Get list of semi major axes lengths
    if isinstance(semimajorlist,list):
        isemimajor=min(semimajorlist)
        fsemimajor=max(semimajorlist)
        nsemimajor=len(semimajorlist)
    else:
        isemimajor=mmlpars.mml_pars(isemimajor,type=float,default=max(xbinsize,ybinsize))
        fsemimajor=mmlpars.mml_pars(fsemimajor,type=float,default=min(imgsiz[0]/2.,imgsiz[1]/2.))
        nsemimajor=mmlpars.mml_pars(nsemimajor,type=int,min=1,default=nsemimajorDEF)
        if   semimajorscl=='log': semimajorlist=list(np.logspace(np.log10(isemimajor),np.log10(fsemimajor),nsemimajor))
        elif semimajorscl=='lin': semimajorlist=list(np.linspace(isemimajor,fsemimajor,nsemimajor))
        elif semimajorscl=='mult':
            semimajorlist=[isemimajor]
            while len(semimajorlist)<nsemimajor:
                semimajorlist.append(semimajorlist[-1]*(1.+semimajormult))
        else: raise Exception('Invalide semimajorscl: {}'.format(semimajorscl))
    # Reverse
    if reverse: semimajorlist=semimajorlist[::-1]
    # Initialize functions
    imgfunc,imgderv=imginterp(imgdat,imgdim=imgdim,interpmeth=interpmeth)
    iparam=copy.deepcopy(params0)
    # Loop over semimajor axis lengths fitting
    paramlist=[]
#    na=0
#    while len(paramlist)<=nsemimajor:
        #ia0=semimajorlist[max(nsemimajor-na-1,0)]
#        na+=1
        # Ensure concentric ellipses
        # if concentric:
        #     if len(paramlist)==0: ia=ia0
        #     else                : ia=min(ia0,paramlist[-1]['b'])
        # else: ia=ia0
        #print ia0,ia
    for ia in semimajorlist:
        # Break if minimum exceeded
        if   ia<isemimajor: continue
        elif ia>fsemimajor: continue
        # iparam['ellip']=(params0['ellip']+iparam['ellip'])/2.
        if verbose: mmlio.verbose('Trying a = {:.2f}'.format(ia))
        iparam=fitellipse(ia,params0=iparam,imgfunc=imgfunc,imgderv=imgderv,imgcen=imgcen,
                          nphi=nphi,errtol=errtol,dertol=dertol,derwid=derwid,
                          miniter=miniter,maxiter=maxiter,pixelsize=(xbinsize+ybinsize)/2.,
                          verbose=verbose)
        # Break for overlap
        if len(paramlist)!=0 and concentric:
            if chkoverlap(iparam,paramlist[-1]): continue
        if verbose: 
            mmlio.verbose('Accepted a = {:.2f}'.format(ia))
            mmlio.verbose('Parameters:')
            pprint.pprint(iparam)
        #print ia0,ia
        paramlist.append(iparam)
    # Return fit parameters
    return paramlist

####################################################################################################################################
# METHOD TO FIT SINGLE ELLIPSE
def fitellipse(semimajor,params0=None,
               imgfunc=None,imgderv=None,imgdat=None,imgdim=None,imgcen=None,
               interpmeth=None,derivmeth=None,fitcenter=None,nphi=None,
               errtol=None,dertol=None,derwid=None,pixelsize=1.0,
               miniter=None,maxiter=None,verbose=None):
    """
    Fits an ellipse to an image
    """
    # Set constants
    derivmethLIST=['interp','rings']
    derivmethDEF='rings'
    fitcenterDEF=True
    nphiDEF=20
    # Pars input
    verbose=mmlpars.mml_pars(verbose,default=False,type=bool)
    semimajor=mmlpars.mml_pars(semimajor,type=float,min=0.)
    params0=mmlpars.mml_pars(params0,type=dict,default={})
    if imgfunc==None or imgderv==None:
        imgdat,imgdim,imgsiz,imgbin,imgcen=pars_imagedata(imgdat,imgdim=imgdim)
        imgfunc,imgderv=imginterp(imgdat,imgdim=imgdim,interpmeth=interpmeth)
    derivmeth=mmlpars.mml_pars(derivmeth,default=derivmethDEF,list=derivmethLIST)
    fitcenter=mmlpars.mml_pars(fitcenter,default=fitcenterDEF,type=bool)
    nphi=mmlpars.mml_pars(nphi,type=int,default=nphiDEF,min=0)
    derwid=mmlpars.mml_pars(derwid,type=float,default=0.05)
    dertol=mmlpars.mml_pars(dertol,type=float,default=0.5)
    errtol=mmlpars.mml_pars(errtol,type=float,default=0.04)
    miniter=mmlpars.mml_pars(miniter,type=int,default=8,min=1)
    maxiter=mmlpars.mml_pars(maxiter,type=int,default=20,min=miniter)
    imgcen=mmlpars.mml_pars(imgcen,type=tuple,default=(0.,0.),nelements=2)
    # Pars parameters
    params0=pars_ellipseparam(params0)
    params0['a']=semimajor
    params0['nphi']=nphi
    nphi0=max(nphi,2.*np.pi/(pixelsize/semimajor))
    mmlio.verbose('nphi {} {}'.format(nphi,nphi0))
    # Initialize stuff
    phiLIST0=np.linspace(0.,2*np.pi,nphi0,endpoint=False)
    phiLIST =np.linspace(0.,2*np.pi,nphi ,endpoint=False)
    params=copy.deepcopy(params0)
    params['err']=1.
    params['niter']=0
    params_min=None
    # Loop until under errtol or max # iterations reached
    fitflag=False
    while not fitflag:
        params['niter']+=1
        if derivmeth=='rings':
            params_out=copy.deepcopy(params) ; params_out['a']*=1.+derwid
            params_inn=copy.deepcopy(params) ; params_inn['a']*=1.-derwid
        # Get positions for provided ellipse
        ix,iy=phi2xy_ellipse(phiLIST,**params)
        if derivmeth=='rings':
            ix_out,iy_out=phi2xy_ellipse(phiLIST,**params_out)
            ix_inn,iy_inn=phi2xy_ellipse(phiLIST,**params_inn)
        # Interpolate
        iimg=imgfunc(ix,iy)
        if derivmeth=='rings':
            iimg_out=imgfunc(ix_out,iy_out)
            iimg_inn=imgfunc(ix_inn,iy_inn)
            iimg_dif=iimg_out-iimg_inn
            ider=np.mean(iimg_dif)/(params_out['a']-params_inn['a'])
        elif derivmeth=='interp':
            ider=float(imgderv(params['a'],params['pa'],params['x0'],params['y0']))
        else: raise Exception('Invalid derivative method: {}'.format(derivmeth))
        params['I0' ]=np.mean(iimg)
        params['der']=ider
        # Fit harmonics
        iharm,irms=fitharmonics(phiLIST,iimg)#,fitcenter=fitcenter)#,verbose=verbose)
        params['rms']=irms
        idxhmax=get_maxharm(params,iharm,imgcen,fitcenter=fitcenter)
        # Get error
        params['err']=abs(iharm[idxhmax]/params['rms'])
        for ikey in iharm.keys(): params[ikey]=iharm[ikey]
        # Get minimum overall
        if params_min==None: params_min=copy.deepcopy(params)
        else:
            if abs(iharm[idxhmax])<(params_min['err']*params_min['rms']):
                params_min=copy.deepcopy(params)
        # Handle maximum number of iterations
        if params['niter'] >= maxiter:
            if verbose: mmlio.verbose('Maximum number of iterations reached [niter={}].'.format(params['niter']))
            params=copy.deepcopy(params)
            fitflag=True
        # Handle tiny derivative
        if abs(params['der']) < dertol*params['rms']:
            if verbose: mmlio.verbose('Slope tolerance reached [der={},dertol={}].'.format(params['der'],dertol))
            fitflag=True
        # Handle error
        if abs(params['err']) < errtol and params['niter'] > miniter:
            if verbose: mmlio.verbose('Error tolerance reached [err={},errtol={}].'.format(params['err'],errtol))
            fitflag=True
        # Break if fit found
        if fitflag:
            if verbose: mmlio.verbose('Breaking...')
        # Adjust appropriate parameter
        else:
            params=adjparams(params,iharm,idxhmax,verbose=verbose)
    # Get final fit to harmonics
    fharm3,frms3=fitharmonics(phiLIST,iimg,ord_beg=3,ord_end=3)
    fharm4,frms4=fitharmonics(phiLIST,iimg,ord_beg=4,ord_end=4)
    params['A3']=fharm3['A3'] ; params['B3']=fharm3['B3']
    params['A4']=fharm4['A4'] ; params['B4']=fharm4['B4']
    # Ensure good behavior of position angle
    if params['pa']<0: params['pa']+=(2.*np.pi)
    params['pa']=params['pa']%np.pi
    # Return parameters
    return params

####################################################################################################################################
# METHOD TO RETURN MAXIMUM HARMONIC
def get_maxharm(params,harm,imgcen,fitcenter=True):
    """
    Returns the key for the maximum harmonic
    """
    # Get maximum harmonic value
    hord=sorted(harm.keys(),key=lambda k: abs(harm[k]),reverse=True)
    # Handle 0 derivative
    if params['der']==0.: return hord[0]
    # Change if invalid value
    imax=-1
    valid=False
    while not valid:
        ivalid=True
        imax+=1
        hmax=hord[imax]
        deltaa=-harm[hmax]/params['der']
        deltab=deltaa*(1.-params['ellip'])
        if hmax=='A1':
            ix0=params['x0']+deltab*np.cos(params['pa']+np.pi/2)
            iy0=params['y0']+deltab*np.sin(params['pa']+np.pi/2)
            idr=np.sqrt((ix0-imgcen[0])**2 + (iy0-imgcen[1])**2)
            if idr > params['b']: ivalid=False
            if not fitcenter: ivalid=False
        if hmax=='B1':
            ix0=params['x0']+deltaa*np.cos(params['pa'])
            iy0=params['y0']+deltaa*np.sin(params['pa'])
            idr=np.sqrt((ix0-imgcen[0])**2 + (iy0-imgcen[1])**2)
            if idr > params['b']: ivalid=False
            if not fitcenter: ivalid=False
        if hmax=='A2' and params['ellip']==0.: ivalid=False
        if hmax=='B2' and params['ellip']+2.*deltab/params['a']<0: ivalid=False
        valid=ivalid
    # Return valid value
    return hmax

####################################################################################################################################
# METHOD TO ADJUST ELLIPSE PARAMETERS
def adjparams(params,harm,hmax,verbose=False):
    """
    Adjusts and return params
    """
    deltaa=-harm[hmax]/params['der']
    deltab=deltaa*(1.-params['ellip'])
    # Adjust parameters based on hmax
    if verbose: mmlio.verbose('Maximum harmonic: {}={}'.format(hmax,harm[hmax]))
    if   hmax=='A1': # A1: minor axis position
        params['x0']+=deltab*np.cos(params['pa']+np.pi/2.)
        params['y0']+=deltab*np.sin(params['pa']+np.pi/2.)
    elif hmax=='B1': # B1: major axis position
        params['x0']+=deltaa*np.cos(params['pa'])
        params['y0']+=deltaa*np.sin(params['pa'])
    elif hmax=='A2': # A2: position angle
        params['pa']+=-2.*deltab/(params['a']*((1.-params['ellip'])**2 - 1.))
    elif hmax=='B2': # B2: ellipticity
        params['ellip']+=2.*deltab/params['a']
    else: raise Exception('Invalid index for maximum harmonic: {}={}'.format(hmax,harm[hmax]))
    params['b']=params['a']*(1.-params['ellip'])
    # Return adjusted params
    return params

####################################################################################################################################
# METHOD TO FIT HARMONICS
def fitharmonics(phi,img,vdict0=None,fitcenter=None,ord_beg=None,ord_end=None,verbose=None):
    """
    Returns the least squares fit to the harmonics
    """
    # Set constants
    minord=1
    maxord=4
    astr=['A{}'.format(iord) for iord in range(minord,maxord+1)]
    bstr=['B{}'.format(iord) for iord in range(minord,maxord+1)]
    vstrall=[]
    for iord in range(len(astr)): vstrall+=[astr[iord],bstr[iord]]
    vdict0DEF={ivstr:0. for ivstr in vstrall}
    # Pars input
    phi=mmlpars.mml_pars(phi,type=np.ndarray,ndim=1)
    img=mmlpars.mml_pars(img,type=np.ndarray,ndim=1,nelements=len(phi))
    vdict0=mmlpars.mml_pars(vdict0,type=dict,keylist=vstrall,default=vdict0DEF) ; vdict0['I0']=np.mean(img)
    fitcenter=mmlpars.mml_pars(fitcenter,default=True,type=bool)
    verbose=mmlpars.mml_pars(verbose,default=False,type=bool)
    if np.isnan(vdict0['I0']): raise Exception('Image data contains a NaN')
    # Set orders
    if fitcenter: ord_begDEF=1
    else        : ord_begDEF=2
    ord_beg=mmlpars.mml_pars(ord_beg,default=ord_begDEF,type=int,min=minord,max=maxord)
    ord_end=mmlpars.mml_pars(ord_end,default=2,type=int,min=ord_beg,max=maxord)
    ordlist=range(ord_beg,ord_end+1)
    # Initalize vector of coefficients
    astrord=[] ; bstrord=[] ; vstrord=[]
    for iord in ordlist:
        astrord.append('A{}'.format(iord))
        bstrord.append('B{}'.format(iord))
        vstrord.append('A{}'.format(iord))
        vstrord.append('B{}'.format(iord))
#    vstrord=astrord+bstrord
    v0=[vdict0[ivstr] for ivstr in vstrord]
    # Fit
    #imgvar=np.sum(img-np.mean(img))/len(img)
    v,success=leastsq(ferr_harm,v0,args=(phi,img,vdict0['I0'],ordlist,1.))
    rms=np.sqrt(np.power(ferr_harm(v,phi,img,vdict0['I0'],ordlist,1.),2).sum()/len(phi))
    if success > 4:
        if verbose: mmlio.verbose('fit success: {}'.format(success))
    if np.allclose(rms,0.):
        if verbose:
            mmlio.verbose('rms=0')
            mmlio.verbose('v={}'.format(v))
            mmlio.verbose('img (min,max)=({},{})'.format(img.min(),img.max()))
            mmlio.verbose('I0={}'.format(vdict0['I0']))
            mmlio.yorn('Continue?')
    # Check for nans
    vnanchk=np.isnan(v)
    if any(vnanchk):
        raise Exception('Detected NaN coefficients: {}'.format(np.array(vstrord)[vnanchk]))
    # Set excluded coefficients to 0
    vdict={vstrord[idx]:v[idx] for idx in range(len(vstrord))}
    for ikey in vstrall:
        if ikey not in vdict:
            pass
#            vdict[ikey]=0.
#            if verbose: mmlio.verbose('{}={} (excluded)'.format(ikey,vdict[ikey]))
        else:
            if verbose: mmlio.verbose('{}={} (included)'.format(ikey,vdict[ikey]))
    # Return fit parameters
    return vdict,rms

####################################################################################################################################
# METHOD TO DETERMINE ERROR IN HARMONICS
def ferr_harm(v,x,y,I0,ordlist,var=1.):
    np.seterr(all='raise')
    nord=len(v)/2
    if len(ordlist) != nord: raise Exception('List of orders must be half the length of coefficients')
    I=I0
    for idx in range(nord):
        I+=v[2*idx  ]*np.sin(float(ordlist[idx])*x)
        I+=v[2*idx+1]*np.cos(float(ordlist[idx])*x)
    err=I-y
    return var*err


####################################################################################################################################
# METHOD TO INTERPOLATE IMAGE ALONG ELLIPSE
def imgellip(imgdat,imgsiz,params,maxpixels=20,window=0.05):
    """
    Returns intensity of an image as a function of phi
    """
             
####################################################################################################################################
# METHOD TO INTERPOLATE IMAGE
def imginterp(imgdat,pos=None,imgdim=None,interpmeth=None):
    """
    Returns the intensity of an image interpolated at specified position
    """
    # Pars input
    imgdat,imgdim,imgsiz,imgbin,imgcen=pars_imagedata(imgdat,imgdim=imgdim)
    xbinsize,ybinsize=imgbin
    nxbin,nybin=imgdat.shape
    if pos==None:
        flagfunc=True
    else:
        flagfunc=False
        pos=mmlpars.mml_pars(pos,type=list,nelements=2)
        pos[0]=mmlpars.mml_pars(pos[0],type=float,range=list(imgdim[0]))
        pos[1]=mmlpars.mml_pars(pos[1],type=float,range=list(imgdim[1]))
    interpmeth=mmlpars.mml_pars(interpmeth,type=str,default='cubic')
    # Get x & y bins for each pixel (at center of pixels)
    ixbin=imgdim[0][0]+xbinsize/2. ; fxbin=imgdim[0][1]-xbinsize/2.
    iybin=imgdim[1][0]+ybinsize/2. ; fybin=imgdim[1][1]-ybinsize/2.
    xbins=np.linspace(ixbin,fxbin,nxbin)
    ybins=np.linspace(iybin,fybin,nybin)
    yybins,xxbins=np.meshgrid(ybins,xbins)
    # Create interpolation function
    ndata=float(len(xxbins.flatten()))
    interfunc = lambda x,y: interpolate.griddata((xxbins.flatten(),yybins.flatten()),imgdat.flatten(),(x,y),
                                                 method=interpmeth,fill_value=0.)
    # Create derivative funciton
    derivfunc = imgderivt(interfunc,min(fxbin,fybin),min(nxbin,nybin),interpmeth=interpmeth)
    # Return functions
    if flagfunc: return interfunc,derivfunc
    # Interpolate
    else:
        interimag=interfunc(pos[0],pos[1])
        derivimag=derivfunc(np.sqrt(pos[0]**2+pos[1]**2),np.arctan(pos[1]/pos[0]))
        return interimag,derivimag

####################################################################################################################################
# METHOD TO RETURN FUNCTION SPECIFYING IMAGE DERIVATIVE
def imgderivt(interfunc,frbin,nrbin,interpmeth=None):
    """
    Returns a function that provides the image derivative along a ray
    """
    # Pars input
    irbin=0.
    frbin=mmlpars.mml_pars(frbin,type=float,min=0.)
    nrbin=mmlpars.mml_pars(nrbin,type=int,min=1)
    interpmeth=mmlpars.mml_pars(interpmeth,type=str,default='cubic')
    # Get bins
    rbins=np.linspace(irbin,frbin,nrbin)
    # Create derivative 
    derivfunc = lambda r,phi,x0,y0: interpolate.griddata((rbins.flatten()[:-1]+np.diff(rbins.flatten())/2.),
                                                         np.diff(interfunc(x0+rbins*np.cos(phi),y0+rbins*np.sin(phi)).flatten())/np.diff(rbins.flatten()),
                                                         (r),method=interpmeth,fill_value=0.)
    return derivfunc

####################################################################################################################################
# METHOD TO PARS DICTIONARY OF ELLIPSE PARAMETERS
def list_ellipseparam():
    """
    Returns a list of keys for ellipse parameters
    """
    listparam=['I0'   , # 1.  Mean intensity
               'B3'   , # 2.  3rd even Fourier coefficient
               'der'  , # 3.  Derivative of intensity
               'a'    , # 4.  Semi-major axis
               'rms'  , # 5.  RMS residual of intensity
               'B4'   , # 6.  4th even Fourier coefficient
               'niter', # 7.  Number of iterations
               'nphi' , # 8.  Number of points used for fit
               'B1'   , # 9.  1st even Fourier coefficient
               'A1'   , # 10. 1st odd Fourier coefficient
               'B2'   , # 11. 2nd even Fourier coefficient
               'A2'   , # 12. 2nd odd Fourier coefficient
               'ellip', # 13. Ellipticity
               'pa'   , # 14. Position angle
               'x0'   , # 15. X center
               'y0'   , # 16. Y center
               'A3'   , # 17. 3rd odd Fourier coefficient
               'A4'   , # 18. 4th odd Fourier coefficient
               'err'  ] # Error
    return listparam
def pars_ellipseparam(param):
    """
    Pars a dictionary of ellipse parameters
    """
    form={
        'A1'   : mmlpars.parsdict(default=0.0 ,type=float       ),
        'A2'   : mmlpars.parsdict(default=0.0 ,type=float       ),
        'A3'   : mmlpars.parsdict(default=0.0 ,type=float       ),
        'A4'   : mmlpars.parsdict(default=0.0 ,type=float       ),
        'B1'   : mmlpars.parsdict(default=0.0 ,type=float       ),
        'B2'   : mmlpars.parsdict(default=0.0 ,type=float       ),
        'B3'   : mmlpars.parsdict(default=0.0 ,type=float       ),
        'B4'   : mmlpars.parsdict(default=0.0 ,type=float       ),
        'I0'   : mmlpars.parsdict(default=0.0 ,type=float,min=0.),
        'a'    : mmlpars.parsdict(default=1.0 ,type=float,min=0.),
        'pa'   : mmlpars.parsdict(default=0.0 ,type=float       ),
        'ellip': mmlpars.parsdict(default=0.01,type=float       ),
        'x0'   : mmlpars.parsdict(default=0.0 ,type=float       ),
        'y0'   : mmlpars.parsdict(default=0.0 ,type=float       ),
        'rms'  : mmlpars.parsdict(default=0.0 ,type=float,min=0.),
        'der'  : mmlpars.parsdict(default=0.0 ,type=float       ),
        'err'  : mmlpars.parsdict(default=1.0 ,type=float,min=0.),
        'niter': mmlpars.parsdict(default=1   ,type=int  ,min=1 ),
        'nphi' : mmlpars.parsdict(default=1   ,type=int  ,min=1 )}
    try:
        param=mmlpars.mml_formpars(param,form)
    except:
        pprint.pprint(param)
        raise
    param['b']=param['a']*(1.-param['ellip'])
    return param

####################################################################################################################################
# METHOD TO PARS IMAGE DATA
def pars_imagedata(imgdat,imgdim=None):
    """
    Parses image data
    """
    # Pars input
    imgdat=mmlpars.mml_pars(imgdat,type=np.ndarray,ndim=2)
    nxbin,nybin=imgdat.shape
    imgdimDEF=[[0.,float(nxbin)],[0.,float(nybin)]]
    imgdim=mmlpars.mml_pars(imgdim,default=imgdimDEF,type=list,nelements=2)
    for idim in range(2):
        imgdim[idim]=mmlpars.mml_pars(imgdim[idim],type=list,nelements=2)
        imgdim[idim][0]=mmlpars.mml_pars(imgdim[idim][0],type=float)
        imgdim[idim][1]=mmlpars.mml_pars(imgdim[idim][1],type=float)
    # Determine image size
    imgsiz=(imgdim[0][1]-imgdim[0][0],imgdim[1][1]-imgdim[1][0])
    # Determine binsizes
    xbinsize=imgsiz[0]/nxbin
    ybinsize=imgsiz[1]/nybin
    imgbin=(xbinsize,ybinsize)
    # Get bins
    xbins=np.linspace(imgdim[0][0]+xbinsize/2.,imgdim[0][1]-xbinsize/2.,nxbin)
    ybins=np.linspace(imgdim[1][0]+ybinsize/2.,imgdim[1][1]-ybinsize/2.,nybin)
    # Determine center
    idxmax=np.argmax(imgdat.flatten())
    xbinmax,ybinmax=np.unravel_index(idxmax,(nxbin,nybin))
    xmax=xbins[xbinmax]
    ymax=ybins[ybinmax]
    imgcen=(xmax,ymax)
    # print 'imgdim={}'.format(imgdim)
    # print 'nbins=({},{})'.format(nxbin,nybin)
    # print 'idx max=({},{},{})'.format(idxmax,xbinmax,ybinmax)
    # print 'img max,mean={},{}'.format(imgdat.flatten()[idxmax],np.mean(imgdat))
#    mmlio.yorn('Continue?')
    # Return output
    return imgdat,imgdim,imgsiz,imgbin,imgcen

####################################################################################################################################
# METHOD TO RETURN XY POSITIONS FOR ANGLES ALONG AN ELLIPSE
def phi2xy_ellipse(phi,**params):
    """
    Returns xy positions for angles along an ellipse
    """
    # Pars input
    params=pars_ellipseparam(params)
    # Determine positions relative to major/minor axis
    xp=params['a']*np.cos(phi)
    yp=params['b']*np.sin(phi)
    rp=np.sqrt((xp**2.)+(yp**2.))
#    rp=params['a']*params['b']/np.sqrt((params['a']*np.sin(phi))**2. + (params['b']*np.cos(phi))**2.)
    # Determine positions
#    x=params['x0']+rp*np.cos(params['pa'])
#    y=params['y0']+rp*np.sin(params['pa'])
#    x=params['x0']+xp*np.cos(params['pa'])-yp*np.sin(params['pa'])
#    y=params['y0']+xp*np.sin(params['pa'])+yp*np.cos(params['pa'])
    x=params['x0']+params['a']*np.cos(params['pa'])*np.cos(phi)-params['b']*np.sin(params['pa'])*np.sin(phi)
    y=params['y0']+params['a']*np.sin(params['pa'])*np.cos(phi)+params['b']*np.cos(params['pa'])*np.sin(phi)
    # Return positions
    return x,y

if __name__ == '__main__':
    main()
