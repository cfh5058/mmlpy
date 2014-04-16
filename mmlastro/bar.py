#!/usr/bin/python
from mmlutils import *
import numpy as np
import os,pprint,copy
import mmlastro.ellipse as efit

FTEST={'ellipse':efit.FTEST,
       'scfAm2' :'/home/langmm/code/python/mine/mmlastro/scfAm2_test.pickle'}

def test(method,**kws):
    """Test bar identification routines"""
    import time,pickle
    # Set defaults
    kws.setdefault('verbose',True)
    mfile=FTEST[method]
    if method=='ellipse':
        mfunc=ellipse
    elif method=='scfAm2':
        mfunc=scfAm2
    else: raise Exception('Invalid method: {}'.format(method))
    # File names
    histfile=kws.pop('histfile',mfile)
    plotfile=kws.pop('plotfile',None)
    if plotfile is None:
        plotfile=os.path.join(os.path.dirname(histfile),'testbarfit_{}.png'.format(method))
    # Load data
    if not kws.get('imgdat',None):
        f=open(histfile,'r')
        kws.update(**parshist(method,pickle.load(f)))
        f.close()
    # Run bar fit
    tic=time.clock()
    out=mfunc(plotfile=plotfile,title='test',**kws)
    toc=time.clock()
    # Print and return output
    print 'Method = {} (time={})'.format(method,toc-tic)
    pprint.pprint(out)
    return out

def parshist(method,hist):
    """Pars histogram to the accepted form"""
    if method=='ellipse':
        import mmlastro.ellipse as efit
        imgdat,imgdim=efit.parshist(hist)
        out={'imgdat':imgdat,'imgdim':imgdim}
    elif method=='scfAm2':
        out={'hist':hist}
    else: raise Exception('Invalid method: {}'.format(method))
    return out

def scfAm2(hist,plotfile=None,ampmin=0.04,runmeth='dphidr',
           title=None,owplot=False,**kws):
    """Determine bar parameters from m=2 part of potential"""
    from mmlutils import mmlmath
    from scipy.interpolate import interp1d
    # Set useful constants
    errval=-999
    phiwid=np.pi
    derwid=1 # bins
    # Pars input
    rad=hist['xs']
    phi=hist['ys']
    nbins=hist['ws'].shape[0]
    amplim=hist['wlim']
    phierr=kws.pop('phierr',10.)*np.pi/180.
    binerr=int(np.ceil(phierr/(phi[1]-phi[0])))
    # Grid & average
    PHI,RAD=np.meshgrid(phi,rad)
    avgwin=nbins/20
    amp=mmlmath.runavg2d(hist['ws'],window=avgwin,weights=hist['as'])
    # Pad in phi direction to allow wrap
    phipad=np.append(phi,phi[0:derwid]+phiwid)
    amppad=np.append(amp,amp[:,0:derwid],axis=1)
    # Take derivative of amplitude in phi
    dA_dphi=np.diff(amppad,n=derwid)/np.diff(phipad,n=derwid)
    # Get phi of max amp in regions of positive and negative dA/dphi for each radius
    idxpos=np.ma.argmin(np.ma.masked_array(abs(amp),mask=dA_dphi< 0),axis=1)
    idxneg=np.ma.argmin(np.ma.masked_array(abs(amp),mask=dA_dphi>=0),axis=1)
    phipos=phi[idxpos] ; fphipos=interp1d(rad,phipos,kind='cubic')
    phineg=phi[idxneg] ; fphineg=interp1d(rad,phineg,kind='cubic')
    # Get derivative of phi in radius
    r1,dphidr1_pos,r2,dphidr2_pos=mmlmath.wrapder(rad,phipos,wrpwin=phiwid,retder2=True)
    r1,dphidr1_neg,r2,dphidr2_neg=mmlmath.wrapder(rad,phineg,wrpwin=phiwid,retder2=True)
    dphidr1=mmlmath.runavg1d((dphidr1_pos+dphidr1_neg)/2.,window=5)
    dphidr2=np.diff(dphidr1)
    # Find run where phi does not change with radius
    if   runmeth=='dphidr': 
        dphidr1_avg,idxmin,idxmax=mmlmath.findrun(dphidr1,phierr)
        rmin_d=r1[idxmin]
        rmax_d=r1[idxmax]
    elif runmeth=='phi'   : 
        dphidr1_avg,idxmin,idxmax=mmlmath.findrun(phipos ,phierr)
        rmin_d=rad[idxmin]
        rmax_d=rad[idxmax]
    else: raise Exception('Invalid method for finding runs: {}'.format(runmeth))
    # Get location of maximum in bar
    idxrmin=np.ma.argmin(np.ma.masked_less(rad,rmin_d))
    idxrmax=np.ma.argmax(np.ma.masked_greater(rad,rmax_d))
    idxbar=np.ma.argmax(np.ma.masked_array(amp,mask=RAD>rmax_d))
    idxramp,idxpamp=np.unravel_index(idxbar,amp.shape)
    Amax=float(amp[idxramp,idxpamp])
    # Find middle of the in radius
    idxravg=(idxrmin+idxrmax)/2
    # Get bar extent in phi at the middle
    # phineg is at leading edge of bar
    # phipos is at trailing edge of bar
    idxpmin=idxpos[idxravg] # Beginning of bar
    idxpmax=idxneg[idxravg] # End of bar
    idxpwid=idxpmax-idxpmin
    if idxpwid<0: idxpwid+=nbins
    idxpavg=(idxpmin+idxpwid/2)%nbins
    # Define values
    rmin=rad[idxrmin]
    ravg=rad[idxravg]
    rmax=rad[idxrmax]
    ramp=rad[idxramp]
    pmin=phi[idxpmin]
    pavg=phi[idxpavg]
    pmax=phi[idxpmax]
    pamp=phi[idxpamp]
    # Determine truth of bar presence
    barFlag=False
    if not barFlag and idxmin==idxmax: barFlag='constant_phi'
    if not barFlag and Amax<ampmin   : barFlag='max_amplitude'
    if not barFlag: barFlag=True
    # Define output
    out={'barFlag':str(barFlag),
         'rmin':float(10**rmin),
         'ramp':float(10**ramp),
         'rmax':float(10**rmax),
         'amp' :Amax,
         'phi' :float(pavg)}
    # Plot
    if mmlfiles.prep_write(plotfile,overwrite=owplot):
        import matplotlib.pyplot as plt
        # Prepare
        plt.clf()
        nplt=2
        ncol=2 ; nrow=int(np.ceil(float(nplt)/float(ncol)))
        fig,axs=plt.subplots(nrow,ncol,squeeze=False,figsize=(1.4*4.*ncol,4.*nrow))
        if title: fig.suptitle(title)
        # Set things
        cbs=[] ; nlvl=50 ; ntic=9
        resclr='k' ; reslin='dashed'
        maxclr='k' ; maxlin='solid'
        cmap=plt.get_cmap('jet')
        ampnorm=mmlplot.SymLogNorm(amplim[0],vmin=-amplim[1],vmax=amplim[1],linscale=0.001)
        ticfmt='% .3f'
        # Define lines
        edat={'rad'   :rad,
              'radder':r1,
              'phi'   :phi,
              'dphidr':dphidr1}
        limits={'rad'   :(rad.min(),rad.max()),
                'phi'   :(phi.min(),phi.max()),
                'dphidr':(-phiwid,phiwid),
                'amp'   :amplim}
        lines={'rad'   :[(rmin,'r','dashed','rmin'),
                         (rmax,'r','dashed','rmax'),
                         (ravg,'r','dotted','ravg'),
                         (ramp,'r','solid' ,'ramp')],
               'phi'   :[(pmin,'g','dashed','phimin'),
                         (pmax,'g','dashed','phimax'),
                         (pavg,'g','dotted','phiavg'),
                         (pamp,'g','solid' ,'phiamp')],
               'dphidr':[],
               'amp'   :[]}
        # Plot contours
        lvls=ampnorm.inverse(np.linspace(0.,1.,nlvl,endpoint=True))
        tics=ampnorm.inverse(np.linspace(0.,1.,ntic,endpoint=True))
        cn=axs[-1][-1].contour(rad,phi,np.transpose(amp),lvls,norm=ampnorm,extend='both')
        cf=axs[-1][-1].contourf(rad,phi,np.transpose(amp),lvls,norm=ampnorm,extend='both')
        plotprofile(axs[-1][-1],rad,fphipos(rad),ls='-r',ywid=phiwid,
                    vlines=lines['rad'],hlines=lines['phi'],
                    xlab='log(R)',ylab='mod phi',
                    xlim=limits['rad'],ylim=limits['phi'])
        plotprofile(axs[-1][-1],rad,fphineg(rad),ls='-b',ywid=phiwid,
                    xlab='log(R)',ylab='mod phi',
                    xlim=limits['rad'],ylim=limits['phi'])
        axs[-1][-1].scatter(ramp,pamp,c='k',marker='x')
        # Plot derivative
        plotprofile(axs[0][0],r1,dphidr1_pos,ls='-r',
                    vlines=lines['rad'],hlines=lines['dphidr'],
                    xlab='log(R)',ylab='dphi/dlogR',
                    xlim=limits['rad'],ylim=limits['dphidr'])
        plotprofile(axs[0][0],r1,dphidr1_neg,ls='-b',
                    xlab='log(R)',ylab='dphi/dlogR',
                    xlim=limits['rad'],ylim=limits['dphidr'])
        plotprofile(axs[0][0],r1,dphidr1,ls='-k',
                    xlab='log(R)',ylab='dphi/dlogR',
                    xlim=limits['rad'],ylim=limits['dphidr'])
        plotprofile(axs[0][0],r2,dphidr2,ls='--k',
                    xlab='log(R)',ylab='dphi/dlogR',
                    xlim=limits['rad'],ylim=limits['dphidr'])
        # Plot profiles
        # xv='a'
        # for i,yv in enumerate(yvars):
        #     iax=axs[i/ncol][i%ncol]
        #     if yv=='phi': ywid=phiwid
        #     else        : ywid=None
        #     plotprofile(iax,edat0[xv],edat0[yv],vlines=lines[xv],hlines=lines[yv],
        #                 xlab=xv,ylab=yv,xlim=limits[xv],ylim=limits[yv],
        #                 xscl=scales[xv],yscl=scales[yv],ywid=ywid)
        # Save plot
        plt.savefig(plotfile,bbox_inches='tight')
        print '    '+plotfile
    return out


def ellipse(imgdat,plotfile=None,owplot=False,title=None,runmeth='phi',**kws):
    """Determine bar parameters from ellipse fits"""
    import mmlastro.ellipse as efit
    from mmlutils import mmlmath
    # Set useful constants
    errval=-1
    phiwid=np.pi
    binerr=3
    # Pars input
    phierr=kws.pop('phierr',10.)*np.pi/180.
    ellerr=kws.pop('ellerr',0.1)
    ellmin=kws.pop('ellmin',0.25)
    rmscut=kws.pop('rmscut',5.)
    errcut=kws.pop('errcut',1.)
    cutmeth=kws.pop('cutmeth','rms')
    # Get ellipse fits
    edat0=efit.plist2dict(efit.fitmult(imgdat,**kws),sortkey='a',array=True)
    edat0['Am4' ]=edat0['A4']/edat0['a']
    edat0['hval']=np.array([edat0[edat0['hmax'][i]][i] for i in range(len(edat0['hmax']))])
    if runmeth=='dphidr':
        edat0['dpdr']=np.append(1,mmlmath.wrapder(edat0['a'],edat0['pa'],wrpwin=phiwid)[1])
    # Make cuts
    idxfit=np.where(edat0['fitMeth']=='Error tolerance')
    idxrms=np.where(abs(edat0['rms'])<rmscut)[0]
    idxerr=np.where(abs(edat0['err'])<errcut)[0]
    if   cutmeth=='fit': idxcut=idxfit
    elif cutmeth=='rms': idxcut=idxrms
    elif cutmeth=='err': idxcut=idxerr
    else: raise Exception('Invalid cut method: {}'.format(cutmeth))
    edat={k:edat0[k][idxcut] for k in edat0.keys()}
    # Extract variables
    a=edat['a']
    b=edat['b']
    p=edat['pa']%phiwid
    e=edat['ellip']
    err=edat['err']
    rms=edat['rms']
    Am4=edat['Am4']
    # Find longest run of constant position angle
    if   runmeth=='dphidr':
        derwid=2
        dpdr1_avg=mmlmath.runavg1d(edat['dpdr'],window=derwid)
        print edat['dpdr'].shape,dpdr1_avg.shape
        print edat0['dpdr'].shape
        print edat0['a'].shape
        pavg,idxmin_p,idxmax_p=mmlmath.findrun(dpdr1_avg,phierr,wrap=phiwid)
    elif runmeth=='phi':
        pavg,idxmin_p,idxmax_p=mmlmath.findrun(p,phierr,wrap=phiwid)
    else: raise Exception('Invalid method for finding runs: {}'.format(runmeth))
    pavg=np.mean(p[idxmin_p:(idxmax_p+1)])
    pmin=(pavg-phierr/2.)
    if pmin<0: pmin+=phiwid
    pmin=pmin%phiwid
    pmax=(pavg+phierr/2.)%phiwid
    # Find inner and outer edge using ellipticity
    ediff=np.append([0],np.diff(e[idxmin_p:(idxmax_p+binerr)]))
    idxerise=np.where(ediff>0       )[0]
    idxefall=np.where(ediff<=-ellerr)[0]
    if len(idxerise)>0: idxmin_e=idxerise[ 0]+idxmin_p
    else              : idxmin_e=errval
    if len(idxefall)>0: idxmax_e=idxefall[-1]+idxmin_p
    else              : idxmax_e=errval
    # Define inner and outer edge
    idxmin=idxmin_p
    if idxmax_e>(idxmax_p+binerr): idxmax=idxmax_p
    elif idxmax_e==errval        : idxmax=idxmax_p
    else                         : idxmax=idxmax_e
    if idxmin==idxmax: idxemax=idxmin
    else             : idxemax=np.argmax(e[idxmin:idxmax])+idxmin
    emax=e[idxemax]
    # Determine truth of bar presence
    barFlag=False
    if not barFlag and idxmin_p==idxmax_p: barFlag='constant_phi'
    if not barFlag and emax<ellmin       : barFlag='max_ellipticity'
    if not barFlag and idxmin_e==errval  : barFlag='rise_in_ellipticity'
    if not barFlag and idxmax_e==errval  : barFlag='fall_in_ellipticity'
    if not barFlag: barFlag=True
    # Generate output
    out={'barFlag':str(barFlag),
         'rmin':a[idxmin],
         'rmax':a[idxmax],
         'ramp':a[idxemax],
         'emax':e[idxemax],
         'phi' :pavg,
         'Am4' :Am4[idxemax]}
    out['amp']=out['emax']
    # Plot
    if mmlfiles.prep_write(plotfile,overwrite=owplot):
        import matplotlib.pyplot as plt
        # List things
        yvars=['pa','ellip','err','rms','Am4','der','hval']
        if runmeth=='dphidr': yvars.append('dpdr')
        limits={'a'    :None,
                'pa'   :(0,phiwid),
                'ellip':(0,1),
                'err'  :(-max(abs(edat0['err'])),max(abs(edat0['err']))),
                'rms'  :(-max(abs(edat0['rms'])),max(abs(edat0['rms']))),
                'Am4'  :(-max(abs(edat0['Am4'])),max(abs(edat0['Am4']))),
                'dpdr' :(-phiwid,phiwid),
                'der'  :(-max(abs(edat0['der'])),max(abs(edat0['der'])))}
        scales={'a'    :'log',
                'pa'   :'linear',
                'ellip':'linear',
                'err'  :'linear',
                'rms'  :'linear',
                'Am4'  :'linear',
                'dpdr' :'linear',
                'der'  :'linear'}
        lines={'a'    :[(a[idxmin  ],'k','dashed','rmin'),
                        (a[idxmax  ],'k','dashed','rmax'),
                        (a[idxmin_p],'g','dashed','rmin (phi)'),
                        (a[idxmax_p],'g','dashed','rmax (phi)'),
                        (a[idxmin_e],'r','dashed','rmin (ellip)'),
                        (a[idxmax_e],'r','dashed','rmax (ellip)'),
                        (a[idxemax ],'r','dotted','r of emax')],
               'pa'   :[(pavg,'g','dotted','avg phi'),
                        (pmin,'g','dashed','min phi'),
                        (pmax,'g','dashed','max phi')],
               'ellip':[(emax,'r','dotted','max ellip')],
               'err'  :[],
               'rms'  :[],
               'Am4'  :[],
               'dpdr' :[],
               'der'  :[]}
        # Prepare
        plt.clf()
        nplt=len(yvars)+1
        ncol=2 ; nrow=int(np.ceil(float(nplt)/float(ncol)))
        fig,axs=plt.subplots(nrow,ncol,squeeze=False,figsize=(1.4*4.*ncol,4.*nrow))
        if title: fig.suptitle(title)
        # Plot ellipses
        efit.plot(edat,imgdat=imgdat,color='m',cmap='bone',axs=axs[-1][-1],
                  linewidth=1,**kws)
        # Plot profiles
        xv='a'
        for i,yv in enumerate(yvars):
            iax=axs[i/ncol][i%ncol]
            if yv=='phi': ywid=phiwid
            else        : ywid=None
            plotprofile(iax,edat0[xv],edat0[yv],xlab=xv,ylab=yv,
                        vlines=lines.get(xv,[]),hlines=lines.get(yv,[]),
                        xlim=limits.get(xv,None),ylim=limits.get(yv,None),
                        xscl=scales.get(xv,None),yscl=scales.get(yv,None),ywid=ywid)
        # Save plot
        plt.savefig(plotfile,bbox_inches='tight')
        print '    '+plotfile
    return out
        
def plotprofile(ax,x,y,xlab=None,ylab=None,xlim=None,ylim=None,xscl=None,yscl=None,
                vlines=[],hlines=[],ywid=None,ls='-b'):
    """Simplifying function for plotting"""
    # Labels
    if xlab: ax.set_xlabel(xlab)
    if ylab: ax.set_ylabel(ylab)
    # Vertical and horizontal lines
    for v in vlines: ax.axvline(v[0],color=v[1],linestyle=v[2],label=v[3])
    for h in hlines: ax.axhline(h[0],color=h[1],linestyle=h[2],label=h[3])
    # Primary
    if ywid:
        brklist=[0]+list(np.where(np.abs(np.diff(y))>ywid/2.)[0])+[len(y)-1]
        for i in range(len(brklist)-1):
            islc=slice(brklist[i]+1,brklist[i+1]+1)
            ax.plot(x[islc],y[islc],ls)
    else: ax.plot(x,y,ls)
    # Limits
    if xlim: ax.set_xlim(xlim)
    if ylim: ax.set_ylim(ylim)
    # Scales
    if xscl: ax.set_xscale(xscl)
    if yscl: ax.set_yscale(yscl)
    # Done

