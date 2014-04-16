"""

bar
====

Functions for dealing with and manipulating bars in halos.


"""

from .. import filt, util, config
from .. import plot as nbplot
from . import mapping, angmom
import numpy as np
import math,os
import matplotlib.pyplot as plt
import matplotlib.cm as mpltcm

def rot(sim0,plot=False,plotfile=None,
        galaxy=None,family=None,cmeth=None,
        nbins=100):
    """Determine bar properties from the rotation rate"""
    # Center
    from ..family import get_family
    from ..galaxy import get_galaxy
    sim1=sim0[get_family(family)] if family else sim0
    sim2=sim1[get_galaxy(galaxy)] if galaxy else sim1
    sim=sim2
    if cmeth: angmom.faceon(sim,mode=cmeth)
    # Set properties
    inpar={}

def Am2(sim0,maxmeth='amp',plot=False,plotfile=None,
        galaxy=None,family=None,cmeth=None,
        rlim=(0.01,30.0),rscl='log',phierr=10.,nbins=100):
    """Determine bar properties from m=2 part of potential"""
    # Convert
    phierr*=np.pi/180.
    # Center
    from ..family import get_family
    from ..galaxy import get_galaxy
    sim1=sim0[get_family(family)] if family else sim0
    sim2=sim1[get_galaxy(galaxy)] if galaxy else sim1
    sim=sim2
    if cmeth: angmom.faceon(sim,mode=cmeth)
    # Set properties
    inpar={'xvar':'r'    ,'xlim':rlim          ,'xscl':rscl    ,
           'yvar':'modaz','ylim':(0.0   ,np.pi),'yscl':'lin'   ,
           'wvar':'Am2'  ,'wlim':(0.0001,1.0  ),'wscl':'symlog',
           'avar':'N','nbins':nbins,'units':['km s**-1','kpc','Msol']}
    # Get histogram
    h=mapping.hist2d(sim,inpar.pop('xvar'),inpar.pop('yvar'),**inpar)
    # Get max amplitude & phi dependence with radius
    rad=h['xs']
    phi=np.zeros(nbins,dtype=float)
    amp=np.zeros(nbins,dtype=float)
    for i in range(nbins):
        imax=np.argmax(h['ws'][i,:])
        phi[i]=h['ys'][imax]
        amp[i]=h['ws'][i,imax]
    # Get limits
    amplim=[0.0,0.1]
    radlim=[np.log10(inpar['xlim'][0]),np.log10(inpar['xlim'][1])]
    philim=inpar['ylim']
    # Get maxidx
    if   maxmeth=='amp': maxidx=np.argmax(amp) ; baridx=np.array([])
    elif maxmeth=='phi': 
        baridx=np.array([]) ; nphibar=0
        for iphi in h['ys']:
            iphires=abs(phi-iphi)/(philim[1]-philim[0])
            idxphi=np.logical_or(iphires<phierr,iphires>(1.-phierr))
            numphi=np.arange(len(idxphi))[idxphi]
            if len(numphi)==0: continue
            nphicur=1 ; idxphicur=[numphi[0]]
            for i in range(len(numphi)-1):
                if (numphi[i+1]-numphi[i])==1:
                    idxphicur.append(numphi[i+1])
                else:
                    if len(idxphicur)>nphibar:
                        nphibar=len(idxphicur)
                        baridx=np.array(idxphicur)
                    idxphicur=[numphi[i+1]]
        maxidx=baridx[np.argmax(amp[baridx])]
            # if np.sum(iphires)>nphibar:
            #     nphibar=np.sum(iphires)
            #     maxidx=np.arange(len(iphires))[iphires][np.argmax(amp[iphires])]
    else: raise Exception('Invalid maxmeth: {}'.format(maxmeth))
    radmax=rad[maxidx]
    # Get phi fit & res
    phifit=phi[maxidx]*np.ones(nbins,dtype=float)
    phires=(phi-phifit)/(philim[1]-philim[0])
    # Determine indices of bar
    if len(baridx)==0:
        barmax=maxidx ; barmin=maxidx
        while True:
            if barmax>=(nbins-1): break
            elif abs(phires[barmax+1])>phierr and abs(phires[barmax+1])<(1.-phierr): break
            else: barmax+=1
        while True:
            if barmin<=0: break
            elif abs(phires[barmin-1])>phierr and abs(phires[barmin-1])<(1.-phierr): break
            else: barmin-=1
        baridx=slice(barmin,barmax+1)
    radlimres=[rad[baridx].min(),rad[baridx].max()]
    philimres=[-phierr,phierr]
    # Set bar prop
    prop={'rmin':float(10**radlimres[0]),
          'rmax':float(10**radmax),#float(10**radlimres[1]),
          'amp' :float(amp[maxidx]),
          'phi' :float(phi[maxidx])}
    # Plot
    if plot:
        plt.clf()
        resclr='g' ; reslin='dashed'
        maxclr='r' ; maxlin='dotted'

        ax1=plt.subplot(2,2,1)
        ax1.plot(rad,phi,'b-',rad,phifit,'r-')
        ax1.vlines(radlimres,philim[0],philim[1],resclr,reslin)
        ax1.vlines(radmax,philim[0],philim[1],maxclr,maxlin)
        ax1.set_xlabel('log(r)')
        ax1.set_ylabel('phi max')
        ax1.set_xlim(radlim)
        ax1.set_ylim(philim)

        ax2=plt.subplot(2,2,2)
        ax2.plot(rad,amp)
        ax2.vlines(radlimres,amplim[0],amplim[1],resclr,reslin)
        ax2.vlines(radmax,amplim[0],amplim[1],maxclr,maxlin)
        ax2.set_xlabel('log(r)')
        ax2.set_ylabel('amp max')
        ax2.set_xlim(radlim)
        ax2.set_ylim(amplim)

        ax3=plt.subplot(2,2,3)
        ax3.plot(rad[baridx],phires[baridx])
        ax3.vlines(radmax,philimres[0],philimres[1],maxclr,maxlin)
        ax3.set_xlabel('log(r)')
        ax3.set_ylabel('phi residual')
        ax3.set_xlim(radlimres)
        ax3.set_ylim(philimres)

        ax4=plt.subplot(2,2,4)
        ax4.plot(rad[baridx],amp[baridx])
        ax4.vlines(radmax,amplim[0],amplim[1],maxclr,maxlin)
        ax4.set_xlabel('log(r)')
        ax4.set_ylabel('amp max')
        ax4.set_xlim(radlimres)
        ax4.set_ylim(amplim)
        
        if plotfile is None:
            plotfile=os.path.join(os.path.expanduser('~'),'Am2testplot.png')
        plt.tight_layout()
        plt.savefig(plotfile)
        print '    '+plotfile
    # Return properties
    return prop

def derAm2(sim0,maxmeth=None,plot=False,plotfile=None,
           galaxy=None,family=None,cmeth=None,ampthresh=0.001,
           rlim=(0.01,30.0),rscl='log',phierr=10.,nbins=100):
    """Determine bar properties from m=2 part of potential"""
    from mmlutils import mmlplot,mmlmath
    from scipy.interpolate import interp1d
    # Get histogram
    h=mapping.hist2d(sim0,'r','modaz',nbins=nbins,
                     xlim=rlim,xscl=rscl,
                     ylim=(0.,np.pi),yscl='lin',
                     wvar='Am2',wlim=(ampthresh/10,1.),wscl='symlog',
                     avar='N',units=['km s**-1','kpc','Msol'],
                     family=family,galaxy=galaxy,cmeth=cmeth)
    # Set constants
    phiwid=np.pi
    derwin=1
    # Pars input
    rad=h['xs']
    phi=h['ys']
    print rlim,h['xlim']
    radlim=h['xlim']
    philim=h['ylim']
    amplim=h['wlim']
    derlim=amplim
    phiderlim=(-phiwid,phiwid)
    phierr*=np.pi/180.
    binerr=int(np.ceil(phierr/(phi[1]-phi[0])))
    # Grid & Average
    PHI,RAD=np.meshgrid(phi,rad)
    avgwin=nbins/20
    amp=h['ws']
    amp=mmlmath.runavg2d(h['ws'],window=avgwin,weights=h['as'])
    # Pad phi & amplitude in phi so thate it wraps
    phipad=np.append(phi,phi[0:derwin]+phiwid)
    amppad=np.append(amp,amp[:,0:derwin],axis=1)
    # Get dA/dphi
    ampderphi=np.diff(amppad,n=derwin)/np.diff(phipad,n=derwin)
    # Get phi of max/min dA/dphi at each radii
    idxpos=np.ma.argmin(np.ma.masked_array(abs(amp),mask=ampderphi< 0),axis=1)
    idxneg=np.ma.argmin(np.ma.masked_array(abs(amp),mask=ampderphi>=0),axis=1)
    phipos=phi[idxpos] ; fphipos=interp1d(rad,phipos,kind='cubic')
    phineg=phi[idxneg] ; fphineg=interp1d(rad,phineg,kind='cubic')
    # Derivatives
    radder1,phidr1_pos,radder2,phidr2_pos=mmlmath.wrapder(rad,phipos,wrpwin=phiwid,retder2=True)
    radder1,phidr1_neg,radder2,phidr2_neg=mmlmath.wrapder(rad,phineg,wrpwin=phiwid,retder2=True)
    phidr1=mmlmath.runavg1d((phidr1_pos+phidr1_neg)/2.,window=5)
    phidr2=np.diff(phidr1)
    # Find run where phi derivative does not change
    derout,idxderbeg,idxderend=mmlmath.findrun(phidr1,phierr)
    radmin=radder1[idxderbeg]
    radmax=radder1[idxderend]
    # derout,idxderbeg,idxderend=mmlmath.findrun(phipos,phierr)
    # radmin=rad[idxderbeg]
    # radmax=rad[idxderend]
    # Get location of maximum in bar
    # phineg is at leading edge of bar
    # phipos is at trailing edge of bar
    idxradmin=np.ma.argmin(np.ma.masked_less(rad,radmin))
    idxradmax=np.ma.argmax(np.ma.masked_greater(rad,radmax))
    idxbar=np.ma.argmax(np.ma.masked_array(amp,mask=RAD>radmax))
    idxrad,idxphi=np.unravel_index(idxbar,amp.shape)
    radamp=rad[idxrad]
    phiamp=phi[idxphi]
    baramp=amp[idxrad,idxphi]
    # Middle
    idxradavg=(idxradmax+idxradmin)/2
    # Bar width
    idxphimin=idxpos[idxradavg] # Beginning of bar
    idxphimax=idxneg[idxradavg] # End of bar
    idxwid=idxphimax-idxphimin
    if idxphimax<idxphimin: idxwid+=nbins
    idxphiavg=idxphimin+idxwid/2
    if idxphiavg>=nbins: idxphiavg-=nbins
    # Define values
    radmin=rad[idxradmin]
    radavg=rad[idxradavg]
    phimin=phi[idxphimin]
    phimax=phi[idxphimax]
    phiavg=phi[idxphiavg]
    # Set bar prop
    prop={'rmin':float(10**radmin),
          'ramp':float(10**radamp),
          'rmax':float(10**radmax),
          'amp' :float(baramp),
          'phi' :float(phiavg)}
    import pprint
    pprint.pprint(prop)
    # Plot
    if plot:
        plt.clf()
        ncol=2 ; nrow=1 ; iplt=-1
        fig,axs=plt.subplots(nrow,ncol,squeeze=False,figsize=(1.4*4.*ncol,4.*nrow))
        fig.suptitle(nbplot.generic.val2str(sim0.properties['time'],label='t',units='Gyr',latex=False))
#                     horizontalalignment='left',x=0.45)
        cbs=[] ; nlvl=50 ; ntic=9
        resclr='k' ; reslin='dashed'
        maxclr='k' ; maxlin='solid'
        cmap=plt.get_cmap('jet')
        ampnorm=mmlplot.SymLogNorm(amplim[0],vmin=-amplim[1],vmax=amplim[1],linscale=0.001)
        dernorm=mmlplot.SymLogNorm(derlim[0],vmin=-derlim[1],vmax=derlim[1],linscale=0.001)
        ticfmt='% .3f'

        # Plot Am2 vs. phi & logr
        ix=rad         ; iy=phi        ; iz=amp
        ixlab='log(R)' ; iylab='phi'   ; iclab='Am2'
        ixlim=radlim   ; iylim=philim  ; iclim=None
        inorm=ampnorm
        if inorm: ilvl=inorm.inverse(np.linspace(0.,1.,nlvl,endpoint=True))
        else    : ilvl=None
        if inorm: itic=inorm.inverse(np.linspace(0.,1.,ntic,endpoint=True))
        else    : itic=None
        iplt+=1 ; irow=iplt/ncol ; icol=iplt%ncol
        iax=axs[irow,icol]
        iax.set_xlabel(ixlab)
        iax.set_ylabel(iylab)
        icn=iax.contour(ix,iy,np.transpose(iz),ilvl,cmap=cmap,norm=inorm,extend='both')
        icf=iax.contourf(ix,iy,np.transpose(iz),ilvl,cmap=cmap,norm=inorm,extend='both')
        iax.vlines(radavg,iylim[0],iylim[1],maxclr,maxlin)
        iax.vlines([radmin,radmax],iylim[0],iylim[1],resclr,reslin)
        iax.hlines(phiavg,ixlim[0],ixlim[1],maxclr,maxlin)
        iax.hlines([phimin,phimax],ixlim[0],ixlim[1],resclr,reslin)

        # iax.scatter(rad,phipos,c='k',marker='^')
        # iax.scatter(rad,phineg,c='k',marker='v')
        iax.scatter(radamp,phiamp,c='k',marker='x')

        def plot_nobreak(axobj,x,y,ywid):
            brklist=[0]+list(np.where(np.abs(np.diff(y))>ywid/2.)[0])+[len(y)-1]
            for i in range(len(brklist)-1):
                islc=slice(brklist[i]+1,brklist[i+1]+1)
                axobj.plot(x[islc],y[islc],'-k')

        plot_nobreak(iax,rad,fphipos(rad),phiwid)
        plot_nobreak(iax,rad,fphineg(rad),phiwid)

        icb=plt.colorbar(icf,ax=iax,format=ticfmt,ticks=itic)
        icb.set_label(iclab)
        cbs.append(icb)
        iax.set_xlim(ixlim)
        iax.set_ylim(iylim)
        iax.set_aspect((ixlim[1]-ixlim[0])/(iylim[1]-iylim[0]))

        # Plot derivative vs. phi & logr
        # ix=rad         ; iy=phi       ; iz=ampderphi
        # ixlab='log(R)' ; iylab='phi'  ; iclab='dAm2/dphi'
        # ixlim=radlim   ; iylim=philim ; iclim=None
        # inorm=dernorm
        # if inorm: ilvl=inorm.inverse(np.linspace(0.,1.,nlvl,endpoint=True))
        # else    : ilvl=None
        # if inorm: itic=inorm.inverse(np.linspace(0.,1.,ntic,endpoint=True))
        # else    : itic=None
        # iplt+=1 ; irow=iplt/ncol ; icol=iplt%ncol
        # iax=axs[irow,icol]
        # icn=iax.contour(ix,iy,np.transpose(iz),ilvl,cmap=cmap,norm=inorm,extend='both')
        # icf=iax.contourf(ix,iy,np.transpose(iz),ilvl,cmap=cmap,norm=inorm,extend='both')
        # iax.scatter(rad,phipos,c='k',marker='^')
        # iax.scatter(rad,phineg,c='k',marker='v')
        # iax.set_xlabel(ixlab)
        # iax.set_ylabel(iylab)
        # iax.vlines([radmin,radmax],iylim[0],iylim[1],maxclr,maxlin)
        # iax.hlines(barphi,ixlim[0],ixlim[1],maxclr,maxlin)
        # iax.hlines([phineg[idxradavg],phipos[idxradavg]],ixlim[0],ixlim[1],resclr,reslin)
        # icb=plt.colorbar(icf,ax=iax,format=ticfmt,ticks=itic)
        # icb.set_label(iclab)
        # cbs.append(icb)
        # iax.set_aspect((ixlim[1]-ixlim[0])/(iylim[1]-iylim[0]))
        # iax.set_xlim(ixlim)
        # iax.set_ylim(iylim)

        # Plot derivative of phi
        ixlab='log(R)' ; iylab='dphi/dlogr'
        ixlim=radlim   ; iylim=phiderlim
        iplt+=1 ; irow=iplt/ncol ; icol=iplt%ncol
        iax=axs[irow,icol]
        iax.set_xlabel(ixlab)
        iax.set_ylabel(iylab)
        iax.vlines(radavg,iylim[0],iylim[1],maxclr,maxlin)
        iax.vlines([radmin,radmax],iylim[0],iylim[1],resclr,reslin)

        iax.plot(radder1,phidr1_pos,'r-')
        iax.plot(radder1,phidr1_neg,'b-')
        iax.plot(radder1,phidr1,'k-')

        # iax.plot(radder2,phidr2_pos,'r--')
        # iax.plot(radder2,phidr2_neg,'b--')
        iax.plot(radder2,phidr2,'k--')

        iax.set_xlim(ixlim)
        iax.set_ylim(iylim)
        iax.set_aspect((ixlim[1]-ixlim[0])/(iylim[1]-iylim[0]))
        
        # Save
        if plotfile is None:
            plotfile=os.path.join(os.path.expanduser('~'),'Am2testplot.png')
            plt.tight_layout()
        plt.savefig(plotfile,bbox_inches='tight')#,bbox_extra_artists=cbs)
        print '    '+plotfile
    # Return properties
    return prop


def ellipse(sim0,plot=False,plotfile=None,
            galaxy=None,family=None,cmeth=None,
            rlim=(0.01,30.0),rscl='log',nbins=100,nellip=100,
            ellmin=0.25,phierr=10.,ellerr=0.1,histmeth='numpy',plotmeth='contourf'):
    """Determine bar properties from ellipse fitting"""
    # Convert
    phierr*=np.pi/180.
    # Center
    from mmlutils import mmlplot,mmlmath
    from ..family import get_family
    from ..galaxy import get_galaxy
    from scipy.interpolate import interp1d
    sim1=sim0[get_family(family)] if family else sim0
    sim2=sim1[get_galaxy(galaxy)] if galaxy else sim1
    sim=sim2
    if cmeth: angmom.faceon(sim,mode=cmeth)
    # Set properties
    xylim=(-rlim[1],rlim[1])
    mlim=(1.0e4,1.0e10)
    radlim=rlim
    philim=(0.,np.pi)
    elllim=(0.,1.)
    phiwid=philim[1]-philim[0]
    inpar={'xvar':'x'    ,'xlim':xylim,'xscl':'lin',
           'yvar':'y'    ,'ylim':xylim,'yscl':'lin',
           'wvar':'mass' ,'wlim':mlim ,'wscl':'log',
           'nbins':nbins,'units':['km s**-1','kpc','Msol'],
           'rethist':True,'histmeth':histmeth,'plotmeth':plotmeth,
           'Na':nellip,'amin':rlim[0],'amax':rlim[1],'ascl':rscl}
    # Get histogram
    edat,hdat=mapping.ellipse(sim,inpar.pop('xvar'),inpar.pop('yvar'),**inpar)
    # Get lists
    elst={k:[] for k in edat[0].keys()}
    for e in edat:
        for k in elst.keys():
            elst[k].append(e[k])
    idxsort=np.array(sorted(range(len(elst['a'])),key=lambda k: elst['a'][k]))
    for k in elst.keys():
        elst[k]=np.array(elst[k])[idxsort]
    rad=elst['a']
    wid=elst['b']
    ell=elst['ellip']
    phi=elst['pa']
    Am4=elst['A4']
    # Get starting index based on increasing ellipticity
    idx0=(list(np.where((ell-ell[0])>ellerr)[0])+[0])[0]
    #idx0=(list(np.where(np.diff(ell)>=0)[0])+[0])[0]
    # Get bar length by constant phi
    phiavg,idxphimin,idxphimax=mmlmath.findrun(phi,phierr,wrap=phiwid)
    # phidiff=phi-phi[idx0]
    # idxwrp=(abs(phidiff)>phiwid/2.)
    # phidiff[idxwrp]-=np.sign(phidiff[idxwrp])*phiwid
    # idxphierr=list(np.where(abs(phidiff[idx0:])>phierr)[0])+[0]
    # idxphimax=max(0,idxphierr[0]-1+idx0)
    # Get bar length by ellipticity
    idxellmax=np.argmax(ell[idxphimin:idxphimax])+idxphimin
    elldiff=ell[idxellmax:]-ell[idxellmax]
    idxellerr=list(np.where(abs(elldiff)>ellerr)[0])+[0]
    idxradmin=(list(np.where(np.diff(ell[idxphimin:idxphimax])>0)[0])+[0])[0]+idxphimin
    idxradmax=max(0,idxellerr[0]-1+idxellmax)
    # Set bar prop
    prop={'rmin':rad[idxradmin],
          'ramp':rad[idxellmax],
          'rmax':rad[idxradmax],
          'amp' :ell[idxellmax],
          'Am4' :Am4[idxradmax],
          'phi' :phiavg}#phi[idxphimax]}
    # Plot
    if plot:
        plt.clf()
        ncol=2 ; nrow=1 ; iplt=-1
        fig,axs=plt.subplots(nrow,ncol,squeeze=False,figsize=(1.4*4.*ncol,4.*nrow))
        fig.suptitle(nbplot.generic.val2str(sim0.properties['time'],label='t',units='Gyr',latex=False))
#                     horizontalalignment='left',x=0.45)
        cbs=[] ; nlvl=50 ; ntic=9
        phiclr='b' ; philin='dashed'
        ellclr='r' ; elllin='dotted'
        radclr='k' ; radlin='solid'
        def plot_nobreak(axobj,x,y,ywid):
            brklist=[0]+list(np.where(np.abs(np.diff(y))>ywid/2.)[0])+[len(y)-1]
            for i in range(len(brklist)-1):
                islc=slice(brklist[i]+1,brklist[i+1]+1)
                axobj.plot(x[islc],y[islc],'-k')

        # Plot phi vs. R
        ix=rad       ; iy=phi
        ixlab='R'    ; iylab='phi'
        ixlim=radlim ; iylim=philim
        ixscl='log'  ; iyscl='linear'
        iplt+=1 ; irow=iplt/ncol ; icol=iplt%ncol
        iax=axs[irow,icol]
        iax.set_xlabel(ixlab)
        iax.set_ylabel(iylab)
        iax.axvline(rad[idxphimax],color=phiclr,linestyle=philin)
        iax.axvline(rad[idxellmax],color=ellclr,linestyle=elllin)
        iax.axvline(rad[idxradmax],color=radclr,linestyle=radlin)
        iax.axhline(phi[idxphimax],color=phiclr,linestyle=philin)

        plot_nobreak(iax,rad,phi,phiwid)

        iax.set_xlim(ixlim)
        iax.set_ylim(iylim)
        iax.set_xscale(ixscl)
        iax.set_yscale(iyscl)
        #iax.set_aspect((ixlim[1]-ixlim[0])/(iylim[1]-iylim[0]))

        # Plot ellipticity vs. rad
        ix=rad        ; iy=ell
        ixlab='R'     ; iylab='ellipticity'
        ixlim=radlim  ; iylim=elllim
        ixscl='log'   ; iyscl='linear'
        iplt+=1 ; irow=iplt/ncol ; icol=iplt%ncol
        iax=axs[irow,icol]
        iax.set_xlabel(ixlab)
        iax.set_ylabel(iylab)
        iax.axvline(rad[idxphimax],color=phiclr,linestyle=philin)
        iax.axvline(rad[idxellmax],color=ellclr,linestyle=elllin)
        iax.axvline(rad[idxradmax],color=radclr,linestyle=radlin)
        iax.axhline(ell[idxellmax],color=ellclr,linestyle=elllin)

        iax.plot(rad,ell,'k-')

        iax.set_xlim(ixlim)
        iax.set_ylim(iylim)
        iax.set_xscale(ixscl)
        iax.set_yscale(iyscl)
        #iax.set_aspect((ixlim[1]-ixlim[0])/(iylim[1]-iylim[0]))
        
        # Save
        if plotfile is None:
            plotfile=os.path.join(os.path.expanduser('~'),'Am2testplot.png')
            plt.tight_layout()
        plt.savefig(plotfile,bbox_inches='tight')#,bbox_extra_artists=cbs)
        print '    '+plotfile
    # Return properties
    return prop
