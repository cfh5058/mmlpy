#!/usr/bin/python
import sys,os,shutil,glob,copy,pprint,scipy,math,re,pickle
import matplotlib.pyplot as plt
import matplotlib
import numpy as np
from mmlutils import *


def snapshot(sim,galaxy='galaxy1',family='disk',fext='000',method='image',
             plotfile=None,plotext='eps',colorbar=False,owhist=False,
             annotate=False,xylabels=False):
    """
    Plot snapshot 
    """
    # Plotting parameters
    textsize=18
    ticksize=15
    # Plot file
    if plotfile==None:
        plotfile=os.path.expanduser('~/{}_{}_{}_{}_{}.{}'.format(sim['runtag'],galaxy,family,method,fext,plotext))
    # Set up default plotting commands
    plotkw0=dict(textsize=textsize,ticksize=ticksize,
                 colorbar=colorbar,annotate=annotate)
    if xylabels: plotkw0.update(xlabel='X [kpc]',ylabel='Y [kpc]')
    else       : plotkw0.update(xlabel=None,ylabel=None)
    inpar=dict(histmeth='numpy',plotmeth='contourf',
               family=family,galaxy=galaxy,
               x='x',xscl='lin',xlim=(0.,0.),
               y='y',yscl='lin',ylim=(0.,0.),
               w='mass',wscl='log',wlim=(0.,0.),
               a='',cmeth='pot')
    inkw0=dict(askuser=True,overwrite=True,owhist=owhist,make_plot=False)
    # Image plotting commands
    if method=='image':
        plotkw=dict(clabel='log(Mass) [M$_{\odot}$]',
                    arwcolor=[0,1,0],arwsize=400,arwtail=0.1,arwwidth=3.0,
                    cbtickloc=[5.0,5.5,6.0,6.5,7.0,7.5,8.0,8.5],
                    cbtickfmt='% .1f',**plotkw0)
        inkw=dict(pmeth='map2d',inpar=inpar,**inkw0)
    # Ellipse plotting commands
    elif method=='ellipse':
        plotkw=dict(clabel='log(Mass) [M$_{\odot}$]',
                    ellipcolor=[1,0,1],ellipalpha=0.5,ellipwidth=3,Na=10,
                    arwcolor=[0,1,0],arwsize=400,arwtail=0.1,arwwidth=3.0,
                    cbtickloc=[5.0,5.5,6.0,6.5,7.0,7.5,8.0,8.5],
                    cbtickfmt='% .1f',**plotkw0)
        inkw=dict(pmeth='ellipse',inpar=inpar,**inkw0)
    # SCF plotting commands
    elif method=='scf':
        plotkw=dict(clabel='Relative m=2 Amplitude',
                    cmap=plt.get_cmap('jet'),fgcolor=[0,0,0],
                    cbtickloc=[-1.0,-0.1,-0.01,-0.001,0,0.001,0.01,0.1,1.0],
                    cbtickfmt="% .3f",**plotkw0)
        inpar.update(w='Am2',wscl='symlog',wlim=(0.0001,1.0),a='N',cmeth='pot')
        inkw=dict(pmeth='map2d',inpar=inpar,**inkw0)
    # Error
    else: raise Exception('Invalid plotting method: {}'.format(method))
    # Plotting things
    inkw.update(plotkw=plotkw,loadkw={'fext':fext,'ftype':'gadget','ptype':family,'pgal':int(galaxy[-1])})
    if method=='scf': inkw['loadkw']['ftype']='scf'
    plt.clf()
    # Get plot
    out=sim.plotsnap(**inkw)
    # Turn of ticks & save file
    plt.setp(out.get_xticklabels(),visible=False)
    plt.setp(out.get_yticklabels(),visible=False)
    plt.savefig(plotfile,bbox_inches='tight')
    print plt.gcf().get_size_inches()
    print '    '+plotfile

def barevolution(sim,galaxy='galaxy1',family='disk',methods=['scfAm2','ellipse'],plotext='eps',
                 colors=None,styles=None,dashes=None,
                 hatchFlag=False       ,hatchColors=None      ,hatchStyles=None      ,
                 shadeFlag=True        ,shadeColors=None      ,shadeAlpha=0.2        ,
                 hatchOverlapFlag=False,hatchOverlapColor=None,hatchOverlapStyle=None,
                 shadeOverlapFlag=False,shadeOverlapColor=None,shadeOverlapAlpha=0.2 ,
                 validFlag=False       ,validColors=None      ,
                 validOverlapFlag=False,validOverlapColor='m' ,
                 tperi=True            ,tperiColor=None       ,tperiStyle='-'      ,
                 tlim=None,tlab='Time (Gyr)',terr=0.11,tavg=2.,tdec=2,rmin=1.,
                 plotfile=None,overwrite=False,owdata=False,owhist=False,sameAxes=True,axs=None,
                 printFlag=True,pbuff='    '):
    """Plot bar properties as a function of time"""
    rgbM=(255,187,255)
    rgbR=(255,181,197)
    rgbB=(202,225,255)
    rgbG=(  0,238,  0)
    # Check number of methods
    if sameAxes:
        if len(methods)>2: raise Exception('Cannot plot more than one variable on the same axes.')
        else: Nax=1
    else:
        Nax=len(methods)
    # Plotting parameters
    Alab={'derAm2' :'Max Am2',
          'ellipse':'Max Ellipticity',
          'scfAm2' :'Max Am2'}
    Alim={'derAm2' :[0.,0.15],
          'ellipse':[0.,1.],
          'scfAm2' :[0.,0.15]}
    # plthgttot=6.
    # labhgt=0.5
    pltwid=18.
    plthgt=6.
    linewidth=4
    textsize=18
    ticksize=15
    labelx_R=1.04
    labelx_L=-0.06
    matplotlib.rcParams['font.size']=ticksize
    ncol=1 ; nrow=int(np.ceil(float(Nax)/float(ncol)))
    # plthgt=(plthgttot-labhgt)/nrow
    # Line properties
    if tperiColor is None: tperiColor=tuple([float(c)/255. for c in rgbG])
    if not isinstance(colors,list): 
        if sameAxes: colors=['b','r']
        else       : colors=len(methods)*['b']
    if not isinstance(styles,list):
        if sameAxes: styles=['-','--']
        else       : styles=len(methods)*['-']
    if not isinstance(dashes,list):
        if sameAxes: dashes=[(20,5) if c=='--' else (None,None) for c in styles]
        else       : dashes=len(methods)*[(None,None)]
    # Shading
    if not isinstance(shadeColors,list):
        if plotext=='eps': shadeColors=[tuple([float(c)/255. for c in rgbB]),
                                        tuple([float(c)/255. for c in rgbR])]
        else             : shadeColors=[c[0] for c in colors]
    if shadeOverlapColor is None:
        if plotext=='eps': shadeOverlapColor=tuple([float(c)/255. for c in rgbM])
        else             : shadeOverlapColor='m'
    # Valid markers
    if not isinstance(validColors,list):
        validColors=[c[0] for c in colors]
    # Hatch
    hatchs=['\\','/','|','-']
    if not isinstance(hatchStyles,list):
        hatchStyles=hatchs[:len(methods)]
    if hatchOverlapStyle is None:
        hatchOverlapStyle='x'
    if not isinstance(hatchColors,list):
        hatchColors=[tuple([float(c)/255. for c in rgbB]),
                     tuple([float(c)/255. for c in rgbR])]
        hatchColors=[c[0] for c in colors]
    if hatchOverlapColor is None:
        hatchOverlapColor='m'
    # Plot file
    if plotfile==None:
        plotfile=os.path.expanduser('~/{}_{}_{}_barevol'.format(sim['runtag'],galaxy,family))
        if len(methods)==1: plotfile+='_'+methods[0]
        plotfile+='.{}'.format(plotext)
    # Get tperi
    if tperi==True:
        parPer,parObs,parTim=sim.orbit_param(cmeth='pot',ftype='gadget',verbose=False,askuser=False)
        tperi=parTim['all']
    # Gather data for each method
    prop={} ; shade={} ; valid={} ; vprint={}
    for idx,m in enumerate(methods):
        # File extension
        if m=='ellipse': fext='*5'
        else           : fext='*5_{}{}'.format(family,galaxy[-1])
        # Load properties (and plot)
        p=sim.stdplots('barprop_{}'.format(m),family,galaxy,fext=fext,
                       overwrite=overwrite,owdata=owdata,owanim=True,owhist=owhist)
        if len(p['time'])==1 and overwrite==False:
            if mmlio.yorn('Only one entry for method {}. Overwrite?'.format(m)):
                p=sim.stdplots('barprop_{}'.format(m),family,galaxy,fext=fext,
                               overwrite=True,owdata=owdata,owanim=True,owhist=owhist)
        # Rotation rate
        p['tder'],p['dphidt']=mmlmath.wrapder(p['time'],p['phi'],wrpwin=np.pi,retder2=False)
        p['tder'  ]=np.append(0,p['tder'  ])
        p['dphidt']=np.append(0,p['dphidt'])
        print p['dphidt']
        # Shading for method
        valid[m]=np.around(p['time'][np.logical_and(p['barFlag']=='True',p['ramp']>rmin)],decimals=tdec)
        shade[m]=mmlmath.findruns(valid[m],terr,val='rel',noOverlap=True,retval=True)
        if 'all' in valid: valid['all']=np.intersect1d(valid['all'],valid[m])
        else             : valid['all']=valid[m]
        # Measurements for method
        vprint[m]={'idxmax':np.array([np.argmax(p['amp'])]),
                   'idxend':np.arange(len(p['time']))[(p['time']>=(p['time'][-1]-tavg))]}
        # Assign properties
        prop[m]=p
    # Shading for all methods
    if 'all' in valid: shade['all']=mmlmath.findruns(valid['all'],terr,val='rel',noOverlap=True,retval=True)
    # Time limits 
    if not tlim:
        tlim=[0.,0.]
        for m,p in prop.iteritems():
            tlim[0]=min(tlim[0],p['time'].min())
            tlim[1]=max(tlim[1],p['time'].max())
    # Set up axes
    if axs is None:
        plt.clf()
        fig,axs=plt.subplots(nrow,ncol,sharex=True,squeeze=False,
                             figsize=(ncol*pltwid,nrow*plthgt))
        save=True
    else: save=False
    # Loop over methods adding shading
    if shadeFlag or shadeOverlapFlag:
        for idx,m in enumerate(methods):
            # Get correct parameters for iteration
            if len(p['time'])==1: 
                mmlio.verbose('Only one entry for method {}. Quitting.'.format(m))
                return
            irow=idx/ncol
            icol=idx%ncol
            if sameAxes: iax=axs[0,0]
            else       : iax=axs[irow,icol]        
            # Shade
            if shadeFlag:
                for ix,iy in shade[m]: iax.axvspan(ix,iy,color=shadeColors[idx],alpha=shadeAlpha)
            if shadeOverlapFlag:# and (not sameAxes or idx==0):
                for ix,iy in shade['all']: iax.axvspan(ix,iy,color=shadeOverlapColor,alpha=shadeOverlapAlpha)
    # Loop over methods adding hatching
    if hatchFlag or hatchOverlapFlag:
        for idx,m in enumerate(methods):
            # Get correct parameters for iteration
            if len(p['time'])==1: 
                mmlio.verbose('Only one entry for method {}. Quitting.'.format(m))
                return
            irow=idx/ncol
            icol=idx%ncol
            if sameAxes: iax=axs[0,0]
            else       : iax=axs[irow,icol]        
            # Hatch
            if hatchFlag:
                for ix,iy in shade[m]: iax.axvspan(ix,iy,color=hatchColors[idx],hatch=hatchStyles[idx],fill=False)
            if hatchOverlapFlag:
                for ix,iy in shade['all']: iax.axvspan(ix,iy,color=hatchOverlapColor,hatch=hatchOverlapStyle,fill=False)
    # Loop over methods, plotting them
    for idx,m in enumerate(methods):
        # Get correct parameters for iteration
        p=prop[m]
        if len(p['time'])==1: 
            mmlio.verbose('Only one entry for method {}. Quitting.'.format(m))
            return
        irow=idx/ncol
        icol=idx%ncol
        if sameAxes: 
            if idx==0: iax=axs[0,0]
            else     : iax=axs[0,0].twinx()
        else       : iax=axs[irow,icol]        
        # Valid times
        if validFlag:
            for ix in valid[m]: iax.axvline(ix,color=validColors[idx])
        if validOverlapFlag and (not sameAxes or idx==0):
            for ix in valid['all']: iax.axvline(ix,color=validOverlapColor)
        # Plot tperi
        if tperi and (not sameAxes or idx==0): 
            iax.axvline(tperi,linestyle=tperiStyle,color=tperiColor,lw=2.*linewidth)
        # Plot amplitude vs. time
        iax.plot(p['time'],p['amp'],c=colors[idx],ls=styles[idx],lw=linewidth,
                 dashes=dashes[idx])
        # Limits
        iax.set_xlim(tlim)
        iax.set_ylim(Alim[m])
        # Labels & ticks
        if not sameAxes and idx<(len(methods)-1):
            plt.setp(iax.get_xticklabels(),visible=False)
        else:
            iax.set_xlabel(tlab,size=textsize)
        iax.set_ylabel(Alab[m])
        if sameAxes and idx==0:
            iax.yaxis.set_label_coords(labelx_L, 0.5)
        else:
            iax.yaxis.set_label_coords(labelx_R, 0.5)
        iax.yaxis.tick_right()
        iax.yaxis.set_label_position("right")
        # Print
        if printFlag:
            parPer,parObs,parTim=sim.orbit_param(cmeth='pot',ftype='gadget',verbose=False,askuser=False)
            print 1*pbuff+'peri  ({}) @ {:3.2f} Gyr, r = {:5.2f} kpc'.format(parObs['fext_peri' ],parObs['tperi' ],parObs['allrperi' ])
            print 1*pbuff+'final ({}) @ {:3.2f} Gyr, r = {:5.2f} kpc'.format(parObs['fext_final'],parObs['tfinal'],parObs['allrfinal'])
            print 1*pbuff+'{}: {}, {}, {}'.format(sim['runtag'],galaxy.capitalize(),family.capitalize(),m)
            for m2 in methods:
                for tloc in ['max','end']:
                    tx=vprint[m2]['idx'+tloc]
                    print 2*pbuff+'@ {} {} (t = {:3.2f} Gyr, {})'.format(tloc.capitalize(),m2,p['time'][tx].max(),p['fext'][tx][-1])
                    print 3*pbuff+'barFlag = {}'.format(p['barFlag'][tx][-1])
                    print 3*pbuff+'amp     = {:3.2f}'.format(np.mean(p['amp'][tx]))
                    print 3*pbuff+'phi     = {:3.2f} rad'.format(np.mean(p['phi'][tx]))
                    print 3*pbuff+'ramp    = {:3.1f} kpc'.format(np.mean(p['ramp'][tx]))
                    print 3*pbuff+'rmax    = {:3.1f} kpc'.format(np.mean(p['rmax'][tx]))
                    print 3*pbuff+'dphidt  = {:3.1f} rad/Gyr'.format(np.mean(p['dphidt'][tx]))
                    if m=='ellipse':
                        print 3*pbuff+'Am4     = {}'.format(np.mean(p['Am4'][tx]))
    # Save
    if save:
        plt.tight_layout()
        plt.savefig(plotfile,bbox_inches='tight')
        print '    '+plotfile
