"""

generic
=======

Flexible and general plotting functions

"""

import numpy as np
import os,pickle
from . import profile, angmom, halo
from .. import config
from ..array import SimArray
from ..units import NoUnit
from ..snapshot import SimSnap

def ComputeMap(sim,xyvar,wvar=0,*arg,**kw):
    """
    Return map
    """
    # Get variable

def ellipse(sim,xvar,yvar,rethist=False,**kwargs):
    """
    Creates an image and fits ellipses to it
    """
    from mmlastro import ellipse
    from mmlutils import mmlplot
    # Check for existing file & return contents if not overwrite
    filename=kwargs.pop('filename',None)
    overwrite=kwargs.get('overwrite',True)
    if filename:
        if not overwrite and os.path.isfile(filename):
            out=ellipse.read(filename,retlist=True)
            return out
    # Get image data
    if 'imgdata' in kwargs: img=kwargs['imgdata']
    else:
        kwargs['filename']=kwargs.get('imgfile',None)
        img=hist2d(sim,xvar,yvar,**kwargs)
    # Get ellipse input
    imgdat,imgdim=ellipse.parshist(img)
    kwargs.update(imgdat=imgdat,imgdim=imgdim)
    tdict={'nsemimajor'  :'Na',
           'isemimajor'  :'amin',
           'fsemimajor'  :'amax',
           'semimajorscl':'ascl'}
    for k in tdict: 
        if k in kwargs: 
            raise Exception('Outdated argument: {}'.format(k))
            # kwargs[tdict[k]]=kwargs.pop(k)
    # Get ellipse data
    out=ellipse.fitmult(filename=filename,overwrite=overwrite,**kwargs)
    # Return data
    if rethist: return out,img
    else      : return out,imgdim

def hist2d(sim0,xvar,yvar,wvar=None,avar=None,method='numpy',**kwargs):
    """
    Returns a dictionary with data for creating a 2d map
    """
    # Identify data type
    pNbody = False if issubclass(sim0.__class__,SimSnap) else True
    if method=='pNbody': pNbody=True
    # Check for existing file & return contents if not overwrite
    filename=kwargs.get('filename',None)
    if filename:
        overwrite=kwargs.get('overwrite',False)
        if not overwrite and os.path.isfile(filename):
            f=open(filename,'rb')
            out=pickle.load(f)
            f.close()
            return out
    # Pars input
    units=kwargs.get('units',['km s**-1','kpc','Msol'])
    xlim=kwargs.get('xlim',(0.,0.)) ; xscl=kwargs.get('xscl','lin')
    ylim=kwargs.get('ylim',(0.,0.)) ; yscl=kwargs.get('yscl','lin')
    wlim=kwargs.get('wlim',(0.,0.)) ; wscl=kwargs.get('wscl','lin')
    if xlim is None: xlim=(0.,0.)
    if ylim is None: ylim=(0.,0.)
    if wlim is None: wlim=(0.,0.)
    gridsize=kwargs.get('gridsize',(100,100))
    nbins=kwargs.get('nbins',None)
    if nbins is not None: gridsize=(nbins,nbins)
    galaxy=kwargs.get('galaxy',None)
    family=kwargs.get('family',None)
    cmeth=kwargs.get('cmeth',None)
    # Select correct particles & orientation
    from . import langmm
    sim=langmm.align(sim0,galaxy=galaxy,family=family,cmeth=cmeth,pNbody=pNbody)
    # Data & units
    if pNbody:
        xo=SimArray(sim.get_var(xvar),sim.get_units(xvar))
        yo=SimArray(sim.get_var(yvar),sim.get_units(yvar))
        wo=SimArray(sim.get_var(wvar),sim.get_units(wvar))
        ao=SimArray(sim.get_var(avar),sim.get_units(avar))
    else:
        xo=sim[xvar]
        yo=sim[yvar]
        wo=sim[wvar] if wvar else SimArray(np.ones(len(sim),dtype=float),units=NoUnit())
        ao=sim[avar] if avar else SimArray(np.ones(len(sim),dtype=float),units=NoUnit())
    # Units
    xo=xo.in_unitsys(units)
    yo=yo.in_unitsys(units)
    #wo=wo.in_unitsys(units)
    #ao=ao.in_unitsys(units)
    # Limits
    if xlim[0]==xlim[1]:
        if xscl=='log': xlim=[np.log10(np.min(xo)),np.log10(np.max(xo))]
        else          : xlim=[np.min(xo),np.max(xo)]
    else:
        if hasattr(xlim,'units'): xlim=xlim.in_units(xo.units)
        else                    : xlim=SimArray(xlim,xo.units)
        if xscl=='log': xlim=[np.log10(xlim[0]),np.log10(xlim[1])]
        else          : xlim=[xlim[0],xlim[1]]
    if ylim[0]==ylim[1]:
        if yscl=='log': ylim=[np.log10(np.min(yo)),np.log10(np.max(yo))]
        else          : ylim=[np.min(yo),np.max(yo)]
    else:
        if hasattr(ylim,'units'): ylim=ylim.in_units(yo.units)
        else                    : ylim=SimArray(ylim,yo.units)
        if yscl=='log': ylim=[np.log10(ylim[0]),np.log10(ylim[1])]
        else          : ylim=[ylim[0],ylim[1]]
    x = np.ma.log10(np.ma.masked_equal(xo,0.)) if xscl=='log' else xo
    y = np.ma.log10(np.ma.masked_equal(yo,0.)) if yscl=='log' else yo
    # Select correct data
    ind = np.ma.where((x>xlim[0]) & (x<xlim[1]) & (y>ylim[0]) & (y<ylim[1]))[0]
    x = x[ind]  ; y = y[ind]
    w = wo[ind] ; a = ao[ind]
    # Histogram
    # numpy
    if method=='numpy':
        hist,xs,ys = np.histogram2d(x,y,weights=a*w,bins=gridsize,range=[xlim,ylim])
        if avar: hist_avg,xs,ys = np.histogram2d(x,y,weights=a,bins=gridsize,range=[xlim,ylim])
        if len(w)==0: 
            hist=np.zeros(gridsize)
            if avar: hist_avg=np.zeros(gridsize)
    # pNbody
    elif method=='pNbody':
        mapkw={}
        if   xvar in [ 'x', 'y', 'z'] and yvar in [ 'x', 'y', 'z']: mapkw['space']='pos'
        elif xvar in ['vx','vy','vz'] and yvar in ['vx','vy','vz']: mapkw['space']='vel'
        else: raise Exception('Invalid xvar/yvar ({}/{}) for pNbody method'.format(xvar,yvar))
        mapkw['shape']=gridsize
        mapkw['view']=xvar[-1]+yvar[-1]
        mapkw['size']=(float(xlim[1]-xlim[0]),float(ylim[1]-ylim[0]))
        # Get histogram
        hist = sim.ComputeMap(mode=wvar,**mapkw)
        if avar: hist_avg=sim.ComputeMap(mode=avar,**mapkw)
        xs=np.linspace(xlim[0],xlim[1],gridsize[1]+1)
        ys=np.linspace(ylim[0],ylim[1],gridsize[0]+1)
    else: raise Exception('Invalid histogram method: {}'.format(method))
    # Average histogram
    if avar:
        good = np.where(hist_avg != 0)
        hist[good] = hist[good]/hist_avg[good]
    # Set units
    ws=SimArray(hist,w.units).in_unitsys(units)
    try: 
        xs = SimArray(.5*(xs[:-1]+xs[1:]), x.units)
        ys = SimArray(.5*(ys[:-1]+ys[1:]), y.units)
    except AttributeError: 
        xs = .5*(xs[:-1]+xs[1:])
        ys = .5*(ys[:-1]+ys[1:])
    # Set histogram limits
    if wlim[0]==wlim[1]: 
        if   wscl=='log'   : wlim=SimArray([np.min(ws[ws>0]),np.max(ws[ws>0])],ws.units)
        elif wscl=='symlog': wlim=SimArray([np.min(np.abs(ws[ws!=0])),np.max(np.abs(ws[ws!=0]))],ws.units)
        else               : wlim=SimArray([np.min(ws),np.max(ws)],ws.units)
    if hasattr(wlim,'units'): wlim=wlim.in_units(ws.units)
    else                    : wlim=SimArray(wlim,ws.units)
    if method=='pNbody': wlim/=10.
    # Assemble output
    out=dict(xvar=xvar,xs=xs,xlim=xlim,xscl=xscl,
             yvar=yvar,ys=ys,ylim=ylim,yscl=yscl,
             wvar=wvar,ws=ws,wlim=wlim,wscl=wscl,
             avar=avar,
             galaxy=galaxy,family=family,cmeth=cmeth)
    if avar: out['as']=hist_avg
    # Save output
    from mmlutils import mmlio,mmlfiles
    #filename=mmlio.askquest('Enter filename to save histogram data:')
    if filename!='None' and mmlfiles.prep_write(filename):
        f=open(filename,'wb')
        pickle.dump(out,f,-1)
        f.close()
    # Return output
    return out

