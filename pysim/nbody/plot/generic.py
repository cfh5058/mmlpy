"""

generic
=======

Flexible and general plotting functions

"""

import numpy as np
from ..analysis import profile, angmom, halo, mapping
from .. import config
from ..snapshot import SimSnap
from ..array import SimArray
from ..units import NoUnit


def plot2d(sim0,xvar,yvar,histmeth='numpy',plotmeth='contourf',axs=None,mksquare=True,
           plotfile=None,overwrite=False,histfile=None,owhist=False,**kwargs):
    """
    Plot a 2D histogram or map
    """
    import ImageDraw,ImageFont,Image,matplotlib
    from matplotlib import ticker, colors, pyplot
    from mmlutils import mmlplot,mmlfiles
    matplotlib.rcParams['path.simplify']=False
    # Check if plot should be overwritten
    if kwargs.get('make_plot',True) and not mmlfiles.prep_write(plotfile,overwrite=overwrite):
        return
    # Identify simulation type
    pNbody = False if issubclass(sim0.__class__,SimSnap) else True
    # Label size
    fontsize=kwargs.get('textsize',None)
    ticksize=kwargs.get('ticksize',None)
    if ticksize:
        matplotlib.rcParams['font.size']=ticksize
    # Pars input
    if plotmeth=='pilimg': 
        axs=None ; clear=None
    else:
        clear=kwargs.pop('clear',True)
        if axs:
            plt=axs
            kwargs.setdefault('xlim',axs.get_xlim())
            kwargs.setdefault('ylim',axs.get_ylim())
        else: 
            if clear: pyplot.clf()
            plt=pyplot.subplot(1,1,1)#,aspect='equal')
    if plotmeth=='pilimg': kwargs.setdefault('gridsize',(512,512))
    else                 : kwargs.setdefault('gridsize',(100,100))
    ret_im=kwargs.pop('ret_im',False)
    # Select proper view of simulation
    sim=sim0
    # galaxy=kwargs.pop('galaxy',None)
    # family=kwargs.pop('family',None)
    # cmeth=kwargs.pop('cmeth',None)
    # sim=langmm.align(sim0,galaxy=galaxy,family=family,cmeth=cmeth,pNbody=pNbody)
    # Get histogram
    if 'imgdata' in kwargs: img=kwargs['imgdata']
    else: img=mapping.hist2d(sim,xvar,yvar,method=histmeth,
                             filename=histfile,overwrite=owhist,**kwargs)
    xvar=img['xvar'] ; xs=img['xs'] ; xlim=img['xlim'] ; xscl=img['xscl']
    yvar=img['yvar'] ; ys=img['ys'] ; ylim=img['ylim'] ; yscl=img['yscl']
    wvar=img['wvar'] ; ws=img['ws'] ; wlim=img['wlim'] ; wscl=img['wscl']
    avar=img['avar']
    galaxy=img['galaxy'] ; family=img['family'] ; cmeth=img['cmeth']
    # Force box
    if plotmeth!='pilimg' and mksquare:
        xaspect=xlim[1]-xlim[0]
        yaspect=ylim[1]-ylim[0]
        plt.set_aspect(xaspect/yaspect)
    # Get ellipse fit
    if plotmeth!='pilimg' and kwargs.get('plot_ellipse',False): 
        ellipse(sim,xvar,yvar,imgdata=img,axs=plt,**kwargs)
    # Labels
    xlabel=kwargs.pop('xlabel',mklabel(units=xs.units,scale=xscl,axes=xvar,latex=(plotmeth!='pilimg')))
    ylabel=kwargs.pop('ylabel',mklabel(units=ys.units,scale=yscl,axes=yvar,latex=(plotmeth!='pilimg')))
    clabel=kwargs.pop('clabel',mklabel(units=ws.units,scale=wscl,axes=wvar,latex=(plotmeth!='pilimg')))
    # Color map
    cmap=kwargs.pop('cmap',None)
    if not cmap:
        famstr=family if family else 'all'
        galstr=galaxy if galaxy else 'all'
        galstr='0' if galstr=='all' else galstr[-1]
        cmaptype='mmlparttype'
        if   histmeth=='numpy' : cmaptype+='_sk0p5'
        elif histmeth=='pNbody': cmaptype+='_sk0p5'
        cmap=mmlplot.get_cmap(famstr+galstr,cmaptype=cmaptype)
    cmap_colors=mmlplot.get_cmaparray(cmap)
    axbgclr=kwargs.pop('bgcolor',cmap_colors[0])
    axfgclr=kwargs.pop('fgcolor',cmap_colors[-1])
    clrmin=1#int(0.01*len(cmap_colors))
    if plotmeth!='pilimg': cmap=colors.ListedColormap(cmap_colors[clrmin:])
    cmap.set_under(color=cmap_colors[0])
    cmap.set_over(color=cmap_colors[-1])
    cmap.set_bad(color=cmap_colors[0]) #no effect for contourf
    # Plot image
    # PIL image
    if plotmeth=='pilimg':
        # Translate to color index
        clrkw=dict(clip=True,scale=wscl,mn=wlim[0],mx=wlim[1],clrmin=clrmin)
        if wscl=='log': clrkw['cd']=10.**((np.log10(clrkw['mx'])+np.log10(clrkw['mn']))/2.-1.)
        matint=mmlplot.map_array(ws,**clrkw)
        if wscl=='log': matint[ws<=0]=0
        # Create PIL image
        matint_T=np.transpose(matint)
        imgpil=Image.fromstring("P",(matint_T.shape[1],matint_T.shape[0]),matint_T.tostring())
        # Color map
        imgpil.putpalette(mmlplot.cmap2palette(cmap))
        # Return image
        if ret_im : return imgpil
    # pyplot method
    else:
        # Make scale linear to allow extend
        if wscl=='log':
            matscl=np.ma.log10(ws)
            matlim=[np.log10(wlim[0]),np.log10(wlim[1])]
            normDEF=None
        else:
            matscl=ws
            matlim=wlim
            if wscl=='symlog': 
                normDEF=mmlplot.SymLogNorm(wlim[0],vmin=-wlim[1],vmax=wlim[1],linscale=0.01)
            else:
                normDEF=None
        norm=kwargs.pop('norm',normDEF)
        matscl_T=np.transpose(matscl)
        # Return image
        if ret_im : return plt.imshow(matscl_T,origin='down',vmin=wlim[0],vmax=wlim[1],norm=norm,
                                      aspect='auto',cmap=cmap,extent=[xlim[0],xlim[1],ylim[0],ylim[1]])
        # Contourf
        if plotmeth=='contourf':
            # Get levels
            nlevels=kwargs.pop('nlevels',51)
            if norm: levelsDEF=norm.inverse(np.linspace(0.,1.,nlevels,endpoint=True))
            else   : levelsDEF=np.linspace(matlim[0],matlim[1],nlevels)
            levels=kwargs.pop('levels',levelsDEF)
            # Plot
            cs = plt.contour(xs,ys,matscl_T,levels,cmap=cmap,norm=norm)
            cf = plt.contourf(xs,ys,matscl_T,levels,cmap=cmap,extend='both',norm=norm)#,**kwargs)
        else: raise Exception('Invalid pyplot plot method: {}'.format(plotmeth))
        # Background color
        #if not axs: 
        plt.set_axis_bgcolor(axbgclr)
        # Axes labels
        #if not axs: 
        if xlabel: plt.set_xlabel(xlabel,size=fontsize)
        if ylabel: plt.set_ylabel(ylabel,size=fontsize)
        # Color bar
        if kwargs.get('colorbar',False):
            cbfmt=kwargs.get('cbtickfmt',None)
            if not cbfmt:
                if   wscl=='log'   : cbfmt="%.2f"
                elif wscl=='symlog': cbfmt="% .2e"
                else               : cbfmt="% .2e"
            cbticks=kwargs.get('cbtickloc',None)
            cb = pyplot.colorbar(cf,format=cbfmt,ticks=cbticks,ax=plt).set_label(r''+clabel,size=fontsize) 
        else: cb=None
        # Legend
        if kwargs.get('legend',False): 
            lg=plt.legend(loc=2)
        else:
            lg=None

    # Annotate
    if kwargs.get('annotate',False):
        kwargs.setdefault('annotate_arrow',True)
        kwargs.setdefault('annotate_text' ,True)
        # Get arrows
        arrows=kwargs.pop('arrows',[])
        if galaxy and kwargs['annotate_arrow'] and not arrows:
            gallist=sim0.galaxies()
            if galaxy in gallist: gallist.remove(galaxy)
            arrows=halo.list_centers_by_var(sim0,gals=gallist,pvars=[xvar,yvar],wvars=[wvar,avar])
            arrows=[(arw[0].in_units(xs.units),arw[1].in_units(ys.units)) for arw in arrows]
        # Get values for annotations
        if pNbody: t_val=SimArray(sim0.atime,sim0.get_units('t'))
        else     : t_val=sim0.properties['time']
        x_val=SimArray(float(xlim[1]-xlim[0]),xs.units)
        y_val=SimArray(float(ylim[1]-ylim[0]),ys.units)
        d_val=(x_val,y_val)
        # Unit labels
        tstr=val2str(t_val,label='t',units='Gyr',latex=False)#(plotmeth!='pilimg'))
        dstr=dim2str(d_val,label='' ,units=None ,latex=False)#(plotmeth!='pilimg'))
        # Text parameters
        text=kwargs.pop('text',[])
        if kwargs['annotate_text'] and not text:
            if plotmeth=='pilimg': text=[['l','u',tstr],['r','u',dstr]]
            else                 : text=[['l','u',tstr]]
        # Annotate
        if plotmeth=='pilimg':
            fntclr = 255
            imgobj = imgpil
        else:
            fntclr = axfgclr
            imgobj = plt
        annotate(imgobj,color=fntclr,xlim=xlim,ylim=ylim,
                 text    =kwargs.get('text'    ,text       ),
                 arrows  =kwargs.get('arrows'  ,arrows     ),
                 arwsize =kwargs.get('arwsize' ,ws.shape[0]),
                 arwtail =kwargs.get('arwtail' ,0.         ),
                 arwcolor=kwargs.get('arwcolor',fntclr     ),
                 arwwidth=kwargs.get('arwwidth',1.         ))
        # Re-adjust limits
        if plotmeth!='pilimg':# and not axs:
            plt.set_xlim((xlim[0],xlim[1])) 
            plt.set_ylim((ylim[0],ylim[1]))
    # Add to existing image
    if plotmeth=='pilimg' and axs: 
        axs[0].paste(imgpil.convert('RGB'),(axs[1],axs[2],axs[1]+matint_T.shape[0]+mathint_T.shape[1]))
        imgpil=axs[0]
    # Save image
    if kwargs.get('make_plot',True) and mmlfiles.prep_write(plotfile,overwrite=overwrite):
        # Rasterize
        # if plotfile.endswith('.eps'):
        #     for d in cf.collections: d.set_rasterized(True)
        # Save
        if config['verbose']: print "Saving "+plotfile
        if   plotmeth=='pilimg'  : imgpil.save(plotfile)
        else                     : pyplot.savefig(plotfile,bbox_inches='tight',
                                                  bbox_extra_artists=[x for x in [lg,cb] if x])
    # Return axes
    if plotmeth=='pilimg': return imgpil
    else                 : return plt

def ellipse(sim,xvar,yvar,axs=None,ellpfile=None,plotfile=None,
            nphi=100,ellipcolor='m',ellipwidth=2.0,ellipalpha=1.0,**kwargs):
    """
    Plot ellipses fitted to histogram
    """
    from mmlastro import ellipse
    # from mmlastro import mmlellipse
    # Create axes/use one provided
    clear=kwargs.pop('clear',True)
    if axs:
        plt=axs
        kwargs.setdefault('xlim',axs.get_xlim())
        kwargs.setdefault('ylim',axs.get_ylim())
    else: 
        import matplotlib.pyplot as pyplot
        if clear: pyplot.clf()
        plt=pyplot.subplot(1,1,1)
    # Get data
    elp,imgdim=mapping.ellipse(sim,xvar,yvar,filename=ellpfile,**kwargs)
    # Plot
    xlim0=plt.get_xlim() ; ylim0=plt.get_ylim()
    phiLIST=np.linspace(0.,2*np.pi,nphi)
    for ie in elp:
        ix,iy=ellipse.phi2xy(phiLIST,**ie)
        # ix,iy=mmlellipse.phi2xy_ellipse(phiLIST,**ie)
        plt.plot(ix,iy,color=ellipcolor,alpha=ellipalpha,lw=ellipwidth,
                 label='a = {}'.format(ie['a']))
    # Reset limits
    if axs: xlim=imgdim[0] ; ylim=imgdim[1]
    else  : xlim=xlim0     ; ylim=ylim0
    plt.set_xlim(xlim)
    plt.set_ylim(ylim)
    # Save image
    if kwargs.get('make_plot',True) and plotfile:
        if config['verbose']: print "Saving "+plotfile
        plt.savefig(plotfile)
    # Return axes
    return plt

def annotate(imgobj,text=[],arrows=[],txtpad=1./20,arwpad=1./40,**kwargs):
    """
    Annotates plots
       text  : List of (x,y,str) text positions and strings
       arrows: List of (x,y) arrow positions
       color : Color of text & arrows
    """
    import ImageDraw,ImageFont,Image
    from matplotlib import ticker, colors, pyplot
    from mmlutils import mmlplot
    from mmlutils.mmlmath import polar_square
    # Identify image object
    pilflag=isinstance(imgobj,Image.Image)
    # Get keywords
    color   =kwargs.get('color'   ,255  )
    arwcolor=kwargs.get('arwcolor',color)
    arwwidth=kwargs.get('arwwidth',1.   )
    arwsize =kwargs.get('arwsize' ,100  )
    arwtail =kwargs.get('arwtail' ,0.   )
    if pilflag:
        draw=kwargs.get('draw',ImageDraw.Draw(imgobj))
        font=kwargs.get('font',ImageFont.load_default())
    # Get limits
    if pilflag:
        ximg=(0,imgobj.size[0])
        yimg=(0,imgobj.size[1])
        xdat=kwargs.get('xlim',ximg)
        ydat=kwargs.get('ylim',yimg)
        xscl=kwargs.get('xscl','lin')
        yscl=kwargs.get('yscl','lin')
    else:
        xdat=kwargs.get('xlim',imgobj.get_xlim())
        ydat=kwargs.get('ylim',imgobj.get_ylim())
        xscl=kwargs.get('xscl',imgobj.get_xscale())
        yscl=kwargs.get('yscl',imgobj.get_yscale())
    if xscl=='log': xdat=(np.log10(xdat[0]),np.log10(xdat[1]))
    if yscl=='log': ydat=(np.log10(ydat[0]),np.log10(ydat[1]))
    if not pilflag:
        ximg=xdat
        yimg=ydat
    # Function to convert text positions
    xtxtpad=txtpad*(ximg[1]-ximg[0])
    ytxtpad=txtpad*(yimg[1]-yimg[0])
    def txtloc(xi,yi):
        # Define positions codes
        if isinstance(xi,str):
            if   xi=='l': xi=xdat[0]
            elif xi=='r': xi=xdat[1]
            else: raise Exception('Invalid x code: {}'.format(xi))
        elif xscl=='log': xi=np.log10(xi)
        if isinstance(yi,str):
            if   yi=='d': yi=ydat[0]
            elif yi=='u': yi=ydat[1]
            else: raise Exception('Invalid y code: {}'.format(yi))
        elif yscl=='log': yi=np.log10(yi)
        # Convert from data to image units
        xo=(xi-xdat[0])*(ximg[1]-ximg[0])/(xdat[1]-xdat[0])+ximg[0]
        yo=(yi-ydat[0])*(yimg[1]-yimg[0])/(ydat[1]-ydat[0])+yimg[0]
        # Ensure positions inside border
        xo=max(ximg[0]+xtxtpad,xo) ; xo=min(ximg[1]-xtxtpad,xo)
        yo=max(yimg[0]+ytxtpad,yo) ; yo=min(yimg[1]-ytxtpad,yo)
        # Set allignments
        if   xo<=(ximg[0]+xtxtpad): xa='left'
        elif xo>=(ximg[1]-xtxtpad): xa='right'
        else                      : xa='center'
        if   yo<=(yimg[0]+ytxtpad): ya='bottom'
        elif yo>=(yimg[1]-xtxtpad): ya='top'
        else                      : ya='center'
        # Rescale log
        if not pilflag:
            if xscl=='log': xo=10**xo
            if yscl=='log': yo=10**yo
        # Return
        return xo,yo,xa,ya
    # Text
    for itxt in text:
        x0,y0,s0=itxt
        x1,y1,ha,va=txtloc(x0,y0)
        if pilflag:
            dwid,dhgt=font.getsize(s0)
            if   ha=='left'  : pass
            elif ha=='center': x1-=dwid/2
            elif ha=='right' : x1-=dwid
            if   va=='bottom': y1+=dhgt
            elif ha=='center': y1+=dhgt/2
            elif ha=='top'   : pass
            draw.text((x1,ximg[1]-y1),s0,fill=color,font=font)
        else: imgobj.text(x1,y1,s0,color=color,ha=ha,va=va)
    # Function to convert arrow positions
    xarwpad=arwpad*(ximg[1]-ximg[0])
    yarwpad=arwpad*(yimg[1]-yimg[0])
    def arwloc(xi,yi):
        iphi=np.arctan((yi/float(ydat[1]-ydat[0]))/(xi/float(xdat[1]-xdat[0])))
        if xi<xdat[0]+xarwpad or xi>xdat[1]-xarwpad or yi<ydat[0]+yarwpad or yi>ydat[1]-yarwpad:
            if   iphi>0 and yi<0: iphi+=np.pi
            elif iphi<0 and xi<0: iphi+=np.pi
            idir=polar_square(iphi,side=(1.-2.*arwpad))
            xo=(idir[0]+0.5)*(ximg[1]-ximg[0])+ximg[0]
            yo=(idir[1]+0.5)*(yimg[1]-yimg[0])+yimg[0]
        else:
            xo=xi
            yo=yi
        # Get position of arrow tail
        xarwtail=(ximg[1]-ximg[0])*arwtail*np.cos(iphi)
        yarwtail=(yimg[1]-yimg[0])*arwtail*np.sin(iphi)
        xe=xo-xarwtail
        ye=yo-yarwtail
        # Rescale log
        if not pilflag:
            if xscl=='log': xo=10**xo ; xe=10**xe
            if yscl=='log': yo=10**yo ; ye=10**ye
        return xo,yo,iphi,xe,ye
    # Arrows
    for iarw in arrows:
        x0,y0=iarw
        x1,y1,phi1,x2,y2=arwloc(x0,y0)
        if pilflag: draw.text([x1,yimg[1]-y1],'@',fill=arwcolor,font=font)
        else      : 
            imgobj.scatter(x1,y1,c=arwcolor,marker=(3,0,-90+phi1*180./np.pi),s=arwsize)
            if arwtail: imgobj.plot([x1,x2],[y1,y2],c=arwcolor,lw=arwwidth)
    # Return
    return

def mklabel(axes=None,units=None,scale=None,latex=False):
    """
    Creates axes label strings
    """
    # Get unit string
    if units:
        if units==NoUnit() or units==1: unistr=''
        else:
            if latex: unistr=units.latex()
            else    : unistr=str(units)
    else:
        unistr=''
    # Get label string
    if axes: axsstr=axes
    else   : axsstr=''
    # Get combined string
    if axsstr and unistr: outstr=axsstr+'/'+unistr
    else                : outstr=axsstr+unistr
    # Add scale string
    if len(outstr)>0 and scale=='log':
        if latex: outstr='log_{10}('+outstr+')'
        else    : outstr='log10('+outstr+')'
    # Return string
    if latex: return r''+'$'+outstr+'$'
    else    : return outstr

def val2str(val,label=None,fmtstr='5.2f',latex=False,units=None):
    """
    Creates a string from a SimArray
    """
    # Convert units
    if units: val=val.in_units(units)
    # Get label string
    if label: 
        labstr=label+' ='
        if latex: labstr=r"$\mathrm{"+labstr+"}$"
        labstr+=" "
    else: labstr=""
    # Get value string
    valstr=str("{:"+fmtstr+"}").format(float(val))
    # Get unit string
    if latex: unistr=' $'+val.units.latex()+'$'
    else    : unistr=' '+str(val.units)
    # Return combined string
    return labstr+valstr+unistr

def dim2str(dim,label=None,fmtstr='3.0f',latex=False,units=None):
    """
    Creates string from SimArray dimension
    """
    # Convert units
    if units: dim=(dim[0].in_units(units[0]),dim[1].in_units(units[1]))
    # Get label string
    if label: 
        labstr=label+' ='
        if latex: labstr=r"$\mathrm{"+labstr+"}$"
        labstr+=" "
    else: labstr=""
    # Get value strings
    xvalstr=str("{:"+fmtstr+"}").format(float(dim[0]))
    yvalstr=str("{:"+fmtstr+"}").format(float(dim[1]))
    # Get unit string
    if latex:
        xunistr=' $'+dim[0].units.latex()+'$'
        yunistr=' $'+dim[1].units.latex()+'$'
    else:
        xunistr=' '+str(dim[0].units)
        yunistr=' '+str(dim[1].units)
    # Combine strings
    if xunistr==yunistr: dimstr=r""+labstr+xvalstr+" x "+yvalstr+xunistr
    else               : dimstr=r""+labstr+xvalstr+xunistr+" x "+yvalstr+yunistr
    return dimstr

def image_pNbody(sim,logscale=False,scalemin=None,scalemax=None,annotate=False,
                 filename=None,ret_im=False,cmap=None,**mapkw0):
    """
    Create an image of the sim using pNbody mapping function
    """
    from mmlutils import mmlplot
    import ImageDraw,ImageFont,Image,copy
    from pNbody import libutil
    # Color map
    if cmap:
        from mmlutils.mmlplot import get_cmaparray
        cmap_colors=mmlplot.get_cmaparray(cmap)
    # Get histogram
    mapkw=libutil.extract_parameters([],mapkw0,sim.defaultparameters)
    hist=sim.CombiMap(mapkw)
    if mapkw['mode'] in ['n','m','mass']:
        nhist=copy.deepcopy(hist)
    else:
        mapkw['mode']='m'
        nhist=sim.CombiMap(mapkw)
    mat=np.ma.masked_where(nhist==0,hist)
    # Map to specified limits
    clrkw=dict(clip=True)
    if cmap: clrkw['clrmin']=int(0.01*len(cmap_colors))
    if logscale: clrkw['scale']='log'
    if scalemin: clrkw['mn']=scalemin
    if scalemax: clrkw['mx']=scalemax/(10.)
    if scalemin and scalemax and logscale:
        clrkw['cd']=10.**((np.log10(scalemin)+np.log10(scalemax))/2.-1.)
    matint=mmlplot.map_array(mat,**clrkw).filled(0)
    print mat.min(),mat.max()
    print matint.min(),matint.max()
    print clrkw
    if ret_im : return matint
    # Create image object from map
    matint_T=np.transpose(matint)
    imgpil=Image.fromstring("P",(matint_T.shape[1],matint_T.shape[0]),matint_T.tostring())
    # Color map
    imgpil.putpalette(mmlplot.cmap2palette(cmap))
    # Annotate
    if annotate:
        imgsiz = imgpil.size
        fntclr = 255
        draw   = ImageDraw.Draw(imgpil)
        font   = ImageFont.load_default()
        fntwid, fnthgt = font.getsize('test')
        locfct = 20
        xlocpad=imgsiz[0]/locfct
        ylocpad=imgsiz[1]/locfct
        # Conversion factors
        time2Gyr=sim.localsystem_of_units.get_UnitTime_in_s()*3.16888e-8/1.0e9
        leng2kpc=sim.localsystem_of_units.get_UnitLength_in_cm()/3.0856780e+21
        # Label specifying time
        tloc = (xlocpad,ylocpad)
        tstr = "t = {:5.2f} {}".format(time2Gyr*sim.atime,'Gyr')
        draw.text(tloc,tstr,fill=fntclr,font=font)
        # Label specifying size
        dstr = "{:3.0f} x {:3.0f} {}".format(leng2kpc*mapkw['size'][0],leng2kpc*mapkw['size'][1],'kpc')
        dwid,dhgt=font.getsize(dstr)
        dloc = (imgsiz[0]-xlocpad-dwid,ylocpad)
        draw.text(dloc,dstr,fill=fntclr,font=font)
    # Save image
    imgpil.save(filename)


def hist2d(xo, yo, weights=None, mass=None, gridsize=(100,100), nbins = None, nlevels=50, make_plot = True, **kwargs):
    """
    Plot 2D histogram for arbitrary arrays that get passed in.

    **Input:**

       *x*: array

       *y*: array

    **Optional keywords:**

       *x_range*: list, array, or tuple
         size(x_range) must be 2. Specifies the X range.

       *y_range*: tuple
         size(y_range) must be 2. Specifies the Y range.

       *gridsize*: (int, int) (default (100,100)) 
         number of bins to use for the 2D histogram

       *nbins*: int
         number of bins for the histogram - if specified, gridsize is set to (nbins,nbins)

       *nlevels*: int
         number of levels to use for the contours

       *logscale*: boolean
         whether to use log or linear spaced contours

       *weights*: numpy array of same length as x and y
         if weights is passed, color corresponds to
         the mean value of weights in each cell

       *mass*: numpy array of masses same length as x andy
         must also have weights passed in.  If you just
         want to weight by mass, pass the masses to weights

       *colorbar*: boolean
         draw a colorbar
         
       *scalemin*: float
         minimum value to use for the color scale

       *scalemax*: float
         maximum value to use for the color scale
    """
    global config
    
    # process keywords
    x_range = kwargs.get('x_range',None)
    y_range = kwargs.get('y_range',None)
    xlogrange = kwargs.get('xlogrange', False)
    ylogrange = kwargs.get('ylogrange', False)
    ret_im = kwargs.get('ret_im',False)

    # Y range
    if y_range is not None :
        if len(y_range) != 2: raise RuntimeError("Range must be a length 2 list or array")
        if ylogrange: y_range = [np.log10(y_range[0]),np.log10(y_range[1])]
    else:
        if ylogrange: y_range = [np.log10(np.min(yo)),np.log10(np.max(yo))]
        else        : y_range = [np.min(yo),np.max(yo)]
    kwargs['y_range'] = y_range
    y = np.log10(yo) if ylogrange else yo
    # X range
    if x_range is not None:
        if len(x_range) != 2: raise RuntimeError("Range must be a length 2 list or array")
        if xlogrange: x_range = [np.log10(x_range[0]),np.log10(x_range[1])]
    else:
        if xlogrange: x_range = [np.log10(np.min(xo)), np.log10(np.max(xo))]
        else        : x_range = [np.min(xo),np.max(xo)]
    kwargs['x_range'] = x_range
    x = np.log10(xo) if xlogrange else xo
    # Grid
    if nbins is not None: gridsize = (nbins,nbins)
    # Select correct data
    ind = np.where((x > x_range[0]) & (x < x_range[1]) &
                   (y > y_range[0]) & (y < y_range[1]))
    x = x[ind[0]]
    y = y[ind[0]]
    
    draw_contours = False
    if weights is not None and mass is not None: 
        draw_contours = True
        weights = weights[ind[0]]
        mass = mass[ind[0]]

        # produce a mass-weighted histogram of average weight values at each bin
        hist, ys, xs = np.histogram2d(y, x, weights=weights*mass, bins=gridsize,range=[y_range,x_range])
        hist_mass, ys, xs = np.histogram2d(y, x, weights=mass,bins=gridsize,range=[y_range,x_range])
        good = np.where(hist_mass > 0)
        hist[good] = hist[good]/hist_mass[good]
            
    else:
        if weights is not None : 
            # produce a weighted histogram
            weights = weights[ind[0]]
        elif mass is not None: 
            # produce a mass histogram
            weights = mass[ind[0]]
      
        hist, ys, xs = np.histogram2d(y, x, weights=weights, bins=gridsize,range=[y_range,x_range])
        
    try: 
        hist = SimArray(hist,weights.units)
    except AttributeError: 
        hist = SimArray(hist)

        
    try: 
        xs = SimArray(.5*(xs[:-1]+xs[1:]), x.units)
        ys = SimArray(.5*(ys[:-1]+ys[1:]), y.units)
    except AttributeError: 
        xs = .5*(xs[:-1]+xs[1:])
        ys = .5*(ys[:-1]+ys[1:])


    if ret_im :
        return make_contour_plot(hist, xs, ys, **kwargs)

    if make_plot : 
        make_contour_plot(hist, xs, ys, nlevels=nlevels, **kwargs)
        if draw_contours:
            make_contour_plot(SimArray(density_mass, mass.units),xs,ys,filled=False,clear=False,
                              nlevels=nlevels/2,colorbar=False,colors='black')

    return hist, xs, ys
    


def gauss_kde(xo, yo, weights=None, mass = None, gridsize = (100,100), nbins = None, nlevels = 10,
              make_plot = True, nmin = None, nmax = None, **kwargs) :

    """
    Plot 2D gaussian kernel density estimate (KDE) given values at points (*x*, *y*). 

    Behavior changes depending on which keywords are passed: 

    If a *weights* array is supplied, produce a weighted KDE.
    
    If a *mass* array is supplied, a mass density is computed. 

    If both *weights* and *mass* are supplied, a mass-averaged KDE of the weights is 
    computed. 

    By default, norm=False is passed to :func:`~pynbody.plot.util.fast_kde` meaning
    that the result returned *is not* normalized such that the integral over the area
    equals one. 

Since this function produces a density estimate, the units of the
    output grid are different than for the output of
    :func:`~pynbody.plot.generic.hist2d`. To get to the same units,
    you must multiply by the size of the cells.

    
    **Input:**

       *xo*: array

       *yo*: array

    **Optional keywords:**

       *mass*: numpy array of same length as x and y 
         particle masses to be used for weighting
    
       *weights*: numpy array of same length as x and y
         if weights is passed, color corresponds to
         the mean value of weights in each cell

       *nmin*: float (default None)
         if *weights* and *mass* are both specified, the mass-weighted
         contours are only drawn where the mass exceeds *nmin*. 
          
       *gridsize*: (int, int) (default: 100,100)
         size of grid for computing the density estimate

       *nbins*: int
         number of bins for the histogram - if specified, gridsize is set to (nbins,nbins)

       *make_plot*: boolean (default: True)
         whether or not to produce a plot
         
    **Keywords passed to** :func:`~pynbody.plot.util.fast_kde`:
       
       *norm*: boolean (default: False) 
         If False, the output is only corrected for the kernel. If True,
         the result is normalized such that the integral over the area 
         yields 1. 
    
       *nocorrelation*: (default: False) If True, the correlation
         between the x and y coords will be ignored when preforming
         the KDE.
         
    **Keywords passed to** :func:`~pynbody.plot.generic.make_contour_plot`:

       *x_range*: list, array, or tuple
         size(x_range) must be 2. Specifies the X range.

       *y_range*: tuple
         size(y_range) must be 2. Specifies the Y range.
       
       *nlevels*: int
         number of levels to use for the contours

       *logscale*: boolean
         whether to use log or linear spaced contours

       *colorbar*: boolean
         draw a colorbar
         
       *scalemin*: float
         minimum value to use for the color scale

       *scalemax*: float
         maximum value to use for the color scale
    """
    from util import fast_kde
    from scipy.stats.kde import gaussian_kde

    global config
    
    # process keywords
    x_range = kwargs.get('x_range',None)
    y_range = kwargs.get('y_range',None)
    xlogrange = kwargs.get('xlogrange', False)
    ylogrange = kwargs.get('ylogrange', False)
    ret_im = kwargs.get('ret_im',False)

    
    if y_range is not None :
        if len(y_range) != 2 : 
            raise RuntimeError("Range must be a length 2 list or array")
    else:
        if ylogrange:
            y_range = [np.log10(np.min(yo)),np.log10(np.max(yo))]
        else:
            y_range = [np.min(yo),np.max(yo)]
            
    if x_range is not None:
        if len(x_range) != 2 :
            raise RuntimeError("Range must be a length 2 list or array")
    else:
        if xlogrange:
            x_range = [np.log10(np.min(xo)), np.log10(np.max(xo))]
        else:
            x_range = [np.min(xo),np.max(xo)]

    if (xlogrange):
        x = np.log10(xo)
    else :
        x = xo
    if (ylogrange):
        y = np.log10(yo)
    else :
        y = yo
        
    if nbins is not None:
        gridsize = (nbins,nbins)

    ind = np.where((x > x_range[0]) & (x < x_range[1]) &
                   (y > y_range[0]) & (y < y_range[1]))

    x = x[ind[0]]
    y = y[ind[0]]

    try:
        xs = SimArray(np.linspace(x_range[0], x_range[1], gridsize[0]+1),x.units)
    except AttributeError:
        xs = np.linspace(x_range[0], x_range[1], gridsize[0]+1)
    xs = .5*(xs[:-1]+xs[1:])
    try:
        ys = SimArray(np.linspace(y_range[0], y_range[1], gridsize[1]+1),y.units)
    except AttributeError:
        ys = np.linspace(y_range[0], y_range[1], gridsize[1]+1)
    ys = .5*(ys[:-1]+ys[1:])
    
    extents = [x_range[0],x_range[1],y_range[0],y_range[1]]
    
    draw_contours = False
    if weights is not None and mass is not None: 
        draw_contours = True
        weights = weights[ind[0]]
        mass = mass[ind[0]]
        
        # produce a mass-weighted gaussian KDE of average weight values at each bin
        density = fast_kde(x, y, weights=weights*mass, gridsize=gridsize,extents=extents, **kwargs)
        density_mass = fast_kde(x, y, weights=mass, gridsize=gridsize,extents=extents, **kwargs)
        good = np.where(density_mass > 0)
        density[good] = density[good]/density_mass[good]

        if nmin is not None : 
            density *= density_mass > nmin

    else:
        # produce a weighted gaussian KDE
        if weights is not None : 
            weights = weights[ind[0]]
        elif mass is not None : 
            weights = mass[ind[0]]

        density = fast_kde(x, y, weights=weights, gridsize=gridsize,
                           extents=extents, **kwargs)

    try: 
        density = SimArray(density,weights.units)
    except AttributeError: 
        density = SimArray(density)

    if ret_im: 
        return make_contour_plot(density,xs,ys,**kwargs)

    if make_plot : 
        make_contour_plot(density,xs,ys,**kwargs)
        if draw_contours:
            make_contour_plot(SimArray(density_mass, mass.units),xs,ys,filled=False,clear=False,colorbar=False,colors='black',
                              scalemin=nmin,scalemax=nmax,nlevels=10)

    return density, xs, ys


def make_contour_plot(arr, xs, ys, x_range=None, y_range=None, nlevels = 20, 
                      axbgclr = 'w', axfgclr = 'k',
                      logscale=True, xlogrange=False, ylogrange=False, 
                      subplot=False, colorbar=False, ret_im=False, cmap=None,
                      clear=True,legend=False, scalemin = None, clip=False,
                      scalemax = None, filename = None, annotate = False, **kwargs) : 
    """
    Plot a contour plot of grid *arr* corresponding to bin centers
    specified by *xs* and *ys*.  Labels the axes and colobar with
    proper units taken from x

    Called by :func:`~pynbody.plot.generic.hist2d` and
    :func:`~pynbody.plot.generic.gauss_density`.
    
    **Input**: 
    
       *arr*: 2D array to plot

       *xs*: x-coordinates of bins
       
       *ys*: y-coordinates of bins

    **Optional Keywords**:
          
       *x_range*: list, array, or tuple (default = None)
         size(x_range) must be 2. Specifies the X range.

       *y_range*: tuple (default = None)
         size(y_range) must be 2. Specifies the Y range.

       *xlogrange*: boolean (default = False)
         whether the x-axis should have a log scale

       *ylogrange*: boolean (default = False)
         whether the y-axis should have a log scale

       *nlevels*: int (default = 20)
         number of levels to use for the contours

       *logscale*: boolean (default = True)
         whether to use log or linear spaced contours

       *colorbar*: boolean (default = False)
         draw a colorbar
         
       *scalemin*: float (default = arr.min())
         minimum value to use for the color scale

       *scalemax*: float (default = arr.max())
         maximum value to use for the color scale
    """


    from matplotlib import ticker, colors
    
    arrows=kwargs.get('arrows',[])

    if not subplot :
        import matplotlib.pyplot as plt
    else :
        plt = subplot
        
    if x_range is None: 
        if subplot: x_range=plt.get_xlim()
        else      : x_range = (np.min(xs),np.max(xs))
    if y_range is None: 
        if subplot: y_range=plt.get_ylim()
        else      : y_range = (np.min(ys),np.max(ys))
    if scalemin is None: scalemin = np.min(arr[arr>0])
    if scalemax is None: scalemax = np.max(arr[arr>0])

    if 'norm' in kwargs: del(kwargs['norm'])

    # Adjust colormap
    if cmap:
        from mmlutils.mmlplot import get_cmaparray
        cmap_colors=get_cmaparray(cmap)
        axbgclr=cmap_colors[0]
        axfgclr=cmap_colors[-1]
        if logscale:
            cmap=colors.ListedColormap(cmap_colors[int(0.01*len(cmap_colors)):])
        cmap.set_under(color=axbgclr)
        cmap.set_over(color=axfgclr)
        cmap.set_bad(color=axbgclr) #no effect 

    levelmin=scalemin/100.
    levelmax=scalemax*100.
    if logscale:
        try:
            levels = np.logspace(np.log10(scalemin),np.log10(scalemax),nlevels)
            cont_color=colors.LogNorm(clip=False)
        except ValueError:
            raise ValueError('crazy contour levels -- try specifying the *levels* keyword or use a linear scale')
            
            return

        if arr.units != NoUnit() and arr.units != 1 :
            cb_label = '$log_{10}('+arr.units.latex()+')$'
        else :
            cb_label = '$log_{10}(N)$'
    else:
        levels = np.linspace(scalemin,scalemax,nlevels)

        cont_color=None
        
        if arr.units != NoUnit() and arr.units != 1 :
            cb_label = '$'+arr.units.latex()+'$'
        else :
            cb_label = '$N$'
    
    if not subplot and clear : plt.clf()

    if logscale: arr=np.ma.masked_less_equal(arr,0).filled(scalemin/100.)

    if ret_im :
        if logscale: arr = np.log10(arr)
        
        return plt.imshow(arr, origin='down', vmin=scalemin, vmax=scalemax,
                          aspect = 'auto',cmap=cmap,
                          #aspect = np.diff(x_range)/np.diff(y_range),cmap=cmap,
                          extent=[x_range[0],x_range[1],y_range[0],y_range[1]])
    # Adjust to log then plot
    if logscale: linarr=np.log10(arr) ; linlvl=np.log10(levels)
    else       : linarr=arr           ; linlvl=levels
    cs = plt.contourf(xs,ys,linarr, linlvl, norm=None,cmap=cmap,extend='both',**kwargs)
    #cs = plt.contourf(xs,ys,arr, levels, norm=cont_color,cmap=cmap,**kwargs)
    plt.gca().set_axis_bgcolor(axbgclr)
    
    if kwargs.has_key('xlabel'):
        xlabel = kwargs['xlabel']
    else :
        try:
            if xlogrange: xlabel=r''+'$log_{10}('+xs.units.latex()+')$'
            else : xlabel = r''+'$x/' + xs.units.latex() +'$'
        except AttributeError:
            xlabel = None

    if xlabel :
        try:
            if subplot :
                plt.set_xlabel(xlabel)
            else:
                plt.xlabel(xlabel)
        except:
            pass

    if kwargs.has_key('ylabel'):
        ylabel = kwargs['ylabel']
    else :
        try:
            if ylogrange: ylabel='$log_{10}('+ys.units.latex()+')$'
            else : ylabel = r''+'$y/' + ys.units.latex() +'$'
        except AttributeError:
            ylabel=None

    if ylabel :
        try:
            if subplot :
                plt.set_ylabel(ylabel)
            else:
                plt.ylabel(ylabel)
        except:
            pass

    
#    if not subplot:
#        plt.xlim((x_range[0],x_range[1]))
#        plt.ylim((y_range[0],y_range[1]))

    if colorbar :
        cbfmt="%.2f" if logscale else "%.2e"
        cb = plt.colorbar(cs, format = cbfmt).set_label(r''+cb_label) 
        #cb = plt.colorbar(cs, format = "%.2e").set_label(r''+cb_label)
        
    if legend : plt.legend(loc=2)


    # Annotate
    if annotate:
        imgsiz = (x_range[1]-x_range[0],y_range[1]-y_range[0])

        # Arrow
        from mmlutils.mmlmath import polar_square
        dx=imgsiz[0]/40.
        dy=imgsiz[1]/40.
        for iarrow in arrows:
            iphi=np.arctan((iarrow[1]/imgsiz[1])/(iarrow[0]/imgsiz[0]))
            idir=polar_square(iphi)
            ix=(idir[0]+0.5)*(imgsiz[0]-dx)+x_range[0]
            iy=(idir[1]+0.5)*(imgsiz[1]-dy)+y_range[0]
            plt.scatter(ix,iy,c=axfgclr,marker=(3,0,-90+iphi*180./np.pi),s=arr.shape[0])
        plt.xlim((x_range[0],x_range[1])) 
        plt.ylim((y_range[0],y_range[1]))

        # Strings
        # imgsiz = (array.SimArray(x_range[1]-x_range[0],units=xs.units),
        #           array.SimArray(y_range[1]-y_range[0],units=ys.units))
        xlocpad=imgsiz[0]/20.
        ylocpad=imgsiz[1]/20.
        for locstr in ['ul','ur','dl','dr']:
            iloc=[0.,0.]
            if   locstr[1]=='l': iloc[0]=x_range[0]+xlocpad ; ha='left'
            elif locstr[1]=='r': iloc[0]=x_range[1]-xlocpad ; ha='right'
            if   locstr[0]=='u': iloc[1]=y_range[1]-ylocpad ; va='top'
            elif locstr[0]=='d': iloc[1]=y_range[0]+ylocpad ; va='bottom'
            plt.text(iloc[0],iloc[1],kwargs.get(locstr+'str',''),color=axfgclr,
                     horizontalalignment=ha, verticalalignment=va)

    # Save
    if (filename): 
        if config['verbose']: print "Saving "+filename
        plt.savefig(filename)



    

                  
def fourier_map(sim, nbins = 100, nmin = 1000, nphi=100, mmin=1, mmax=7, rmax=10, 
                levels = [.01,.05,.1,.2], subplot = None, ret = False, **kwargs) : 
    """

    Plot an overdensity map generated from a Fourier expansion of the
    particle distribution. A :func:`~pynbody.analysis.profile.Profile`
    is made and passed to :func:`~pynbody.plot.util.inv_fourier` to
    obtain an overdensity map. The map is plotted using the usual
    matplotlib.contour. 
    
    **Input**:

    *sim* :  a :func:`~pynbody.snapshot.SimSnap` object

    **Optional Keywords**:
    
    *nbins* (100) : number of radial bins to use for the profile

    *nmin* (1000) : minimum number of particles required per bin 

    *nphi* (100)  : number of azimuthal bins to use for the map

    *mmin* (1)    : lowest multiplicity Fourier component

    *mmax* (7)    : highest multiplicity Fourier component

    *rmax* (10)   : maximum radius to use when generating the profile

    *levels* [0.01,0.05,0.1,0.2] : tuple of levels for plotting contours
    
    *subplot* (None) : Axes object on which to plot the contours
    
    """
    from . import util

    if subplot is None : 
        import matplotlib.pylab as plt

    else : 
        plt = subplot

    p = profile.Profile(sim,max=rmax,nbins=nbins)
    phi,phi_inv = util.inv_fourier(p,nmin,mmin,mmax,nphi)

    rr,pp = np.meshgrid(p['rbins'],phi)

    xx = (rr*np.cos(pp)).T
    yy = (rr*np.sin(pp)).T

    plt.contour(xx,yy,phi_inv,levels,**kwargs)
    
    if ret: 
        return xx,yy,phi_inv


def prob_plot(x,y,weight,nbins=(100,100),extent=None,axes=None,**kwargs) : 
    """ 

    Make a plot of the probability of y given x, i.e. p(y|x). The
    values are normalized such that the integral along each column is
    one.

    **Input**: 

    *x*: primary binning axis

    *y*: secondary binning axis

    *weight*: weights array

    *nbins*: tuple of length 2 specifying the number of bins in each direction

    *extent*: tuple of length 4 speciphysical extent of the axes
     (xmin,xmax,ymin,ymax)

    **Optional Keywords**:

    all optional keywords are passed on to the imshow() command

    """

    import matplotlib.pylab as plt

    assert(len(nbins)==2)
    grid = np.zeros(nbins)

    if extent is None : 
        extent = (min(x),max(x),min(y),max(y))

    xbinedges = np.linspace(extent[0],extent[1],nbins[0]+1)
    ybinedges = np.linspace(extent[2],extent[3],nbins[1]+1)
    

    for i in xrange(nbins[0]) : 
        
        ind = np.where((x > xbinedges[i])&(x < xbinedges[i+1]))[0]
        h, bins = np.histogram(y[ind],weights=weight[ind], bins = ybinedges, density = True)
        grid[:,i] = h

    if axes is None: 
        im = plt.imshow(grid,extent=extent,origin='lower',**kwargs)

    else : 
        im = axes.imshow(grid,extent=extent,origin='lower',**kwargs)

    cb = plt.colorbar(im,format='%.2f')
    cb.set_label(r'$P(y|x)$')
    
    return grid, xbinedges, ybinedges
