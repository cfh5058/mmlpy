#!/usr/bin/python
import matplotlib
#matplotlib.use('Agg')
import matplotlib.pyplot as plt
import mmlio,mmlstring,mmlpars,mmlfiles,mmlmath
import os,glob,copy
import numpy as np
#from mmlnbody import mmlgadget

LIST_SCALES=['lin','log']
LIST_PTYPS=['all','visi','gas','halo','disk','bulge','stars','bndry','star','dm']
LIST_PGALS=[0,1,2]

####################################################################################################################################
# METHOD TO ANIMATE FRAMES
####################################################################################################################################
def ffmpeg(frameTemp,animFile,overwrite=False,run=True,ndigits=None,rmtemp=True,verbose=False,tempDir=None,**input_kw):
    # Get frames
    if isinstance(frameTemp,list):
        frameMeth='list'
        frameList=frameTemp
        frameTemp=None
    else:
        frameMeth='temp'
        frameTemp=mmlpars.mml_pars(frameTemp,type=str)
        frameList=mmlfiles.ls(frameTemp)
    # Get number of digits required
    if not isinstance(ndigits,int): ndigits=int(np.ceil(np.log10(float(len(frameList)))))
    # Pars input and handle file overwrite
    anim_form={
        'r'       : mmlpars.parsdict(default=10,type=int,min=1),
        'sameq'   : mmlpars.parsdict(default=True,type=bool   ),
        'loglevel': mmlpars.parsdict(default='quiet',type=str )
        }
#        'b'       : mmlpars.parsdict(default='200k',type=str  ),
    anim_kw=mmlpars.mml_formpars(input_kw,anim_form)
    for ianimkey in anim_kw.keys():
        if ianimkey not in anim_form: del anim_kw[ianimkey]
    # Handle file overwrite
    if not mmlfiles.prep_write(animFile,overwrite): return None
    # Gen file extentsion
    frameExt=os.path.splitext(frameList[0])[1]
    # Get temporary file template
    countform='%0{}d'.format(ndigits)
    if not isinstance(tempDir,str): tempDir=os.path.join(os.path.dirname(animFile),'temp{}'.format(np.random.randint(1000)))
    mmlfiles.mkdirs(tempDir)
    tempFrame=os.path.join(tempDir,'fram_'+countform+frameExt)
    # Create symbolic links
    if   frameMeth=='list':
        for i,f in enumerate(frameList,start=1): os.symlink(f,tempFrame % i)
    elif frameMeth=='temp':
        cmdsym='x=1; for i in {}; do itmpFrame=$(printf {} $x); ln -s "$i" "$itmpFrame"; x=$(($x+1)); done'.format(frameTemp,tempFrame)
    # Create keyword string
    kwstr=''
    for ikey in anim_kw.keys():
        if isinstance(anim_kw[ikey],bool):
            if anim_kw[ikey]:
                kwstr+='-{} '.format(ikey)
        else:
            kwstr+='-{} {} '.format(ikey,anim_kw[ikey])
    # Create command string
    cmdffmpeg='ffmpeg '+kwstr+' -i {} {}'.format(tempFrame,animFile)
    # Execute command
    if run:
        # Create symbolic links
        if frameMeth=='temp':
            symout=os.system(cmdsym)
            if symout!=0: mmlio.verbose('WARNING: There may have been a problem during the creation of symbolic links to frames. [out={}]'.format(symout))
        # Create animation
        aniout=os.system(cmdffmpeg)
        if aniout!=0: mmlio.verbose('WARNING: There may have been a problem running ffmpeg. [out={}]'.format(aniout))
        if rmtemp:
            for tf in mmlfiles.ls(os.path.join(tempDir,'*')): os.unlink(tf)
            os.rmdir(tempDir)
        # Check file
        if not os.path.isfile(animFile): raise Exception('Animation failed.')
        else:
            if verbose: mmlio.verbose('Created animation:')
            if verbose: print '    '+animFile
    # Return command string
    return cmdffmpeg


####################################################################################################################################
####################################################################################################################################
# METHODS TO PARS & SET PLOT OPTIONS

####################################################################################################################################
# METHOD TO PARS PLOT OPTIONS
def parsopt(inopt=None,outextra=False):
    inopt=mmlpars.mml_pars(inopt,default={},type=dict)
    xboxDEF=dict(pad=105, alpha=0.0)
    yboxDEF=dict(pad=105, alpha=0.0)
    form={
        'tickfontsize': mmlpars.parsdict(default=9.    ,type=float,min=0.             ),
        'textfontsize': mmlpars.parsdict(default=10.   ,type=float,min=0.             ),
        'figsize'     : mmlpars.parsdict(default=(5.,5.)    ,type=tuple,nelements=2   ),
        'axsize'      : mmlpars.parsdict(default=(2.,2.)    ,type=tuple,nelements=2   ),
        'limpad'      : mmlpars.parsdict(default=0.1        ,type=float               ),
        'bgclr'       : mmlpars.parsdict(default=(0,0,0)    ,type=tuple,nelements=3   ),
        'fgclr'       : mmlpars.parsdict(default=(1,1,1)    ,type=tuple,nelements=3   ),
        'pad'         : mmlpars.parsdict(default=1.0        ,type=float,min=0.0       ),
        'wpad'        : mmlpars.parsdict(default=0.5        ,type=float,min=0.0       ),
        'hpad'        : mmlpars.parsdict(default=0.5        ,type=float,min=0.0       ),
        'plotdir'     : mmlpars.parsdict(default='~/utils/python/plotdir',type=str    ),
        'plotext'     : mmlpars.parsdict(default='png'      ,type=str                 ),
        'animext'     : mmlpars.parsdict(default='mp4'      ,type=str                 ),
        'xlabpos'     : mmlpars.parsdict(default=(0.5,-0.06),type=tuple,nelements=2   ),
        'ylabpos'     : mmlpars.parsdict(default=(-0.06,0.5),type=tuple,nelements=2   ),
        'xlabbox'     : mmlpars.parsdict(default=xboxDEF    ,type=dict                ),
        'ylabbox'     : mmlpars.parsdict(default=yboxDEF    ,type=dict                ),
        'xunits'      : mmlpars.parsdict(default=''         ,type=str                 ),
        'yunits'      : mmlpars.parsdict(default=''         ,type=str                 ),
        'title'       : mmlpars.parsdict(default=''         ,type=str                 ),
        'outplot'     : mmlpars.parsdict(default=False      ,type=bool                ),
        'showplot'    : mmlpars.parsdict(default=False      ,type=bool                ),
        'label'       : mmlpars.parsdict(default=''         ,type=str                 )
        }
    outopt0=mmlpars.mml_formpars(inopt,form)
    if outextra:
        outopt={} ; extrakw={}
        for ikey in outopt0:
            if ikey in form:
                outopt[ikey]=outopt0[ikey]
            else:
                extrakw[ikey]=outopt0[ikey]
        return outopt,extrakw
    else:
        return outopt0

####################################################################################################################################
# METHOD TO RETURN LIST OF TEXT OBJECTS
def get_textobj(axesObj):
    """
    Returns a list of text objects associated with an axes object
    """
    textList=[]
    textList+=[axesObj.get_title()]
    textList+=[axesObj.get_xlabel(),axesObj.get_ylabel()]
    textList+=axesObj.get_xticklabels(minor=True)
    textList+=axesObj.get_yticklabels(minor=True)
    for ichild in axesObj.get_children():
        if isinstance(ichild,matplotlib.text.Text):
            textList.append(ichild)
    return textList

####################################################################################################################################
# METHOD TO RETURN LIST OF OBJECTS
def get_objlist(axesObj):
    """
    Returns a list of objects associated with an axes object
    """
    objList={'text':[],'spines':[]}
    for ichild in axesObj.get_children():
        if isinstance(ichild,matplotlib.text.Text):
            objList['text'].append(ichild)
        elif isinstance(ichild,matplotlib.spines.Spine):
            objList['spines'].append(ichild)
    return objList

####################################################################################################################################
# METHOD TO SET FIGURE PROPERTIES
def set_figprop(figObj,options=None):
    """
    Sets the properties of a figure object
    """
    options=parsopt(options)
    axesList=[]
    for ichild in figObj.get_children():
        if isinstance(ichild,matplotlib.axes.Axes):
            axesList.append(ichild)
        if isinstance(ichild,matplotlib.text.Text):
            ichild.set_color(options['fgclr'])
    figObj.set_edgecolor(options['fgclr'])
    figObj.set_facecolor(options['bgclr'])
    figObj.set_size_inches(options['figsize'])
    set_axesprop(axesList,options=options)
#    print len(axesList)
#    figObj.tight_layout(pad=options['pad'],h_pad=options['hpad'],w_pad=options['wpad'])

####################################################################################################################################
# METHOD TO SET AXES PROPERTIES
def set_axesprop(axesList,options=None):
    """
    Sets the properties of a list of axes objects
    """
    options=parsopt(options)
    for iax in axesList:
        iax.set_axis_bgcolor(options['bgclr'])
        iax.tick_params(colors=options['fgclr'],labelsize=options['tickfontsize'])
        objList=get_objlist(iax)
        for ispine in objList['spines']:
            ispine.set_color(options['fgclr'])
        for itext in objList['text']:
            itext.set_color(options['fgclr'])
            itext.set_fontsize(options['textfontsize'])
        xlab=iax.get_xlabel()
        ylab=iax.get_ylabel()
        iax.set_xlabel(xlab,color=options['fgclr'],size=options['textfontsize'],
                       bbox=options['xlabbox'],transform=iax.transAxes)
        iax.set_ylabel(ylab,color=options['fgclr'],size=options['textfontsize'],
                       bbox=options['ylabbox'],transform=iax.transAxes)
        iax.xaxis.set_label_coords(options['xlabpos'][0],options['xlabpos'][1],transform=iax.transAxes)
        iax.yaxis.set_label_coords(options['ylabpos'][0],options['ylabpos'][1],transform=iax.transAxes)
    return
        
####################################################################################################################################
# METHOD TO SET TEXT PROPERTIES
def set_textprop(textList,options=None):
    """
    Sets the properties of a list of text objects
    """
    options=parsopt(options)
    for itx in textList:
        itx.set_color(options['fgclr'])
    return

####################################################################################################################################
####################################################################################################################################
# COLORMAP METHODS & CLASSES

####################################################################################################################################
# METHOD TO VIEW COLORMAPS
def show_cmaps(cmaptype=None,verbose=None,inclrev=None):
    """
    Plots colormaps
    """
    # Set constants
    fnamebase='/home/langmm/utils/python/mine/files/colormaps_'
    cmaptypeLIST=list_cmaptypes()
    # Pars input
    cmaptype=mmlpars.mml_pars(cmaptype,default='all',list=['all']+cmaptypeLIST)
    verbose=mmlpars.mml_pars(verbose,default=True,type=bool)
    inclrev=mmlpars.mml_pars(inclrev,default=False,type=bool)
    fname=fnamebase+cmaptype+'.png'
    if verbose: print '[mmlplot.show_cmaps] Beginning to plot colormaps...'
    # Create type list
    if cmaptype=='all': cmaplist=cmaptypeLIST
    else              : cmaplist=[cmaptype]
    # Get list of colormaps
    mapdict={}
    ntyps=len(cmaplist)
    nmaps=0
    nmaplist=[]
    for icmaptype in cmaplist:
        mapdict[icmaptype]=list_cmaps(cmaptype=icmaptype,inclrev=inclrev)
        nmaplist.append(len(mapdict[icmaptype]))
        nmaps+=len(mapdict[icmaptype])
    nmapsmax=max(nmaplist)
    # Set plot size
    cbhgt=0.15
    opad=cbhgt
    cbwid=6.
    lbwid=3.
    figsiz=(cbhgt*nmaps+3*opad,(cbwid+lbwid)*ntyps+2*opad)
    opad_relx=opad/figsiz[0]
    opad_rely=opad/figsiz[1]
    lbwid_rely=lbwid/figsiz[1]
    lbwid_relcb=lbwid/cbwid
    # Initialize plot
    matplotlib.rc('text', usetex=False)
    a=np.outer(np.arange(0,1,0.01),np.ones(10))
    fig=plt.figure(figsize=figsiz)
    fig.subplots_adjust(hspace=lbwid_relcb,
                        bottom=opad_rely,top=1.-opad_rely-lbwid_rely,
                        left=2*opad_relx,right=1.-opad_relx)
    # Plot colormaps
#    print ntyps,nmapsmax,ntyps*nmapsmax
    for ntypprev,icmaptype in enumerate(cmaplist):
#        print icmaptype
        nmapprev=nmapsmax*ntypprev
        # Get list of colormaps
        imaplist=sorted(mapdict[icmaptype])
        inmap=len(imaplist)
        # Loop over colormaps plotting
        for i, m in enumerate(imaplist):
#            print m
            iax=plt.subplot(ntyps,nmapsmax,1+i+nmapprev)
            iax.spines['top'   ].set_color('none')
            iax.spines['bottom'].set_color('none')
            iax.spines['left'  ].set_color('none')
            iax.spines['right' ].set_color('none')
            iax.axes.get_xaxis().set_ticks([])
            iax.axes.get_yaxis().set_ticks([])
#            plt.axis("off")
            plt.imshow(a,aspect='auto',cmap=get_cmap(m,cmaptype=icmaptype),origin="lower")
            plt.title(m,rotation=90,fontsize=10,va='bottom')
            if i==0: plt.ylabel(icmaptype,fontsize=10)
    # Save figure
#    fig.tight_layout()
    plt.savefig(fname,dpi=100,facecolor='white')
    if verbose: print '    {}'.format(fname)
    return

####################################################################################################################################
# METHOD TO RETURN LIST OF VALID COLORMAP TYPES
def list_cmaptypes():
    """
    Returns a list of supported colormap types
    """
    typeLIST=['matplotlib','pnbody','mmlmonograd','mmldualgrad','mmlparttype','mmlpgal','mmlptyp']
    typeLIST+=['mmlparttype_sk'+skew for skew in ['0p25','0p5','0p75']]
    return typeLIST

####################################################################################################################################
# METHOD TO RETURN LIST OF MMLLIGHT CMAPS
def list_mmllight():
    """
    Returns a list of supported mmllight cmaps
    """
    nlgt=6
    typeLIST=['mmllight{}'.format(ilgt) for ilgt in range(nlgt)]
    return typeLIST

####################################################################################################################################
# METHOD TO RETURN LIST MML CMAPS
def list_mmlcmapskew():
    """
    Returns a list of supported cmap skews
    """
    skewLIST=[0.25,0.5,0.75]
    return skewLIST
def list_mmlmonograd():
    """
    Returns a list of supported 1 color mml gradient cmaps
    """
    cmapLIST0=list_mmlcolors()
    cmapLISTs=[]
    for iskew in list_mmlcmapskew():
        cmapLISTs+=[icmap+'_sk{}'.format(mmlstring.dec2str(iskew,trimzero=False)) for icmap in cmapLIST0]
    cmapLISTr=[icmap+'_r' for icmap in cmapLISTs]
    cmapLIST=cmapLISTs+cmapLISTr
    return cmapLIST
def list_mmldualgrad():
    """
    Returns a list of supported 2 color mmlgradient cmaps
    """
    cmapLIST1=list_mmlcolors()
    cmapLIST2=[]
    for idx,icmap1 in enumerate(cmapLIST1):
        for icmap2 in cmapLIST1:
#            if icmap1 == icmap2: continue
#            if icmap1 in icmap2: continue
#            if icmap2 in icmap1: continue
            cmapLIST2.append(icmap1+'_'+icmap2)
    cmapLISTs=[]
    for iskew in list_mmlcmapskew():
        cmapLISTs+=[icmap+'_sk{}'.format(mmlstring.dec2str(iskew,trimzero=False)) for icmap in cmapLIST2]
    cmapLISTr=[icmap+'_r' for icmap in cmapLISTs]
    cmapLIST=cmapLISTs+cmapLISTr
    return cmapLIST
def list_mmlparttype():
    """
    Returns a list of default color maps for SPH simulations
    """
    ngal=3
    cmapLIST0=LIST_PTYPS+['star']+['dm']
    cmapLIST=[]
    for igal in range(ngal):
        cmapLIST+=[icmap+str(igal) for icmap in cmapLIST0]
    return cmapLIST
def list_mmlpgal():
    """
    Returns a list of default color maps for galaxy colors
    """
    cmapLIST=LIST_PGALS
    return cmapLIST
def list_mmlptyp():
    """
    Returns a list of default color maps for particle type colors
    """
    cmapLIST=LIST_PTYPS
    return cmapLIST
def list_mmlcolors():
    """
    Returns a list of defined mml color strings
    """
    colorLIST=['cyan' ,'green' ,'blue' ,'magenta' ,'red' ,'yellow' ,
               'cyan2','green2','blue2','magenta2','red2','yellow2',
               'grey']
    return colorLIST
def get_mmlcolor(cstr):
    """
    Returns an rgb tuple for a corresponding color string
    """
    cstr=mmlpars.mml_pars(cstr,list=list_mmlcolors())
    if   cstr=='cyan'    : rgbclr=(0.0,1.0,1.0)
    elif cstr=='cyan2'   : rgbclr=(0.0,0.6,1.0)
    elif cstr=='green'   : rgbclr=(0.0,1.0,0.0)
    elif cstr=='green2'  : rgbclr=(0.6,1.0,0.0)
    elif cstr=='blue'    : rgbclr=(0.0,0.0,1.0)
    elif cstr=='blue2'   : rgbclr=(0.0,0.6,1.0)
    elif cstr=='magenta' : rgbclr=(1.0,0.0,1.0)
    elif cstr=='magenta2': rgbclr=(1.0,0.0,0.6)
    elif cstr=='red'     : rgbclr=(1.0,0.0,0.0)
    elif cstr=='red2'    : rgbclr=(1.0,0.6,0.0)
    elif cstr=='yellow'  : rgbclr=(1.0,1.0,0.0)
    elif cstr=='yellow2' : rgbclr=(0.8,1.0,0.0)
    elif cstr=='grey'    : rgbclr=(0.5,0.5,0.5)
    else: raise Exception('[mmlplot.get_mmlcolor] Invalid color string: {}'.format(cstr))
    return rgbclr

####################################################################################################################################
# METHOD TO RETURN LIST OF VALID COLORMAPS
def list_cmaps(cmaptype=None,inclrev=None):
    """
    Returns a list of colormaps
    """
    # Set constants
    cmaptypeLIST=list_cmaptypes()
    # Pars input
    cmaptype=mmlpars.mml_pars(cmaptype,default='all',list=['all']+cmaptypeLIST)
    inclrev=mmlpars.mml_pars(inclrev,default=True,type=bool)
    # Get list of colormaps
    if cmaptype=='all': cmaptypes=cmaptypeLIST
    else              : cmaptypes=[cmaptype]
    maps=[]
    for icmaptype in cmaptypes:
        if cmaptype=='all': pfix=icmaptype+'_'
        else              : pfix=''
        if icmaptype=='matplotlib':
            maps+=[pfix+m for m in matplotlib.cm.datad]
            maps+=[m for m in matplotlib.cm.datad]
        elif icmaptype=='pnbody':
            from pNbody import parameters
#            import pNbody.parameters as pNparm
            flist=glob.glob(os.path.join(parameters.PALETTEDIR,'*'))
            maps+=[pfix+os.path.basename(ifile) for ifile in flist]
            maps+=[pfix+ilight for ilight in list_mmllight()]
        elif icmaptype=='mmlmonograd':
            maps+=[pfix+imono for imono in list_mmlmonograd()]
        elif icmaptype=='mmldualgrad':
            maps+=[pfix+idual for idual in list_mmldualgrad()]
        elif icmaptype.startswith('mmlparttype'):
            maps+=[pfix+ipart for ipart in list_mmlparttype()]
        elif icmaptype=='mmlptyp':
            maps+=[pfix+ipart for ipart in list_mmlptyp()]
        elif icmaptype=='mmlpgal':
            maps+=[pfix+ipart for ipart in list_mmlptyp()]
        else: raise Exception('[mmlplot.list_cmaps] Invalid colormap type: {}'.format(icmaptype))
    # Remove reverses
    if not inclrev:
        maps_nor=copy.deepcopy(maps)
        for imap in maps:
            if imap.endswith("_r"): maps_nor.remove(imap)
        maps=maps_nor
    # Check length
    if len(maps)==0: print '[mmlplot.list_cmaps] WARNING: no colormaps returned for {}'.format(cmaptype)
    # Return list
    return maps

####################################################################################################################################
# RETURN DICTIONARY SPECIFYING TYPE TO ALPHA MAPPING
def get_typalpha(outstr=False):
    typelist=LIST_PTYPS
    typedict={typelist[ityp]:ityp for ityp in range(len(typelist))}
    alphadict_str={'all'  :0.01  ,
                   'visi' :0.75  ,
                   'star' :0.75  ,
                   'gas'  :0.5   ,
                   'halo' :0.25  ,
                   'disk' :1.0   ,
                   'bulge':0.5   ,
                   'stars':0.25  ,
                   'bndry':0.25  }
    alphadict_int={typedict[ikey]:alphadict_str[ikey] for ikey in typelist}
    if outstr: return alphadict_str
    else     : return alphadict_int

####################################################################################################################################
# RETURN DICTIONARY SPECIFYING TYPE TO COLOR MAPPING
def get_typcolor(outstr=False):
    typelist=LIST_PTYPS
    typedict={typelist[ityp]:ityp for ityp in range(len(typelist))}
    colordict_str={'all'  :'grey'  ,
                   'visi' :'blue'  ,
                   'star' :'blue'  ,
                   'gas'  :'red2'  ,
                   'halo' :'green' ,
                   'dm'   :'green' ,
                   'disk' :'cyan'  ,
                   'bulge':'red'   ,
                   'stars':'yellow',
                   'bndry':'magenta'}
    colordict_int={typedict[ikey]:colordict_str[ikey] for ikey in typelist}
    if outstr: return colordict_str
    else     : return colordict_int

####################################################################################################################################
# RETURN DICTIONARY SPECIFYING GALAXY TO ALPHA MAPPING
def get_galalpha():
    alphadict={0:1.0,1:1.0,2:0.5}
    return alphadict

####################################################################################################################################
# RETURN DICTIONARY SPECIFYING GALAXY TO COLOR MAPPING
def get_galcolor():
    colordict={0:'green',1:'blue',2:'magenta'}
    return colordict

####################################################################################################################################
# METHOD TO RETURN COLORMAP
def get_cmap(m0,cmaptype=None,reverse=None):
    """
    Returns a given colormap
    """
    # Set constants
    fnamebase='/home/langmm/utils/python/mine/files/colormaps_'
    cmaptypeLIST=list_cmaptypes()
    # Pars input
    m0=mmlpars.mml_pars(m0,type=str)
    cmaptype=mmlpars.mml_pars(cmaptype,default=cmaptypeLIST[0],list=cmaptypeLIST)
    reverse=mmlpars.mml_pars(reverse,default=False,type=bool)
    # Handle type in name
    if m0[-2:]=='_r':
        reverse=not reverse
        m=m0[:-2]
    else:
        m=m0
    for ictype in cmaptypeLIST:
        if ictype+'_' in m0:
            cmaptype=ictype
            m=m0.split(ictype+'_')[-1]
            break
    m=mmlpars.mml_pars(m,list=list_cmaps(cmaptype))
    # Get cmap
    if cmaptype=='matplotlib':
        if reverse:
            cmap=matplotlib.cm.get_cmap(m+'_r')
        else:
            cmap=matplotlib.cm.get_cmap(m)
    elif cmaptype=='pnbody':
        import pNbody.palette as palette
        if m in list_mmllight():
            pobj=palette.Palette(name='light')
        else:
            pobj=palette.Palette(name=m)
#        print '[mmlplot.get_cmap] max rgba: {}'.format((pobj.r.max(),pobj.g.max(),pobj.b.max(),pobj.t.max()))
        if reverse:
            rlist=pobj.r[::-1]/256.
            glist=pobj.g[::-1]/256.
            blist=pobj.b[::-1]/256.
            tlist=pobj.t[::-1]/256.
        else:
            rlist=pobj.r/256.
            glist=pobj.g/256.
            blist=pobj.b/256.
            tlist=pobj.t/256.
        if m in list_mmllight():
            if   m=='mmllight0': pass
            elif m=='mmllight1':
                rlist=pobj.b/256.
                glist=pobj.r/256.
                blist=pobj.g/256.
            elif m=='mmllight2':
                rlist=pobj.g/256.
                glist=pobj.b/256.
                blist=pobj.r/256.
            elif m=='mmllight3':
                rlist=pobj.g/256.
                glist=pobj.r/256.
            elif m=='mmllight4':
                rlist=pobj.b/256.
                blist=pobj.r/256.
            elif m=='mmllight5':
                glist=pobj.b/256.
                blist=pobj.g/256.
        if reverse:
            cmap=matplotlib.colors.ListedColormap(name='pnbody_'+m+'_r',colors=np.column_stack((rlist,glist,blist,tlist)))
        else:
            cmap=matplotlib.colors.ListedColormap(name='pnbody_'+m,colors=np.column_stack((rlist,glist,blist,tlist)))
    elif cmaptype.startswith('mmlparttype'):
        if '_sk' in cmaptype: skew=cmaptype.split('_')[-1]
        else                : skew='sk0p25'
        colorgalDICT=get_galcolor()
        colortypDICT=get_typcolor(outstr=True)
        colortyp=colortypDICT[m[:-1]]
        colorgal=colorgalDICT[int(float(m[-1]))]
        try:
            cmap=get_cmap('{}_{}_{}'.format(colorgal,colortyp,skew),cmaptype='mmldualgrad')
        except:
            mmlio.yorn('{}_{}_{}'.format(colorgal,colortyp,skew))
            raise
    elif cmaptype=='mmlptyp': cmap=get_cmap('{}_sk0p25'.format(get_typcolor(outstr=True)[m]),cmaptype='mmlmonograd')
    elif cmaptype=='mmlpgal': cmap=get_cmap('{}_sk0p25'.format(get_galcolor()[int(float(m))]),cmaptype='mmlmonograd')
    elif cmaptype in ['mmlmonograd','mmldualgrad']:
        # Set skew
        if '_sk' in m:
            msplt=m.split('_sk')
            mb=msplt[0]
            fmid=mmlstring.str2dec(msplt[-1])
        else:
            mb=m
            fmid=list_mmlcmapskew()[0]
        # Reverse color order if necessary
        if reverse:
            cbeg=(1.0,1.0,1.0)
            cend=(0.0,0.0,0.0)
            fmid=1.0-fmid
        else:
            cbeg=(0.0,0.0,0.0)
            cend=(1.0,1.0,1.0)
        # Get color range
        if cmaptype=='mmlmonograd':
            cmid=get_mmlcolor(mb)
            cdict={'red'  :[(0.0 ,cbeg[0],cbeg[0]),
                            (fmid,cmid[0],cmid[0]),
                            (1.0 ,cend[0],cend[0])],
                   'green':[(0.0 ,cbeg[1],cbeg[1]),
                            (fmid,cmid[1],cmid[1]),
                            (1.0 ,cend[1],cend[1])],
                   'blue' :[(0.0 ,cbeg[2],cbeg[2]),
                            (fmid,cmid[2],cmid[2]),
                            (1.0 ,cend[2],cend[2])]}
        elif cmaptype=='mmldualgrad':
            cstr1,cstr2=mb.split('_')
            if reverse:
                cmid10=get_mmlcolor(cstr2)
                cmid20=get_mmlcolor(cstr1)
            else:
                cmid10=get_mmlcolor(cstr1)
                cmid20=get_mmlcolor(cstr2)
            cmid1=tuple((np.array(cmid10)+np.array(cbeg))/2.)
            cmid2=tuple((np.array(cmid20)+np.array(cend))/2.)
            cmid=tuple((np.array(cmid1)+np.array(cmid2))/2.)
            fbin=min(fmid,1.-fmid)/2.
            fmid1=fmid-fbin
            fmid2=fmid+fbin
            cdict={'red'  :[( 0.0  , cbeg[0] , cbeg[0]  ),
                            ( fmid1, cmid1[0], cmid1[0] ),
                            ( fmid , cmid[0] , cmid[0]  ),
                            ( fmid2, cmid2[0], cmid2[0] ),
                            ( 1.0  , cend[0] , cend[0]  )],
                   'green':[( 0.0  , cbeg[1] , cbeg[1]  ),
                            ( fmid1, cmid1[1], cmid1[1] ),
                            ( fmid , cmid[1] , cmid[1]  ),
                            ( fmid2, cmid2[1], cmid2[1] ),
                            ( 1.0  , cend[1] , cend[1]  )],
                   'blue' :[( 0.0  , cbeg[2] , cbeg[2]  ),
                            ( fmid1, cmid1[2], cmid1[2] ),
                            ( fmid , cmid[2] , cmid[2]  ),
                            ( fmid2, cmid2[2], cmid2[2] ),
                            ( 1.0  , cend[2] , cend[2]  )]}
        if reverse:
            cmap=matplotlib.colors.LinearSegmentedColormap('mml_'+m+'_r',cdict)
        else:
            cmap=matplotlib.colors.LinearSegmentedColormap('mml_'+m,cdict)
    else: raise Exception('[mmlplot.get_cmap] Invalid colormap type: {}'.format(cmaptype))
    # Return colormap
    return cmap

####################################################################################################################################
# METHOD TO RETURN COLORMAP ARRAY
def get_cmaparray(cmap):
    """
    Returns an array of cmap rgb values
    """
    nmin=0. ; nmax=256.
    x=np.arange(nmin,nmax)
    norm=matplotlib.colors.Normalize()
    norm.autoscale(x)
    smap=matplotlib.cm.ScalarMappable(norm=norm,cmap=cmap)
    rgba=smap.to_rgba(x)
    return rgba

####################################################################################################################################
# CONVERT CMAP TO PALETTE
def cmap2palette(cmap):
    """
    Converts a colormap to a PIL compatible palette
    """
    rgba=get_cmaparray(cmap)
    rgb=255.*rgba[:,:3]
    pal=[]
    for idx in range(rgb.shape[0]):
        pal.extend(tuple(rgb[idx,:].flatten().astype('int')))
##     print 'rgb min,max: {},{}'.format(rgb.min(),rgb.max())
##     pal=''
##     for idx in range(rgb.shape[0]):
##         ir=int(rgb[idx,0])
##         ig=int(rgb[idx,1])
##         ib=int(rgb[idx,2])
##         pal=pal+chr(ir)+chr(ig)+chr(ib)
##     pal=(255.*np.array(map(lambda x: cmap(x)[0:3], np.linspace(0., 1.,256))).ravel()).astype('int')
    return pal

####################################################################################################################################
# CONVERT ARRAY TO RGBA IMAGE
def array2rgb(mat,cmap=None,cmapname=None,**mapkw):
    """
    Maps an array to RGB image
    """
    if cmap==None:
        cmapname=mmlpars.mml_pars(cmapname,default='matplotlib_jet',list=list_cmaps())
        cmap=get_cmap(cmapname)
    x=array2map(mat,**mapkw)
    norm=matplotlib.colors.Normalize()
    norm.autoscale(x)
    smap=matplotlib.cm.ScalarMappable(norm=norm,cmap=cmap)
    rgba=smap.to_rgba(x)
    return rgba
    
####################################################################################################################################
# CONVERT ARRAY TO INTEGER MAPPING
def array2map(mat,**input):
    """
    Maps an array onto an integer array corresponding to a colormap
    """
    return map_array(mat,**input)

####################################################################################################################################
# METHOD TO MAP ARRAY TO A COLORMAP
def map_array(mat,scale=None,mn=None,mx=None,cd=None,clip=None,clrmin=None,clrmax=None):
    """
    Maps an array onto an integer array corresponding to a colormap
    """
    # Pars input
    rm=np.ravel(mat).astype('float')
    scale=mmlpars.mml_pars(scale,default='log',list=['linear','log','lin'])
    if scale=='lin': scale='linear'
    mn=mmlpars.mml_pars(mn,default=float(min(rm)),type=float)
    mx=mmlpars.mml_pars(mx,default=float(max(rm)),type=float,min=mn)
    if mn==mx:
        mn=min(rm)
        mx=max(rm)
    cd=mmlpars.mml_pars(cd,default=float(rm.mean()),type=float)
    if cd==0: cd=float(rm.mean())
    if scale=='linear': cd=0.
    clip=mmlpars.mml_pars(clip,default=False,type=bool)
    clrmin=float(mmlpars.mml_pars(clrmin,default=0,type=int,range=[0,254]))
    clrmax=float(mmlpars.mml_pars(clrmax,default=255,type=int,range=[clrmin,255]))
    # Handle negative on logscale
    if scale=='log':
        mat=np.ma.masked_less_equal(mat,0.)
        if mn<0.: mn=mat.min()
        if mx<mn: mx=mat.max()
        mat=mat.filled(mn)
    # Clip array
    if clip: mat=mat.clip(mn,mx)
    # Scale array
    if scale=='log':
        matscl = clrmin+(clrmax-clrmin)*np.log(1.+(mat-mn)/(cd)) / np.log(1.+(mx-mn)/(cd))
    elif scale=='linear':
        matscl = clrmin+(clrmax-clrmin)*(mat-mn)/(mx-mn)
    else: raise Exception('[mmlplot.map_array] Invalid scale: {}'.format(scale))
    # Clip array
#    if clip: matscl=matscl.clip(0,255)
    # Return scaled array
    return matscl.astype('uint8')

####################################################################################################################################
# METHOD TO PLOT IMAGE
def hist2imshow(hist,xbins=None,ybins=None,extent=None,axs=None,**map_kw):
    """
    Method to plot histogram data as image
    """
    verbose=False
    # Pars input
    nx,ny=hist.shape
    if not isinstance(extent,tuple):
        if isinstance(xbins,np.ndarray) and (len(xbins)==nx or len(xbins)==(nx+1)):
            xmin=xbins.min() ; xmax=xbins.max()
        else:
            xmin=float(0.)   ; xmax=float(nx)
        if isinstance(ybins,np.ndarray) and (len(ybins)==ny or len(ybins)==(ny+1)):
            ymin=ybins.min() ; ymax=ybins.max()
        else:
            ymin=float(0.)   ; ymax=float(ny)
        extent=(xmin,xmax,ymin,ymax)
    # Transpose array
    imgdat=np.transpose(hist)
    if verbose:
        print imgdat.min()
        print np.mean(imgdat)
        print imgdat.max()
        mmlio.yorn('Continue?')
    # Scale array
    imgrgb=array2rgb(imgdat,**map_kw)
    # Plot array with first element in bottom left
    if axs is None: axs=plt
    img=axs.imshow(imgrgb,extent=extent,origin='lower')
    # Return image object
    return img

####################################################################################################################################
# Method TO PLOT CONTOUR
def hist2contour(hist,xbins=None,ybins=None,extent=None,ncont=None,levels=None,
                 mn=None,mx=None,scale=None,axs=None,**map_kw):
    """
    Method to plot histogram data as image
    """
    # Pars input
    nx,ny=hist.shape
    if not isinstance(extent,tuple):
        if isinstance(xbins,np.ndarray) and (len(xbins)==nx or len(xbins)==(nx+1)):
            xmin=xbins.min() ; xmax=xbins.max()
        else:
            xmin=float(0.)   ; xmax=float(nx)
        if isinstance(ybins,np.ndarray) and (len(ybins)==ny or len(ybins)==(ny+1)):
            ymin=ybins.min() ; ymax=ybins.max()
        else:
            ymin=float(0.)   ; ymax=float(ny)
        extent=(xmin,xmax,ymin,ymax)
    ncont=mmlpars.mml_pars(ncont,default=11,type=int,min=1)
    rm=np.ravel(hist).astype('float')
    scale=mmlpars.mml_pars(scale,default='log',list=['linear','log'])
    mn=mmlpars.mml_pars(mn,default=float(min(rm)),type=float)
    mx=mmlpars.mml_pars(mx,default=float(max(rm)),type=float,min=mn)
    if mn==mx:
        mn=min(rm)
        mx=max(rm)
    map_kw=dict(map_kw,mn=mn,mx=mx,scale=scale)
    # Transpose array
#    imgdat=hist
    imgdat=np.transpose(hist)
    # Get levels
    if levels==None:
        if   scale=='linear': levels=np.linspace(mn,mx,ncont)
        elif scale=='log'   : levels=mmlmath.logspace(mn,mx,ncont,minord=None,addzero=(mn*mx<0))
        else: raise Exception('Invalid scale method: {}'.format(scale))
    # Get colors
#    mmlio.yorn(str(levels))
    lvlclrs=array2rgb(np.arange(len(levels)))
    # Plot array with first element in bottom left
    if axs is None: axs=plt
    ctr=axs.contour(imgdat,levels=levels,colors=lvlclrs,extent=extent,origin='lower')
    plt.clabel(ctr,inline=1,fontsize=10)
    # Return image object
    return ctr

####################################################################################################################################
# CLASSES
####################################################################################################################################
import matplotlib.cm as mplibcm


class SymLogNorm(matplotlib.colors.Normalize):
    """
    The symmetrical logarithmic scale is logarithmic in both the
    positive and negative directions from the origin.

    Since the values close to zero tend toward infinity, there is a
    need to have a range around zero that is linear.  The parameter
    *linthresh* allows the user to specify the size of this range
    (-*linthresh*, *linthresh*).
    """
    def __init__(self,  linthresh, linscale=1.0,
                 vmin=None, vmax=None, clip=False):
        """
        *linthresh*:
        The range within which the plot is linear (to
        avoid having the plot go to infinity around zero).

        *linscale*:
        This allows the linear range (-*linthresh* to *linthresh*)
        to be stretched relative to the logarithmic range.  Its
        value is the number of decades to use for each half of the
        linear range.  For example, when *linscale* == 1.0 (the
        default), the space used for the positive and negative
        halves of the linear range will be equal to one decade in
        the logarithmic range. Defaults to 1.
        """
        matplotlib.colors.Normalize.__init__(self, vmin, vmax, clip)
        self.linthresh = float(linthresh)
        self._linscale_adj = (linscale / (1.0 - np.e ** -1))

    def __call__(self, value, clip=None):
        if clip is None:
            clip = self.clip

        result, is_scalar = self.process_value(value)
        self.autoscale_None(result)
        vmin, vmax = self.vmin, self.vmax

        if vmin > vmax:
            raise ValueError("minvalue must be less than or equal to maxvalue")
        elif vmin == vmax:
            result.fill(0)
        else:
            if clip:
                mask = np.ma.getmask(result)
                result = np.ma.array(np.clip(result.filled(vmax), vmin, vmax),
                                     mask=mask)
            # in-place equivalent of above can be much faster
            resdat = self._transform(result.data)
            resdat -= self._lower
            resdat /= (self._upper - self._lower)

        if is_scalar:
            result = result[0]
        return result

    def _transform(self, a):
        """
        Inplace transformation.
        """
        masked = np.abs(a) > self.linthresh
        sign = np.sign(a[masked])
        log = (self._linscale_adj + np.log(np.abs(a[masked]) / self.linthresh))
        log *= sign * self.linthresh
        a[masked] = log
        a[~masked] *= self._linscale_adj
        return a

    def _inv_transform(self, a):
        """
        Inverse inplace Transformation.
        """
        masked = np.abs(a) > (self.linthresh * self._linscale_adj)
        sign = np.sign(a[masked])
        exp = np.exp(sign * a[masked] / self.linthresh - self._linscale_adj)
        exp *= sign * self.linthresh
        a[masked] = exp
        a[~masked] /= self._linscale_adj
        return a

    def _transform_vmin_vmax(self):
        """
        Calculates vmin and vmax in the transformed system.
        """
        vmin, vmax = self.vmin, self.vmax
        arr = np.array([vmax, vmin]).astype(np.float)
        self._upper, self._lower = self._transform(arr)

    def inverse(self, value):
        if not self.scaled():
            raise ValueError("Not invertible until scaled")
        else:
            self._transform_vmin_vmax()
        val = np.ma.asarray(value)
        val = val * (self._upper - self._lower) + self._lower
        return self._inv_transform(val)

    def autoscale(self, A):
        """
        Set *vmin*, *vmax* to min, max of *A*.
        """
        self.vmin = np.ma.min(A)
        self.vmax = np.ma.max(A)
        self._transform_vmin_vmax()

    def autoscale_None(self, A):
        """ autoscale only None-valued vmin or vmax """
        if self.vmin is not None and self.vmax is not None:
            pass
        if self.vmin is None:
            self.vmin = np.ma.min(A)
        if self.vmax is None:
            self.vmax = np.ma.max(A)
        self._transform_vmin_vmax()

class VectorMappable(mplibcm.ScalarMappable):
    """
    This is a mixin class to support vector -> RGBA/HSLA mapping.  Handles
    normalization and colormapping
    """
    
    def __init__(self, norm=None, cmap=None):
        """
        *norm* is an instance of :class:`mmlplot.Normalize` or one of
        its subclasses, used to map luminance vector to unitvectors. *cmap* is a
        :mod:`cm` colormap instance, for example :data:`cm.jet`
        """

        self.callbacksSM = cbook.CallbackRegistry()

        if cmap is None: cmap = get_cmap()
        if norm is None: norm = Normalize()

        self._A = None
        self.norm = norm
        self.cmap = get_cmap(cmap)
        self.colorbar = None
        self.update_dict = {'array':False}

    def set_colorbar(self, im, ax):
        'set the colorbar image and axes associated with mappable'
        self.colorbar = im, ax

    def to_rgba(self, x, alpha=None, bytes=False):
        '''Return a normalized rgba array corresponding to *x*. If *x*
        is already an rgb array, insert *alpha*; if it is already
        rgba, return it unchanged. If *bytes* is True, return rgba as
        4 uint8s instead of 4 floats.
        '''
        if alpha is None:
            _alpha = 1.0    
        else:
            _alpha = alpha
        try:
            if x.ndim == 3:
                if x.shape[2] == 3:
                    if x.dtype == np.uint8:
                        _alpha = np.array(_alpha*255, np.uint8)
                    m, n = x.shape[:2]
                    xx = np.empty(shape=(m,n,4), dtype = x.dtype)
                    xx[:,:,:3] = x
                    xx[:,:,3] = _alpha  
                elif x.shape[2] == 4:
                    xx = x
                else:
                    raise ValueError("third dimension must be 3 or 4")
                if bytes and xx.dtype != np.uint8:
                    xx = (xx * 255).astype(np.uint8)
                return xx
        except AttributeError:
            pass    
        x = ma.asarray(x)
        x = self.norm(x)
        x = self.cmap(x, alpha=alpha, bytes=bytes)
        return x  

#class get_cmap():
#
#import matplotlib.colors as mplibcolors
#class Normalize(mplibcolors.Normalize):
