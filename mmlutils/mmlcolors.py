#!/usr/bin/python

####################################################################################################################################
# METHOD TO CREATE A SIMPLE LINEAR COLORMAP
####################################################################################################################################
def cmap(clr1,clr2,clrType=None,name=None):
    '''
    NAME:
        mmlcolors.cmap
    PURPOSE:
        To create a simple linear colormap from two colors.
    CALLING:
        cmapout=cmap(clr1,clr2[,clrType=None,name=None])
    ARGUMENTS:
        clr1:    3x1 list of intial color values
        clr2:    3x1 list of final color values
    KEYWORDS:
        clrType: String specifying what type of colors clr1 & clr2 are
                 "rgb" [DEFAULT], "hsv", or "hls"
        name:    String to use to name colormap
    OUTPUT:
        cmapout: New instance of Colormap
    '''
    import mmlpars,numpy
    import matplotlib.colors as colors
    import matplotlib.cm as cm
    # Pars input
    clr1=mmlpars.mml_pars(clr1,type=list)
    clr2=mmlpars.mml_pars(clr2,type=list)
    clrType=mmlpars.mml_pars(clrType,default='rgb',type=str,list=['rgb','hsv','hls'])
    # Handle rgb colors
    if clrType=='rgb':
        # Create dictionary for segmented color map
        cdict={'red':   [(0.0,clr1[0],clr1[0]),
                         (1.0,clr2[0],clr2[0])],
               'green': [(0.0,clr1[1],clr1[1]),
                         (1.0,clr2[1],clr2[1])],
               'blue':  [(0.0,clr1[2],clr1[2]),
                         (1.0,clr2[2],clr2[2])]}
        # Create linear segmented colormap
        cmapout=colors.LinearSegmentedColormap(name,cdict)
    # Handle other color types
    else:
        # Create array in original colorspace
        clrlist_0=numpy.linspace(clr1[0],clr2[0],256)
        clrlist_1=numpy.linspace(clr1[1],clr2[1],256)
        clrlist_2=numpy.linspace(clr1[2],clr2[2],256)
        clrarr=numpy.concatentate(([clrlist_0],[clrlist_1],[clrlist_2]),axis=0).T
        # Convert to rgb space
        if clrType=='hsv':
            clrarr=colors.hsv_to_rgb(clrarr)
        elif clrType=='hls':
            clrarr=colors.hls_to_rgb(clrarr)
        # Create color map from array
        cmapout=colors.ListedColormap(clrarr,name=name)
    # Return colormap
    return cmapout
