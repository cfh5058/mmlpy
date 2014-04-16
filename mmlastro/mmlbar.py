#!/usr/bin/python
from mmlutils import *
import mmlphot
import numpy as np

####################################################################################################################################
# METHOD TO GET BAR PROPERTIES
####################################################################################################################################
def get_barstat(pos,includehist=None,deltaEllip=None,
                range=None,bins=100,normed=None,weights=None,
                isolevels=None,nisolevels=10,
                plotflag=False,plotname=None):
    """
    Returns dictionary of bar properties
    """
    # Pars input
    includehist=mmlpars.mml_pars(includehist,default=False,type=bool)
    deltaEllip=mmlpars.mml_pars(deltaEllip,default=0.1,type=float,min=0.)
    # Create histograms
    H_xy,xedges_xy,yedges_xy=np.histogram2d(pos[:,0],pos[:,1],bins=bins,range=range,
                                            normed=normed,weights=weights)
    H_xz,xedges_xz,yedges_xz=np.histogram2d(pos[:,0],pos[:,2],bins=bins,range=range,
                                            normed=normed,weights=weights)
    # Get conversion factor from image to physical
    imag2phys_xy=(xedges_xy.max()-xedges_xy.min())/float(bins-1)
    imag2phys_xz=(xedges_xz.max()-xedges_xz.min())/float(bins-1)
    # Create isophot objects
    iso_xy=mmlphot.Isophotes(H_xy,isolevels=isolevels,nisolevels=nisolevels)
    iso_xz=mmlphot.Isophotes(H_xz,isolevels=isolevels,nisolevels=nisolevels)
    # Get ellipticity
    deltaEllip_xy=iso_xy.ellip[1:]-iso_xy.ellip[:-1]
    deltaEllip_xz=iso_xz.ellip[1:]-iso_xz.ellip[:-1]
    # Find turning point in ellipticity
    print type(deltaEllip_xy)
    idxout_xy=np.where(deltaEllip_xy >= deltaEllip)[0]
    idxout_xz=np.where(deltaEllip_xz >= deltaEllip)[0]
#    idxout_xy=np.array(np.where(deltaEllip_xy >= deltaEllip))
#    idxout_xz=np.array(np.where(deltaEllip_xz >= deltaEllip))
    if len(idxout_xy)==0: idxend_xy=0
    else: idxend_xy=idxout_xy.min()
    if len(idxout_xz)==0: idxend_xz=0
    else: idxend_xz=idxout_xz.min()
#    idxend_xy is 0 if len(idxout_xy)==0 else idxout_xy.min()
#    idxend_xz is 0 if len(idxout_xz)==0 else idxout_xz.min()
    # Create dictionary of properties
    statdict={
        'ellipmax':iso_xy.ellip.max(),
        'pa':iso_xy.phi[idxend_xy],
        'x' :iso_xy.major[idxend_xy]*imag2phys_xy,
        'y' :iso_xy.minor[idxend_xy]*imag2phys_xy,
        'z' :iso_xz.minor[idxend_xz]*imag2phys_xz
        }
    statdict['r']=np.array([statdict['x'],statdict['y'],statdict['z']]).max()
    if includehist:
        statdict['hist_xy']=H_xy
        statdict['xbin_xy']=xedges_xy
        statdict['ybin_xz']=yedges_xy
        statdict['hist_xz']=H_xz
        statdict['xbin_xz']=xedges_xz
        statdict['ybin_xz']=yedges_xz
    # Plot
    if plotflag: iso_xy.plot(fname=plotname)
    # Return dictionary
    return statdict

if __name__ == '__main__':
    main()
