#!/usr/bin/python
import numpy as np
import matplotlib.pyplot as plt
from mmlutils import *

####################################################################################################################################
# CLASS FOR ELLIPSE FITTING
####################################################################################################################################
class IsophotalEllipse:

    def __init__(self,imdata,isolevel=None):
        # Set values
        self._setImdata(imdata)
        if isolevel is None: isolevel = np.mean(self.imdata)
        self._setLevel(isolevel)
        self._fitted = False
        self._a = 1 #major axis
        self._b = 1 #minor axis
        self._phi = 0 #rotation
        self._fitpoints = 100
        self._x0 = self._imdata.shape[0]/2.0
        self._y0 = self._imdata.shape[1]/2.0
        self._fixcen = False
        self._keeppix = .01 #this could be auto-tuned but generally seems to give a good answer

    def _getImdata(self):
        return self._imdata
    def _setImdata(self,val):
        val = np.array(val)
        if len(val.shape) != 2: raise ValueError('2D image data not provided')
        self._imdata = val
        self._fitted = False
    imdata = property(_getImdata,_setImdata)

    def _getFitpoints(self):
        return self._fitpoints
    def _setFitpoints(self):
        self._fitted = False
    nfitpoints = property(_getFitpoints,_setFitpoints)

    def _getLevel(self):
        return self._level
    def _setLevel(self,val):
        self._level = val
        self._fitted = False
    isolevel = property(_getLevel,_setLevel)

    def _getFixcen(self):
        return self._fixcen
    def _setFixcen(self,val):
        self._fixcen = val
        self._fitted = False
    fixcenter = property(_getFixcen,_setFixcen)

    def _getKeeppix(self):
        return self._keeppix
    def _setKeeppix(self,val):
        self._keeppix = val
        self._fitted = False
    keeppix=property(_getKeeppix,_setKeeppix,doc="""
    The pixels to keep for fitting, either as a fraction of the total
    (if less than 1) or as a fixed number of pixels (if greater than 1)
    """)

    def polar(self,npoints):
        """
        returns the full ellipse as an (rho,theta) tuple
        """
        th = np.linspace(0,2*pi,npoints)
        return self._th2r(th),th

    def _th2r(self,th):
        a,b = self._a,self._b
        aterm = a*np.sin(th-self._phi)
        bterm = b*np.cos(th-self._phi)
        r = a*b*(aterm*aterm+bterm*bterm)**-0.5
        return r

    def cartesian(self,npoints):
        """
        returns the full ellipse as an (x,y) tuple
        """
        r,th = self.polar(npoints)
        x = r*np.cos(th)
        y = r*np.sin(th)
        return x+self._x0,y+self._y0

    def _fitEllipse(self):
        """
        fits an ellipse to the specified isophot
        """
        from scipy.optimize import leastsq
        from scipy.ndimage import map_coordinates
        # Get the difference between image and isophot
        diff = self._imdata - self._level
        maxdiff = max(np.max(diff)**2,np.min(diff)**2)
        maxsz = np.min(diff.shape)/2
        # Create grid of image coordinates
#        ximg=np.arange(diff.shape[0])
#        yimg=np.arange(diff.shape[1])
#        Ximg,Yimg = np.meshgrid(ximg,yimg)
        Ximg,Yimg = np.meshgrid(np.arange(diff.shape[0]),np.arange(diff.shape[1]))
#        mask=np.logical_and(self._imdata==0)
        mask_img=self._imdata==0.
        Ximg=np.ma.masked_array(Ximg,mask=mask_img)
        Yimg=np.ma.masked_array(Yimg,mask=mask_img)
        Zimg=np.ma.masked_array(diff,mask=mask_img)

        mask_lin=np.reshape(mask_img,-1)
        Xlin=np.reshape(Ximg,-1)
        Ylin=np.reshape(Yimg,-1)
        Zlin=np.reshape(Zimg,-1)
        sorti = (Zlin**2).argsort(fill_value=100.*Zlin.max())
        nkeep = self._keeppix*Zlin.size if self._keeppix <= 1 else self._keeppix
        nkeep = 0.0001*Zlin.size
        loseidx=sorti[nkeep:]
        keepidx=sorti[:nkeep]

        mask_lin[loseidx]=True
        mask_img=np.reshape(mask_lin,diff.shape)
        Ximg.mask=mask_img
        Yimg.mask=mask_img
        Zimg.mask=mask_img

        xi = Xlin[keepidx]
        yi = Ylin[keepidx]

##         xi = xi[sorti][:nkeep]
##         yi = yi[sorti][:nkeep]
##         keepdiffs = diff.ravel()[sorti][:nkeep]
        
##         ims=plt.imshow(diff)
##         cs=plt.contour(Ximg,Yimg,Zimg,[0.0],colors='m')
##         p = cs.collections[0].get_paths()[0]
##         v = p.vertices
##         xi = v[:,0]
##         yi = v[:,1]


#        plt.imshow(diff)
#        plt.plot(xi,yi,'o')
#        plt.show()

        def fmin(vals,self):
            """
            subfunction to determine fit
            """
            # Set object values based on input
            self._a = np.abs(vals[0])
            self._b = np.abs(vals[1])
            self._phi = vals[2]
            if not self._fixcen:
                self._x0 = vals[3]
                self._y0 = vals[4]
            # Get difference
            xo,yo = xi-self._x0,yi-self._y0
            sep = (xo*xo+yo*yo)**0.5 - self._th2r(np.arctan2(yo,xo))
            # Return difference
            return list(sep)
#            return list(sep*keepdiffs)

        # Perform least squares fit
        v0 = np.array([self._a,self._b,self._phi]) if self._fixcen else np.array([self._a,self._b,self._phi,self._x0,self._y0])
        diag = [1./diff.shape[0],1./diff.shape[1],1.] if self._fixcen else [1./diff.shape[0],1./diff.shape[1],1.,1.,1.]
        soln,cov,infodict,mesg,ier = leastsq(fmin,v0,args=(self,),full_output=True,diag = diag)
        print soln
        # Sort solution variables and handle error
        if not ier in range(1,5): raise Exception('Possible fit problem: [{}] {}'.format(ier,mesg))
        else:
            if not ier==1: print 'Possible fit problem: [{}] {}'.format(ier,mesg)
        if self._fixcen:
            self._a,self._b,self._phi = soln
        else:
            self._a,self._b,self._phi,self._x0,self._y0 = soln
        self.lastier = ier
        self.lastmesg = mesg
        # Relabel major vs. minor axis as necessary
        if self._a < self._b:
            self._a,self._b = self._b,self._a
            self._phi += np.pi/2
        # Handle wrapping of angle (forces range [0,2*pi))
        if self._phi >= 2*np.pi:
            self.phi -= 2*np.pi*np.floor(self.phi/2/np.pi)
        elif self._phi < 0:
            self._phi += 2*np.pi*np.floor(-self.phi/2/np.pi)
        # Change fitted flag to true
        self._fitted = True

    def getDiff(self,fractional=False):
        """
        returns the difference between the fitted ellipse and the flux level
        """
        from scipy.ndimage import map_coordinates
        # Fit if necessary
        if not self._fitted: self._fitEllipse()
        # Calculate the residual
        if fractional:
            res = map_coordinates(self._imdata-self._level,self.cartesian(self._fitpoints),mode='nearest')
        else:
            res = map_coordinates(1-self._level/self._imdata,self.cartesian(self._fitpoints),mode='nearest')
        return res

    @property
    def ecc(self):
        """
        returns eccentricity
        """
        # Fit if necessary
        if not self._fitted: self._fitEllipse()
        # Calculate eccentricity
        ratio = self._b/self._a
        ecc = (1-ratio*ratio)**0.5
        return ecc

    @property
    def ellip(self):
        """
        return ellipticity
        """
        # Get eccentricity
        ecc = self.ecc
        # Calculate ellipticity
        # ellip = (ratio*ratio-1)**(-0.5)
        # ratio = sqrt(1+1/ellip*ellip)
        # b=a/ratio
        ellip = (1-ecc*ecc)**(-0.5)
        return ellip

    @property
    def major(self):
        """
        returns major axis length in image coordinates 
        """
        # Fit if necessary
        if not self._fitted: self._fitEllipse()
        # Return major axis
        return self._a
    a = major

    @property
    def minor(self):
        """
        returns major axis length in image coordinates
        """
        # Fit if necessary
        if not self._fitted: self._fitEllipse()
        # Return minor axis
        return self._b
    b = minor

    @property
    def phi(self):
        """
        returns rotation angle clockwise from x-axis
        """
        # Fit if necessary
        if not self._fitted: self._fitEllipse()
        # Return phi
        return self._phi

    def plot(self,clf=True):
        """
        plots the ellipse
        """
        from matplotlib import pyplot as plt
        # Fit if necessary
        if not self._fitted: self._fitEllipse()
        # CLF?
        if clf: plt.clf()
        # Plot image
        plt.imshow(self.imdata.T)
        # Plot ellipse
        plt.plot(*self.cartesian(self._fitpoints),**{'c':'k'})

        

####################################################################################################################################
# METHOD TO FIT MULTIPLE ELLIPSE TO IMAGE
####################################################################################################################################
class Isophotes:

    def __init__(self,imdata,isolevels=None,nisolevels=None):
        # Pars input
        minlvlDEF=imdata[imdata!=0].min()*100.
        maxlvlDEF=imdata[imdata!=0].max()
        if minlvlDEF == maxlvlDEF:
            raise Exception('Image data is singular, isophotes impossible.')
        if isolevels == None:
            nisolevels = mmlpars.mml_pars(nisolevels,default=5,type=int)
            if minlvlDEF > 0:
                isolevels = np.logspace(np.log10(minlvlDEF),np.log10(maxlvlDEF),nisolevels)
            else:
                print '[mmlphot.isophotes] WARNING: linear spacing forced by data range.'
                isolevels = np.linspace(minlvlDEF,maxlvlDEF,nisolevels)
        # Set variables
        self._setImdata(imdata)
        self._setLevels(isolevels)

    def _getImdata(self):
        return self._imdata
    def _setImdata(self,val):
        val = np.array(val)
        if len(val.shape) != 2: raise ValueError('2D image data not provided')
        self._imdata = val
        self._fitted = False
    imdata = property(_getImdata,_setImdata)

    def _getLevels(self):
        return self._levels
    def _setLevels(self,val):
        self._levels = val
        self._nlevels = len(val)
        self._fitted = False
        self._isophots = [IsophotalEllipse(self.imdata,isolevel=iiso) for iiso in self.isolevels]
    isolevels = property(_getLevels,_setLevels)

    def _getIsophots(self):
        return self._isophots
    isophotes = property(_getIsophots)

    def _fit(self):
        for iphot in self.isophotes:
            print iphot._level
            iphot._fitEllipse()
        self._fitted = True

    @property
    def ecc(self):
        """
        Returns list of eccentricities
        """
        # Fit if neccessary
        if not self._fitted: self._fit()
        # Return eccentricities
        return np.array([iphot.ecc for iphot in self.isophotes])

    @property
    def ellip(self):
        """
        Returns list of ellipticities
        """
        # Fit if neccessary
        if not self._fitted: self._fit()
        # Return ellipticities
        return np.array([iphot.ellip for iphot in self.isophotes])

    @property
    def major(self):
        """
        Returns list of major axes
        """
        # Fit if neccessary
        if not self._fitted: self._fit()
        # Return major axes
        return np.array([iphot.major for iphot in self.isophotes])
    a = major

    @property
    def minor(self):
        """
        Returns list of minor axes
        """
        # Fit if neccessary
        if not self._fitted: self._fit()
        # Return minor axes
        return np.array([iphot.minor for iphot in self.isophotes])
    b = minor

    @property
    def phi(self):
        """
        Returns list of angles clockwise from x-axes
        """
        # Fit if neccessary
        if not self._fitted: self._fit()
        # Return angles
        return np.array([iphot.phi for iphot in self.isophotes])

    def plot(self,clf=True,fname='Isophotes.png'):
        """
        Plots ellipses over the image
        """
        from matplotlib import pyplot as plt
        # Fit if neccessary
        if not self._fitted: self._fit()
        # CLF?
        if clf: plt.clf()
        # Plot image
        plt.imshow(self.imdata.T)
        # Plot ellipses
        for iphot in self._isophots:
            plt.plot(*iphot.cartesian(iphot._fitpoints),**{'c':'k'})
        # Save figure
        plt.savefig(fname)

def IsophotesFromFile(fname,**extra_kw):
    """
    Returns Isophotes class for an image file
    """
    imdata=mmlimage.read_imagdata(fname)
    return Isophotes(imdata,**extra_kw)
    
def test_Isophotes(fname=None,plotname=None):
    """
    Tests the Isophotes class and methods
    """
    # Set constants
    fnameDEF='/home/langmm/utils/python/mine/files/testisophotes_img.png'
    plotnameDEF='/home/langmm/utils/python/mine/files/testisophotes_sol.png'
    # Pars input
    fname=mmlpars.mml_pars(fname,default=fnameDEF,type=str)
    plotname=mmlpars.mml_pars(plotname,default=plotnameDEF,type=str)
    # Load image data
    isoObj=IsophotesFromFile(fname)
    # Plot
    isoObj.plot(fname=plotname)
    # Print data
    print plotname
    print 'ellip: {}'.format(self.ellip)
    print 'major: {}',format(self.major)
    print 'phi:   {}',format(self.phi)
    # Return control
    return isoObj
        

    

