#!/usr/bin/python
from PIL import Image
import numpy as np

####################################################################################################################################
# METHOD TO CONVERT IMAGE OBJECT TO DATA
def image2array(im):
    """
    Converts PIL image object to array
    """
#    a = numpy.array(im.getdata()).reshape(im.size[0],im.size[1],3)
    a = np.array(im)
    return a

####################################################################################################################################
# METHOD TO CONVERT DATA TO IMAGE OBJECT
def array2image(a):
    """
    Converts array to PIL image object
    """
    return Image.fromarray(a)

####################################################################################################################################
# METHOD TO READ IMAGE DATA FROM FILE
def read_imagdata(fname):
    """
    Returns ndarray of image data
    """
    from PIL import ImageFilter
    # Create image object
    im=Image.open(fname)
    im_grey=im.convert('L')
    im_grey=im_grey.filter(ImageFilter.SMOOTH_MORE)
#    im_grey=im_grey.resize((100,100),ImageFilter.BILINEAR)
    # Convert to array
    imdata=image2array(im_grey)
    # Return array
    return imdata

