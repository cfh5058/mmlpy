"""

projection
==========

Routines for projecting a simulation onto a 2D map

"""

def map(sim,r=201732.223771,phi=,theta):

def ComputeMap(sim,*arg,**kw):
    """
    Return an image in form of a matrix (nx x ny float array) 
    """
    params =    extract_parameters(arg,kw,self.defaultparameters)
    obs         = params['obs']
    x0          = params['x0']
    xp          = params['xp']
    alpha       = params['alpha']
    mode        = params['mode']
    view        = params['view']
    r_obs       = params['r_obs']
    eye         = params['eye']
    dist_eye    = params['dist_eye']
    foc         = params['foc']
    space       = params['space']
    persp       = params['persp']
    clip        = params['clip']
    size        = params['size']
    shape       = params['shape']
    cut         = params['cut']
    frsp        = params['frsp']
    filter_name = params['filter_name']
    filter_opts = params['filter_opts']
