#!/usr/bin/python 

import simulo
import nbody
import status

def loadsnap(*args,**kwargs):
    """
    Loads a snapshot using nbody
    """
    return nbody.load(*args,**kwargs)

def loadsim(*args,**kwargs):
    """
    Loads simulation info using simulo
    """
    return simulo.loadsimobj(*args,**kwargs)

def statsim(*args,**kwargs):
    """
    Reports on the status of a run
    """
    return status.writestat(*args,**kwargs)

def allsim(*args,**kwargs):
    """
    Performs method for all simulations
    """
    return simulo.allsim(*args,**kwargs)
        

__all__=['loadsnap','loadsim','statsim','allsim']
