#!/usr/bin/env python
import main
import simlist
import simfile
import mmlgadget
import snapshot

def load(fname,ftype=None,*args,**kwargs):
    """
    Loads a file using the appropriate class
    """
    if ftype not in simfile.LIST_SNAPFORM:
        for c in simfile.LIST_SNAPFORM:
            cmod=simfile.ftype2snapobj(c)
            if cmod._can_load(fname): ftype=c
    if ftype in simfile.LIST_SNAPFORM:
        fmod=simfile.ftype2snapobj(ftype)
        return fmod._load(fname,*args,**kwargs)
    else:
        raise IOError('Format not understood or does not exist: %r'%fname)

from snapshot import _new as new

__all__=['load','new']

