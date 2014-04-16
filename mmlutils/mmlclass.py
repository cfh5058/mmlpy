####################################################################################################################################
#
# MEAGAN LANG'S CLASSES
#
####################################################################################################################################
import numpy as np
import copy

class mydict(dict):
    '''
    A modification of the dictionary class.
    '''

    def __init__(self, *args, **kwargs):
        '''
        Method for initializing the mydict class.
        '''
        self.update(*args, **kwargs)

    def update(self, *args, **kwargs):
        '''
        Method for updating the mydict class with arbitrary keywords.
        '''
        for k, v in dict(*args, **kwargs).iteritems():
            self[k] = v

class mydict_nd(dict):
    '''
    Multidimensional dictionary class
    '''

    def __init__(self,val,*args):
        ndim=len(args)
        ival=copy.deepcopy(val)
        idim=0
        for iarg in reversed(args):
            idim+=1
            if not isinstance(iarg,list): raise Exception('Arguments must be lists ({}).'.format(iarg))
            idict={ikey:copy.deepcopy(ival) for ikey in iarg}
#            idict=dict.fromkeys(iarg,ival)
            if idim==ndim:
                self.update(**idict)
            else:
                ival=copy.deepcopy(idict) ; del idict

        


