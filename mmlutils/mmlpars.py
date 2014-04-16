####################################################################################################################################
#
# MEAGAN LANG'S PARSING METHODS
#
####################################################################################################################################
import mmlclass,copy,os

class parsdict(mmlclass.mydict):

    def __init__(self,default=None,type=None,nelements=None,ndim=None,shape=None,
                 min=None,max=None,range=None,list=None,ispath=None,isdir=None,isfile=None,
                 keylist=None):
        '''
        Method to initialize the parsdict class.
        '''
        mmlclass.mydict.__init__(self,default=default,type=type,nelements=nelements,ndim=ndim,shape=shape,
                                 min=min,max=max,range=range,list=list,ispath=ispath,isdir=isdir,isfile=isfile,
                                 keylist=None)

    def pars(self,inval):
        '''
        Method to use a parsdict object to pars an input variable.
        '''
        import numpy,collections

        # Set default if inval undefined
        if inval == None:
            if self['default'] == None:
                raise Exception('Variable not set default undefined.')
            else:
                outval=copy.deepcopy(self['default'])
        else:
            outval=inval

        # Check variable class
        if self['type'] != None:
            if isinstance(self['type'],(list,tuple)): typList=self['type']
            else                                    : typList=[self['type']]
            validtype=False
            for itype in typList:
                if isinstance(outval,itype): validtype=True
            if not validtype:
                raise Exception('Variable is not correct type {}, but is type {}.'.format(self['type'],type(outval)))
                
        # Check # of elements
        if self['nelements'] != None:
            if not isinstance(self['nelements'],int):
                raise Exception('NELEMENT fields must be an integer.')
            if len(numpy.array(outval)) != self['nelements']:
                raise Exception('Variable does not have correct number of elements.')

        # Check number of dimensions of variable
        if self['ndim'] != None:
            if not isinstance(self['ndim'],int):
                raise Exception('NDIM field must be an integer.')
            if not isinstance(outval,numpy.ndarray):
                raise Exception('NDIM field is only valid for arrays.')
            if len(outval.shape) != self['ndim']:
                raise Exception('Variable does not have the specified number of dimensions.')

        # Check shape of variable
        if self['shape'] != None:
            if not isinstance(self['shape'],tuple):
                raise Exception('SHAPE field must be a tuple.')
            if numpy.array(outval).shape != self['shape']:
                raise Exception('Variable needs to have shape {}, not {}.'.format(self['shape'],numpy.array(outval).shape))

        # Check range
        if self['range'] != None:
            if not isinstance(self['range'],list):
                raise Exception('RANGE field must be a list.')
            self['min']=self['range'][0]
            self['max']=self['range'][1]
        if self['min'] != None:
            if not isinstance(self['min'],(float,int,long)):
                raise Exception('MIN field must be a float or int.')
            if isinstance(outval,collections.Iterable):
                if not all(outval >= self['min']):
                    raise Exception('One or more elements in variable is smaller than specified minimum.')
            else:
                if not outval >= self['min']:
                    raise Exception('One or more elements in variable is smaller than specified minimum.')
        if self['max'] != None:
            if not isinstance(self['max'],(float,int,long)):
                raise Exception('MAX field must be a float or int')
            if isinstance(outval,collections.Iterable):
                if not all(outval <= self['max']):
                    raise Exception('One or more elements in variable is larger than specified maximum.')
            else:
                if not outval <= self['max']:
                    raise Exception('One or more elements in variable is larger than specified maximum.')

        # Check list
        if self['list'] != None:
            if not isinstance(self['list'],list):
                raise Exception('LIST field must be a list')
            if isinstance(outval,collections.Iterable) and not isinstance(outval,str):
                for ix in outval:
                    if not ix in self['list']:
                        raise Exception('One or more elements in variable is not in specified list.')
            else:
                if not outval in self['list']:
                    raise Exception('One or more elements in variable is not in specified list.')

        # Check if variable is a path
        if self['ispath'] != None:
            if not isinstance(self['ispath'],bool):
                raise Exception('ISPATH field must be a bool')
            if not isinstance(outval,str):
                raise Exception('Paths must be strings')
            if not self['ispath']==isinstance(outval,str):
                raise Exception('Variable test for path is {}'.format(isinstance(outval,str)))

        # Check if variable is a file
        if self['isfile'] != None:
            if not isinstance(self['isfile'],bool):
                raise Exception('ISFILE field must be a bool')
            if not isinstance(outval,str):
                raise Exception('Files must be strings')
            if not self['isfile']==os.path.isfile(outval):
                raise Exception('Variable test for file is {}'.format(os.path.isfile(outval)))

        # Check if variable is a directory
        if self['isdir'] != None:
            if not isinstance(self['isdir'],bool):
                raise Exception('ISDIR field must be a bool')
            if not isinstance(outval,str):
                raise Exception('Directories must be strings')
            if not self['isdir']==os.path.isdir(outval):
                raise Exception('Variable test for directory is {}'.format(os.path.isdir(outval)))

        # Check if variable has correct keys
        if self['keylist'] != None:
            if not isinstance(self['keylist'],list):
                raise Exception('KEYLIST field must be a list')
            if not isinstance(outval,dict):
                raise Exception('KEYLIST field only valid for dictionaries, not type {}'.format(type(outval)))
            for ikey in self['keylist']:
                if ikey not in outval:
                    raise Exception('Variable does not have required key: {}'.format(ikey))
            
        return outval

def mml_pars(inval,*args, **kwargs):
    pd=parsdict(*args, **kwargs)
    outval=pd.pars(inval)

    return outval

def mml_formpars(indict,formdict):
    import mmlclass
    if indict == None:
        outdict=mmlclass.mydict()
    else:
        outdict=indict
    for k,v in formdict.iteritems():
        if k in outdict:
            outdict[k]=v.pars(outdict[k])
        else:
            outdict[k]=v.pars(None)

    return outdict
