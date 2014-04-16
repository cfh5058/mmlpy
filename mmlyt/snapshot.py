#!/usr/bin/env python

class NbodySnap(object):
    def __init__(self):
        """
        Initialize standard parameters
        """
        self._arrays={}
        self._num_particles=0

    @property
    def filename(self): return self._fname

    def __len__(self): return self._num_particles

    def __repr__(self): 
        if self._fname!="":
            return "<NbodySnap \""+self._filename+"\" len="+str(len(self))+">"
        else :
            return "<NbodySnap len="+str(len(self))+">"

    ################################################################################################################################
    # METHODS FOR GETTING AND SETTING THINGS
    def __getitem__(self,i):
        """
        Returns an array or subset of particles
        """
        # Array
        if   isinstance(i,str): return self._get_array(i)
        # Slice
        elif isinstance(i,slice): pass
        # Family
        # Mask array
        elif isinstance(i,np.ndarray) and np.issubdtype(np.bool,i.dtype):
            if len(i.shape)>1 or i.shape[0]>len(self):
                raise ValueError, "Incorrect shape for masking array"
            else: return self[np.where(i)]
        # List of indexes
        # Single index
        else: raise TypeError

    def __setitem__(self,name,i):
        """
        Assigns an array
        """
        # Get index
        if isinstance(i,(list,tuple)):
            index=name[1]
            name=name[0]
        else: index=None
        # Select array
        if isinstance(i,array.SimArray):
            ax=item
        else:
            ax=np.asanyarray(i).view(array.SimArray)
        # Add to list
        if name not in self.keys():
            # Get dimensions
            try: ndim = len(ax[0])
            except TypeError : ndim = 1
            except IndexError: ndim = ax.shape[-1] if len(ax.shape) > 1 else 1
            # Get data type
            dtype = self._get_preferred_dtype(name)
            if dtype is None : dtype = getattr(item,'dtype',None)
            # Create the array
            self._create_array(name,ndim,dtype=dtype)
        # Add the array
        self._set_array(name,ax,index)

    def __delitem__(self,name):
        """
        Method for removing arrays
        """
        # Only remove family array if it is one
        if name in self._family_arrays :
            assert name not in self._arrays # mustn't have simulation-level array of this name
            del self._family_arrays[name]
        # Otherwise remove it completely
        else: del self._arrays[name]

    def __getattr__(self,name):
        """
        Method for referencing families
        """
        try: return self[family.get_family(name)]
        except ValueError: pass
        raise AttributeError("%r object has no attribute %r"%(type(self).__name__,name))

    def __setattr__(self,name,val):
        """
        Prevent writing to families
        """
        if name in family.family_names(): raise AttributeError, "Cannot assign family name "+name
        return object.__setattr__(self,name,val)

    def __delattr__(self,name):
        """
        Allow for 'persistence' later
        """
        object.__delattr__(self,name)

    ################################################################################################################################
    # METHODS CONTROLING STANDARD DICTIONARY BEHAVIOR
    def keys(self): return self._arrays.keys()
    def has_key(self): return name in self.keys()
    def values(self): return [self[k] for k in self.keys()]
    def items(self): return [(k,self[k]) for k in self.keys()]
    def get(self,key,alternative=None) :
        """Returns self[key] if in self else alternative"""
        try: return self[key]
        except KeyError: return alternative
    def iterkeys(self): #yield k for k in self.keys()
        for k in self.keys(): yield k
    __iter__ = iterkeys
    def itervalues(self): #yield self[k] for k in self:
        for k in self: yield self[k]
    def iteritems(self): #yield (k,self[k]) for k in self
        for k in self: yield (k,self[k])

    ################################################################################################################################
    # METHODS CONTROLING NON-STANDARD DICTIONARY BEHAVIOR
    def has_family_key(self,name): return name in self.family_keys()
    def loadable_keys(self,fam=None): return []
    def familty_keys(self,fam=None):
        """
        Arrays available from a sub-family view
        """
        if fam is not None:
            return [x for x in self._family_arrays if fam in self._family_arrays[x]]
        else:
            return self._family_arrays.keys()

    ################################################################################################################################
    # UNIT MANIPULATION
    def conversion_context(self):
        """
        Returns a dictionary of scale factors
        """
        d = {}
        wanted = ['a','h']
        for x in wanted:
            if x in self.properties:
                d[x] = self.properties[x]
        return d
    def original_units(self):
        """
        Converts arrays to original units
        """
        self.physical_units(distance=self.infer_original_units('km'),
                            velocity=self.infer_original_units('km s^-1'),
                            mass=self.infer_original_units('Msol'), persistent=False)
    def physical_units(self,distance='kpc',velocity='km s^-1',mass='Msol',persistent=True):
        """
        Converts arrays to be consistent with specified units
        """
        dims = [units.Unit(x) for x in distance, velocity, mass, 'a', 'h']
        # Get arrays
        all = self._arrays.values()
        for x in self._family_arrays: all+=self._family_arrays[x].values()
        # Convert arrays
        for ar in all: self._autoconvert_array_unit(ar.ancestor, dims)

        for k in self.properties.keys() :
            v = self.properties[k]
            if isinstance(v, units.UnitBase) :
                try :
                    new_unit = v.dimensional_project(dims)
                except units.UnitsException :
                    continue
                new_unit = reduce(lambda x,y: x*y, [a**b for a,b in zip(dims, new_unit[:3])])
                new_unit*=v.ratio(new_unit, **self.conversion_context())
                self.properties[k] = new_unit

        if persistent :
            self._autoconvert =dims
        else :
            self._autoconvert=None

    def infer_original_units(self, dimensions) :
        """
        Return a units with the same physical dimensions but in the current file units
        """
        dimensions = units.Unit(dimensions)
        d = dimensions.dimensional_project(self._file_units_system+["a","h"])
        new_unit = reduce(lambda x,y: x*y, [a**b for a,b in zip(self._file_units_system, d)])
        return new_unit

    def _default_units_for(self,array_name):
        """
        Construct and return units based on what is known about the array
        """
        array_name = self._array_name_1D_to_ND(array_name) or array_name
        u = units._default_units.get(array_name,None)
        if u is not None :
            u = self.infer_original_units(u)
        return u

    ################################################################################################################################
    # LOADING METHODS
    def _load_array(self,array_name,fam=None):
        """Overwrite this"""
        raise IOError, "No lazy-loading implemented"

    ################################################################################################################################
    # WRITING METHODS
    def write(self,fmt=None,filename=None):
        """
        Method to write files
        """
        if filename is None and "<" in self.filename :
            raise RuntimeError, 'Cannot infer a filename; please provide one (use obj.write(filename="filename"))'
        if fmt is None :
            if not hasattr(self, "_write") :
                raise RuntimeError, 'Cannot infer a file format; please provide one (e.g. use obj.write(filename="filename", fmt=pynbody.tipsy.TipsySnap)'
            self._write(self, filename)
        else :
            fmt._write(self, filename)

def _new(n_particles=0,**ptypes):
    """
    Initializes a blank snapshot
    """
    if len(families)==0: families = {'dm':n_particles}

    tot_particles=0
    for k,v in families.items():
        assert isinstance(v,int)
        tot_particles+=v

    x = NbodySnap()
    x._num_particles = tot_particles
    x._fname = "<created>"

#    x._create_arrays(["pos","vel"],3)
#    x._create_arrays(["mass"],1)

#    rt = 0
#    for k,v in t_fam :
#        x._family_slice[k] = slice(rt,rt+v)
#        rt+=v

#    x._decorate()
    return x
