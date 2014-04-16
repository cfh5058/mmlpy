"""

gadget
======

Implements classes and functions for handling gadget files; you rarely
need to access this module directly as it will be invoked
automatically via pynbody.load.

"""


from . import snapshot,array, units
from . import family, galaxy
from . import config
from . import config_parser
from . import util

import ConfigParser

import numpy as np
#Needed to unpack things
import struct
import sys
import copy
import os
import warnings
import errno

#This is set here and not in a config file because too many things break 
#if it is not 6
N_TYPE = 6

_type_map = {}
for x in family.family_names() :
    try :
        pp =  [int(q) for q in config_parser.get('gadget-type-mapping',x).split(",")]
        qq=np.array(pp)
        if (qq >= N_TYPE).any() or (qq < 0).any() :
            raise ValueError,"Type specified for family "+x+" is out of bounds ("+pp+")." 
        _type_map[family.get_family(x)] = pp
    except ConfigParser.NoOptionError :
        pass
_rev_type_map = {}
for f,t in _type_map.iteritems() :
    for it in t:
        if it in _rev_type_map: _rev_type_map[it].append(f)
        else                  : _rev_type_map[it]=[f]

_name_map, _rev_name_map = util.setup_name_maps('gadget-name-mapping', gadget_blocks=True)
_translate_array_name = util.name_map_function(_name_map, _rev_name_map)

_sphblock_names = config_parser.get('gadget-1-blocks','sph-blocks').split(',')
_sphblock_names = [q.upper().ljust(4) for q in _sphblock_names]


def gadget_type(fam) :
    if fam == None:
        return list(np.arange(0,N_TYPE))
    else :
        return _type_map[fam]

def gadget_typelen(sim,typ):
    Ntyp=gadget_typeidx(sim,typ).sum()
    return Ntyp

def gadget_typeidx(sim,typ):
    famlist = set(_rev_type_map[typ])
    tidx=np.ones(sim._num_particles,dtype=bool)
    typlist=set([])
    # Loop over families associate with the type
    for f in famlist:
        # Determine intersection of all these families
        itidx=np.zeros(sim._num_particles,dtype=bool)
        if f in sim._family_slice:
            itidx[sim._family_slice[f]]=True
        tidx=np.logical_and(tidx,itidx)
        # Find intersecting types that could still be counted
        if len(typlist)==0: typlist=set(_type_map[f])
        else              : typlist=typlist.intersection(_type_map[f])
    # Remove the indices of particles of intervening types
    typlist.remove(typ)
    for t in typlist:
        if len(famlist.symmetric_difference(_rev_type_map[t]))==0:
            raise Exception('types {} and {} have identical families!'.format(typ,t))
        else: tidx[gadget_typeidx(sim,t)]=False
    # Return indices of particles in the type
    return tidx

def test_typeidx(sim):
    for t in _rev_type_map.keys():
        idx=gadget_typeidx(sim,t)
        print t,sim[idx].families(),sim.header.npart[t],idx.sum()

class GadgetBlock(object) :
    """Class to describe each block.
    Each block has a start, a length, and a length-per-particle"""
    def __init__(self, gfile, name='    ', **exkw):
    # def __init__(self, start=0, length=0,partlen=0,dtype=np.float32,p_types=np.zeros(N_TYPE,bool),name=None) :
        keylist=['start','partlen','dtype','length','npart']
        #File containing block
        self.file=gfile
        #Name of block
        self._name=name[0:4]
        #Other attributes that can be set at read in
        for k in keylist:
            if k in exkw:
                if k=='dtype': setattr(self,'_data_type',exkw[k])
                else         : setattr(self,'_'+k       ,exkw[k])

    def get_start(self,p_type=-1):
        """Returns the starting postion for a given particle type in the block"""
        start=self.start
        if not self.name=='HEAD': start+=(self.get_spart(p_type=p_type)*self.partlen)
        return start

    @property
    def name(self): return self._name
    @property
    def start(self): return self._start
    @property
    def npart(self):
        """Number of particles of each type in this block"""
        if not hasattr(self,'_npart'):
            #Set particle number based on header
            fnpart=self.file.header.npart
            self._npart=np.zeros(N_TYPE,dtype=fnpart.dtype)
            # Handle mass and gas differently
            if self.name=='MASS': 
                self._npart=fnpart*(self.file.header.mass==0)
            elif self.name in _sphblock_names: 
                sphtyp=_type_map[family.get_family('gas')]
                self._npart[sphtyp]=fnpart[sphtyp]
            else:
                self._npart=fnpart
        return self._npart
    @property
    def partlen(self):
        """Bytes per particle"""
        if self.name=='HEAD': raise Exception('HEAD block does not contain per particle info')
        if not hasattr(self,'_partlen'):
            self._partlen=0 if self.tpart==0 else int(self.length/self.tpart)
        return self._partlen
    @property
    def data_type(self):
        """Block data type"""
        if self.name=='HEAD': raise Exception('HEAD block does not contain per particle info')
        if not hasattr(self,'_data_type'):
            #Set the partlen, using better heuristic
            #3D float/double array
            if self.name in ['POS ','VEL ','ACCE']:
                if   self.partlen==24: self._data_type=np.float64
                elif self.partlen==12: self._data_type=np.float32
                elif self.partlen==0 : self._data_type=np.float32
            #1D int/long array
            elif self.name == 'ID  ':
                if   self.partlen==8 : self._data_type=np.int64
                elif self.partlen==4 : self._data_type=np.int32
                elif self.partlen==0 : self._data_type=np.int32
            #1D float/double array
            else :
                if   self.partlen==8 : self._data_type=np.float64
                elif self.partlen==4 : self._data_type=np.float32
                elif self.partlen==0 : self._data_type=np.float32
            if not hasattr(self,'_data_type'):
                raise Exception('No data type found for {} block. (plen={})'.format(self.name,self.partlen))
        return self._data_type
    @property
    def end(self):
        if not hasattr(self,'_end'): self._end=self.start+self.length
        return self._end
    @property
    def tpart(self): return self.get_tpart()
    def get_tpart(self,p_type=-1):
        """Total number of particles of a given type in the block"""
        if self.name=='HEAD': raise Exception('HEAD block does not contain per particle info')
        if p_type==-1: _tpart=self.npart.sum()
        else         : _tpart=self.npart[p_type]
        return _tpart
    @property
    def spart(self): return self.get_spart()
    def get_spart(self,p_type=-1):
        """Number of particle before a particle type in the block"""
        if self.name=='HEAD': raise Exception('HEAD block does not contain per particle info')
        if p_type==-1: _spart=0
        else         : _spart=self.npart[:p_type],sum()
        return _spart
    @property
    def length(self):
        """Length of block in bytes"""
        if not hasattr(self,'_length'): 
            if self.name=='HEAD': self._length=256
            else                : self._length=self.tpart*self.partlen
        return self._length
    @length.setter
    def length(self,value): setattr(self,'_length',value)
    @property
    def p_types(self):
        """Boolean array specifying which particle types are in this block"""
        if self.name=='HEAD': raise Exception('HEAD block does not contain per particle info')
        if not hasattr(self,'_p_types'): self._p_types=(self.npart>0)
        return self._p_types
    @property
    def dim(self):
        """Block dimension"""
        if self.name=='HEAD': raise Exception('HEAD block does not contain per particle info')
        return self.partlen/(np.dtype(self.data_type).itemsize)

    def rw_block(self, rwid, fd=None, filename=None, seek=False):
        """Reads/writes blocks"""
        # Open file
        closefile=False
        if fd is None:
            closefile=True
            seek=True
            if filename is None: filename=self.file._filename
            # Try for existing, if not creat it
            try: fd = open(filename,"r+")
            except IOError as (err, strerror):
                if err == errno.ENOENT: fd = open(filename, "w+")
        # Read/write header
        self.rw_head(rwid,fd,seek=seek)
        # Read
        if rwid in ['R','r','read']:
            # Record position
            self._start = fd.tell()
            # Stop if empty block
            if self.length==0: return
            # Handle data
            if self.name=='HEAD': self.data=self.rw_data(rwid,fd,seek=seek)
            else                : fd.seek(self.length,1)
        # Write
        elif rwid in ['W','w','write']: self.rw_data(rwid,fd,seek=seek)
        # Error
        else: raise Exception('Invalid rwid: {}'.format(rwid))
        # Read/write footer
        self.rw_foot(rwid,fd,seek=seek)
        if closefile: fd.close()

    def rw_data(self,rwid,fd,p_type=-1,data=None,seek=False):
        """
        Read/Write a full block of data
            p_type=-1 for all types
            KeyError if particle type not present
            ValueError if too many particles
        """
        if seek: fd.seek(self.get_start(p_type=p_type),0)
        if p_type==-1: plist=range(N_TYPE)
        else         : plist=[p_type]
        # Read
        if rwid in ['R','r','read']:
            # Header
            if self.name=='HEAD':
                rawdata=fd.read(self.length)
                if len(rawdata)!=self.length:
                    raise IOError, "Could not read {} block in {}".format(self.name,self.file._filename)
                data=_construct_gadget_header(rawdata,self.file.endian)
            # Arrays
            else:
                data=[None]*N_TYPE
                # Loop over types
                for p in plist:
                    n_type=self.get_tpart(p_type=p)*self.dim
                    if n_type==0: continue
                    data[p]=np.fromfile(fd,dtype=self.data_type,count=n_type,sep='')
                    if self.file.endian!='=': data[p]=data[p].byteswap(True)
            # Return data
            return data
        # Write
        elif rwid in ['W','w','write']:
            if not data: data=self.data
            # Header
            if self.name=='HEAD': 
                #Construct new header with the passed header and overwrite npart with the file header.
                #This has ref. semantics so use copy  
                head=copy.deepcopy(data)
                head.npart=np.array(self.file.header.npart)
                #Write header body
                fd.write(head.serialize())
                #Seek 48 bytes forward, to skip the padding (which may contain extra data) 
                fd.seek(48,1)
            # Arrays
            else:
                # Check for correct amount of data
                if np.size(data) > self.get_tpart(p_type)*self.dim:
                    raise ValueError, "{}: Space for {} particles of type {} in file {}, but {} requested".format(self.name,self.tpart,p_type,self.file._filename,np.shape(data)[0])
                # Check type
                dt = np.dtype(self.data_type)
                bt = data.dtype
                if dt.kind!=bt.kind: 
                    raise ValueError, "{}: Data of incorrect type ({}) passed for block of type {}".format(self.name,bt,dt)
                # Check if bytes should be swapped
                if self.file.endian!='=': data=data.byteswap(False)
                # Write data
                fd.write(np.ravel(data.astype(dt)).tostring())
        # Error
        else: raise ValueError, "Invalid rwid {}".format(rwid)

    def rw_head(self,rwid,fd,retdata=False,seek=False):
        """Reads/writes block header"""
        if self.file.format2: fmt=self.file.endian+'I4sIII' ; nbytes=5*4
        else                : fmt=self.file.endian+'I'      ; nbytes=4
        blkheadsize = 4 + 4*1;#1 int and 4 chars
        headroom = 2 * 4;#Relative location of next block; the extra 2 uints are for storing the headers in format2
        # Read
        if rwid in ['R','r','read']:
            data=fd.read(nbytes)
            #If we have run out of file, we don't want an exception,
            #we just want a zero length empty block
            if len(data)!=nbytes:
                #print '[rw_head] Empty block ',len(data),nbytes
                self._name='    '
                self._length=0
                return
            head=struct.unpack(fmt,data)
            if self.file.format2:
                if head[0]!=blkheadsize or head[3]!=blkheadsize or head[4]!=head[2]-headroom :
                    raise IOError, "Corrupt header record. Possibly incorrect file format"
                #Don't include the two "record_size" indicators in the total length count
                self._name=head[1]
                self._length=head[4]
            else:
                self._length=head[0]
                try:
                    self._name = self.file.block_names[0]
                    self.file.block_names = self.file.block_names[1:]
                except IndexError:
                    if self.file.extra == 0:
                        warnings.warn("Run out of block names in the config file. Using fallbacks: UNK*",RuntimeWarning)
                    self._name = "UNK"+str(self.file.extra)
                    self.file.extra+=1
            # Ensure header is the right size
            if self.name=='HEAD' and self.length != 256:
                raise IOError, "Mis-sized HEAD block in "+self.file._filename
        # Write
        elif rwid in ['W','w','write']:
            if self.file.format2: head=[blkheadsize,self.name,self.length+headroom,blkheadsize,self.length]
            else            : head=[self.length]
            data=struct.pack(fmt,*head)
            if retdata: return data
            if seek: fd.seek(self.start-len(data),0)
            fd.write(data)
        # Error
        else: raise Exception('Invalid rwid: {}'.format(rwid))

    def rw_foot(self,rwid,fd,retdata=False,seek=False):
        """Reads/writes block footer"""
        fmt=self.file.endian+'I' ; nbytes=4
        # Read
        if rwid in ['R','r','read']:
            data = fd.read(nbytes)
            if len(data) != nbytes: raise IOError, "Could not read {} block footer".format(self.name)
            (blocksize,)=struct.unpack(fmt,data)
            if blocksize!=self.length:
                raise IOError, "{} block footer ({}) does not match header ({})".format(self.name,blocksize,self.length)
        # Write
        elif rwid in ['W','w','write']:
            data=struct.pack(fmt,self.length)
            if retdata: return data
            else      : fd.write(data)
        # Error
        else: raise Exception('Invalid rwid: {}'.format(rwid))


def _output_order_gadget(all_keys,incldregs=True) :

    out = []
    out_dregs = copy.copy(all_keys)
    for X in map(str.strip,config_parser.get('gadget-default-output', 'field-ordering').split(',')) :
        if X in out_dregs :
            del out_dregs[out_dregs.index(X)]
            out.append(X)
        else: print 'Skipping gadget field: {}'.format(X)
    if incldregs:
        return out+out_dregs
    else:
        return out
    
def _construct_gadget_header(data,endian='=') :
    """Create a GadgetHeader from a byte range read from a file."""
    npart = np.zeros(N_TYPE, dtype=np.uint32)
    mass = np.zeros(N_TYPE)
    time = 0.
    redshift = 0.
    npartTotal=np.zeros(N_TYPE,dtype=np.int32)
    num_files=0
    BoxSize=0.
    Omega0=0.
    OmegaLambda=0.
    HubbleParam=0.
    NallHW=np.zeros(N_TYPE,dtype=np.int32)
    if data == '':
        return
    fmt= endian+"IIIIIIddddddddiiIIIIIIiiddddiiIIIIIIiiif48s"
    if struct.calcsize(fmt) != 256:
        raise Exception, "There is a bug in gadget.py; the header format string is not 256 bytes"
    (npart[0], npart[1],npart[2],npart[3],npart[4],npart[5],
    mass[0], mass[1],mass[2],mass[3],mass[4],mass[5],
    time, redshift,  flag_sfr, flag_feedback,
    npartTotal[0], npartTotal[1],npartTotal[2],npartTotal[3],npartTotal[4],npartTotal[5],
    flag_cooling, num_files, BoxSize, Omega0, OmegaLambda, HubbleParam,flag_stellarage, flag_metals,
    NallHW[0], NallHW[1],NallHW[2],NallHW[3],NallHW[4],NallHW[5],
    flag_entropy_instead_u, flag_doubleprecision, flag_ic_info, lpt_scalingfactor,fill) = struct.unpack(fmt, data)

    header=GadgetHeader(npart,mass,time,redshift,BoxSize,Omega0,OmegaLambda,HubbleParam,num_files)
    header.flag_sfr=flag_sfr
    header.flag_feedback=flag_feedback
    header.npartTotal=npartTotal
    header.flag_cooling=flag_cooling
    header.flag_stellarage=flag_stellarage
    header.flag_metals=flag_metals
    header.NallHW=NallHW
    header.flag_entropy_instead_u=flag_entropy_instead_u       
    header.flag_doubleprecision=flag_doubleprecision
    header.flag_ic_info=flag_ic_info
    header.lpt_scalingfactor=lpt_scalingfactor
    header.endian=endian

    return header
    
class GadgetHeader(object) :
    """Describes the header of gadget class files; this is all our metadata, so we are going to store it inline"""
    def __init__ (self,npart, mass, time, redshift, BoxSize,Omega0, OmegaLambda, HubbleParam, num_files=1 ) :
        "Construct a header from values, instead of a datastring."""
        assert(len(mass) == 6)
        assert(len(npart) == 6)
        # Mass of each particle type in this file. If zero,
        # particle mass stored in snapshot.
        self.mass = mass
        # Time of snapshot
        self.time = time
        # Redshift of snapshot
        self.redshift = redshift
        # Boolean to test the presence of star formation
        self.flag_sfr=False
        # Boolean to test the presence of feedback
        self.flag_feedback=False
        # Boolean to test the presence of cooling
        self.flag_cooling=False
        # Number of files expected in this snapshot
        self.num_files=num_files
        # Box size of the simulation
        self.BoxSize=BoxSize
        # Omega_Matter. Note this is Omega_DM + Omega_Baryons
        self.Omega0=Omega0
        # Dark energy density
        self.OmegaLambda=OmegaLambda
        # Hubble parameter, in units where it is around 70.
        self.HubbleParam=HubbleParam
        # Boolean to test whether stars have an age
        self.flag_stellarage=False
        # Boolean to test the presence of metals
        self.flag_metals=False
        self.flag_entropy_instead_u=False      # flags that IC-file contains entropy instead of u
        self.flag_doubleprecision=False  # flags that snapshot contains double-precision instead of single precision
        self.flag_ic_info=False
        # flag to inform whether IC files are generated with Zeldovich approximation,
        # or whether they contain 2nd order lagrangian perturbation theory ICs.
        #    FLAG_ZELDOVICH_ICS     (1)   - IC file based on Zeldovich
        #    FLAG_SECOND_ORDER_ICS  (2)   - Special IC-file containing 2lpt masses
        #    FLAG_EVOLVED_ZELDOVICH (3)   - snapshot evolved from Zeldovich ICs
        #    FLAG_EVOLVED_2LPT      (4)   - snapshot evolved from 2lpt ICs
        #    FLAG_NORMALICS_2LPT    (5)   - standard gadget file format with 2lpt ICs
        # All other values, including 0 are interpreted as "don't know" for backwards compatability.
        self.lpt_scalingfactor=0.    # scaling factor for 2lpt initial conditions  
        self.endian=""
        #Number of particles
        self.npart = np.array(npart,dtype=np.uint32)
        if (npart < 2**31).all() :
            # First 32-bits of total number of particles in the simulation
            self.npartTotal=np.array(npart,dtype=np.int32)
            # Long word of the total number of particles in the simulation.
            # At least one version of N-GenICs sets this to something entirely different.
            self.NallHW=np.zeros(N_TYPE,dtype=np.int32)
        else :
            self.header.NallHW = np.array(npart/2**32,dtype=np.int32)
            self.header.npartTotal = np.array(npart - 2**32*self.header.NallHW,dtype=np.int32)

    
    def serialize(self) :
        """This takes the header structure and returns it as a packed string"""
        fmt= self.endian+"IIIIIIddddddddiiIIIIIIiiddddiiIIIIIIiiif"
        #Do not attempt to include padding in the serialised data; the most common use of serialise 
        #is to write to a file and we don't want to overwrite extra data that might be present
        if struct.calcsize(fmt) != 256-48:
            raise Exception, "There is a bug in gadget.py; the header format string is not 256 bytes"
        #WARNING: On at least python 2.6.3 and numpy 1.3.0 on windows, castless code fails with:
        #SystemError: ..\Objects\longobject.c:336: bad argument to internal function
        #This is because self.npart, etc, has type np.uint32 and not int.
        #This is I think a problem with python/numpy, but cast things to ints until I can determine how widespread it is. 
        data=struct.pack(fmt,int(self.npart[0]), int(self.npart[1]),int(self.npart[2]),int(self.npart[3]),int(self.npart[4]),int(self.npart[5]),
        self.mass[0], self.mass[1],self.mass[2],self.mass[3],self.mass[4],self.mass[5],
        self.time, self.redshift,  self.flag_sfr, self.flag_feedback,
        int(self.npartTotal[0]), int(self.npartTotal[1]),int(self.npartTotal[2]),int(self.npartTotal[3]),int(self.npartTotal[4]),int(self.npartTotal[5]),
        self.flag_cooling, self.num_files, self.BoxSize, self.Omega0, self.OmegaLambda, self.HubbleParam,self.flag_stellarage, self.flag_metals,
        int(self.NallHW[0]), int(self.NallHW[1]),int(self.NallHW[2]),int(self.NallHW[3]),int(self.NallHW[4]),int(self.NallHW[5]),
        self.flag_entropy_instead_u, self.flag_doubleprecision, self.flag_ic_info, self.lpt_scalingfactor)
        return data


class GadgetFile(object):
    """Class for all gadget snapshots"""
    def __init__(self, filename=None, sim=None, **kwargs):
        if sim is None: self.readfile(filename,**kwargs)
        else:
            if filename is None: filename=sim.filename
            if not sim._num_particles:
                self.readfile(filename,**kwargs)
            else:
                self.initfile(filename,sim,**kwargs)
        # if sim is None: self=GadgetReadFile(filename,**kwargs)
        # else          : 
        #     if filename is None: filename=sim.filename
        #     if not sim._num_particles:
        #         self=GadgetReadFile(filename,**kwargs)
        #     else:
        #         self=GadgetWriteFile(filename,sim,**kwargs)

    def get_block_types(self,block, npart):
        """ Set up the particle types in the block, with a heuristic,
        which assumes that blocks are either fully present or not for a given particle type"""
        #This function is horrible.
        p_types = np.zeros(N_TYPE,bool)
        if block.length == npart.sum()*block.partlen:
            p_types= np.ones(N_TYPE, bool)
            return p_types
        #Blocks which contain a single particle type
        for n in np.arange(0,N_TYPE) :
            if block.length == npart[n]*block.partlen :
                p_types[n] = True
                return p_types
        #Blocks which contain two particle types
        for n in np.arange(0,N_TYPE) :
            for m in np.arange(0,N_TYPE) :
                if block.length == (npart[n]+npart[m])*block.partlen :
                    p_types[n] = True
                    p_types[m] = True
                    return p_types
        #Blocks which contain three particle types
        for n in np.arange(0,N_TYPE) :
            for m in np.arange(0,N_TYPE) :
                for l in np.arange(0,N_TYPE) :
                    if block.length == (npart[n]+npart[m]+npart[l])*block.partlen :
                        p_types[n] = True
                        p_types[m] = True
                        p_types[l] = True
                        return p_types
        #Blocks which contain four particle types
        for n in np.arange(0,N_TYPE) :
            for m in np.arange(0,N_TYPE) :
                if block.length == (npart.sum() - npart[n]-npart[m])*block.partlen :
                    p_types = np.ones(N_TYPE, bool)
                    p_types[n] = False
                    p_types[m] = False
                    return p_types
        #Blocks which contain five particle type
        for n in np.arange(0,N_TYPE) :
            if block.length == (npart.sum() -npart[n])*block.partlen :
                p_types = np.ones(N_TYPE, bool)
                p_types[n] = False
                return p_types
        print block.length,npart.sum()*block.partlen
        raise ValueError, "Could not determine particle types for block"

    def check_format(self, fd):
        """This function reads the first character of a file and, depending on its value, determines
        whether we have a format 1 or 2 file, and whether the endianness is swapped. For the endianness,
        it then determines the correct byteorder string to pass to struct.unpack. There is not string
        for 'not native', so this is more complex than it needs to be"""
        fd.seek(0,0)
        (r,) = struct.unpack('=I',fd.read(4))
        if r == 8 :
            self.endian = '='
            self.format2 = True
        elif r == 134217728 :
            if sys.byteorder == 'little':
                self.endian = '>'
            else :
                self.endian = '<'
            self.format2 = True
        elif r == 65536 :
            if sys.byteorder == 'little':
                self.endian = '>'
            else :
                self.endian = '<'
            self.format2 = False
        elif r == 256 :
            self.endian = '='
            self.format2 = False
        else :
            raise IOError, "File corrupt. First integer is: "+str(r)
        fd.seek(0,0)
        return

    def get_block(self, name, p_type, p_toread) :
        """Get a particle range from this file, starting at p_start,
        and reading a maximum of p_toread particles"""
        p_read = 0
        cur_block = self.blocks[name]
        parts = self.get_block_parts(name, p_type)
        p_start = self.get_start_part(name, p_type)
        if p_toread > parts :
            p_toread = parts
        fd=open(self._filename, 'rb')
        fd.seek(cur_block.start+int(cur_block.partlen*p_start),0)
        #This is just so that we can get a size for the type
        dt = np.dtype(cur_block.data_type)
        n_type = p_toread*cur_block.partlen/dt.itemsize
        data=np.fromfile(fd, dtype=cur_block.data_type, count=n_type, sep = '')
        fd.close()
        if self.endian != '=' :
            data=data.byteswap(True)
        return (p_toread, data)

    def get_block_parts(self, name, p_type):
        """Get the number of particles present in a block in this file"""
        if not self.blocks.has_key(name) :
            return 0
        cur_block = self.blocks[name]
        if p_type == -1 :
            return cur_block.length/cur_block.partlen
        else :
            return self.header.npart[p_type]*cur_block.p_types[p_type]

    def get_start_part(self, name, p_type) :
        """Find particle to skip to before starting, if reading particular type"""
        if p_type == -1:
            return 0
        else :
            if not self.blocks.has_key(name) :
                return 0
            cur_block = self.blocks[name]
            return (cur_block.p_types*self.header.npart)[0:p_type].sum().astype(long)

    def get_block_dims(self, name):
        """Get the dimensionality of the block, eg, 3 for POS, 1 for most other things"""
        if not self.blocks.has_key(name) :
            return 0
        return self.blocks[name].dim
    
    def add_file_block(self, name, **inkw):
        """Add a block to the block table at the end of the file. Do not actually write anything"""
        if self.blocks.has_key(name) :
            raise KeyError,"Block "+name+" already present in file. Not adding"
        #Get last block
        if self.format2: headsize=5*4
        else           : headsize=4
        footsize=4
        if len(self.blocks.values())==0:
            inkw['start']=headsize
        else:
            lb=max(self.blocks.values(), key=lambda val: val.start)
            inkw['start']=lb.start+lb.length+headsize+footsize #For the block header, and footer of the previous block
        #Make new block
        block=GadgetBlock(self,name=name,**inkw)
        #Add name and block
        self.blocks[name]=block
        self.block_names.append(name)
 
    #The following functions are for writing blocks back to the file
    def write_block(self, name, p_type, big_data, filename=None) :
        """Write a full block of data in this file. Any particle type can be written. If the particle type is not present in this file, 
        an exception KeyError is thrown. If there are too many particles, ValueError is thrown. 
        big_data contains a reference to the data to be written. Type -1 is all types"""
        try:
            cur_block=self.blocks[name]
        except KeyError:
            raise KeyError, "Block "+name+" not in file "+self._filename
        parts = self.get_block_parts(name, p_type)
        p_start = self.get_start_part(name, p_type)
        MinType=np.ravel(np.where(cur_block.p_types * self.header.npart))[0]
        MaxType=np.ravel(np.where(cur_block.p_types * self.header.npart))[-1]
        #Have we been given the right number of particles?
        if np.size(big_data) > parts*self.get_block_dims(name):
            raise ValueError, "Space for "+str(parts)+" particles of type "+str(p_type)+" in file "+self._filename+", "+str(np.shape(big_data)[0])+" requested."
        #Do we have the right type?
        dt = np.dtype(cur_block.data_type)
        bt=big_data.dtype
        if bt.kind != dt.kind : 
            raise ValueError, "Data of incorrect type passed to write_block"
        #Open the file
        if filename == None : 
            fd = open(self._filename, "r+b")
        else :
            fd = open(filename, "r+b")
        #Seek to the start of the block
        fd.seek(cur_block.start+cur_block.partlen*p_start,0)
        #Add the block header if we are at the start of a block
        if p_type == MinType  or p_type < 0:
            data=cur_block.rw_head('W',fd,retdata=True)
            #Better seek back a bit first. (cur_block is the start of the block, after the header)
            fd.seek(-len(data),1)
            fd.write(data)

        # Swap bytes
        if self.endian != '=': big_data=big_data.byteswap(False)

        #Actually write the data
        #Make sure to ravel it, otherwise the wrong amount will be written, 
        #because it will also write nulls every time the first array dimension changes.
        d=np.ravel(big_data.astype(dt)).tostring()
        fd.write(d)
        if p_type == MaxType or p_type < 0: cur_block.rw_foot('w',fd)

        fd.close()

    def write_header(self, head_in, filename=None) :
        """Write a file header. Overwrites npart in the argument with the npart of the file, so a consistent file is always written."""
        #Construct new header with the passed header and overwrite npart with the file header. 
        #This has ref. semantics so use copy
        head=copy.deepcopy(head_in)
        head.npart=np.array(self.header.npart)
        if filename == None: filename = self._filename
        #a mode will ignore the file position, and w truncates the file.
        try :
            fd = open(filename, "r+")
        except IOError as (err, strerror):
            #If we couldn't open it because it doesn't exist open it for writing.
            if err == errno.ENOENT :
                fd = open(filename, "w+")
            #If we couldn't open it for any other reason, reraise exception
            else :
                raise IOError(err,strerror)
        fd.seek(0) #Header always at start of file
        #Write block header
        self.blocks['HEAD'].data=head
        self.blocks['HEAD'].rw_block('w',fd)
        #self.rw_block_head('W',fd,"HEAD",blocksize=256)
        #Write header body
        # fd.write(head.serialize())
        # #Write block footer
        # #Seek 48 bytes forward, to skip the padding (which may contain extra data)
        # fd.seek(48,1)
        # self.rw_block_foot('W',fd,blocksize=256)
        fd.close()

# class GadgetReadFile(GadgetFile) :
#     """Gadget file management class. Users should access gadget files through
#     :class:`~pynbody.gadget.GadgetSnap`."""

#     def __init__(self, filename, makeopt={}, **exkw) :
        

    def readfile(self, filename, makeopt={}, **exkw) : 
        global config
        self._filename=filename
        if not os.path.isfile(self._filename): self._filename+='.0'
        self.blocks = {}
        self.endian=''
        self.format2=True
        t_part = 0
        fd=open(filename, "rb")
        self.check_format(fd)
        #If format 1, load the block definitions.
        if not self.format2 :
            self.block_names = config_parser.get('gadget-1-blocks',"blocks").split(",")
	    self.block_names = [q.upper().ljust(4) for q in self.block_names]
            #Adjust based on makeopt
            mkopt2block=dict(OUTPUTPOTENTIAL      ='POT ',
                             OUTPUTACCELERATION   ='ACCE',
                             OUTPUTCHANGEOFENTROPY='ENDT',
                             OUTPUTTIMESTEP       ='TSTP')
            for opt in mkopt2block:
                if opt not in makeopt: makeopt[opt]=False
                if not makeopt[opt]: self.block_names.remove(mkopt2block[opt])
	    #This is a counter for the fallback
	    self.extra = 0
            if config['verbose']: print self.block_names
        while True:
            block=GadgetBlock(self)
            block.rw_block('r',fd)
            if block.length==0: break
            #Do special things for the HEAD block
            if block.name == "HEAD" :
                self.header=block.data
                #Handle blocks that header rules out
                if not self.format2:
                    #Remove MASS block from list if not included
                    if self.header.npart[self.header.mass==0].sum()==0:
                        self.block_names.remove('MASS')
                    #Remove SPH blcoks if not included
                    if self.header.npart[_type_map[family.get_family('gas')]].sum()==0:
                        for isphblock in _sphblock_names: 
                            if isphblock in self.block_names: self.block_names.remove(isphblock)
                    if config['verbose']: print self.block_names
            #Add block to list of blocks
            self.blocks[block.name] = block
        #and we're done.
        fd.close()

        # Make a mass block if one isn't found.
        # In the header, mass is a double
        if 'MASS' not in self.blocks: self.blocks['MASS'] = GadgetBlock(self,name='MASS',partlen=8,dtype=np.float64)

    def initfile(self, filename, sim, format2=False, flag_ic=False, idx_file=0, num_files=None, **exkw) : 
        """Initializes files based on simulation info"""
        self.endian='=' # write with default endian of this system
        self.format2=format2
        self.blocks={}
        self.block_names=[]
        # Initialize properties
        sim._decorate()
        # Get npartTotal & mass
        npartTotal = np.zeros(N_TYPE, int)
        massarr = np.zeros(N_TYPE, np.float64)
        softarr = np.zeros(N_TYPE, float)
        for t in _rev_type_map.keys():
            npartTotal[t] = gadget_typelen(sim,t)
            if npartTotal[t]==0: continue
            tidx=gadget_typeidx(sim,t)
            imass=np.double(sim[tidx]['mass'][0])
            isoft=sim[tidx]['eps' ][0]
            if np.all(sim[tidx]['mass']==imass): massarr[t]=imass
            if np.all(sim[tidx]['eps' ]==isoft): softarr[t]=isoft
            #Note that if we have families with identical type maps, we cannot
            #determine which type each individual particle is
        # How to split files
        if num_files is None:
            if len(sim)*12. > 2**31-1: 
                num_files=int(np.ceil(len(sim)*12./(2**31-1)))
                warnings.warn("Data too large to fit into a single gadget file. Splitting into {}".format(num_files))
            else                     : num_files=0
        if num_files==0:
            self._filename=filename
            npart=npartTotal
        else:
            warnings.warn("File splitting untested. Double check output.")
            self._filename=filename+'.0'
            npart=npartTotal/num_files
            # Extras go in last file
            if idx_file==(num_files-1): npart+=npartTotal%num_files
        # Header
        # npart, mass, time, redshift, BoxSize,Omega0, OmegaLambda, HubbleParam, num_files=1 
        self.header=GadgetHeader(npart,massarr,sim.properties["time"],sim.properties["z"],
                                 sim.properties["boxsize"].in_units(sim['pos'].units, **sim.conversion_context()),
                                 sim.properties["omegaM0"],sim.properties["omegaL0"],sim.properties["h"],num_files)
        self.header.npartTotal=npartTotal
        self.add_file_block('HEAD')
        # Get keys
        all_keys=set(sim.loadable_keys()).union(sim.keys()).union(sim.family_keys())
        all_keys = [ k for k in all_keys if not sim.is_derived_array(k) and not k in ["x","y","z","vx","vy","vz"] ] 
        # Add IDs if missing (ensuring correct type order)
        if 'iord' not in all_keys:
            sim['iord']=array.SimArray(np.zeros(len(sim),dtype=np.int32))
            iiord=0L ; fiord=0L
            for t in _rev_type_map.keys():
                tidx=gadget_typeidx(sim,t)
                fiord+=tidx.sum()
                sim[tidx]['iord']=np.arange(iiord,fiord,dtype=np.int32)
                iiord=fiord
            all_keys.append('iord')
        # Remove mass & softenings from keylist if can use mass arr
        if np.all(massarr[npart!=0]!=0) and 'mass' in all_keys: all_keys.remove('mass')
        if np.all(softarr[npart!=0]!=0) and 'eps'  in all_keys: all_keys.remove('eps' )
        # Remove potential if this is an IC file
        if flag_ic and 'phi' in all_keys: all_keys.remove('phi')
        # Select only gadget keys
        allkeys = _output_order_gadget(all_keys,incldregs=False)
        for k in allkeys:
            p_types = np.zeros(N_TYPE,bool)
            for t in _rev_type_map.keys() :
                # Skip types w/o particles or mass (if on mass key)
                if npart[t]==0: continue
                if k=='mass' and massarr[t]!=0: continue
                tidx=gadget_typeidx(sim,t)
                # Things can be derived for some families but not others
                try:
                    if sim[tidx].is_derived_array(k): continue
                    dtype = sim[tidx][k].dtype
                    p_types[t] += True
                    try:
                        partlen = np.shape(sim[tidx][k])[1]*dtype.itemsize
                    except IndexError:
                        partlen = dtype.itemsize
                except KeyError: pass
            # Create block
            name=_translate_array_name(k).upper().ljust(4)[0:4]
            if p_types.sum():
                self.add_file_block(name,npart=npart*p_types,partlen=partlen,dtype=dtype)
        if config['verbose']: print self.block_names


class oldGadgetWriteFile(GadgetFile) :
    """Class for write-only snapshots, as when we are creating a new set of files from, eg, a TipsySnap.
        Should not be used directly. block_names is a list so we can specify an on-disc ordering."""
    def __init__(self, filename, npart, block_names, header, format2=True) :
        self.header=header
        self._filename = filename
        self.endian='=' # write with default endian of this system
        self.format2=format2
        self.blocks={}
        self.header.npart = np.array(npart)
        #Set up the positions
        header_size = 4
        if format2 :
            header_size += 3*4 + 4
        footer_size = 4
        #First block is just past the header. 
        cur_pos = 256 + header_size + footer_size
        for block in block_names :
            #Add block if present for some types
            if block.types.sum() :
                b_part = npart * block.types
                b=GadgetBlock(self,start=cur_pos+header_size, partlen=block.partlen, length=block.partlen*b_part.sum(), dtype=block.dtype,p_types=block.types)
                cur_pos += b.length+header_size+footer_size
                self.blocks[block.name] = b

class WriteBlock :
    """Internal structure for passing data around between file and snapshot"""
    def __init__(self, partlen=4, dtype=np.float32, types = np.zeros(N_TYPE,bool), name = "    ") :
        #Bytes per particle in file
        self.partlen=partlen
        #Data type of block
        self.dtype = dtype
        #Types of particle this block contains
        self.types = types
        self.name = name

class GadgetSnap(snapshot.SimSnap):
    """Main class for reading Gadget-2 snapshots. The constructor makes a map of the locations
    of the blocks, which are then read by _load_array"""
    _class_inherited = ['_makefile_name','_unitfile_name']

    def __init__(self, filename=None, **kwargs) :
        global config
        super(GadgetSnap,self).__init__()
        if filename: self._setup(filename=filename,**kwargs)

    def _setup(self, **kwargs) :
        # File names
        self._setup_filenames(**kwargs)
        # Files
        self._setup_files(**kwargs)
        # Header
        self._setup_header(**kwargs)
        # Family slices
        self._setup_family_slice(**kwargs)
        # Loadable keys
        self._setup_loadable_keys(**kwargs)
        # Decorate
        self._decorate()
        # Galaxy slices
        self._setup_galaxy_slice(**kwargs)

    @property
    def filename(self):
        return self._filename
    @filename.setter
    def filename(self,val):
        if isinstance(val,str):
            if not self._num_particles:
                try:
                    fd=open(val)
                except IOError:
                    fd=open(val+".0")
                    val+='.0'
                fd.close()
            if val[-2:]=='.0': self._filename=val[:-2]
            else             : self._filename=val

    def _setup_filenames(self, **kwargs):
        # Initialize filenames
        self.filename=kwargs.get('filename',None)
        self.unitfile=kwargs.get('unitfile',None)
        self.makefile=kwargs.get('makefile',None)

    def _setup_files(self, **kwargs):
        #Initialize stuff
        if not hasattr(self,'_files'): self._files=[]
        self._files=kwargs.get('files',self._files)
        #Fill in list of files
        if len(self._files)==0:
            if self.__len__()*12. > 2**31-1:
                print len(self),'>',(2**31-1)/12.
                #raise IOError,"Data too large to fit into a single gadget file, and splitting not implemented. Cannot write."
            # Read/create files
            first_file = GadgetFile(sim=self,makeopt=self.makefile,**kwargs)
            self._files.append(first_file)
            files_expected = self._files[0].header.num_files
            for i in np.arange(1, files_expected):
                kwargs['filename'] = first_file._filename[:-1]+str(i)
                tmp_file=GadgetFile(sim=self,makeopt=self.makefile,**kwargs)
                if not self.check_headers(tmp_file.header, self._files[0].header) :
                    warnings.warn("file "+str(i)+" is not part of this snapshot set!",RuntimeWarning)
                    continue
                self._files.append(tmp_file)

    def _setup_header(self,**kwargs):
        if not hasattr(self,'header'): self.header=copy.deepcopy(kwargs.get('header',None))
        if not self.header:
            #Set up global header
            self.header=copy.deepcopy(self._files[0].header)
            #Get number of particles
            npart = np.zeros_like(self._files[0].header.npart)
            for f in self._files: npart=npart+f.header.npart
            self.header.npart=npart
            #Check and fix npartTotal and NallHW if they are wrong.
            if npart is not self.header.npartTotal+2**32*self.header.NallHW :
                self.header.NallHW = npart/2**32
                self.header.npartTotal = npart - 2**32*self.header.NallHW
                for f in self._files :
                    f.header.npartTotal = self.header.npartTotal
                    f.header.NallHW = self.header.NallHW
        self._num_particles = self.header.npart.sum()
        
    def _setup_family_slice(self, setslice=True, **kwargs):
        npart = kwargs.get('npart',self.header.npart)
        _family_slice = {}
        for x in _type_map :
            max_t=_type_map[x]
            _family_slice[x] = slice(npart[0:np.min(max_t)].sum(),npart[0:np.max(max_t)+1].sum())
        if setslice: self._family_slice=_family_slice
        else: return _family_slice

    def _setup_loadable_keys(self, **kwargs):
        self._loadable_keys = set([])
        for f in self._files :
            self._loadable_keys = self._loadable_keys.union(set(f.blocks.keys()))
        if 'HEAD' in self._loadable_keys: self._loadable_keys.remove('HEAD')
        #Add default mapping to unpadded lower case if not in config file.
        for nn in self._loadable_keys : 
            mm = nn.lower().strip()
            if not nn in _rev_name_map :
                _rev_name_map[nn] = mm
            if not mm in _name_map :
                _name_map[mm] = nn
        #Use translated keys only
        self._loadable_keys = [_translate_array_name(x, reverse=True) for x in self._loadable_keys]
        #Set up block list, with attached families, as a caching mechanism
        self._block_list = self.get_block_list()

    def _setup_galaxy_slice(self, setslice=True, **kwargs):
        if not hasattr(self,'_galids'): self._galids=[]
        self._galids=kwargs.get('galids',self._galids)
        _galaxy_slice = {}
        #Fix galids
        galids=sorted(self._galids)
        if len(galids)==0: galids=[long(self['iord'].min()),long(self['iord'].max())+1]
        if galids[ 0]  > self['iord'].min(): galids=[long(self['iord'].min())]+galids
        if galids[-1]  < self['iord'].max(): galids=galids+[long(self['iord'].max())]
        if galids[-1] == self['iord'].max(): galids[-1]+=1
        self._galids = galids
        #Find galaxies
        igal0=1
        for igal in range(1,len(galids)):
            iid = galids[igal-1]
            fid = galids[igal]
            idx = np.logical_and(self['iord']>=iid,self['iord']<fid)
            if np.any(idx): 
                _galaxy_slice[galaxy.get_galaxy('galaxy{}'.format(igal0))] = idx
                igal0+=1
        #Set slice
        if setslice: self._galaxy_slice=_galaxy_slice
        else: return _galaxy_slice

    def _append_family_slice(self, solf, **kwargs):
        if solf.__class__ is not GadgetSnap:
            raise Exception('Incompatible class {}'.format(solf.__class__))
        # Record things
        oldnpart1=self.header.npart
        oldnpart2=solf.header.npart
        oldn1=oldnpart1.sum()
        oldn2=oldnpart2.sum()
        oldfamslice1=self._family_slice
        oldfamslice2=solf._family_slice
        newnpart=oldnpart1+oldnpart2
        newn=oldn1+oldn2
        return self._setup_family_slice(npart=newnpart,setslice=False)

    def _append(self, solf, **kwargs) :
        if solf.__class__ is not GadgetSnap:
            raise Exception('Incompatible class {}'.format(solf.__class__))
        # Check header
        chkattr=['BoxSize','Omega0','OmegaLambda','HubbleParam']
        for ichk in chkattr:
            if getattr(self.header,ichk)!=getattr(solf.header,ichk):
                raise Exception('Incompatible header attribute: {}={},{}'.format(ichk,getattr(self.header,ichk),getattr(solf.header,ichk)))
#        chkmass=np.intersect1d(np.where(self.header.mass>0)[0],np.where(solf.header.mass>0)[0])
#        if not np.array_equal(self.header.mass[chkmass],solf.header.mass[chkmass]):
        if not np.array_equal(self.header.mass,solf.header.mass):
            raise Exception('Incompatible header attribute: mass={},{}'.format(self.header.mass,solf.header.mass))
        # Append header
        self.header.npart=self.header.npart+solf.header.npart

    def _downsample(self, ntot, **kwargs):
        """Downsample by types"""
        ntot0=len(self)
        # Allocate
        idx=np.array([],dtype=int)
        npart=np.zeros(len(self.header.npart),dtype=int)
        massarr=np.zeros(len(self.header.npart),dtype=np.float64)
        # Down sample by type
        for i in range(len(self.header.npart)):
            npart[i]=int(self.header.npart[i]*(float(ntot)/float(ntot0)))
            massarr[i]=self.header.mass[i]*(np.float64(ntot0)/np.float64(ntot))
            typidx=np.arange(ntot0,dtype=int)[gadget_typeidx(self,i)]
            np.random.shuffle(typidx)
            idx=np.append(idx,typidx[:npart[i]])
        # Create new snapshot
        pnb=copy.deepcopy(self[np.sort(idx)])
        pnb.header.npart=npart
        pnb.header.mass=massarr
        # Return
        return pnb

    def _init_file_units_system(self,unitfile=None):
        """Gets unit system from Gadget unitfile"""
        #cosmo = (sim._hdf['Parameters']['NumericalParameters'].attrs['ComovingIntegrationOn'])!=0
        # Load units from file
        if unitfile and os.path.isfile(unitfile):
            del self.unitfile
            self.unitfile=unitfile
        if self.unitfile:
            param=self.unitfile
            if param['ComovingIntegrationOn']:
                vel_unit = units.Unit(str(param['UnitVelocity_in_cm_per_s'])+' cm s**-1 a**1/2')
                dist_unit = units.Unit(str(param['UnitLength_in_cm'])+' cm h**-1')
                mass_unit = units.Unit(str(param['UnitMass_in_g'])+' g h**-1')
            else:
                vel_unit = units.Unit(str(param['UnitVelocity_in_cm_per_s'])+' cm s**-1')
                dist_unit = units.Unit(str(param['UnitLength_in_cm'])+' cm')
                mass_unit = units.Unit(str(param['UnitMass_in_g'])+' g')
        # Set to default units
        else:
            print 'Setting default units'
            vel_unit = config_parser.get('gadget-units', 'vel')
            dist_unit = config_parser.get('gadget-units', 'pos')
            mass_unit = config_parser.get('gadget-units', 'mass')
        # Set file units system
        return [vel_unit,dist_unit,mass_unit,"K"]
        

    @property
    def unitfile(self):
        extlist=['param','gdpar']
        if not hasattr(self,'_unitfile_dict'):
            if not hasattr(self,'_unitfile_name'):
                from mmlutils.mmlfiles import search as search_file
                self._unitfile_name=search_file(self.filename,extlist,nlvlup=1)
            if self._unitfile_name:
                from pysim.files.gadget import rw_param
                self._unitfile_dict=rw_param('R',self._unitfile_name)
            else:
                self._unitfile_dict={}
                raise Exception('Unitfile invalid: {}'.format(self._unitfile_name))
        return self._unitfile_dict
    @unitfile.setter
    def unitfile(self,val): 
        if val and os.path.isfile(val):
            del self.unitfile
            self._unitfile_name=val
    @unitfile.deleter
    def unitfile(self):
        del self._unitfile_name
        del self._unitfile_dict

    @property
    def makefile(self):
        extlist=['Makefile','makefile','gdmk']
        if not hasattr(self,'_makefile_dict'):
            if not hasattr(self,'_makefile_name'):
                from mmlutils.mmlfiles import search as search_file
                self._makefile_name=search_file(self.filename,extlist,nlvlup=1)
            if self._makefile_name:
                from pysim.files.gadget import rw_makefile
                self._makefile_dict=rw_makefile('R',self._makefile_name)
            else:
                self._makefile_dict={}
                print self.filename
                raise Exception('Makefile invalid: {}'.format(self._makefile_name))
        return self._makefile_dict
    @makefile.setter
    def makefile(self,val): 
        if val and os.path.isfile(val):
            del self.makefile
            self._makefile_name=val
        # else: 
        #     raise Exception('Invalid makefile: {}'.format(val))
    @makefile.deleter
    def makefile(self):
        del self._makefile_name
        del self._makefile_dict

    @staticmethod
    def _get_npart(self):
        if self.__class__ is GadgetSnap:
            npart = self.header.npart
        else:
            npart = np.zeros(N_TYPE, int)
            for t in _rev_type_map.keys() : 
                npart[t] = gadget_typelen(self,t)
        return npart

    def loadable_family_keys(self, fam=None) :
        """Return list of arrays which are loadable for specific families, 
        but not for all families."""
        warnings.warn("loadable_family_keys functionality has now been incorporated into loadable_keys", warnings.DeprecationWarning)
        return self.loadable_keys(fam)


    def loadable_keys(self, fam=None) :
        if hasattr(self,'_loadable_keys'): 
            if fam is not None : 
                return [x for x in self._loadable_keys if self._family_has_loadable_array(fam, x)]
            else :
                return [x for x in self._loadable_keys if self._family_has_loadable_array(None, x)]
        else: return []


    def _family_has_loadable_array(self, fam, name) :
        """Returns True if the array can be loaded for the specified family.
        If fam is None, returns True if the array can be loaded for all families."""
        if name in self._block_list:
            if fam is not None :
                return fam in self._block_list[name]
            else :
                return set(self.families()) <= set(self._block_list[name])
        else:
            return False

    def get_block_list(self):
        """Get list of unique blocks in snapshot, with the types they refer to"""
        b_list = {}
        for f in self._files :
            for (n,b) in f.blocks.iteritems() :
                if n=='HEAD': continue
                if b_list.has_key(n) :
                    b_list[n] += b.p_types
                else :
                    b_list[n] = np.array(b.p_types, dtype=bool)
        #Special case mass. Note b_list has reference semantics.
        if b_list.has_key("MASS") :
            b_list["MASS"] += np.array(self.header.mass,dtype=bool)
        #Translate this array into families and external names
        out_list={}
        for k,b in b_list.iteritems() :
            b_name = _translate_array_name(k,reverse=True)
            #Make this be only if there are actually particles of that type in the snap
            b_types = [ f for f in self.families() if b[np.intersect1d(gadget_type(f), np.ravel(np.where(self.header.npart != 0)))].all() ]
            out_list[b_name] = b_types
        return out_list

    def get_block_parts(self, name, family) :
        """Get the number of particles present in a block, of a given type"""
        total=0
        for f in self._files:
            total+= sum([ f.get_block_parts(name, gfam) for gfam in gadget_type(family)])
        #Special-case MASS
        if name == "MASS" :
            total+= sum([ self.header.npart[p]*np.array(self.header.mass[p],dtype=bool) for p in gadget_type(family)])
        return total

    def check_headers(self, head1, head2) :
        """Check two headers for consistency"""
        if ( head1.time != head2.time or head1.redshift!= head2.redshift or
           head1.flag_sfr != head2.flag_sfr or
           head1.flag_feedback != head2.flag_feedback or
           head1.num_files != head2.num_files or
           head1.BoxSize != head2.BoxSize or
           head1.Omega0 != head2.Omega0 or
           head1.OmegaLambda != head2.OmegaLambda or
           head1.HubbleParam != head2.HubbleParam  or
           head1.flag_stellarage != head2.flag_stellarage or
           head1.flag_metals != head2.flag_metals) :
            return False
        #Check array quantities
        if (((head1.mass - head2.mass) > 1e-5*head1.mass).any()  or 
                (head1.npartTotal != head2.npartTotal).any()) :
            return False
        #  At least one version of N-GenICs writes a header file which
        #  ignores everything past flag_metals (!), leaving it uninitialised.
        #  Therefore, we can't check them.
        return True
    def _get_array_type(self, name) :
        """Get the type for the array given in name"""
        g_name = _translate_array_name(name)
        return self._get_array_type_g(g_name)
    
    def _get_array_type_g(self, name) :
        """Get the type for the array given in name"""
        return self._files[0].blocks[name].data_type


    def _get_array_dims(self, name) :
        """Get the dimensions of an array; ie, is it 3d or 1d"""
        g_name = _translate_array_name(name)
        return self._files[0].get_block_dims(g_name)

    def _load_array(self, name, fam=None) :
        """Read in data from a Gadget file.
        If fam != None, loads only data for that particle family"""
        #g_name is the internal name
        g_name = _translate_array_name(name)

        if not self._family_has_loadable_array( fam, name) :
            if fam is None and name in self._block_list:
                raise KeyError,"Block "+name+" is not available for all families"
            else :
                raise IOError, "No such array on disk"

        ndim = self._get_array_dims(name)

        if ndim == 1:
            dims = [self.get_block_parts(g_name, fam),]
        else:
            dims = [self.get_block_parts(g_name, fam), ndim]

        p_types = gadget_type(fam)

        #Get the data. Get one type at a time and then concatenate. 
        #A possible optimisation is to special-case loading all particles.
        data = np.array([], dtype = self._get_array_type(name))
        for p in p_types :
            #Special-case mass
            if g_name == "MASS" and self.header.mass[p] != 0. :
                data = np.append(data, self.header.mass[p]*np.ones(self.header.npart[p],dtype=data.dtype))
            else :
                data = np.append(data, self.__load_array(g_name, p))

        if fam is None :
            self[name] = data.reshape(dims,order='C').view(array.SimArray)
            self[name].set_default_units(quiet=True)
        else :
            self[fam][name] = data.reshape(dims,order='C').view(array.SimArray)
            self[fam][name].set_default_units(quiet=True)


    def __load_array(self, g_name, p_type) :
        """Internal helper function for _load_array that takes a g_name and a gadget type, 
        gets the data from each file and returns it as one long array."""
        data=np.array([],dtype=self._get_array_type_g(g_name))
        #Get a type from each file
        for f in self._files:
            f_parts = f.get_block_parts(g_name, p_type)
            if f_parts == 0:
                continue
            (f_read, f_data) = f.get_block(g_name, p_type, f_parts)
            if f_read != f_parts :
                raise IOError,"Read of "+f._filename+" asked for "+str(f_parts)+" particles but got "+str(f_read)
            data = np.append(data, f_data)
        return data

    @staticmethod
    def _can_load(f) :
        """Check whether we can load the file as Gadget format by reading
        the first 4 bytes"""
        try:
            fd=open(f)
        except IOError:
            try:
                fd=open(f+".0")
            except:
                return False
                #If we can't open the file, we certainly can't load it...
        (r,) = struct.unpack('=I',fd.read(4))
        fd.close()
        #print '[nbody.gadget.py] r='+str(r)
        #First int32 is 8 for a Gadget 2 file, or 256 for Gadget 1, or the byte swapped equivalent.
        if r == 8 or r == 134217728 or r == 65536 or r == 256 :
            return True
        else :
            return False
    
    @staticmethod
    def _write(self, filename=None, unitfile=None, format2=False, overwrite=False, flag_ic=False, Nfiles=None) :
        """Write an entire Gadget file (actually an entire set of snapshots)."""
        # Pars filename
        if not filename:
            if self.__class__ is not GadgetSnap :
                raise Exception,"Please specify a filename to write a new file."
            else:
                filename=self.filename
            if os.path.isfile(filename) and not overwrite:
                print 'File already exists and overwrite not set'
                print '    '+filename
        # Get object to write
        if self.__class__ is GadgetSnap: gdsim=self
        else                           : gdsim=self.convert(GadgetSnap,unitfile=unitfile)
        #Convert units
        self.original_units()
        #Write headers
        if filename != None :
            if np.size(gdsim._files) > 1 :
                for i in np.arange(0, np.size(gdsim._files)) :
                    ffile = filename+"."+str(i)
                    if os.path.isfile(ffile) and not overwrite:
                        print 'File already exists and overwrite not set'
                        print '    '+ffile
                        return
                    gdsim._files[i].write_header(gdsim.header,ffile) 
            else :
                if os.path.isfile(filename) and not overwrite:
                    print 'File already exists and overwrite not set'
                    print '    '+filename
                else:
                    if os.path.isfile(filename): os.remove(filename)
                gdsim._files[0].write_header(gdsim.header,filename) 
        else :
            #Call write_header for every file. 
            [ f.write_header(gdsim.header) for f in gdsim._files ]
        #Call _write_array for every array.
        for x in gdsim._files[0].block_names:
            if x=='HEAD': continue
            print '[_write]',x
            GadgetSnap._write_array(self, _translate_array_name(x,reverse=True), filename=filename)

    @staticmethod
    def _oldwrite(self, filename=None, unitfile=None, format2=False, overwrite=False, flag_ic=False) :
        """Write an entire Gadget file (actually an entire set of snapshots)."""
        # Pars filename
        if not filename:
            if self.__class__ is not GadgetSnap :
                raise Exception,"Please specify a filename to write a new file."
            else:
                filename=self.filename
            if os.path.isfile(filename) and not overwrite:
                print 'File already exists and overwrite not set'
                print '    '+filename

        # Parameter file
        if unitfile:
            from pysim.files.gadget import rw_param
            param = rw_param('R',unitfile)
            # Units
            if param['ComovingIntegrationOn']:
                vel_unit = units.Unit(str(param['UnitVelocity_in_cm_per_s'])+' cm s**-1 a**1/2')
                dist_unit = units.Unit(str(param['UnitLength_in_cm'])+' cm h**-1')
                mass_unit = units.Unit(str(param['UnitMass_in_g'])+' g h**-1')
            else:
                vel_unit = units.Unit(str(param['UnitVelocity_in_cm_per_s'])+' cm s**-1')
                dist_unit = units.Unit(str(param['UnitLength_in_cm'])+' cm')
                mass_unit = units.Unit(str(param['UnitMass_in_g'])+' g')
            self.physical_units(distance=dist_unit, velocity=vel_unit, mass=mass_unit)
            # Header parameters
            self.properties.setdefault('boxsize',array.SimArray(param['BoxSize'],dist_unit))
            self.properties.setdefault('omegaM0',array.SimArray(param['Omega0']))
            self.properties.setdefault('omegaL0',array.SimArray(param['OmegaLambda']))
            self.properties.setdefault('h'      ,array.SimArray(param['HubbleParam']))
            self.properties.setdefault('a'      ,self._time)
            self.properties.setdefault('z'      ,1./(1.+self.properties['a']))
        
        
        with self.lazy_derive_off :
            #If caller is not a GadgetSnap, construct the GadgetFiles, 
            all_keys=set(self.loadable_keys()).union(self.keys()).union(self.family_keys())
            all_keys = [ k for k in all_keys if not self.is_derived_array(k) and not k in ["x","y","z","vx","vy","vz"] ] 
            #This code supports (limited) format conversions
            if self.__class__ is not GadgetSnap :
                #Splitting the files correctly is hard; the particles need to be reordered, and
                #we need to know which families correspond to which gadget types.
                #So don't do it.

                #Make sure the data fits into one files. The magic numbers are:
                #12 - the largest block is likely to  be POS with 12 bytes per particle. 
                #2**31 is the largest size a gadget block can safely have
                if self.__len__()*12. > 2**31-1:
                    raise IOError,"Data too large to fit into a single gadget file, and splitting not implemented. Cannot write."              
                #Make npart
                npart = np.zeros(N_TYPE, int)
                arr_name=(self.keys()+self.loadable_keys())[0]
                for t in _rev_type_map.keys() : 
                    npart[t] = gadget_typelen(self,t)
                    #Note that if we have families with identical type maps, we cannot
                    #determine which type each individual particle is
                #Make mass array
                massarr = np.zeros(N_TYPE, float)
                for t in _rev_type_map.keys():
                    if npart[t]==0: continue
                    tidx=gadget_typeidx(self,t)
                    im=self[tidx]['mass'][0]
                    if np.all(self[tidx]['mass']==im): massarr[t]=im
                #Construct a header
                # npart, mass, time, redshift, BoxSize,Omega0, OmegaLambda, HubbleParam, num_files=1 
                gheader=GadgetHeader(npart,massarr,self.properties["a"],self.properties["z"],
                                     self.properties["boxsize"].in_units(self['pos'].units, **self.conversion_context()),
                                     self.properties["omegaM0"],self.properties["omegaL0"],self.properties["h"],1)
                # Add IDs if missing (ensuring correct type order)
                if 'iord' not in all_keys:
                    self['iord']=array.SimArray(np.zeros(len(self),dtype=np.int32))
                    iiord=0L ; fiord=0L
                    for t in _rev_type_map.keys():
                        tidx=gadget_typeidx(self,t)
                        fiord+=tidx.sum()
                        self[tidx]['iord']=np.arange(iiord,fiord,dtype=np.int32)
                        iiord=fiord
                    all_keys.append('iord')
                # Select keys
                all_keys = _output_order_gadget(all_keys,incldregs=False)
                if np.all(massarr[npart!=0]!=0): 
                    if 'mass' in all_keys: all_keys.remove('mass')
                if flag_ic and 'phi' in all_keys: all_keys.remove('phi')
                #Construct the block_names; each block_name needs partlen, data_type, and p_types, 
                #as well as a name. Blocks will hit the disc in the order they are in all_keys.
                #First, make pos the first block and vel the second.
                block_names = []
                for k in all_keys :
                    types = np.zeros(N_TYPE,bool)
                    for t in _rev_type_map.keys() : 
                        if npart[t]==0: continue
                        if k=='mass' and massarr[t]!=0: continue
                        tidx=gadget_typeidx(self,t)
                        try :
                            #Things can be derived for some families but not others
                            if self[tidx].is_derived_array(k): continue
                            dtype = self[tidx][k].dtype
                            types[t] += True
                            try :
                                partlen = np.shape(self[tidx][k])[1]*dtype.itemsize
                            except IndexError:
                                partlen = dtype.itemsize
                        except KeyError:
                            pass
                    bb=WriteBlock(partlen, dtype=dtype, types = types, name = _translate_array_name(k).upper().ljust(4)[0:4])
                    block_names.append(bb)
                    print _translate_array_name(k).upper().ljust(4)[0:4],partlen,dtype,types
                #Create an output file
                out_file = GadgetWriteFile(filename, npart, block_names, gheader, format2=format2)
                #Write the header
                out_file.write_header(gheader,filename) 
                #Write all the arrays    
                for x in all_keys :
                    g_name = _translate_array_name(x).upper().ljust(4)[0:4]
                    #for fam in self.families() :
                    for t in _rev_type_map.keys() : 
                        if npart[t]==0: continue
                        tidx=gadget_typeidx(self,t)
                        if tidx.sum()==0: continue
                        try:
                            #Things can be derived for some families but not others
                            if self[tidx].is_derived_array(x): 
                                print g_name,t,'derived'
                                continue
                            data = self[tidx][x]
                            out_file.write_block(g_name, t, data, filename=filename)
                            print g_name,t,np.sum(data),data.units,filename
                        except KeyError :
                            pass
                return

            #Write headers
            if filename != None :
                if np.size(self._files) > 1 :
                    for i in np.arange(0, np.size(self._files)) :
                        ffile = filename+"."+str(i)
                        self._files[i].write_header(self.header,ffile) 
                else :
                    if os.path.isfile(self._filename) and not overwrite:
                        print 'File already exists and overwrite not set'
                        print '    '+self._filename
                    self._files[0].write_header(self.header,filename) 
            else :
                #Call write_header for every file. 
                [ f.write_header(self.header) for f in self._files ]
            #Call _write_array for every array.
            for x in all_keys :
                GadgetSnap._write_array(self, x, filename=filename)

    @staticmethod
    def _write_array(self, array_name, fam=None, filename=None) :
        """Write a data array back to a Gadget snapshot, splitting it across files."""
        write_fam = fam or self.families()
        
        #Make the name a four-character upper case name, possibly with trailing spaces
        g_name = _translate_array_name(array_name).upper().ljust(4)[0:4]
        nfiles=np.size(self._files)
        #Find where each particle goes
        f_parts = [ f.blocks[g_name].tpart for f in self._files ]
        #If there is no block corresponding to this name in the file, 
        # add it (so we can write derived arrays).
        if np.sum(f_parts) == 0:
            #Get p_type
            p_types=np.zeros(N_TYPE,bool)
            npart=np.zeros(N_TYPE,np.int32)
            for gfam in _rev_type_map.keys(): 
                gidx=gadget_typeidx(self,gfam)
                #We get the particle types we want by trying to load all 
                #particle types (from memory) and seeing which ones work
                p_types[gfam]=self[gidx].has_key(array_name)
                if p_types[gfam] :
                    ashape = np.shape(self[gidx][array_name])
                    #If the partlen is 1, append so the shape array has the right shape.
                    if np.size(ashape) < 2: ashape = (ashape[0],1)
                    npart[gfam]=ashape[0]
            if p_types.sum() :
                dtype=self[array_name].dtype
                per_file = npart/nfiles
                blk_kw={'name':g_name,'npart':per_file,'partlen':ashape[1]*dtype.itemsize,'dtype':dtype}
                for f in self._files[:-1]: f.add_file_block(**blk_kw)
                blk_kw['npart']=npart-(nfiles-1)*per_file
                self._files[-1].add_file_block(**blk_kw)

        #Find where each particle goes
        npart_prev=np.zeros(N_TYPE,dtype=np.int32)
        for i in np.arange(0,nfiles) :
            #Set up filename
            if filename != None:
                ffile = filename + '.'+str(i)
                if nfiles==1: ffile = filename
            else: ffile = None

            #Add data from each type
            data=np.array([],dtype=self[array_name].dtype)
            for t in _rev_type_map.keys():
                s_part=npart_prev[t]
                f_part=self._files[i].header.npart[t]
                tidx=gadget_typeidx(self,t)
                print '[_write_array]',array_name,t,self[array_name][tidx].shape,len(tidx),f_part
                idata=self[array_name][tidx][s_part:(s_part+f_part)]
                # s_part=self._files[i].blocks[g_name].spart
                # f_part=self._files[i].blocks[g_name].npart[t]
                # idata=self[tidx][array_name] #[s_part:(s_part+f_part-1)]
                #Special-case MASS. 
                if g_name == "MASS" and self.header.mass[t] != 0.:
                    #Warn if there are now different masses for this particle type, 
                    # as this information cannot be represented in this snapshot.
                    minmass = np.min(idata)
                    maxmass = np.max(idata)
                    if minmass!=maxmass: 
                        warnings.warn("Cannot write variable masses for type "+str(t)+", as masses are stored in the header.",RuntimeWarning)
                    elif self.header.mass[t] != minmass :
                        self.header.mass[t] = minmass
                        self._files[i].write_header(self.header, filename=ffile)
                else:
                    data=np.append(data,idata)
            #Add data to block & write
            self._files[i].blocks[g_name].data=data
            self._files[i].blocks[g_name].rw_block('w',seek=True,filename=ffile)
            #Advance npart counter
            npart_prev+=self._files[i].header.npart

@GadgetSnap.decorator
def do_units(sim) :
    #cosmo = (sim._hdf['Parameters']['NumericalParameters'].attrs['ComovingIntegrationOn'])!=0
    sim.init_file_units_system()

@GadgetSnap.decorator
def do_properties(sim) :
    list_prop=['time','a','omegaM0','omegaL0','boxsize','z','h']
    prop0={}
    # Units
    vel_unit,dist_unit,mass_unit,temp_unit=sim._file_units_system[:]
    time_unit=dist_unit/vel_unit
    # If header initialize
    if hasattr(sim,'header'):
        h = sim.header
        prop0['time'] = array.SimArray(h.time,time_unit)
        prop0['a'] = h.time
        prop0['omegaM0'] = h.Omega0
        #prop0['omegaB0'] = ... This one is non-trivial to calculate
        prop0['omegaL0'] = h.OmegaLambda
        prop0['boxsize'] = array.SimArray(h.BoxSize,dist_unit)
        prop0['z'] = h.redshift
        prop0['h'] = h.HubbleParam
    # If parameter file exists
    elif sim.unitfile:
        param=sim.unitfile
        # print param
        # print param['ComovingIntegrationOn']
        if param['ComovingIntegrationOn']:
            vel_unit = units.Unit(str(param['UnitVelocity_in_cm_per_s'])+' cm s**-1 a**1/2')
            dist_unit = units.Unit(str(param['UnitLength_in_cm'])+' cm h**-1')
            mass_unit = units.Unit(str(param['UnitMass_in_g'])+' g h**-1')
        else:
            vel_unit = units.Unit(str(param['UnitVelocity_in_cm_per_s'])+' cm s**-1')
            dist_unit = units.Unit(str(param['UnitLength_in_cm'])+' cm')
            mass_unit = units.Unit(str(param['UnitMass_in_g'])+' g')
        prop0['omegaM0'] = param['Omega0']
        prop0['omegaL0'] = param['OmegaLambda']
        prop0['boxsize'] = array.SimArray(param['BoxSize'],dist_unit)
        prop0['h'] = param['HubbleParam']
        if param['ComovingIntegrationOn']:
            prop0['time'] = array.SimArray(param['TimeBegin'])
            prop0['a'] = sim.properties.get('time',prop0['time'])
            prop0['z']=(1./sim.properties.get('a',prop0['a']))-1.
        else:
            prop0['time'] = array.SimArray(param['TimeBegin'],dist_unit/vel_unit)
            prop0['a'] = 1.0
            prop0['z']=0.
    # Otherwise
    else: prop0={p:None for p in list_prop}
    # Add parameters, but don't overwrite existing ones
    for p in list_prop: sim.properties.setdefault(p,prop0[p])

@GadgetSnap.decorator
def do_softenings(sim) :
    #Add array of softenings
    dist_unit=sim._file_units_system[1]
    sim['eps']=array.SimArray(np.zeros(len(sim),dtype=np.float32),dist_unit)
    for i,ityp in enumerate(config_parser.get('gadget-type-names','types').split(',')):
        epspar='Softening'+ityp.strip().capitalize()
        if epspar in sim.unitfile: 
            sim['eps'][gadget_typeidx(sim,i)]=sim.unitfile[epspar]



    # def __init__(self, filename, unitfile=None, makefile=None, galids=[]) :
    #     from mmlutils.mmlfiles import search as search_file
    #     global config
    #     super(GadgetSnap,self).__init__()
    #     self._files = []
    #     self._filename=filename
    #     self._unitfile=unitfile ; self._unitfile_opt={}
    #     self._makefile=makefile ; self._makefile_opt={}
    #     self._galids=galids
    #     npart = np.empty(N_TYPE)
    #     #Check whether the file exists, and get the ".0" right
    #     try:
    #         fd=open(filename)
    #     except IOError:
    #         fd=open(filename+".0")
    #         #The second time if there is an exception we let it go through
    #         filename = filename+".0"
    #     fd.close()
    #     if filename[-2:] == ".0" :
    #         self._filename = filename [:-2]
    #     #Initialize unitfile and makefile
    #     if not self._unitfile: self._unitfile=search_file(filename,['param','gdpar'],nlvlup=1)
    #     if not self._makefile: self._makefile=search_file(filename,['Makefile','makefile','gdmk'],nlvlup=1)
    #     #Read in unitfile & makefile options
    #     if self._unitfile:
    #         from pysim.files.gadget import rw_param
    #         self._unitfile_opt=rw_param('R',self._unitfile)
    #     if self._makefile:
    #         from pysim.files.gadget import rw_makefile
    #         self._makefile_opt=rw_makefile('R',self._makefile)
    #     #Read the first file and use it to get an idea of how many files we are expecting.
    #     first_file = GadgetFile(filename,makeopt=self._makefile_opt)
    #     self._files.append(first_file)
    #     files_expected = self._files[0].header.num_files
    #     npart = np.array(self._files[0].header.npart)
    #     for i in np.arange(1, files_expected):
    #         filename = filename[:-1]+str(i)
    #         tmp_file=GadgetFile(filename,makeopt=self._makefile_opt)
    #         if not self.check_headers(tmp_file.header, self._files[0].header) :
    #             warnings.warn("file "+str(i)+" is not part of this snapshot set!",RuntimeWarning)
    #             continue
    #         self._files.append(tmp_file)
    #         npart=npart+tmp_file.header.npart
    #     #Set up things from the parent class
    #     self._num_particles = npart.sum()
    #     #Set up global header
    #     self.header=copy.deepcopy(self._files[0].header)
    #     self.header.npart = npart
    #     #Check and fix npartTotal and NallHW if they are wrong.
    #     if npart is not self.header.npartTotal+2**32*self.header.NallHW :
    #         self.header.NallHW = npart/2**32
    #         self.header.npartTotal = npart - 2**32*self.header.NallHW
    #         for f in self._files :
    #             f.header.npartTotal = self.header.npartTotal
    #             f.header.NallHW = self.header.NallHW


    #     self._family_slice = {}
    #     self._galaxy_slice = {}

    #     self._loadable_keys = set([])
    #     self._family_keys=set([])
    #     self._family_arrays = {}
    #     self._arrays = {}
    #     self.properties = {}
        
    #     #Set up _family_slice
    #     for x in _type_map :
    #         max_t=_type_map[x]
    #         self._family_slice[x] = slice(npart[0:np.min(max_t)].sum(),npart[0:np.max(max_t)+1].sum())

    #     #Set up _loadable_keys
    #     for f in self._files :
    #         self._loadable_keys = self._loadable_keys.union(set(f.blocks.keys()))

    #     #Add default mapping to unpadded lower case if not in config file.
    #     for nn in self._loadable_keys : 
    #         mm = nn.lower().strip()
    #         if not nn in _rev_name_map :
    #             _rev_name_map[nn] = mm
    #         if not mm in _name_map :
    #             _name_map[mm] = nn

    #     #Use translated keys only
    #     self._loadable_keys = [_translate_array_name(x, reverse=True) for x in self._loadable_keys]
    #     #Set up block list, with attached families, as a caching mechanism
    #     self._block_list = self.get_block_list()
        
    #     self._decorate()
    
    #     #Galaxy slice
    #     #print 'iord',self['iord'].min(),self['iord'].max()
    #     orddiff = np.diff(np.sort(self['iord']))
    #     galids.sort()
    #     if galids[ 0]  > self['iord'].min(): galids=[long(self['iord'].min())]+galids
    #     if galids[-1]  < self['iord'].max(): galids=galids+[long(self['iord'].max())]
    #     if galids[-1] == self['iord'].max(): galids[-1]+=1
    #     for igal in range(1,len(galids)):
    #         iid = 0 if igal==0 else galids[igal-1]
    #         fid = galids[igal]
    #         idx = np.logical_and(self['iord']>=iid,self['iord']<fid)
    #         if np.any(idx): 
    #             #print idx.sum(),iid,fid,'galaxy{}'.format(igal)
    #             self._galaxy_slice[galaxy.get_galaxy('galaxy{}'.format(igal))] = idx
        
    #     #Add array of softenings
    #     self['eps']=array.SimArray(np.zeros(npart.sum(),dtype=np.float32),self['pos'].units)
    #     for i,ityp in enumerate(config_parser.get('gadget-type-names','types').split(',')):
    #         epspar='Softening'+ityp.strip().capitalize()
    #         if epspar in self._unitfile_opt: 
    #             self['eps'][npart[0:i].sum():npart[0:i+1].sum()]=self._unitfile_opt[epspar]

