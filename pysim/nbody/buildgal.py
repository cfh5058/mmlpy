"""

buildgal
========

Implements classes and functions for handling buildgal files.

"""

from __future__ import with_statement # for py2.5                                                                                                                                    
from __future__ import division

from . import snapshot, array, util
from . import family
from . import units
from . import config, config_parser
from . import chunk

import ConfigParser

import struct, os
import numpy as np
import gzip
import sys
import warnings
import copy
import types
import math

N_TYPE=12

_type_map = {}
for x in family.family_names() :
    try :
        pp =  [int(q) for q in config_parser.get('buildgal-type-mapping',x).split(",")]
        qq=np.array(pp)
        _type_map[family.get_family(x)] = pp
    except ConfigParser.NoOptionError :
        pass
_rev_type_map = {}
for f,t in _type_map.iteritems() : 
    for it in t:
        if it in _rev_type_map: _rev_type_map[it].append(f)
        else                  : _rev_type_map[it]=[f]

_name_map, _rev_name_map = util.setup_name_maps('buildgal-name-mapping')
_translate_array_name = util.name_map_function(_name_map, _rev_name_map)

def buildgal_type(fam):
    if fam == None:
        return list(np.arange(0,N_TYPE))
    else :
        return _type_map[fam]

def buildgal_typelen(sim,typ):
    Ntyp=buildgal_typeidx(sim,typ).sum()
    return Ntyp

def buildgal_typeidx(sim,typ):
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
    # Remove the number of particles of intervening types
    typlist.remove(typ)
    for t in typlist:
        if len(famlist.symmetric_difference(_rev_type_map[t]))==0:
            raise Exception('types {} and {} have identical families!'.format(typ,t))
        else: tidx[buildgal_typeidx(sim,t)]=False
    # Return index of particles in the type
    return tidx

class BuildgalSnap(snapshot.SimSnap):
    _class_inherited = ['_unitfile_name']
    _basic_loadable_keys = set(['mass','pos','vel','phi','tag'])

    def __init__(self, filename, **kwargs):
        super(BuildgalSnap,self).__init__()
        if filename:
            self._setup(filename=filename,**kwargs)
    def _setup(self, **kwargs):
        global config
        from mmlutils.mmlfiles import search as search_file
        # Setup file names
        self.filename=kwargs.get('filename',None)
        self.unitfile=kwargs.get('unitfile',None)
        self._time=0.0
        # Initialize unitfile & set units
        if not self._num_particles:
            if self.filename and os.path.isfile(self.filename):
                self._decorate()
                #Read in file
                if config['verbose'] : print>>sys.stderr, "BuildgalSnap: loading ",self.filename
                from pysim.files.buildgal import rw_snap
                fdict=rw_snap('R',self.filename)
                self._num_particles=fdict['header']['nbody']
                # Pars to correct variables
                for nm in fdict.keys():
                    if nm=='header': self.header=fdict['header']
                    else:
                        if nm=='pot': iarr=np.array(self._G)*fdict[nm]
                        else        : iarr=fdict[nm]
                        tnm=_translate_array_name(nm,reverse=True)
                        self[tnm]=iarr.view(array.SimArray)
                        self[tnm].set_default_units(quiet=True)
        #Set up things from the parent class
        self._family_slice = {}

        self._loadable_keys = set([])
        self._family_keys=set([])
        self._family_arrays = {}
        self.properties = {}

        # Set up family slice
        self._setup_family_slice(**kwargs)
        
    def _setup_family_slice(self, setslice=True, **kwargs):
        tagarr=kwargs.get('tag',self['tag'])
        _family_slice={}
        for x in _type_map :
            tags_t=_type_map[x] ; isl=np.zeros(len(tagarr),dtype=bool)
            for itag_t in tags_t: isl=np.logical_or(isl,tagarr==itag_t)
            _family_slice[x] = isl
        if setslice: self._family_slice=_family_slice
        else: return _family_slice

    def _append_family_slice(self, solf, **kwargs):
        if solf.__class__ is not BuildgalSnap:
            raise Exception('Incompatible class {}'.format(solf.__class__))
        dtype=getattr(self['tag'],'dtype',None)
        newtag=array.SimArray(np.zeros(len(self),dtype=dtype),self['tag'].units)
        newtag=np.hstack((self['tag'],solf['tag']))
        return self._setup_family_slice(setslice=False,tag=newtag)

    def _init_file_units_system(self,unitfile=None):
        """Get unit system from unit file"""
        if unitfile and os.path.isfile(unitfile):
            del self.unitfile
            self.unitfile=unitfile
        if self.unitfile:
            param=self.unitfile
            vel_unit = units.Unit(str(param['UnitVelocity_in_cm_per_s'])+' cm s**-1')
            dist_unit = units.Unit(str(param['UnitLength_in_cm'])+' cm')
            mass_unit = units.Unit(str(param['UnitMass_in_g'])+' g')
            self._G=array.SimArray(6.672e-8,'cm cm**2 s**-2 g**-1').in_unitsys([vel_unit,dist_unit,mass_unit])
            print 'G={} {}'.format(np.array(self._G),self._G.units)
        else:
            print '[_init_file_units_system] Incorrect units'
            vel_unit = config_parser.get('buildgal-units', 'vel')
            dist_unit = config_parser.get('buildgal-units', 'pos')
            mass_unit = config_parser.get('buildgal-units', 'mass')
            self._G=array.SimArray(1.0,units='0')
        return [units.Unit(x) for x in [vel_unit,dist_unit,mass_unit,"K"]]

    @property
    def filename(self):
        return self._filename
    @filename.setter
    def filename(self,val):
        if isinstance(val,str):
            if not self._num_particles and not os.path.isfile(val):
                raise Exception('Invalid file name: {}'.format(val))
            self._filename=val

    @property
    def unitfile(self):
        extlist=['units','bgunits']
        if not hasattr(self,'_unitfile_dict'):
            if not hasattr(self,'_unitfile_name'):
                from mmlutils.mmlfiles import search as search_file
                self._unitfile_name=search_file(self.filename,extlist,nlvlup=2)
            if self._unitfile_name:
                from pysim.files.buildgal import rw_units
                self._unitfile_dict=rw_units('R',self._unitfile_name)
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


    @staticmethod
    def _can_load(f) :
        try:
            fid=open(f,'r')
            headline=fid.readline()
            fid.close()
            if len(headline.split())==9: return True
            # check = BuildgalSnap(f, verbose=False)
            # del check
        except :
            raise
            return False
        return True

    @staticmethod
    def _write(self, filename=None, overwrite=False) :
        """
        Write a buildgal snapshot file
        """
        # Format conversion
        if self.__class__ is not BuildgalSnap:
            if filename == None :
                raise Exception,"Please specify a filename to write a new file."

            # Construct tag
            if 'tag' not in self:
                tagarr=np.zeros(self._num_particles,dtype=int)
                for t in _rev_type_map.keys():
                    tidx=buildgal_typeidx(self,t)
                    tagarr[tidx]=t
                self['tag']=tagarr.view(array.SimArray)
        else:
            if not filename: filename=self._filename

        # Prevent overwrite
        if not overwrite and os.path.isfile(filename):
            print 'File already exists and overwrite not set'
            print '    '+filename

        # Create dictionary of arrays for write out
        fdict={}
        for nm in BuildgalSnap._basic_loadable_keys:
            tnm=_translate_array_name(nm)
            fdict[tnm]=np.array(self[nm])

        # Write files
        from pysim.files.buildgal import rw_snap
        rw_snap('W',filename,fdict,overwrite=overwrite)
                

@BuildgalSnap.decorator
def do_units(sim) :
    sim.init_file_units_system()
