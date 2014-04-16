"""

buildgal
========

Implements classes and functions for handling buildgal files.

"""

from __future__ import with_statement # for py2.5                                                                                                                                    
from __future__ import division

from . import snapshot, array, util
from . import family, galaxy
from . import units
from . import config, config_parser

import ConfigParser

import struct, os
import numpy as np
import gzip
import sys
import warnings
import copy
import types
import math

_type_map = {}
for x in family.family_names() :
    try :
        pp =  [q.strip() for q in config_parser.get('scf-type-mapping',x).split(",")]
        qq=np.array(pp)
        _type_map[family.get_family(x)] = pp
    except ConfigParser.NoOptionError :
        pass
_rev_type_map = {}
for f,t in _type_map.iteritems() : 
    for it in t:
        if it in _rev_type_map: _rev_type_map[it].append(f)
        else                  : _rev_type_map[it]=[f]

_name_map, _rev_name_map = util.setup_name_maps('scf-name-mapping')
_translate_array_name = util.name_map_function(_name_map, _rev_name_map)

def scf_type(fam):
    if fam == None:
        tlist=set([])
        for f in _type_map.keys():
            tlist=tlist.union(_type_map[f])
        return list(tlist)
    else :
        return _type_map[fam]

def scf_typelen(sim,typ):
    Ntyp=scf_typeidx(sim,typ).sum()
    return Ntyp

def scf_typeidx(sim,typ):
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
        else: tidx[scf_typeidx(sim,t)]=False
    # Return index of particles in the type
    return tidx

class ScfSnap(snapshot.SimSnap):
    _basic_loadable_keys = set(['mass','pos','vel'])
    _class_inherited = ['_unitfile_name']
    # _basic_loadable_keys = set(['mass','x','y','z','vx','vy','vz'])

    def __init__(self, filename, **kwargs):
        super(ScfSnap,self).__init__()
        if filename:
            self._setup(filename=filename,**kwargs)

    def _setup(self, **kwargs):
        global config

        self.filename=kwargs.get('filename',None)
        self.unitfile=kwargs.get('unitfile',None)

        self._family_slice = {}

        self._loadable_keys = set([])
        self._family_keys=set([])
        self._family_arrays = {}
        self._arrays = {}
        self.properties = {}
        
        #Initialize unitfile & set units
        if self.filename and os.path.isfile(self.filename):
            #Read in file
            if config['verbose'] : print>>sys.stderr, "ScfSnap: loading ",filename
            from pysim.files.scf import rw_snap
            fdict=rw_snap('R',self._filename)
            self.header=fdict['header']
            self._num_particles=self.header['nbody']
            self._decorate()
            # Pars to correct variables
            for nm in fdict.keys():
                if nm=='header': self.header=fdict['header']
                else:
                    if nm in ['pot','pm2']: iarr=np.array(self.properties['G'])*fdict[nm]
                    else                  : iarr=fdict[nm]
                    tnm=_translate_array_name(nm,reverse=True)
                    #print nm,tnm
                    self[tnm]=iarr.view(array.SimArray)
                    self[tnm].set_default_units(quiet=True)
                    #print nm,tnm,self[tnm].units
                    
        if self.filename:
            #Set up _family_slice
            for x in _type_map :
                tags_t=_type_map[x]
                for itag_t in tags_t: 
                    if itag_t in self.filename:
                        self._family_slice[x] = slice(0,self._num_particles)

            #Set up _galaxy_slice
            for x in galaxy.galaxy_names():
                if self.filename.endswith(x[-1]):
                    self._galaxy_slice[galaxy.get_galaxy(x)]=slice(0,self._num_particles)

    def _init_file_units_system(self,unitfile=None):
        """Gets unit system from SCF unit file"""
        # Load units from file
        if unitfile and os.path.isfile(unitfile):
            del self.unitfile
            self.unitfile=unitfile
        if self.unitfile:
            param=self.unitfile
            vel_unit = units.Unit(str(param['UnitVelocity_in_cm_per_s'])+' cm s**-1')
            dist_unit = units.Unit(str(param['UnitLength_in_cm'])+' cm')
            mass_unit = units.Unit(str(param['UnitMass_in_g'])+' g')
        # Set to default units
        else:
            vel_unit = config_parser.get('scf-units', 'vel')
            dist_unit = config_parser.get('scf-units', 'pos')
            mass_unit = config_parser.get('scf-units', 'mass')
        # Set file units system
        return [vel_unit,dist_unit,mass_unit,"K"]

    # FILENAME
    @property
    def filename(self): return getattr(self,'_filename',None)
    @filename.setter
    def filename(self,val):
        if isinstance(val,str): self._filename=val

    # UNITFILE
    @property
    def unitfile(self):
        extlist=['units','scfunits']
        if not hasattr(self,'_unitfile_dict'):
            if not hasattr(self,'_unitfile_name'):
                from mmlutils.mmlfiles import search as search_file
                self._unitfile_name=search_file(self.filename,extlist,nlvlup=2)
            if self._unitfile_name:
                from pysim.files.scf import rw_units
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
            if len(headline.split()) in [2,3]: return True
            # check = ScfSnap(f, verbose=False)
            # del check
        except :
            raise
            return False
        return True

    @staticmethod
    def _write(self, filename=None, overwrite=False, outformat=False, **exkw) :
        """
        Write a scf snapshot file
        """
        # Format conversion
        if self.__class__ is not ScfSnap:
            if filename == None :
                raise Exception,"Please specify a filename to write a new file."

            # Construct header
            header=dict(nbody=getattr(self,'_num_particles'),
                        time=self.properties.get('time',0.))

        else:
            if not filename: filename=self._filename
            header=self.header

        if outformat: header.setdefault('bhmass',0.)

        # Prevent overwrite
        if not overwrite and os.path.isfile(filename):
            print 'File already exists and overwrite not set'
            print '    '+filename

        # Create dictionary of arrays for write out
        fdict={}
        for nm in ScfSnap._basic_loadable_keys:
            tnm=_translate_array_name(nm,reverse=True)
            fdict[nm]=np.array(self[tnm])
        fdict['header']=header

        # Write files
        from pysim.files.scf import rw_snap
        rw_snap('W',filename,fdict,overwrite=overwrite,flag_out=outformat,**exkw)



@ScfSnap.decorator
def do_units(sim): sim.init_file_units_system()

@ScfSnap.decorator
def do_properties(sim):
    list_prop=['time','G']
    prop0={}
    # Units
    vel_unit,dist_unit,mass_unit,temp_unit=sim._file_units_system[:]
    time_unit=array.SimArray(1.0,'s').in_unitsys([vel_unit,dist_unit,mass_unit]).units
    # If header initialize
    if hasattr(sim,'header'):
        h = sim.header
        prop0['time'] = array.SimArray(h['time'],time_unit)
        prop0['G'   ] = array.SimArray(6.672e-8,'cm cm**2 s**-2 g**-1').in_unitsys([vel_unit,dist_unit,mass_unit])
    else:
        prop0={p:None for p in list_prop}
    # Add parameters, but preserve existing ones
    for p in list_prop: sim.properties.setdefault(p,prop0[p])
