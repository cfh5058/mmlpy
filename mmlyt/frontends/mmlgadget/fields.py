"""
Fields specific to Streaming data

Author: Matthew Turk <matthewturk@gmail.com>
Affiliation: KIPAC/SLAC/Stanford
Homepage: http://yt-project.org/
License:
  Copyright (C) 2008-2011 Matthew Turk.  All Rights Reserved.

  This file is part of yt.

  yt is free software; you can redistribute it and/or modify
  it under the terms of the GNU General Public License as published by
  the Free Software Foundation; either version 3 of the License, or
  (at your option) any later version.

  This program is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  GNU General Public License for more details.

  You should have received a copy of the GNU General Public License
  along with this program.  If not, see <http://www.gnu.org/licenses/>.
"""

from yt.data_objects.field_info_container import \
    FieldInfoContainer, \
    NullFunc, \
    TranslationFunc, \
    FieldInfo, \
    ValidateParameter, \
    ValidateDataField, \
    ValidateProperty, \
    ValidateSpatial, \
    ValidateGridType
import yt.data_objects.universal_fields

KnownMMLGadgetFields = FieldInfoContainer()
add_mmlgadget_field = KnownMMLGadgetFields.add_field

MMLGadgetFieldInfo = FieldInfoContainer.create_with_fallback(FieldInfo)
add_field = MMLGadgetFieldInfo.add_field

# add_mmlgadget_field("mass",function=NullFunc,take_log=True,
#                     validators = [ValidateDataField("mass")],
#                     units=r"\rm{g}")

# KnownMMLGadgetFields["mass"]._projected_units =r'\rm{g}'

# def _Mass(field,data): return data["mass"]
# add_field("Mass",function=_Mass,take_log=True,
#           units=r'\rm{g}')

_particle_field_list = ["mass",
                        "position_x",
                        "position_y",
                        "position_z",
                        "momentum_x",
                        "momentum_y",
                        "momentum_z",
                        "angmomen_x",
                        "angmomen_y",
                        "angmomen_z",
                        "r",
                        "n",
                        "id"]        

for pf in _particle_field_list:
    pfunc = particle_func("particle_%s" % (pf))
    add_field("particle_%s" % pf, function=pfunc,
              validators = [ValidateSpatial(0)],
              particle_type=True)          

def _ParticleMassMsun(field, data):
    particles = data["particle_mass"].astype('float64')
    return particles/1.989e33

add_field("ParticleMass",
          function=_ParticleMass, validators=[ValidateSpatial(0)],
          particle_type=True)
add_field("ParticleMassMsun",
          function=_ParticleMassMsun, validators=[ValidateSpatial(0)],
          particle_type=True)
