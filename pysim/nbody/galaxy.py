"""

galaxy
======

This module defines the Galaxy class which represents
galaxies of particles (e.g. primary, secondar).
New Galaxy objects are automatically registered so that
snapshots can use them in the normal syntax (snap.primary,
snap.secondar, etc).

In practice the easiest way to make use of the flexibility
this module provides is through adding more galaxies 
in your config.ini.

"""

from . import config_parser


_registry = []

def galaxy_names(with_aliases=False) :
    """Returns a list of the names of all particle galaxies.
    If with_aliases is True, include aliases in the list."""
    
    global _registry
    l = []
    for o in _registry :
        l.append(o.name)
        if with_aliases :
            for a in o.aliases : l.append(a)
    return l

def get_galaxy(name, create=False) :
    """Returns a galaxy corresponding to the specified string.  If the
    galaxy does not exist and create is False, raises ValueError. If
    the galaxy does not exist and create is True, an appropriate
    object is instantiated, registered and returned."""
    
    if isinstance(name, Galaxy) :
        return name

    name = name.lower() # or should it check and raise rather than just convert? Not sure.
    for n in _registry :
        if n.name==name or name in n.aliases :
            return n

    if create :
        return Galaxy(name)
    else :
        raise ValueError, name+" not a galaxy" # is ValueError the right thing here?

class Galaxy(object) :
    def __init__(self, name, aliases=[]) :
        if name!=name.lower() :
            raise ValueError, "Galaxy names must be lower case"
        if name in galaxy_names(with_aliases=True) :
            raise ValueError, "Galaxy name "+name+" is not unique"
        for a in aliases :
            if a!=a.lower() :
                raise ValueError, "Aliases must be lower case"

        self.name = name
        self.aliases = aliases
        _registry.append(self)

    def __repr__(self) :
        return "<Galaxy "+self.name+">"

    def __reduce__(self) :
        return get_galaxy, (self.name, True), {"aliases": self.aliases}

    def __iter__(self) :
        # Provided so a single galaxy can be treated as a list of galaxies
        yield self

    def __str__(self) :
	return self.name

    def __cmp__(self, other) :
	return cmp(str(self),str(other))


# Instantiate the default galaxies as specified
# by the configuration file

g = globals()
for f in config_parser.options('galaxies') :
    aliases = config_parser.get('galaxies', f)
    g[f] = Galaxy(f,map(str.strip,aliases.split(","))) 
