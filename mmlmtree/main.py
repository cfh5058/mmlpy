#!/usr/bin/python
import sys,os,shutil,glob,copy,pprint,scipy,math,importlib
import numpy as np
from mmlutils import *
import simlist,simfile,simutil,simstat

####################################################################################################################################
# METHOD TO GET COSSIM OBJECT FROM STRING
def runtag2simobj(runtag=None,**exkw):
    """
    Returns the cossim object for a given runtag
    """
    runlist=mmlparam.par_taglist(simlist,'cossim')
    if runtag in runlist:
        simdict=mmlparam.loadpar(simlist,'cossim',runtag)
        simobj=cossim(**simdict)
    else:
        simobj=asksimobj()
    return simobj
def asksimobj(**simdict):
    """
    Returns an cossim object provided user input
    """
    simdict=mmlparam.parspar(simlist,'cossim',inpar=simdict,askuser=True)
    simobj=cossim(**simdict)
    return simobj

####################################################################################################################################
####################################################################################################################################
# COSSIM CLASS
class cossim(dict):
    """
    A dictionary subclass with simulation keys.
    """
    
    def __init__(self,**inkw):
        """
        Initializes mmlsim objects
        """
        # Get cossim dict
        self.update(mmlparam.parspar(simlist,'cossim',inpar=inkw))
        if 'runtag' not in self: self['runtag']=mmlparam.par2tag(simlist,'cossim',self)
        # Add file dictionary
        fdict=simfile.files(simstr=self)
        setattr(self,'fdict',fdict)

    def run(self,mtype,method=None,**method_kw):
        """
        Handles running different types of methods for runs
        """
        # Pars input
        mtype=mmlpars.mml_pars(mtype,list=simlist.LIST_METHTYP)
        # Pars general options
        if mtype=='general':
            if method not in simlist.DICT_METHODS[mtype]: method=mmlio.askselect('Select a valid {} method:'.format(mtype),simlist.DICT_METHODS[mtype])
            if   method=='get_simobj'  : 
                out=self
            elif method=='get_snapshot': out=self.get_snapshot(**method_kw)
            elif method=='get_halosnap': out=self.get_halosnap(**method_kw)
            elif method=='get_groups'  : out=self.get_groups(**method_kw)
            elif method=='get_intlist' : out=self.get_intlist(**method_kw)
            elif method=='plot_intlist': out=simstat.plotcuts(self,**method_kw)
            else: raise Exception('Invalid method: {}'.format(method))
        # Pars other options
        else:
            runmod=importlib.import_module('mmlmtree.'+mtype)
            out=runmod.run(self,method=method,**method_kw)
        # Return output
        return out

    def get_snapshot(self,snapnum=None):
        """
        Returns GADGET output data for a snapshot
        """
        import pNbody
        if not isinstance(snapnum,int): snapnum=mmlio.askquest('Enter a valid snapshot number:',default=0,dtype='int')
        fsnap=self.fdict['snapbase'].format(snapnum)
        nb=pNbody.Nbody(p_name=fsnap,ftype='gadget',cosmosim=True,unitsfile=self.fdict['gdpm'])
        return nb

    def get_halosnap(self,snapnum=None,nb=None,haloid=None):
        """
        Retrieves pNbody object containing only halo particle data
        """
        if not isinstance(snapnum,int): snapnum=mmlio.askquest('Enter a valid snapshot number:',default=0,dtype='int')
        if not isinstance(haloid,long): haloid=mmlio.askquest('Enter a valid halo ID:',default=1L,dtype='long')
        # Get group number
        groupnum=simutil.haloid2groupnum(self,snapnum=snapnum,haloid=haloid)
        if groupnum==-1: 
            mmlio.verbose('Haloid {} could not be found in groups data for snapshot {}.'.format(haloid,snapnum))
            return
        # Load groups
        grps=self.get_groups(snapnum=snapnum,inclparticles=True)
        # Get pNbody object if not provided
        if nb==None: nb=self.get_snapshot(snapnum=snapnum)
        halonb=nb.selectp(lst=grps['pids'][groupnum],from_num=True)
        return halonb

    def convert_pos(self,pos_gdt,z):
        """
        Converts GADGET positions to physical units
        """
        a=1./(1.+z)
        pos_phy=a*pos_gdt
        return pos_phy
    def convert_vel(self,vel_gdt,z):
        """
        Converts GADGET velocities to physical units
        """
        gdparam=mmlio.rwdict('R',self.fdict['gdpm'])
        a=1./(1.+z)
        omega_v=gdparam['OmegaLambda']
        omega_m=gdparam['Omega0']
        H_0=gdparam['HubbleParam']
        H_a=H_0*np.sqrt(omega_m/(a**3)+omega_v+(1.-omega_m-omega_v)/(a**2))
        # Velocities already peculiar
        vel_phy=a*vel_gdt#+H_a*r
#        vel_phy=a*np.sqrt(a)*vel_gdt#+H_a*r
        return vel_phy

    def get_groups(self,snapnum=None,inclparticles=None):
        """
        Returns the group dictionary for a snapshot
        """
        if not isinstance(snapnum,int): snapnum=mmlio.askquest('Enter a valid snapshot number:',default=0,dtype='int')
        inclparticles=mmlpars.mml_pars(inclparticles,default=False,type=bool)
        # Get file names
        ffofcat =self.fdict['grpfofcat' ].format(snapnum)
        fsubcat =self.fdict['grpsubcat' ].format(snapnum)
        ffofprop=self.fdict['grpfofprop'].format(snapnum)
        fsubprop=self.fdict['grpsubprop'].format(snapnum)
        if inclparticles:
            ftypes  =self.fdict['grptypes'  ].format(snapnum)
            fids    =self.fdict['grpids'    ].format(snapnum)
            fpos    =self.fdict['grppos'    ].format(snapnum)
            fvel    =self.fdict['grpvel'    ].format(snapnum)
        # Get data
        dfofcat=simfile.readgrpcat(ffofcat,grptyp='fof')
        dsubcat=simfile.readgrpcat(fsubcat,grptyp='sub')
        dfofprop=simfile.readgrpprop(ffofprop,grptyp='fof')
        dsubprop=simfile.readgrpprop(fsubprop,grptyp='sub')
        if inclparticles:
            dtypes=simfile.readgrptype(ftypes)
            dids=simfile.readgrpids(fids)
            dpos=simfile.readgrppos(fpos)
            if os.path.isfile(fvel): dvel=simfile.readgrpvel(fvel)
            else: dvel=None
#                nb=self.get_snapshot(snapnum=snapnum)
#                idx=np.zeros(dids['ids'].shape,dtype=np.int64)
#                for ip in range(len(dids['ids'])): idx[ip]=nb.getindex(dids['ids'][ip])
#                dvel={'Npart':dids['Npart'],'vel':nb.vel[idx,:]}
#                with open(fvel,'wb') as file:
#                    np.array(dvel['Npart'],dtype=np.int32  ).tofile(fvel)
#                    np.array(dvel['vel'  ],dtype=np.float32).tofile(fvel)
        # Create group dictionary
        g=dict(N    =dsubcat['len'],
               cm   =dsubprop['cm'  ],
               cmv  =dsubprop['cmv' ],
               mtot =dsubprop['mtot'],
               mgas =dsubprop['mgas'],
               Parent=np.zeros(dsubcat['Ngrp'],dtype=int),
               ppos=[],pvel=[],pids=[],ptype=[],
               ppos_phys=[],pvel_phys=[],pmass=[])
        for ihpar in simlist.LIST_HALOREAD:
            if ihpar not in g: g[ihpar]=[]
        ind=np.zeros(dfofcat['Ngrp'],dtype=np.int32)
        suma=0L
        for igrp in range(dfofcat['Ngrp']):
            ind[igrp]=suma
            suma=suma+dfofcat['subs'][igrp]
        g['Parent'][ind[ind<dsubcat['Ngrp']]]=1
        # Add particle data
        if inclparticles:
            g['Npart']=dtypes['Npart']
            for isub in range(dsubcat['Ngrp']):
                idx=np.arange(dsubcat['offset'][isub],dsubcat['offset'][isub]+dsubcat['len'][isub],dtype=np.int32)
                g['pids' ].append(dids['ids'][idx])
                g['ppos' ].append(dpos['pos'][idx,:])
                g['ptype'].append(dtypes['type'][idx])
                if dvel!=None: 
                    g['pvel' ].append(dvel['vel'][idx,:])
                    # Converted mass & position
                    z=self.fdict['redlist'][snapnum]
                    m=np.ones(g['pids'][isub].shape,dtype=np.float32)*g['mtot'][isub]/g['N'][isub]
                    r=self.convert_pos(g['ppos'][isub],z)
                    v=self.convert_vel(g['pvel'][isub],z)
#                    rcom=self.convert_pos(g['cm'][isub],z)
#                    vcom=self.convert_vel(g['cmv'][isub],z)
                    rcom,vcom=mmlmath.center_mass(m,r,vel=v,vflag=True)
                    r=(r-rcom).astype(np.float32)
                    v=(v-vcom).astype(np.float32)
                    # Wrap positions
                    rsign=np.ma.masked_equal(r,0)/np.ma.masked_equal(np.abs(r),0)
                    rsign=rsign.filled(1.)
                    rwrap=(np.abs(r)>(1000.*self['boxsize']/2.))
                    r[rwrap]=rsign[rwrap,:]*(np.abs(r[rwrap,:])-1000.*self['boxsize'])
                    # Add mass, pos, & vel to dictionary
                    g['pmass'].append(m)
                    g['ppos_phys'].append(r)
                    g['pvel_phys'].append(v)
                    # Angular momentum
                    Lxyz=np.sum(mmlmath.angmom(m,r,v),axis=0)
                    g['Lx'].append(Lxyz[0])
                    g['Ly'].append(Lxyz[1])
                    g['Lz'].append(Lxyz[2])
                    # Spin
                    M=g['mtot'][isub]
                    G=43007.1 # GADGET units
                    R=mmlmath.absvect(r).max()
                    h=np.sqrt(np.sum(Lxyz**2))/M
                    V=np.sqrt(G*M/R)
                    g['spin_tot'].append(mmlmath.spin(h,V,R))
                    g['htot'].append(h)
                    g['Vtot'].append(V)
                    g['G'].append(G)
        # Return dictionary
        return g

    def get_intlist(self,askuser=True,**exkw): return simstat.mkcuts_added(self,askuser=askuser,**exkw)

####################################################################################################################################
# WALK THROUGH METHOD
def walk(mtype=None,runtag=None,verbose=True,**extra_kw):
    """
    Method to walk users through performing different run methods 
    """
    # Get method type if not provided
    mtypeLIST=simlist.LIST_METHTYP
    if mtype not in mtypeLIST: mtype=mmlio.askselect('Select a type of method:',mtypeLIST)
    # Get run info
    if verbose: mmlio.verbose('Getting info on the run...')
    simobj=runtag2simobj(runtag,verbose=verbose,**extra_kw)
    # Call run
    if verbose: mmlio.verbose('Beginning {} for {} simulation...'.format(mtype.upper(),simobj['runtag']),border=True)
    output=simobj.run(mtype,verbose=verbose,**extra_kw)
    if verbose: print [mtype,simobj['runtag']]
    if output==None: return
    else           : return output


if __name__ == '__main__':
    walk()

        
