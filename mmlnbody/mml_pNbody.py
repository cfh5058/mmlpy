####################################################################################################################################
#
# MEAGAN LANG'S PNBODY CLASS
#
####################################################################################################################################
import math,os,shutil,copy,cPickle,sys,pNbody,types,pprint
import numpy as np
import matplotlib as mplib
#import astropysics.models as apymods
#import astropysics.phot as apyphot
from mpl_toolkits.axes_grid1 import make_axes_locatable
import ImageDraw,ImageFont
import pNbody
#from pNbody import libutil
#print dir(pNbody)
import pNbody.libutil as libutil
import pNbody.mapping as mapping
from pNbody.libutil import getval
from mmlutils import *
from mmlastro import *
import simlist,simplot,simcalc,simprof

class mmlNbodyDefault:

  def __init__(self,p_name=None,pos=None,vel=None,mass=None,num=None,tpe=None,ftype=None,status='old',
               byteorder=sys.byteorder,pio='no',local=False,log=None,unitsfile=None,cosmosim=None,
               makefile=None,scffiledict=None,snapfile0=None,idgal2=None,ptype=None,verbose=None,**nbody_kw):

    # Set constants
    self.typList=simlist.LIST_PTYPBASE
    self.typDict={self.typList[ityp]:ityp for ityp in range(len(self.typList))}
    self.cosmosim=mmlpars.mml_pars(cosmosim,default=False,type=bool)
    self.verbose=mmlpars.mml_pars(verbose,default=False,type=bool)
    self.scffiledict=scffiledict
    self.snapfile0=snapfile0
    if self.verbose: print '[mml_pNbody.mmlNbodyDefault.__init__]'
    # Files
    unitexts={'gadget'  :['param','gdpar'],
              'scf'     :['units','scfunits'],
              'buildgal':['units','bgunits']}
    makeexts={'gadget'  :['Makefile','makefile','gdmk']}
    if p_name:
      if not unitsfile and ftype in unitexts: unitsfile=mmlfiles.search(p_name,unitexts[ftype],nlvlup=1)
      if not makefile  and ftype in makeexts: makefile =mmlfiles.search(p_name,makeexts[ftype],nlvlup=1)
    # GADGET files
    if ftype=='gadget':
      if unitsfile: self.read_gparams(unitsfile)
      if makefile:  self.read_makefile(makefile)
    # SCF files
    elif ftype=='scf':
      ptype=mmlpars.mml_pars(ptype,default='disk',type=str,list=self.typList)
    # Set flags
    self.update_flags(flag_init=True,**nbody_kw)
    # Initialize base class
    pNbody.main.NbodyDefault.__init__(self,p_name=p_name,pos=pos,vel=vel,mass=mass,num=num,tpe=tpe,
                                      ftype=ftype,status=status,byteorder=byteorder,pio=pio,local=local,
                                      log=log,unitsfile=unitsfile)
    # Get SCF data for GADGET sims
#    if   ftype=='gadget': self.read_scffiles(scffiledict)
#    elif ftype=='scf'   : pass
#    else                : pass
    # Array dependent things
    if len(self.num)>0:
      # Add galaxy assignment array
      if idgal2==None:
        if hasattr(self,'idgal2'): idgal2=self.idgal2
        else                     : idgal2=self.get_idgal2()
      self.set_galaxy(idgal2=idgal2)
      # Galaxy stuff
      if not self.cosmosim:
        # Center on primary galaxy
        com,vcom=self.get_center(method='mass',centergal=1,centertyp='all',move=True,vcomflag=True)

  def families(self): return ['star']+self.get_typlist()
  def galaxies(self): return ['galaxy{}'.format(g) for g in list(np.unique(self.gal))]
  def __len__(self): return self.nbody

  def listpar(self,partyp,**exkw):
    from pysim.simlist import listpar as uplistpar
    exkw['famlist']=self.families()
    exkw['gallist']=self.galaxies()
    exkw['comlist']=simlist.LIST_COMMETHS
    return uplistpar(partyp,**exkw)
   
  def par2tag(self,partyp,inpar,**exkw):
    from pysim.simlist import par2tag as uppar2tag
    exkw['famlist']=self.families()
    exkw['gallist']=self.galaxies()
    exkw['comlist']=simlist.LIST_COMMETHS
    return uppar2tag(partyp,inpar,**exkw)

  def update_flags(self,**flag_kw):
    pass

  def init_units(self):
    """
    Reads info from parameter file
    """
    from pNbody import io,param,parameters
    self.unitsparameters = param.Params(parameters.UNITSPARAMETERFILE,None)
    if self.unitsfile!=None:
      try:
        gparams = io.read_params(self.unitsfile)
        self.unitsparameters.set('UnitLength_in_cm',        gparams['UnitLength_in_cm'])
        self.unitsparameters.set('UnitMass_in_g',           gparams['UnitMass_in_g'])
        self.unitsparameters.set('UnitVelocity_in_cm_per_s',gparams['UnitVelocity_in_cm_per_s'])
        self.unitsparameters.set('filename',self.unitsfile)
      except:
        try:
          self.unitsparameters = param.Params(self.unitsfile,None)
        except:
          raise IOError(015,'format of unitsfile %s unknown ! Pease check.'%(self.unitsfile))

    # define local system of units it it does not exists
    self.set_local_system_of_units()

  def open_and_read(self,name,readfct):
    """
    Open and read file
    """
    from pNbody import io,mpi

    # check p_name
    if self.pio=='yes' or mpi.mpi_IsMaster():
      io.checkfile(name)

    # get size
    if self.pio=='yes' or mpi.mpi_IsMaster():
      isize = os.path.getsize(name)

    # open file
    if self.pio=='yes' or mpi.mpi_IsMaster():
      f = open(name,'r')
    else:
      f = None


    # read the file
    readfct(f)

    if self.pio=='yes' or mpi.mpi_IsMaster():
      fsize = f.tell()
    else:
      fsize = None

    if self.pio=='yes' or mpi.mpi_IsMaster():
      if fsize != isize and self.verbose: mmlio.verbose('File is not read completely.')

    # close file
    if self.pio=='yes' or mpi.mpi_IsMaster():
      f.close()

  def index_norepeats(self,cut=0.):
    """
    Returns the indicies of repeated particle positions
    """
    a=np.lexsort(self.pos.T)
    diff=np.ones(self.nbody)
    diff[1:]=np.sqrt((np.diff(self.pos[a,:],axis=0)**2).sum(axis=1))
    idx=np.sort(a[diff>cut])
    return idx

  def countrepeats(self,poscut=0.,velcut=0.):
    """
    Checks if there are repeated particle positions
    """
    dtype=np.float64
    pos=self.pos.astype(dtype)
    vel=self.vel.astype(dtype)
    apos=np.lexsort(pos.T)
    avel=np.lexsort(vel.T)
    diff_pos=np.ones(self.nbody,dtype=dtype) ; diff_vel=np.ones(self.nbody,dtype=dtype)
    diff_pos[:-1]=np.sqrt((np.diff(pos[apos,:],axis=0)**2).sum(axis=1))
    diff_vel[:-1]=np.sqrt((np.diff(vel[avel,:],axis=0)**2).sum(axis=1))
    diff=np.ones(diff_pos.shape,dtype=dtype)
    diff[diff_pos<=poscut]=0
#    diff[np.logical_and(diff_pos<=poscut,diff_vel<=velcut)]=0
    nr=(diff==0).sum()
    print 'nr={}'.format(nr)
    diff2=np.diff(np.hstack(([1],diff,[1])))
    idx_beg,=np.where(diff2<0)
    idx_end,=np.where(diff2>0)
    rpos=self.pos[apos,:][diff_pos<=poscut,:]
    rvel=self.vel[avel,:][diff_vel<=velcut,:]
    print 'nr pos={}'.format(rpos.shape[0])
    print 'nr vel={}'.format(rvel.shape[0])
    print rpos
    print idx_end-idx_beg
    print len(idx_beg),len(idx_end)
    return nr

  def fitprofile(self,ptyp=None,pgal=None,cmeth=None,nbin=None,coord=None,**exkw):
    """
    Fits profile
    """
    # Set constants
    typList=self.get_typlist()
    typeprof={'halo':'NFW','disk':'FLATDISK','bulge':'HERNQUIST'}
    for ityp in typList: 
      if ityp not in typeprof: typeprof[ityp]='HERNQUIST'
    # Pars input
    ptyp=mmlpars.mml_pars(ptyp,default='disk',list=typList)
    pgal=mmlpars.mml_pars(pgal,default=1,list=self.get_gallist())
    nbin=mmlpars.mml_pars(nbin,default=100,type=int)
    coord=mmlpars.mml_pars(coord,default='r',list=['r','z'])
    profid=exkw.pop('profid',typeprof[ptyp])
    # Get index
    pidx=self.get_index(ptyp=ptyp,pgal=pgal)
    com=self.get_center(method=cmeth,idxcent=pidx)
    # Get histogram
    if   coord=='r':
      if ptyp=='disk': x=mmlmath.xyz2rxy(self.pos[pidx,:]-com)
      else           : x=mmlmath.xyz2rxyz(self.pos[pidx,:]-com)
    elif coord=='z': x=np.abs((self.pos[pidx,:]-com)[:,2])
    m=self.mass[pidx,:]
    mhist,xbins=np.histogram(x,bins=nbin,weights=m)
    # Get data to fit
    menc=np.cumsum(mhist)
    xenc=xbins[1:]
    # Fit
    exkw[coord]=xenc
    profpar=mmlprofiles.fit(profid,menc,**exkw)
    # Return parameters
    return profpar

  def printinfo(self):
    """
    Prints basic info on simulation to the screen
    """
    # Set constants
    typList=simlist.LIST_PTYPBASE
    mtotarr=np.array(self.massarr)*np.array(self.npart)
    rsoft=mmlio.askquest('At what radius (in kpc) should the mean interparticle separation be taken?',dtype='float',default=4.0)
#    rsoft=4.0
    rbuff=0.1
    zbuff=0.1
    massfact0=(10.**10)
    typwid=len(max(typList,key=len))
#    R1_kpc=apymods.NFWModel().Mvir_to_Rvir(M1_msol)
    massunit='10^12 Msol' ; massfact=massfact0/(10.**12)
    # File name and time
    print self.p_name[0]
    print 't = {} Gyr'.format(self.atime)
    # Center
    com0=self.get_center(centergal=1,centertyp='disk',method='dens',move=True)
    com=self.get_center(centergal=0,centertyp='all',method='mass',move=False)
    print 'COM = {} kpc'.format(com)
    for idx in range(len(typList)):
      if self.npart[idx]>0:
        icom=self.get_center(centergal=0,centertyp=typList[idx],method='mass',move=False)
        print '    {}: {}'.format(typList[idx].ljust(typwid),icom)
    # Print number of galaxies
    ngal=self.get_ngal()
    print 'NGAL = {}'.format(ngal)
    # SCF STUFF
    try:
      idxbox=np.abs(self.get_var('rxyz'))<15.
      scfvec=np.ma.masked_invalid(self.get_var('scfpotm2rel',idxtot=idxbox))
      scfidx=np.argmax(scfvec)
      posvec=self.pos[idxbox,:]
      scfpos=posvec[scfidx,:]
      scfmax=scfvec[scfidx]
      print 'SCFSTATS:'
      print '    SCFMAX={}'.format(scfmax)
      print '    SCFPOS={}'.format(scfpos)
      print '    SCFRAD={}'.format(np.sqrt(np.sum(scfpos*scfpos)))
    except:
      raise
      print 'NO SCF DATA AVAILABLE'
    # Particle number
    print 'PARTICLE NUMBER ({}):'.format(self.nbody)
    for idx in range(len(typList)):
      if self.npart[idx]>0:
        print '    {}: {}'.format(typList[idx].ljust(typwid),self.npart[idx]/ngal)
    # Particle mass
    print 'PARTICLE MASSES (IN {}):'.format(massunit)
    for idx in range(len(typList)):
      if self.npart[idx]>0:
        print '    {}: {}'.format(typList[idx].ljust(typwid),self.massarr[idx]*massfact)
    # Component mass
    print 'COMPONENT MASSES ({}):'.format(mtotarr.sum()*massfact/ngal)
    for idx in range(len(typList)):
      if self.npart[idx]>0:
        print '    {}: {}'.format(typList[idx].ljust(typwid),mtotarr[idx]*massfact/ngal)
    # Interparticle separations
    print 'MEAN INTER-PARTICLE SEPARATIONS (@ {} kpc):'.format(rsoft)
    for idx in range(len(typList)):
      if self.npart[idx]>0:
        cylindric=(typList[idx]=='disk')
        intersep=self.get_intersep(rsoft,ptyp=typList[idx],rbuff=rbuff,zbuff=zbuff,cylindric=cylindric)
        print '    {}: {}'.format(typList[idx].ljust(typwid),intersep)
    # Particles encloed
    print 'NUMBER OF PARTICLES ENCLOSED (@ {} kpc):'.format(rsoft)
    for idx in range(len(typList)):
      if self.npart[idx]>0:
        cylindric=(typList[idx]=='disk')
        nin,vin=self.get_ninside(rsoft,ptyp=typList[idx],cylindric=cylindric)
        print '    {}: {} (vol={})'.format(typList[idx].ljust(typwid),nin,vin)
    # Softenings
    print 'PARTICLE SOFTENINGS:'
    sftlist=self.get_softenings()
    for idx in range(len(typList)):
      if self.npart[idx]>0:
        print '    {}: {} kpc'.format(typList[idx].ljust(typwid),sftlist[idx])
    # Component spin
    print 'COMPONENT SPIN (LAMBDA):'
    for idx in range(len(typList)):
      if self.npart[idx]>0:
        lspin=self.get_spin(ptyp=typList[idx],pgal=1)
        print '    {}: {} '.format(typList[idx].ljust(typwid),lspin)
    # Component stability
#    stabmeth='toomre81'
    stabmeth='efstat82'
#    stabmeth='op73'
    print 'COMPONENT STABILITY ({})'.format(stabmeth)
    for idx in range(len(typList)):
      if self.npart[idx]>0:
        stabcrit=self.get_stability(ptyp=typList[idx],pgal=1,method=stabmeth)
        print '    {}: {} '.format(typList[idx].ljust(typwid),stabcrit)
    return

  ##################################################################################################################################
  # METHOD TO RETURN A DICTIONARY OF NBODY SIMULATION STATISTICS
  def get_statlist(self,pgal=1,centergal=1,centertyp='all'):
    try:
      slist=self.statlist
    except:
      # Center galaxy
      oldcm=self.get_center(centergal=centergal,centertyp=centertyp,move=True)
      # Calculate quatitites
      cm=self.cm() ; cv=self.cv()
      haloModel={'rho0':0.,'rc':0.}
      diskModel={'A':0.,'l':0.,'h':0.}
      bulgeModel={'A':0.,'r0':0.}
#      haloModel=self.get_haloprop('nfw',pgal=pgal)
#      diskModel=self.get_diskprop('sechdisk',pgal=pgal)
#      bulgeModel=self.get_bulgeprop('hernquist',pgal=pgal)
      barModel=self.get_barprop(pgal=pgal,method='dens')
      if hasattr(self,'scfpot') and hasattr(self,'scfpotm2'):
        barModel_scf=self.get_barprop(pgal=pgal,method='scfm2')
        barModel['m2']=barModel_scf['A']
      else:
        barModel['m2']=0.
      eps=max(self.get_softenings())
      # Form dictionary
      slist={
        'fname':self.p_name[0],'time':self.atime,
        'comx' :cm[0],'comy' :cm[1],'comz' :cm[2],
        'vcomx':cv[0],'vcomy':cv[1],'vcomz':cv[2],
        'clause':self.get_clause(),'sigma':self.get_sigmar(),
        'halo_A':haloModel['rho0'],'halo_r':haloModel['rc'],
        'disk_A':diskModel['A'],'disk_r':diskModel['l'],'disk_z':diskModel['h'],
        'bulge_A':bulgeModel['A'],'bulge_r':bulgeModel['r0'],
        'bar_phi':barModel['pa'],'bar_amp':barModel['A'],'bar_x':barModel['x'],'bar_y':barModel['y'],'bar_z':barModel['z'],
        'bar_m2':barModel['m2'],
        'barstat':barModel,
        }
#        'ekin':self.Ekin(),'epot':self.Epot(eps),
      # Add dictionary to pNbody object for later use
      self.statlist=slist
      
    return slist

####################################################################################################################################
# METHODS THAT MODIFY THE PNBODY OBJECT
####################################################################################################################################

  ##################################################################################################################################
  # METHOD TO DETERMINE GAL ARRAY
  def set_galaxy(self,idgal2=None):
    # Preallocate
    if idgal2==None:
      if hasattr(self,'idgal2'): idgal2=self.idgal2
      else: idgal2=self.get_idgal2()
    self.idgal2=idgal2
    self.gal=np.ones(self.num.shape,dtype=np.int)
    # Fill with galaxy 2
    if isinstance(idgal2,long) or isinstance(idgal2,int):
      self.gal[self.num >= idgal2]*=2
    return

  ##################################################################################################################################
  # METHOD TO PUT GALAXY ON AN ORBIT
  def add_orbit(self,M1_msol,M2_msol,R1_pc=None,R2_pc=None,rinit=1.,rperi=0.5,ecc=1.5,incl=0.,**exkw):
    # Set constants
#      R1_kpc=apymods.NFWModel().Mvir_to_Rvir(M1_msol) ; R1_pc=1000.*R1_kpc
#      R2_kpc=apymods.NFWModel().Mvir_to_Rvir(M2_msol) ; R2_pc=1000.*R2_kpc
    gconst=4.302e-3 # pc (km/s)**2 / Msol
    exkw['units']='galaxy' # pc/Msol
    # Get virial radii
    if R1_pc == None: R1_pc=mmlprofiles.mvir2rvir(M1_msol,**exkw) ; R1_kpc=R1_pc/1000.
    else: R1_kpc=R1_pc/1000.
    if R2_pc == None: R2_pc=mmlprofiles.mvir2rvir(M2_msol,**exkw) ; R2_kpc=R2_pc/1000.
    else: R2_kpc=R2_pc/1000.
    # Give input units based on virial radii
    rtot_pc=R1_pc+R2_pc
    rini_pc=rinit*rtot_pc
    rper_pc=rperi*rtot_pc
    # Determine intial position/velocity of orbit
    xorb_pc=np.array([rini_pc,0.,0.])
    vorb_kms,phiorb=mmlmath.orbit(M1_msol,M2_msol,rini_pc,rper_pc,ecc,gconst)
    # Shift particle positions and velocities to be on orbit
    self.pos[:,:]+=xorb_pc/self.localsystem_of_units.convertionFactorTo(pNbody.units.Unit_pc)
    self.vel[:,:]+=vorb_kms/self.localsystem_of_units.convertionFactorTo(pNbody.units.Unit_kms)
    # Incline the orbit
    self.rotate(angle=incl,axis=[0,-1,0])
    # Force pericenter to negative x-axis
    self.rotate(angle=phiorb,axis=[0,0,1])
    # Print some info
    verbose=True
#    verbose=self.verbose
    if verbose:
        if ecc == 0:              print 'CIRCULAR orbit'
        elif ecc > 0 and ecc < 1: print 'ELLIPTICAL orbit' 
        elif ecc == 1:            print 'PARABOLIC orbit'
        elif ecc > 1:             print 'HYPERBOLIC orbit'
        print '    M1    = {:4.2e} Msol'.format(M1_msol)
        print '    M2    = {:4.2e} Msol'.format(M2_msol)
        print '    R1    = {:5f} kpc '.format(R1_kpc)
        print '    R2    = {:5f} kpc '.format(R2_kpc)
        print '    rinit = {:5f} kpc '.format(rini_pc/1000.)
        print '    rperi = {:5f} kpc '.format(rper_pc/1000.)
        print '    ecc   = {:2.2f}     '.format(ecc)
        print '    pos   = [{:5f},{:5f},{:5f}] kpc '.format(xorb_pc[0]/1000.,xorb_pc[1]/1000.,xorb_pc[2]/1000.)
        print '    vel   = [{:5f},{:5f},{:5f}] km/s'.format(vorb_kms[0],vorb_kms[1],vorb_kms[2])
        print '    phi   = {:4.1f} deg '.format(phiorb*(180./math.pi))
        print '    theta = {:4.1f} deg '.format(incl*(180./math.pi))
    # Return control
    return

  ##################################################################################################################################
  # METHOD TO KINEMATICALLY HEAT THE GALAXY'S DISK BY SOME FACTOR
  def heatdisk(self,heatfact):
    diskindex=2
    print '[pNbody.heatdisk] Check that the disk is being heated.'
    print 'before: {}'.format(self.vel[self.tpe==diskindex,2][0])
    self.vel[self.tpe==diskindex,2]*=heatfact
    print 'after:  {}'.format(self.vel[self.tpe==diskindex,2][0])
#    nbdisk=self.select('disk')
#    for inum in nbdisk.num:
#      self.pos[self.getindex(inum),2]*=heatfact
#    nbdisk.pos[:,2]*=heatfact

  ##################################################################################################################################
  # METHOD THAT CONVERTS PNBODY UNITS
  def convert_units(self,unitDict=None):
    # Set constants
    unitList=[
      'UnitDensity',
      'UnitEnergy',
      'UnitLength',
      'UnitMass',
      'UnitSpecEnergy',
      'UnitSurfaceDensity',
      'UnitTime',
      'UnitVelocity']
    pureUnitList=[
      'UnitLength',
      'UnitMass',
      'UnitTime']
    unitForm={iunit:mmlpars.parsdict(default=getattr(self.localsystem_of_units,iunit)) for iunit in unitList}
    unitDict=mmlpars.mml_formpars(unitDict,unitForm)
#    newsystem_of_units=pNbody.units.UnitSystem('new',unitDict.values())
    unitVars=dict.fromkeys(unitList,[])
#    unitVars['UnitDensity'   ]=['Density','rho']
    unitVars['UnitLength'    ]=['pos','hsml','boxsize']#,'Hsml']
    unitVars['UnitMass'      ]=['mass','mass_tot','massarr']
    unitVars['UnitSpecEnergy']=['pot','u']
    unitVars['UnitTime'      ]=['atime','tstp']
    unitVars['UnitVelocity'  ]=['vel']
#    unitVars['UnitVelocity_Time']=['acce']
#    unitVars['UnitEnergy_Temp_Time']=['endt']
    # Get a list of arrays
#    list_arr=self.get_list_of_array() ; print '[pNbody.gadget.convert_units] List of arrays: {}'.format(list_arr)
#    list_var=self.get_list_of_vars()  ; print '[pNbody.gadget.convert_units] List of vars: {}'.format(list_var)
    # Loop over units
    for iunit in unitList:
      convert_factor=self.localsystem_of_units.convertionFactorTo(unitDict[iunit])
      for iattr in unitVars[iunit]:
        if hasattr(self,iattr):
          ivar=getattr(self,iattr)
          if ivar != None:
            if isinstance(ivar,list):
              for iele in range(len(ivar)): ivar[iele]*=convert_factor
            else:
              ivar*=convert_factor
            setattr(self,iattr,ivar)
      setattr(self.localsystem_of_units,iunit,unitDict[iunit])
    # Assign new system of units
#    setattr(self,'localsystem_of_units',newsystem_of_units)
    return unitDict
    
  ##################################################################################################################################
  # METHOD THAT ADDS TWO PNBODY OBJECTS
  def append(self,solf,preserve_num=False,newgal=True,do_not_sort=False):
    # Ensure that objects match
    if solf.ftype != self.ftype:
      raise Exception('Append Error: file types do not match. ({},{})').format(self.ftype,solf.ftype)
    if solf.get_list_of_array() != self.get_list_of_array():
      raise Exception('Append Error: arrays do not match.')
    if len(solf.npart) != len(self.npart):
      raise Exception('Append Error: # of particle types do not match. ({},{})').format(len(self.npart),len(solf.npart))
    # Loop over types
    self_npart=self.npart
    solf_npart=solf.npart
    # Adjust IDs
    if not preserve_num:
      self.get_num()
      solf.get_num()
    if newgal: idgal2=self.get_idgal2(preserve_num)
    else:      idgal2=self.num.max()+1
    solf.num+=idgal2
    self.idgal2=idgal2
    # Adjust particle galaxy tags
    if newgal: solf.gal+=1
    # Loop over arrays combining
    names=self.get_list_of_array()
    for name in names:
      vec1=getattr(self,name)
      vec2=getattr(solf,name)
      vec=np.concatenate((vec1,vec2))
      setattr(self,name,vec)
    # Sort by particle type
    if do_not_sort:
      pass
    else:
      #self.sort_type()
      sequence = self.tpe.argsort()
      for name in names:
        vec=getattr(self,name)
        vec=np.take(vec,sequence,axis=0)
        setattr(self,name,vec)
    # Add header info
    self.nbody=self.nbody+solf.nbody
    self.npart=self.npart+solf.npart
    self.npart_tot = self.get_npart_tot()
    self.nbody_tot = self.get_nbody_tot()
    self.init()
    # Set galaxy ID
    self.idgal2=idgal2
    self.set_galaxy(idgal2=self.idgal2)

####################################################################################################################################
# OVERLOADED PNBODY METHODS
####################################################################################################################################

  ##################################################################################################################################
  # METHOD FROM PNBODY GADGET CLASS FOR SELECTING PARTICLES BY TYPE
  def select(self,*arg,**kw):
    """ 
    Return an N-body object that contain only particles of a 
    certain type, defined by in gadget:

    gas		: 	gas particles
    halo	: 	halo particles
    disk	:	disk particles
    bulge	:	bulge particles
    stars	:	stars particles
    bndry	:	bndry particles
    
    sph         : 	gas with u > u_c
    sticky      : 	gas with u < u_c
    """   
    index = {'gas':0,'halo':1,'disk':2,'bulge':3,'stars':4,'bndry':5,'stars1':1,'halo1':2}
    # this allows to write nb.select(('gas','disk'))
    if len(arg)==1:
      if type(arg[0])==types.TupleType:
        arg = arg[0]
    tpes = arg    
    # create the selection vector
    c = np.zeros(self.nbody)
    # add different types
    for tpe in tpes:
      # Handle string type codes
      if type(tpe) == types.StringType:
        # Handle specialized type codes
      	if   (tpe=='sph'):
	  c = c+(self.u>self.critical_energy_spec)*(self.tpe==0)
      	elif (tpe=='sticky'):
	  c = c+(self.u<self.critical_energy_spec)*(self.tpe==0)
      	elif (tpe=='diskbulge' or tpe=='star'):
	  c = c+(self.tpe==2)+(self.tpe==3)
      	elif (tpe=='wg'):
	  c = c+(self.u>0)*(self.tpe==0)
      	elif (tpe=='cg'):
	  c = c+(self.u<0)*(self.tpe==0)
      	elif (tpe=='all'):
	  return self
        # Handle galaxies
        elif 'gal' in tpe:
          c = c+(self.gal==int(float(tpe[-1])))
        # Handle non-existant type
        elif not index.has_key(tpe):
          print "unknown type, do nothing %s"%(tpe)
	  return self
        # Handle standard type
	else:
          i = index[tpe]
          c = c+(self.tpe==i)
      # Handle integer type codes
      elif type(tpe) == types.IntType:	
	c = c+(self.tpe==tpe)
    # Return particles selected by conditional array
    return self.selectc(c)
      
  ##################################################################################################################################
  # METHOD TO OVERLOAD SELECTC FOR PNBODY OBJECTS (COPYS FILES)
  def selectc(self,c,local=False):
    """
    Return an N-body object that contain only particles where the
    corresponding value in c is not zero.
    c is a nx1 Nbody array.

    c      : the condition vector
    local  : local selection (True) or global selection (False)
    """
    new = pNbody.Nbody(status='new',ftype=self.ftype[6:],local=local,unitsfile=self.unitsfile,
                       idgal2=self.idgal2,cosmosim=self.cosmosim)
    # now, copy all var linked to the model
    for name in self.get_list_of_vars(): setattr(new, name, getattr(self,name))
    # now, copy and compress all array linked to the model
    for name in self.get_list_of_array():
      vec = getattr(self,name)
      setattr(new, name, np.compress(c,vec,axis=0))
    # other vars
    new.init()
    return new

  ##################################################################################################################################
  # METHOD TO OVERLOAD SELECTI FOR PNBODY OBJECTS (COPYS FILES)
  def selecti(self,i,local=False):
    """
    Return an N-body object that contain only particles having
    their index (not id) in i.

    i      : vector containing indexes
    local  : local selection (True) or global selection (False)
    """
    new = pNbody.Nbody(status='new',ftype=self.ftype[6:],local=local,unitsfile=self.unitsfile)
    # now, copy all var linked to the model
    for name in self.get_list_of_vars(): setattr(new, name, getattr(self,name))
    # now, copy and compress all array linked to the model
    for name in self.get_list_of_array():
      vec = getattr(self,name)
      setattr(new, name, vec[i])
    # other vars
    new.init()
    return new

  ##################################################################################################################################
  # METHOD TO OVERLOAD set_local_system_of_units FOR PNBODY OBJECTS (ALLOWS SILENT RUN MODE)
  def set_local_system_of_units(self,params=None,UnitLength_in_cm=None,UnitVelocity_in_cm_per_s=None,UnitMass_in_g=None,
                                unitparameterfile=None,gadgetparameterfile=None):
    """
    Set local system of units using UnitLength_in_cm,UnitVelocity_in_cm_per_s,UnitMass_in_g
      1)  if nothing is given, we use self.unitsparameters to obtain these values
      2)  if UnitLength_in_cm,UnitVelocity_in_cm_per_s,UnitMass_in_g are given, we use them
      2b) if UnitLength_in_cm,UnitVelocity_in_cm_per_s,UnitMass_in_g are given in a dictionary
      3)  if unitparameterfile   is given we read the parameters from the file (units parameter format)
      4)  if gadgetparameterfile is given we read the parameters from the file (gadget param format)
    """
    # GADGETPARAMETERFILE
    if gadgetparameterfile!=None:
      params = io.read_params(gadgetparameterfile)
      if self.verbose: print "Units Set From %s"%gadgetparameterfile
    # UNITPARAMETERFILE
    elif unitparameterfile!=None:
      unitsparameters = pNbody.param.Params(unitparameterfile,None)
      params = {}
      params['UnitLength_in_cm']         = unitsparameters.get('UnitLength_in_cm')
      params['UnitVelocity_in_cm_per_s'] = unitsparameters.get('UnitVelocity_in_cm_per_s')
      params['UnitMass_in_g']            = unitsparameters.get('UnitMass_in_g')
      if self.verbose: print "Units Set From %s"%unitparameterfile
    # PARAMS DICTIONARY
    elif params!=None:
      if self.verbose: print "Units Set From %s"%params
    # INDIVIDUAL UNITS
    elif UnitLength_in_cm!=None and UnitVelocity_in_cm_per_s!=None and UnitMass_in_g!=None:
      params = {}
      params['UnitLength_in_cm']         = UnitLength_in_cm
      params['UnitVelocity_in_cm_per_s'] = UnitVelocity_in_cm_per_s
      params['UnitMass_in_g']            = UnitMass_in_g
      if self.verbose: print "Units Set From UnitLength_in_cm,UnitVelocity_in_cm_per_s,UnitMass_in_g"
    # DEFAULT TO STANDARD VALUES
    else:
      params = {}
      params['UnitLength_in_cm']         = self.unitsparameters.get('UnitLength_in_cm')
      params['UnitVelocity_in_cm_per_s'] = self.unitsparameters.get('UnitVelocity_in_cm_per_s')
      params['UnitMass_in_g']            = self.unitsparameters.get('UnitMass_in_g')
      #if self.verbose: print "Units Set From %s (%s)"%("self.unitsparameters",self.unitsparameters.filename)
    # now, create the system of units
    self.localsystem_of_units = pNbody.units.Set_SystemUnits_From_Params(params)
    return
  

####################################################################################################################################
# METHODS THAT RETURN INFO ON THE PNBODY OBJECT
####################################################################################################################################

  ##################################################################################################################################
  # METHOD FOR RETURNING THE MEAN NUMBER DENSITY AT A GIVEN RADIUS
  def get_numdens(self,r,z=None,ptyp=None,pgal=None,cylindric=None,rbuff=None,zbuff=None):
    """
    Returns the number density at a specified radius or list of radii
    """
    # Pars input
    if not isinstance(r,list):
      flagscalar=True
      r=[mmlpars.mml_pars(r,type=float,min=0.)]
      z=[mmlpars.mml_pars(z,type=float,default=0.)]
    else:
      flagscalar=False
      z=mmlpars.mml_pars(z,type=list,nelements=len(r),default=[0. for ir in r])
    cylindric=mmlpars.mml_pars(cylindric,default=False,type=bool)
    rbuff=mmlpars.mml_pars(rbuff,default=0.1,min=0.,type=float)
    zbuff=mmlpars.mml_pars(zbuff,default=0.1,min=0.,type=float)
    # Get indeces for selected particles & galaxy
    paridx=self.get_index(ptyp=ptyp,pgal=pgal)
    # Get radial positions of selected particles
    if cylindric:
      rpos=self.get_var('rxy',idxtot=paridx)
      zpos=self.get_var('z',idxtot=paridx)
    else:
      rpos=self.get_var('rxyz',idxtot=paridx)
    # Loop over radii
    dens=[]
    for iridx in range(len(r)):
      ir=r[iridx]
      irmin = max(ir-rbuff,0.)
      irmax = ir+rbuff
      if cylindric:
        iz=z[iridx]
        izmin = iz-zbuff
        izmax = iz+zbuff
        ivol = np.pi*(irmax**2 - irmin**2)*(izmax-izmin)
        inumr = np.logical_and(rpos>=irmin,rpos<=irmax)
        inumz = np.logical_and(zpos>=izmin,zpos<=izmax)
        inum = np.logical_and(inumr,inumz).sum()
      else:
        ivol = (4.0/3.0)*np.pi*(irmax**3 - irmin**3)
        inum = np.logical_and(rpos>=irmin,rpos<=irmax).sum()
      dens.append(inum/ivol)
    # Return correct format
    if flagscalar:
      return dens[0]
    else:
      return dens
      
  ##################################################################################################################################
  # METHOD FOR COUNTING PARTICLES ENCLOSED WITHIN A GIVEN RADIUS
  def get_ninside(self,r,z=None,ptyp=None,pgal=None,cylindric=None):
    # Pars input
    if not isinstance(r,list):
      flagscalar=True
      r=[mmlpars.mml_pars(r,type=float,min=0.)]
      z=[mmlpars.mml_pars(z,type=float,default=r[0],min=0.)]
    else:
      flagscalar=False
      z=mmlpars.mml_pars(z,type=list,nelements=len(r),default=[ir for ir in r])
    cylindric=mmlpars.mml_pars(cylindric,default=False,type=bool)
    # Get indeces for selected particles & galaxy
    paridx=self.get_index(ptyp=ptyp,pgal=pgal)
    # Get radial positions of selected particles
    if cylindric:
      rpos=self.get_var('rxy',idxtot=paridx)
      zpos=self.get_var('z',idxtot=paridx)
    else:
      rpos=self.get_var('rxyz',idxtot=paridx)
    # Loop over radii
    num=[] ; vol=[]
    for iridx in range(len(r)):
      ir=r[iridx]
      if cylindric:
        iz=z[iridx]
        ivol = np.pi*(ir**2)*(2.*iz)
        inumr = np.logical_and(rpos>=0.,rpos<=ir)
        inumz = np.logical_and(zpos>=-iz,zpos<=iz)
        inum = np.logical_and(inumr,inumz).sum()
      else:
        ivol = (4.0/3.0)*np.pi*(ir**3)
        inum = np.logical_and(rpos>=0.,rpos<=ir).sum()
      num.append(inum)
      vol.append(ivol)
    # Return correct format
    if flagscalar:
      return num[0],vol[0]
    else:
      return num,vol
      
  ##################################################################################################################################
  # METHOD FOR RETURNING THE MEAN INTERPARTICLE SEPARATION AT A GIVEN RADIUS
  def get_intersep(self,r,z=None,ptyp=None,pgal=None,cylindric=None,rbuff=None,zbuff=None):
    """
    Returns the mean interparticle separation at specified radius or list of radii
    """
    # Pars input
    if not isinstance(r,list):
      flagscalar=True
      r=[mmlpars.mml_pars(r,type=float,min=0.)]
      z=[mmlpars.mml_pars(z,type=float,default=0.)]
    else:
      flagscalar=False
      z=mmlpars.mml_pars(z,type=list,nelements=len(r),default=[0. for ir in r])
    # Get number density
    numdens=np.array(self.get_numdens(r,z=z,ptyp=ptyp,pgal=pgal,cylindric=cylindric,rbuff=rbuff,zbuff=zbuff))
    # Mask zeros to avoid divide by error
    numdens=np.ma.masked_equal(numdens,0.)
    # Calculate interparticles separation
    intersep=np.array(numdens)**(-1./3.)
    # Return correct format
    if flagscalar:
      return intersep[0]
    else:
      return list(intersep)

  ##################################################################################################################################
  # METHOD FOR RETURNING THE LIST OF VALID PARTICLE TYPES FOR THE SIMULATION
  def get_typlist(self):
    typListFull=simlist.LIST_PTYPBASE
    typlist=[ityp for ityp in typListFull if self.npart[typListFull.index(ityp)]>0]
    return typlist

  ##################################################################################################################################
  # METHOD FOR RETURN THE LIST OF VALID GALAXIES FOR THE SIMULATION
  def get_gallist(self):
    try:
      gallist=[0]+list(np.unique(self.gal))
    except:
      gallist=[0,1,2]
    return gallist

  ##################################################################################################################################
  # METHOD FOR RETURNING THE NUMBER OF GALAXIES IN A SIMULATION
  def get_ngal(self):
    gallist=self.get_gallist()
    ngal=max(gallist)
    return ngal
  
  ##################################################################################################################################
  # METHOD FOR DETERMINING VPEAK IN THE ROTATION CURVE
  def get_Vpeak(self,ptyp=None,pgal=None,idxtot=None):
    if idxtot==None: idxtot=self.get_index(ptyp=ptyp,pgal=pgal)
    vphi=self.get_var('vphi',idxtot=idxtot)
    r=self.get_var('rxy',idxtot=idxtot)
    histn,rbins=np.histogram(r,bins=100)
    histv,rbins=np.histogram(r,bins=rbins,weights=vphi)
    histv=np.ma.masked_where(histn==0,histv) 
    histn=np.ma.masked_where(histn==0,histn) 
    Vpeak=np.max(histv/histn)
    mmlio.verbose('Vpeak={}'.format(Vpeak))
    return Vpeak
    
  ##################################################################################################################################
  def get_stability(self,ptyp=None,pgal=None,idxtot=None,method=None):
    from mmlastro import mmldisk
    if idxtot==None: idxtot=self.get_index(ptyp=ptyp,pgal=pgal)
    method=mmlpars.mml_pars(method,default='op73',list=mmldisk.list_stabmeths())
    if  method=='op73':
      Trot=np.sum(self.get_var('Erot',idxtot=idxtot))
      W=np.sum(self.get_var('Ep',idxtot=idxtot))
      crit=mmldisk.stability_op73(Trot,W)
    elif method=='efstat82':
      G=self.get_gconst()
      Vpeak=self.get_Vpeak(ptyp=ptyp,pgal=pgal,idxtot=idxtot)
      Md=np.sum(self.mass[self.get_index(ptyp='disk',pgal=pgal)])
      Rd=4. # CHANGE THIS!!!
      crit=mmldisk.stability_efstat82(Vpeak,Md,Rd,G)
    elif method=='toomre81':
      G=self.get_gconst()
      R=4.
      sigr=self.get_sigmaR(ptyp=ptyp,pgal=pgal,idxtot=idxtot)
      crit=None
#      kappa
#      crit=mmldisk.stability_toomre81(R,sigr,kappa,rhoS,G)
    else: raise Exception('Invalid stability criterion: {}'.format(method))
    return crit
      
  ##################################################################################################################################
  # METHOD FOR RETURNING SPIN PARAMETER
  def get_spin(self,ptyp=None,pgal=None,idxtot=None,R=None):
    if pgal==None or pgal==0:
      pgal=1
      mmlio.verbose('Set pgal=1 by default')
    if idxtot==None: idxtot=self.get_index(ptyp=ptyp,pgal=pgal)
    # Select particles by radius
    rvec=self.get_var('rxyz')
    if not isinstance(R,float): R=rvec[idxtot].max()
    idxR=(rvec<=R)
    idxRtot=np.logical_and(idxtot,idxR)
    # Calculate quantities
    J=self.get_var('Ltot',idxtot=idxRtot)
    M=np.sum(self.mass[idxRtot]) # Just selected component(s)
#    M=np.sum(self.mass[idxR])    # All components
    G=self.get_gconst()
    V=np.sqrt(G*M/R)
    # Calc lambda spin param (Bullock et al. 2001)
    lspin=mmlmath.spin(J,M,V,R)
    return lspin

  ##################################################################################################################################
  def get_gconst(self):
    G_cgs=6.674e-8 # cm^3/g/s
    uR=self.localsystem_of_units.convertionFactorTo(pNbody.units.Unit_cm)
    uV=self.localsystem_of_units.convertionFactorTo(pNbody.units.Unit_cm/pNbody.units.Unit_s)
    uT=self.localsystem_of_units.convertionFactorTo(pNbody.units.Unit_s )
    uM=self.localsystem_of_units.convertionFactorTo(pNbody.units.Unit_g )
#    print 'uR={}'.format(uR)
#    print 'uM={}'.format(uM)
#    print 'uT={}'.format(uT)
#    print 'uR/uV={}'.format(uR/uV)
    G=G_cgs*(uR**-3)*uM*(uT**2)
#    mmlio.verbose('G={}'.format(G))
    return G
    
  ##################################################################################################################################
  # METHOD FOR RETURNING THE CIRCULAR ANGULAR MOMENTUM
  def get_Lcirc(self,ptyp=None,pgal=None,idxtot=None,mag=False):
    if idxtot==None: idxtot=self.get_index(ptyp=ptyp,pgal=pgal)
    mass=self.get_var('m')[idxtot]
    pos=self.get_var('pos')[idxtot,:]
    vcirc=self.get_vcirc(ptyp=ptyp,pgal=pgal,idxtot=idxtot,mag=False)
    angmom=mmlmath.angmom(mass,pos,vcirc)
    if mag: angmom=mmlmath.absvect(angmom)
    return angmom
    
  ##################################################################################################################################
  # METHOD FOR RETURNING THE CIRCULAR VELOCITY
  def get_vcirc(self,ptyp=None,pgal=None,idxtot=None,mag=False):
    if idxtot==None: idxtot=self.get_index(ptyp=ptyp,pgal=pgal)
    vel=self.get_var('vel')[idxtot]
    vtot=self.get_var('vrxyz')[idxtot]
    try:    pot=self.get_var('pot')[idxtot]
    except: pot=np.zeros(vtot.shape)
    vcirc=np.sqrt(np.abs(pot))
    if not mag:
      vx=vcirc*vel[:,0]/vtot
      vy=vcirc*vel[:,1]/vtot
      vz=vcirc*vel[:,2]/vtot
      vcirc=np.vstack((vx,vy,vz)).T
    return vcirc
  
  ##################################################################################################################################
  # METHOD FOR RETURNING THE RADIAL VELOCITY DISPERSION
  def get_sigmaR(self,ptyp=None,pgal=None,idxtot=None):
    if idxtot==None: idxtot=self.get_index(ptyp=ptyp,pgal=pgal)
    vrxy=self.get_var('vrxy')[idxtot]
    sigmaR=np.sqrt(np.mean((vrxy-np.mean(vrxy))**2))
    return sigmaR

  ##################################################################################################################################
  # METHOD FOR RETURNING THE RADIAL VELOCITY DISPERSION
  def get_sigmar(self):
    return np.sqrt(np.mean((self.vrxyz()-np.mean(self.vrxyz()))**2.))

  ##################################################################################################################################
  # METHOD FOR RETURN THE CLAUSIUS ENERGY
  def get_clause(self):
    try:
      clause=np.sum(self.mass*np.sum(self.acce*self.pos,axis=0))
    except:
      clause=0.0
    return clause

  ##################################################################################################################################
  # METHOD TO DETERMINE MINIMUM ID OF A SECOND GALAXY
  def get_idgal2(self,preserve_num=True):
    if self.cosmosim: idgal2=self.num.max()+1
    else:
      if preserve_num:
        idgal2=long(9*mmlmath.oom(float(self.num.max()),nsig=1,method='CEIL'))
      else:
        idgal2=long(9*mmlmath.oom(float(self.nbody_tot),nsig=1,method='CEIL'))
    return idgal2
  
  def center(self,method=None,ret_cen=False,**kwargs):
    kwargs['method']=method
    kwargs.setdefault('centergal',0)
    kwargs.setdefault('centertyp','all')
    kwargs.setdefault('move',True)
    kwargs.setdefault('vcomflag',True)
    if isinstance(kwargs['centergal'],str): kwargs['centergal']=int(float(kwargs['centergal'][-1]))
    x=self.get_center(**kwargs)
    if ret_cen: return x

  ##################################################################################################################################
  # METHOD TO CENTER PARTICLES
  def get_center(self,method=None,centergal=1,centertyp='all',idxcent=None,move=False,vcomflag=False,navg=None,rmax=None):
    # Pars input
    if method not in simlist.LIST_COMMETHS+[None]:
      print methLIST
      mmlio.yorn('{}'.format(method))
    method=mmlpars.mml_pars(method,list=simlist.LIST_COMMETHS,default=simlist.LIST_COMMETHS[0])
    if method=='com': method='mass'
    # Handle individual galaxies/particle types
    if centergal=='v' or centertyp=='v':
        return self.get_center_by_part(method=method,centergal=centergal,centertyp=centertyp,idxcent=idxcent,
                                       vcomflag=vcomflag,move=move,navg=navg,rmax=rmax)
    # Handle average over N particles
    if method.endswith('_avg'):
      imeth=method.split('_avg')[0]
      navgDEF=100
    else:
      imeth=method
      navg=0
    # Get indices
    if idxcent==None: idxcent=self.get_index(ptyp=centertyp,pgal=centergal,rmax=rmax)
    if not np.any(idxcent):
      com=np.array([0.,0.,0.])
      vcom=np.array([0.,0.,0.])
      if vcomflag: return com,vcom
      else       : return com
    # Get data
    m=self.mass[idxcent].astype(float)
    pos=self.pos[idxcent].astype(np.float64)
    vel=self.vel[idxcent].astype(np.float64)
    if imeth=='pot':
      if hasattr(self,'pot'): pot=self.pot
      else                  : pot=None
      if pot==None: imeth='none'
      else        : pot=pot[idxcent].astype(np.float64)
    # Get center
    if   imeth=='mass'    : out=mmlmath.center_mass(m,pos,vel=vel,vflag=vcomflag)
    elif imeth=='pot'     : out=mmlmath.center_potential(pot,m,pos,vel=vel,vflag=vcomflag,navg=navg)
    elif imeth=='dens'    : out=mmlmath.center_density(m,pos,vel=vel,vflag=vcomflag,navg=navg)
    elif imeth=='none'    : 
      com=np.array([0.,0.,0.])
      vcom=np.array([0.,0.,0.])
      if vcomflag: out=(com,vcom)
      else       : out=com
    else: raise Exception('Invalid centering method: {}'.format(method))
    # Set output
    if vcomflag: com,vcom=out
    else       : com=out
    # Adjust positions by center
    if move:
      if self.verbose: mmlio.verbose('Recentering... Old center of mass: {} ({},{})'.format(com,centertyp,centergal))
      self.pos=(self.pos-com).astype(np.float32)
      if hasattr(self,'scfpos'): self.scfpos=(self.scfpos-com).astype(np.float32)
      if vcomflag: 
        self.vel=(self.vel-vcom).astype(np.float32)
        if hasattr(self,'scfvel'): self.scfvel=(self.scfvel-vcom).astype(np.float32)
    # Return output
    if vcomflag: return com,vcom
    else       : return com

  ##################################################################################################################################
  # METHOD TO CENTER PARTICLES BY TYPE
  def get_center_by_part(self,method=None,centergal=0,centertyp='all',move=False,vcomflag=False,**inkw):
    if centergal=='v': galList=self.get_gallist()
    else             : galList=[centergal]
    if centertyp=='v': typList=self.get_typlist()
    else             : typList=[centertyp]
    # Loop over galaxies
    centerList={}
    for ityp in typList:
      for igal in galList:
        # Get center for ith galaxy and type
        icenter=self.get_center(method=method,centergal=igal,centertyp=ityp,move=False,vcomflag=vcomflag,**inkw)
        centerList['{}{}'.format(ityp,igal)]=icenter
        # Move if necessary
        if move:
          # Get index for ith galaxy and type
          iidx=self.get_index(ptyp=ityp,pgal=igal)
          if not np.any(iidx): continue
          # Recover results
          if vcomflag: com,vcom=icenter 
          else       : com=icenter
          # Move positions
          self.pos[iidx,:]=(self.pos[iidx,:]-com).astype(np.float32)
          if hasattr(self,'scfpos'): self.scfpos[iidx,:]=(self.scfpos[iidx,:]-com).astype(np.float32)
          # Move velocities
          if vcomflag:
            self.vel[iidx,:]=(self.vel[iidx,:]-vcom).astype(np.float32)
            if hasattr(self,'scfvel'): self.scfvel[iidx,:]=(self.scfvel[iidx,:]-vcom).astype(np.float32)
    return centerList

  ##################################################################################################################################
  # METHOD TO CENTER PARTICLES BY TYPE
  def get_center_by_type(self,centergal=0,typList=None,move=False,method=None):
    typListDEF=simlist.LIST_PTYPBASE
    typdict={typListDEF[ityp]:ityp for ityp in range(len(typListDEF))}
    typList=mmlpars.mml_pars(typList,default=typListDEF,type=list)
    # Center on primary's center of mass by type
    centerList=[]
    for ityp in typList:
      if not self.npart[typdict[ityp]]>0: continue
      icenter=self.get_center(method=method,centergal=centergal,centertyp=ityp,move=False)
      centerList.append(icenter)
      if move:
        idxboth=self.get_index(ptyp=ityp,pgal=0)
        for icoord in range(3): self.pos[idxboth,icoord]+=(-icenter[icoord])
        if hasattr(self,'scfpos'):
          for icoord in range(3): self.scfpos[idxboth,icoord]+=(-icenter[icoord])
    return centerList

  ##################################################################################################################################
  # METHOD FOR GETTING INDICES FOR PARTICLES BY GALAXY AND TYPE
  def get_index(self,ptyp=None,pgal=None,rmax=None,rbox=None):
    typelist=simlist.LIST_PTYPBASE
    typedict={typelist[ityp]:ityp for ityp in range(len(typelist))}
    # Pars input
    if isinstance(ptyp,str):
      ptyp=[ptyp]
    else:
      ptyp=mmlpars.mml_pars(ptyp,default=['all'],type=list)
    pgal=mmlpars.mml_pars(pgal,default=0,list=self.get_gallist())
    # Find type indices
    typeindex=np.array(self.nbody_tot*[False])
    for iptyp in ptyp:
#      if not self.npart_tot[typedict[iptyp.lower()]] > 0: continue
      if iptyp.lower() == 'all':
        itypeindex=np.array(self.nbody_tot*[True])
      elif iptyp.lower() in ['visi','star']:
        itypeindex=np.array(self.nbody_tot*[False])
        for iiptyp in simlist.LIST_PTYPVISI:
          itypeindex=np.logical_or(itypeindex,(self.tpe==typedict[iiptyp.lower()]))
      elif iptyp in simlist.LIST_PTYPBASE:
        itypeindex=(self.tpe==typedict[iptyp.lower()])
      else: raise Exception('Invalid particle type: {}'.format(iptyp))
      typeindex=np.logical_or(typeindex,itypeindex)
    # Find galaxy indicies
    if pgal == 0:
      galindex=np.array(self.nbody_tot*[True])
    else:
      galindex=(self.gal==pgal)
    # Find rmax indices
    if not isinstance(rmax,float):
      rmaxindex=np.array(self.nbody_tot*[True])
    else:
      rmax=mmlpars.mml_pars(rmax,type=float,min=0.)
      rxyz=self.get_var('rxyz')
      rmaxindex=(rxyz<=rmax)
    # Find rbox indicies
    if not isinstance(rbox,float):
      rboxindex=np.array(self.nbody_tot*[True])
    else:
      xmaxindex=(np.abs(self.pos[:,0])<=rbox)
      ymaxindex=(np.abs(self.pos[:,1])<=rbox)
      zmaxindex=(np.abs(self.pos[:,2])<=rbox)
      rboxindex=np.logical_and(np.logical_and(xmaxindex,ymaxindex),zmaxindex)
    # Find total selection indices
    rindex=np.logical_and(rmaxindex,rboxindex)
    totindex=np.logical_and(np.logical_and(typeindex,galindex),rindex)
    # Return indices
    return totindex

  ##################################################################################################################################
  # METHOD TO GET UNITS FOR DIFFERENT ARRAYS
  def get_units(self,mode):
#['UnitDensity', 'UnitDic', 'UnitEnergy', 'UnitLength', 'UnitMass', 'UnitSpecEnergy', 'UnitSurfaceDensity', 'UnitTime', 'UnitVelocity', '__doc__', '__init__', '__module__', 'convertionFactorTo', 'dic_of_factors', 'dic_of_powers', 'get_UnitDensity_in_cgs', 'get_UnitEnergy_in_cgs', 'get_UnitLength_in_cm', 'get_UnitMass_in_g', 'get_UnitTime_in_s', 'get_UnitVelocity_in_cm_per_s', 'info', 'into']
    # Pars mode
    unitsys=self.localsystem_of_units
    # print dir(unitsys)
    # mmlio.yorn('?')
    if not mode: unit=1
    elif mode in ['n','none','gal','potm2rel']: 
      unit=1
    elif mode in ['t','time']:
      unit=str(unitsys.get_UnitTime_in_s())+' s'
    elif mode in ['m','mass']: 
      unit=str(unitsys.get_UnitMass_in_g())+' g'
    elif mode in ['pos','x','y','z','r','rxyz','R','rxy']:
      unit=str(unitsys.get_UnitLength_in_cm())+' cm'
    elif mode in ['vel','vx','vy','vz','vr','vxyz','vR','vxy','vr_sph','vr_cyl','vphi','vcir','vtheta','vphi_sph','vphi_cyl','vdyn']:
      unit=str(unitsys.get_UnitVelocity_in_cm_per_s())+' cm s**-1'
    elif mode in ['pot','potm2']:
      unit=None
    elif mode in ['Ek','Ep','Etot','Erot']:
      unit=str(0.5*unitsys.get_UnitMass_in_g()*(unitsys.get_UnitVelocity_in_cm_per_s()**2))+' g cm**2 s**-2'
    elif mode in ['Lxyz','Lx','Ly','Lz','Ltot','Lmag']:
      unit=str(unitsys.get_UnitMass_in_g()*
               unitsys.get_UnitLength_in_cm()*
               unitsys.get_UnitVelocity_in_cm_per_s())+' g cm cm**2 s**-2'
    elif mode.startswith('scf'): unit=self.get_units(mode.split('scf')[-1])
    elif mode.endswith('_enc'): unit=self.get_units(mode.split('_enc')[0])
    elif mode in ['spin','spintot']: unit=1
    elif mode in ['phi','theta']: unit=1
    elif mode in ['wphi','wcir','wphi_sph','wphi_cyl']:
      unit=str(unitsys.get_UnitVelocity_in_cm_per_s()/unitsys.get_UnitLength_in_cm())+' s**-1'
    elif mode == 'R2w':
      unit=str(unitsys.get_UnitVelocity_in_cm_per_s()*unitsys.get_UnitLength_in_cm())+' cm**2 s**-1'
    elif mode in ['oortA','oortB']: unit=1
    elif mode in ['epicycle','epicycle_noaprx']: unit=None
    elif mode == 'surfdens' :
      unit=str(unitsys.get_UnitMass_in_g()/(unitsys.get_UnitLength_in_cm()**2))+' g cm**-2'
    elif mode == 'toomreQ'  : unit=1
    else: unit=None
    return unit

  def has_var(self,name):
    '''
    Return true if the object pNbody has
    a variable called self.name
    '''
    get_list_of_vars = self.get_list_of_vars()
    try:
      getattr(self,name)
      return True
    except AttributeError:
      return False

  def has_array(self,name):
    '''
    Return true if the object pNbody has
    an array called self.name
    '''
    list_of_array = self.get_list_of_array()
    try:
      list_of_array.index(name)
      return True
    except ValueError:
      try:
        self.get_var(name)
      except AttributeError:
        return False

  ##################################################################################################################################
  # METHOD TO RETURN DIFFERENT ARRAYS
  def get_var(self,mode,ptyp=None,pgal=None,idxtot=None,rmax=None,rbox=None):
    from pNbody.libutil import getval
    if idxtot==None: idxtot=self.get_index(ptyp=ptyp,pgal=pgal,rmax=rmax,rbox=rbox)
    nbods=np.sum(self.npart)
    # Pars mode
    if not mode: var=np.ones(len(self.mass[idxtot]),dtype=float)
    elif mode in ['n','none']: var=np.ones(len(self.mass[idxtot]),dtype=float)
    elif mode in ['m','mass']: var=self.mass[idxtot]
    elif mode == 'pos'       : var=self.pos[idxtot,:]
    elif mode == 'vel'       : var=self.vel[idxtot,:]
    elif mode == 'pot':
      try:    var=self.pot[idxtot]
      except: var=np.zeros(np.sum(idxtot))
    elif mode == 'scfpot':
      try:    var=self.scfpot[idxtot]
      except: 
        self.read_scffiles(self.scffiledict)
        var=self.scfpot[idxtot]
    elif mode == 'scfpotm2':
      try:    var=self.scfpotm2[idxtot]
      except: 
        self.read_scffiles(self.scffiledict)
        var=self.scfpotm2[idxtot]
    elif mode == 'scfpotm2rel':
      try:    var=self.scfpotm2[idxtot]/self.scfpot[idxtot]
      except: 
        self.read_scffiles(self.scffiledict)
        var=self.scfpotm2[idxtot]/self.scfpot[idxtot]
    elif mode == 'gal':
      try:    var=self.gal[idxtot]
      except: var=None
    elif mode == 'Ek'  : 
      vtot=self.get_var('vxyz',idxtot=idxtot)
      mass=self.get_var('m',idxtot=idxtot)
      var=0.5*mass*vtot*vtot
    elif mode == 'Ep'  : 
      mass=self.get_var('m',idxtot=idxtot)
      try:    pot=self.get_var('pot',idxtot=idxtot)
      except: pot=np.zeros(mass.shape)
      var=mass*pot
    elif mode == 'Etot': var=self.get_var('Ek',idxtot=idxtot)+self.get_var('Ep',idxtot=idxtot)
    elif mode == 'Erot': var=0.5*self.mass[idxtot]*(self.get_var('vcir',idxtot=idxtot)**2)
    elif mode.startswith('sigma_'):
      varstr=mode.split('_')
      if len(varstr)!=3: raise Exception('Variable must be of the form sigma_x_y, not {}'.format(mode))
      x=self.get_var(varstr[1],idxtot=idxtot)
      y=self.get_var(varstr[2],idxtot=idxtot)
      var=mmlmath.dispersion(x,y)
    elif mode == 'vsig':
      vtot=self.get_var('vxyz',idxtot=idxtot)
      var=np.sqrt(np.mean((vtot-np.mean(vtot))**2))
    elif mode in ['Lxyz','Lx','Ly','Lz','Ltot','Lmag']:
      m=self.get_var('m',idxtot=idxtot)
      r=self.get_var('pos',idxtot=idxtot)
      v=self.get_var('vel',idxtot=idxtot)
      Lxyz=mmlmath.angmom(m,r,v)
      if   mode == 'Lxyz': var=Lxyz
      elif mode == 'Lmag': var=mmlmath.absvect(Lxyz)
      elif mode == 'Ltot': var=np.sqrt(np.sum(np.sum(Lxyz,axis=0)**2))
      elif mode == 'Lx'  : var=Lxyz[:,0]
      elif mode == 'Ly'  : var=Lxyz[:,1]
      elif mode == 'Lz'  : var=Lxyz[:,2]
    elif mode.endswith('_enc'):
      mdbase=mode.split('_enc')[0]
      R=self.get_var('rxyz',idxtot=idxtot)
      X=self.get_var(mdbase,idxtot=idxtot)
      var=mmlmath.xenc(R,X)
    elif mode == 'spin': 
      J=self.get_var('Lxyz',idxtot=idxtot)
      M=self.get_var('mass',idxtot=idxtot)
      R=self.get_var('rxyz',idxtot=idxtot)
      G=self.get_gconst()
      var=mmlmath.spinvec(J,M,R,G)
    elif mode == 'spintot':
      J=self.get_var('Ltot',idxtot=idxtot)
      M=self.get_var('mass',idxtot=idxtot).sum()
      R=self.get_var('rxyz',idxtot=idxtot).max()
      G=self.get_gconst()
      V=np.sqrt(G*M/R)
      var=mmlmath.spin(J,M,V,R)
    # 1D position
    elif mode == 'x'        : var=self.pos[idxtot,0]
    elif mode == 'y'        : var=self.pos[idxtot,1]
    elif mode == 'z'        : var=self.pos[idxtot,2]
    elif mode == 'rxy'      : var=mmlmath.xyz2rxy(self.pos[idxtot,:])
    elif mode == 'rxyz'     : var=mmlmath.xyz2rxyz(self.pos[idxtot,:])
    elif mode == 'phi'      : var=mmlmath.xyz2phi(self.pos[idxtot,:])
    elif mode == 'theta'    : var=mmlmath.xyz2theta(self.pos[idxtot,:])
    elif mode == 'scfx'     : var=self.scfpos[idxtot,0]
    elif mode == 'scfy'     : var=self.scfpos[idxtot,1]
    elif mode == 'scfz'     : var=self.scfpos[idxtot,2]
    elif mode == 'scfrxy'   : var=mmlmath.xyz2rxy(self.scfpos[idxtot,:])
    elif mode == 'scfrxyz'  : var=mmlmath.xyz2rxyz(self.scfpos[idxtot,:])
    elif mode == 'scfphi'   : var=mmlmath.xyz2phi(self.scfpos[idxtot,:])
    elif mode == 'scftheta' : var=mmlmath.xyz2theta(self.scfpos[idxtot,:])
    # 1D velocity
    elif mode == 'vx'       : var=self.vel[idxtot,0]
    elif mode == 'vy'       : var=self.vel[idxtot,1]
    elif mode == 'vz'       : var=self.vel[idxtot,2]
    elif mode == 'vxy'      : var=mmlmath.xyz2rxy(self.vel[idxtot,:])
    elif mode == 'vxyz'     : var=mmlmath.xyz2rxyz(self.vel[idxtot,:]) 
#    elif mode == 'vrxy'     : var=mmlmath.xyz2rxy(self.vel[idxtot,:])
#    elif mode == 'vrxyz'    : var=mmlmath.xyz2rxyz(self.vel[idxtot,:])
    elif mode == 'vr'       : var=self.get_var('vr_sph',idxtot=idxtot)
    elif mode == 'vR'       : var=self.get_var('vr_cyl',idxtot=idxtot)
    elif mode == 'vr_sph'   : var=mmlmath.vrad(self.pos[idxtot,:],self.vel[idxtot,:],coord='sph')
    elif mode == 'vr_cyl'   : var=mmlmath.vrad(self.pos[idxtot,:],self.vel[idxtot,:],coord='cyl')
    elif mode == 'vphi'     : var=self.get_var('vphi_sph',idxtot=idxtot)
    elif mode == 'vcir'     : var=self.get_var('vphi_cyl',idxtot=idxtot)
    elif mode == 'vphi_sph' : var=mmlmath.vphi(self.pos[idxtot,:],self.vel[idxtot,:],coord='sph')
    elif mode == 'vphi_cyl' : var=mmlmath.vphi(self.pos[idxtot,:],self.vel[idxtot,:],coord='cyl')
    elif mode == 'vtheta'   : var=mmlmath.vtheta(self.pos[idxtot,:],self.vel[idxtot,:])
    elif mode == 'wphi'     : var=self.get_var('wphi_sph',idxtot=idxtot)
    elif mode == 'wcir'     : var=self.get_var('wphi_cyl',idxtot=idxtot)
    elif mode == 'wphi_sph' : var=self.get_var('vphi_sph',idxtot=idxtot)/np.ma.masked_equal(self.get_var('rxyz',idxtot=idxtot),0)
#      v=self.get_var('vphi_sph',idxtot=idxtot)
#      r=self.get_var('rxyz',idxtot=idxtot)*self.localsystem_of_units.convertionFactorTo(pNbody.units.Unit_km)
#      var=v/np.ma.masked_equal(r,0)
    elif mode == 'wphi_cyl' : var=self.get_var('vphi_cyl',idxtot=idxtot)/np.ma.masked_equal(self.get_var('rxy' ,idxtot=idxtot),0)
#      v=self.get_var('vphi_cyl',idxtot=idxtot)
#      r=self.get_var('rxy',idxtot=idxtot)*self.localsystem_of_units.convertionFactorTo(pNbody.units.Unit_km)
#      var=v/np.ma.masked_equal(r,0)
    elif mode == 'R2w'      : var=(self.get_var('rxy',idxtot=idxtot)**2)*self.get_var('wphi_cyl',idxtot=idxtot)
    elif mode == 'scfvx'    : var=self.scfvel[idxtot,0]
    elif mode == 'scfvy'    : var=self.scfvel[idxtot,1]
    elif mode == 'scfvz'    : var=self.scfvel[idxtot,2]
    elif mode == 'scfvxy'   : var=mmlmath.xyz2rxy(self.scfvel[idxtot,:])
    elif mode == 'scfvxyz'  : var=mmlmath.xyz2rxyz(self.scfvel[idxtot,:])
    elif mode == 'vdyn'     :
      renc=self.get_var('rxyz' ,idxtot=idxtot)
      menc=self.get_var('m_enc',idxtot=idxtot)
      G=self.get_gconst()
      var=np.sqrt(G*menc/renc) 
    elif mode == 'oortA'    :
      R=self.get_var('rxy'     ,idxtot=idxtot)
      V=self.get_var('vcir'    ,idxtot=idxtot)
      var=mmlmath.oorts_const(R,V,method='A')
    elif mode == 'oortB'    :
      R=self.get_var('rxy'     ,idxtot=idxtot)
      V=self.get_var('vcir'    ,idxtot=idxtot)
      var=mmlmath.oorts_const(R,V,method='B')
    elif mode == 'epicycle' :
      R=self.get_var('rxy'     ,idxtot=idxtot)
      V=self.get_var('vcir'    ,idxtot=idxtot)
      var=mmlmath.epicycle(R,V)
    elif mode == 'epicycle_noaprx':
      R=self.get_var('rxy'     ,idxtot=idxtot)
      V=self.get_var('vcir'    ,idxtot=idxtot)
      var=mmlmath.epicycle(R,V,approx=False)
    elif mode == 'surfdens' :
      M=self.get_var('mass',idxtot=idxtot)
      R=self.get_var('rxy' ,idxtot=idxtot)
      var=mmlmath.surfdens(M,R)
    elif mode == 'toomreQ'  :
      vdisp=self.get_var('sigma_rxy_vR',idxtot=idxtot)
      sdens=self.get_var('surfdens',idxtot=idxtot)
      kappa=self.get_var('epicycle',idxtot=idxtot)
      G=self.get_gconst()
      var=mmlmath.toomreq(vdisp,kappa,sdens,G)
      if var.min()==var.max(): mmlio.yorn('toomreq null: G={}'.format(G))
    else:
      try:
        var=getval(self,mode)
        vshp=var.shape
        if   len(vshp)==1: var=var[idxtot]
        elif len(vshp)==2: var=var[idxtot,:]
        else: mmlio.verbose('Variable has more than 2 dimensions (shape={}). Unsure which should be sampled...'.format(vshp))
      except:
        var=None
    if var==None: raise AttributeError('No such sim attribute: {}'.format(mode))
    return var
    
  ##################################################################################################################################
  # METHOD TO RETURN LIMITS OF VARIABLES
  def get_limits(self,modeList=None,limpad=None,mirror=None,scale=None,
                 ptyp=None,pgal=None,idxtot=None,lim=None):
    # Set constants
    modeListDEF=['n','m','rho','theta','rxy','rxyz','vxy','vxyz','nan','Rz','scfpot','scfpotm2',
                 'phi','x','y','z','pos','vx','vy','vz','vr','vt','vcir','vrad','pos','vel']
    # Pars input
    lim=mmlpars.mml_pars(lim,default={},type=dict)
    limpad=mmlpars.mml_pars(limpad,default=0.0,type=float)
    if isinstance(modeList,str): modeList=[modeList]
    modeList=mmlpars.mml_pars(modeList,default=modeListDEF,type=list)
    if idxtot == None: idxtot=self.get_index(ptyp,pgal)
    # Create limits for each variable
    for ivar in modeList:
      if not lim.has_key(ivar) or (lim.has_key(ivar) and lim[ivar][0]==lim[ivar][1]):
        imethodDEF,iscaleDEF,imirrorDEF=simplot.get_histdef(ivar)
        imirror=mmlpars.mml_pars(mirror,default=imirrorDEF,type=bool)
        iscale=mmlpars.mml_pars(scale,default=iscaleDEF,type=str)
#        print '[get_limits] ivar={}'.format(ivar)
#        print '[get_limits] imirror={} ({})'.format(imirror,imirrorDEF)
#        print '[get_limits] iscale={} ({})'.format(iscale,iscaleDEF)
        if   ivar.endswith('_enc'): imin,imax=self.get_limits(ivar.split('_enc')[0],idxtot=idxtot,lim=lim)
        elif ivar == 'n'          : imin,imax=(1.,float(len(self.mass[idxtot])))
        elif ivar == 'm'          : imin,imax=(self.mass[idxtot].min(),self.mass[idxtot].sum())
        elif 'spin' in ivar       : imin,imax=(0.,1.)
        elif ivar == 'phi'        : imin,imax=(0.,math.pi)
        elif ivar == 'theta'      : imin,imax=(0.,math.pi)
        elif ivar == 'rho':
          rxyz=self.get_var('rxyz',idxtot=idxtot)
          mass=self.get_var('mass',idxtot=idxtot)
          rlim=(rxyz.min(),rxyz.max())
          mmin,mmax=(mass.max(),mass.sum())
          vmin,vmax=tuple((4./3.)*math.pi*(np.array(rlim)**3))
          imin,imax=(mmax/vmax,10.*mmin/vmin)
        elif ivar in ['pos','Rz','x','y','z','scfpos']:
          varArr=self.get_var('rxyz',idxtot=idxtot)
          imin,imax=(varArr.min(),varArr.max())
        elif ivar in ['vel','scfvel','vr','vt']:
          varArr=self.get_var('vrxyz',idxtot=idxtot)
          imin,imax=(varArr.min(),varArr.max())
        elif ivar in ['scfpot','scfpotm2']:
          if hasattr(self,'scfpot'):
            imin=abs(self.scfpot[idxtot]).min()
            imax=abs(self.scfpot[idxtot]).max()
            imean=abs(self.scfpot[idxtot]).mean()
            imax/=10.
            imax=0.1
#            mmlio.yorn('[get_limits] {}: {}'.format(ivar,(imin,imean,imax)))
          else:
            lim[ivar]=(0.,0.)
            continue
        else:
          varArr=self.get_var(ivar,idxtot=idxtot)
          if varArr == None:
            lim[ivar]=(0.,0.)
            continue
          else:
            imin=abs(varArr).min()
            imax=abs(varArr).max()
        # Ensure that masked constant is not returned
        if isinstance(imin,np.ma.masked_array) or isinstance(imax,np.ma.masked_array):
          lim[ivar]=(0.,0.)
          continue
#        if imin==imax:
#          lim[ivar]=(imin,imax)
#          continue
        # Add dependence on scale
        if   iscale == 'linear':
          if imirror:
            minmir=-imax
            maxmir=+imax
          else:
            minmir=imin
            maxmir=imax
          abspad=limpad*(maxmir-minmir)/2.
          minpad=-abspad
          maxpad=+abspad
        elif iscale == 'log':
          if imirror:
            minmir=-imax
            maxmir=+imax
            abspad=limpad*np.log10(maxmir)
            maxpad=(10**(+abspad))*maxmir - maxmir
            minpad=(10**(+abspad))*minmir - minmir
#            if ivar in ['scfpot','scfpotm2']:
#              maxpad=0.
#              minpad=0.
#            ilim=(minmir+minpad,maxmir+maxpad)
#            mmlio.verbose('mode={} abspad={} minpad={} maxpad={} lim={}'.format(ivar,abspad,minpad,maxpad,ilim))
#            minmir=10**(-np.log10(imax))
#            maxmir=10**(+np.log10(imax))
          else:
            minmir=imin
            maxmir=imax
            abspad=limpad*(np.log10(maxmir)-np.log10(minmir))/2.
            minpad=10**(-abspad)*minmir - minmir
            maxpad=10**(+abspad)*maxmir - maxmir
        else: raise Exception('Unidentified scale {}'.format(iscale))
#        if imirror:
#          lim[ivar]=(-(1.+limpad)*imax,
#                     +(1.+limpad)*imax)
#        else:
#          lim[ivar]=((1.-limpad)*imin,
#                     (1.+limpad)*imax)
        lim[ivar]=(minmir+minpad,maxmir+maxpad)
#        if ivar in ['scfpot','scfpotm2']:
#            print '[get_limits] {}: {}'.format(ivar,lim[ivar])
    # Output limits
    if len(modeList) == 1:
      outlim=lim[modeList[0]]
#      print '[get_limits] {} {} {}: {}'.format(modeList[0],iscale,imirror,outlim)
    else:
      outlim=lim
    return outlim

  ##################################################################################################################################
  # METHOD FOR RETURN DICTIONARIES OF NBODY STATISTICS
  def get_statdict(self,method,smeth=None,pgal=None,ptyp=None,rmax=None,
                   cmeth=None,cgal=None,ctyp=None,cmeth0=None,cgal0=None,ctyp0=None,**addkeys):
    # Pars input
    method=mmlpars.mml_pars(method,type=str)
    method=method.lower()
    # Determine what statistics to return
    if   method=='general': statdict=self.get_genstat(pgal=pgal,centergal=cgal,centertyp=ctyp,centermeth=cmeth)
    elif method=='var'    : statdict=self.get_varstat(smeth,centermeth=cmeth,ctyp=ctyp,cgal=cgal,rmax=rmax)
    elif method=='profile': statdict=self.get_profstat(pgal=pgal,centergal=cgal,centertyp=ctyp,centermeth=cmeth)
    elif method=='com'    : statdict=self.get_comstat(centermeth=cmeth,centermeth0=cmeth0,centergal0=cgal0,centertyp0=ctyp0,rmax=rmax)
    elif method=='bar'    : statdict=self.get_barstat(method=smeth,pgal=pgal,ptyp=ptyp)
    elif method=='bulk'   : statdict=self.get_bulkstat(ptyp=ptyp,pgal=pgal,centermeth=cmeth,rmax=rmax)
    else: raise Exception('Invalid method: {}'.format(method))
    # Add standard entries
    statdict['fname']=self.p_name[0]
    statdict['time']=self.atime
    # Return statistics
    return statdict
  
  ##################################################################################################################################
  # METHOD TO RETURN DICTIONARY OF GENERAL SIM PROPERTIES
  def get_genstat(self,pgal=1,centergal=1,centertyp='all',centermeth=None):
    # Center galaxy
    oldcm=self.get_center(method=centermeth,centergal=centergal,centertyp=centertyp,move=True)
    # Calculate quatitites
    cm=self.cm() ; cv=self.cv()
    barModel=self.get_barstat(pgal=pgal,method='dens')
    eps=max(self.get_softenings())
    # Form dictionary
    slist={
      'comx' :cm[0],'comy' :cm[1],'comz' :cm[2],
      'vcomx':cv[0],'vcomy':cv[1],'vcomz':cv[2],
      'clause':self.get_clause(),'sigma':self.get_sigmar(),
      }
    # Add dictionary to pNbody object for later use
    self.genstat=slist
      
    return slist

  ##################################################################################################################################
  # METHOD TO RETURN A DICTIONARY OF NBODY SIMULATION PROFILE INFO
  def get_profstat(self,pgal=1,centergal=1,centertyp='all',centermeth=None):
    # Center galaxy
    oldcm=self.get_center(method=centermeth,centergal=centergal,centertyp=centertyp,move=True)
    # Calculate quatitites
    haloModel=self.get_haloprop('nfw',pgal=pgal)
    diskModel=self.get_diskprop('sechdisk',pgal=pgal)
    bulgeModel=self.get_bulgeprop('hernquist',pgal=pgal)
    # Form dictionary
    slist={
      'halo_A':haloModel['rho0'],'halo_r':haloModel['rc'],
      'disk_A':diskModel['A'],'disk_r':diskModel['l'],'disk_z':diskModel['h'],
      'bulge_A':bulgeModel['A'],'bulge_r':bulgeModel['r0'],
      }
    # Add dictionary to pNbody object for later use
    self.profstat=slist
    
    return slist

  ##################################################################################################################################
  # METHOD TO RETURN DICTIONARY OF BULK SIMULATION STATISTICS
  def get_bulkstat(self,ptyp=None,pgal=None,rmax=None,centermeth=None):
    """
    Returns a dictionary of bulk simulation statistics
    """
    # Set constants
    inb=copy.deepcopy(self)
    # Pars input
    ptyp=mmlpars.mml_pars(ptyp,default='all',list=['all']+simlist.LIST_PTYPBASE)
    pgal=mmlpars.mml_pars(pgal,default=1,list=self.get_gallist())
    # Recenter
    oldcm=inb.get_center(method=centermeth,centergal=pgal,centertyp=ptyp,move=True)
    rmaxDEF=float(np.abs(inb.pos[inb.get_index(ptyp='disk',pgal=pgal),:]).max())
    rmax=mmlpars.mml_pars(rmax,default=rmaxDEF,type=float,min=0.)
    bulkstat={}
    # Get index
    idxtot=inb.get_index(ptyp=ptyp,pgal=pgal,rmax=rmax)
    # Calculate energies
    bulkstat['Ep']=inb.get_potenergy(idxtot=idxtot).sum()
    bulkstat['Ek']=inb.get_kinenergy(idxtot=idxtot).sum()
    # Velocity dispersion
    bulkstat['vsig']=inb.get_veldisp(idxtot=idxtot)
    # Angular momentum
    Ltot=inb.get_var('Lxyz',idxtot=idxtot)
    bulkstat['Lx']=Ltot[:,0].sum()
    bulkstat['Ly']=Ltot[:,1].sum()
    bulkstat['Lz']=Ltot[:,2].sum()
    # Circular angular momentum
    Ltot_circ=inb.get_Lcirc(idxtot=idxtot)
    bulkstat['Lx_circ']=Ltot_circ[:,0].sum()
    bulkstat['Ly_circ']=Ltot_circ[:,1].sum()
    bulkstat['Lz_circ']=Ltot_circ[:,2].sum()
    return bulkstat

  ##################################################################################################################################
  # METHOD TO RETURN DICTIONARY OF CENTER OF MASS PROPERTIES
  def get_comstat(self,centermeth=None,centermeth0=None,centergal0=None,centertyp0=None,rmax=None):
    # Set constants
    inb=copy.deepcopy(self)
    typlistDEF=simlist.LIST_PTYPS
    gallistDEF=simlist.LIST_PGALS
    typlistINC=['all','visi']+inb.get_typlist()
    gallistINC=inb.get_gallist()
    posstr=['x','y','z']
    if len(gallistINC)>2:
      galflag=False
    else:
      galflag=True
      gallistINC.remove(1)
    # Center galaxy
    oldcm=inb.get_center(method=centermeth0,centergal=centergal0,centertyp=centertyp0,rmax=rmax,move=True)
#    newcm=inb.get_center(centergal=1,centertyp='all',move=False,method=centermeth0,rmax=rmax)
    # Loop over particle types
    comstat={}
    for ipgal in gallistDEF:
      for iptyp in typlistDEF:
        # Get com
        if iptyp in typlistINC and ipgal in gallistINC:
          icom=inb.get_center(method=centermeth,centergal=ipgal,centertyp=iptyp,move=False,rmax=rmax)
        # Skip if particle not in snapshot 
        else: icom=[0.,0.,0.]
        # Add to dictionary
        for posidx in range(len(posstr)):
          comstat['{}{}_{}'.format(iptyp,ipgal,posstr[posidx])]=icom[posidx]
          if galflag: comstat['{}{}_{}'.format(iptyp,1,posstr[posidx])]=icom[posidx]
    # Return dictionary
    return comstat

  ##################################################################################################################################
  # METHOD TO RETURN DICTIONARY OF VARIABLE PROPERTIES
  def get_varstat(self,method,centermeth=None,ctyp='all',cgal=1,rmax=None):
    import simstat
    # Set constants
    inb=copy.deepcopy(self)
    typlistDEF=simlist.LIST_PTYPS
    gallistDEF=simlist.LIST_PGALS
    typlistINC=['all','visi']+inb.get_typlist()
    gallistINC=inb.get_gallist()
    if len(gallistINC)>2:
      galflag=False
    else:
      galflag=True
      gallistINC.remove(1)
    # Pars input
    varlist=simstat.list_varstat(method)
    # Center galaxy
    oldcm=inb.get_center(method=centermeth,centergal=cgal,centertyp=ctyp,move=True,rmax=rmax)
    # Loop over particle types
    varstat={}
    for ipgal in gallistDEF:
      for iptyp in typlistDEF:
        # Get index
        if iptyp in typlistINC and ipgal in gallistINC:
          idx=inb.get_index(ptyp=iptyp,pgal=ipgal,rmax=rmax)
        else:
          idx=np.array([False]*len(inb.mass))
        # Loop over variables
        for imeth in varlist:
            # Get variable
            if np.any(idx): ivararr=inb.get_var(imeth,idxtot=idx)
            else          : ivararr=0.0
            # Sum variable if necessary
            if isinstance(ivararr,np.ndarray): ivar=np.sum(ivararr)
            else                             : ivar=ivararr
            # Add to dictionary
            varstat['{}_{}{}'.format(imeth,iptyp,ipgal)]=ivar
            if galflag: varstat['{}_{}{}'.format(imeth,iptyp,1)]=ivar
    # Return dictionary
    return varstat

  ##################################################################################################################################
  # METHOD TO RETURN DICTIONARY OF PROFILES
  def oldget_profile(self,method,Rlist,space=None,galaxy=None,centermeth=None):
    """
    Returns a dictionary of profiles
    """
    import main as mmlnbody
    # Set constants
    inb=copy.deepcopy(self)
    typlistDEF=simlist.LIST_PTYPS
    typlistINC=['all','visi']+inb.get_typlist()
    # Pars input
    method=mmlpars.mml_pars(method,default='M',list=simprof.list_fpar('profile')['mode'])
    Rlist=np.array(mmlpars.mml_pars(Rlist,type=list))
    space=mmlpars.mml_pars(space,default='r',list=simprof.list_fpar('profile')['space'])
    galaxy=mmlpars.mml_pars(galaxy,default=1,list=inb.get_gallist())
    Rlist0=np.array([0.]+list(Rlist))
    # Determine recorded variable
    profdict={}
    if   method=='M'    : profdict[space]=list(Rlist0[:-1]+np.diff(Rlist0)/2.)
    elif method=='Menc' : profdict[space]=list(Rlist)
    elif method=='spin' : profdict[space]=list(Rlist)
    elif method=='omega': profdict[space]=list(Rlist0[:-1]+np.diff(Rlist0)/2.)
    elif method=='Vcirc': profdict[space]=list(Rlist)
    else: raise Exception('Invalid profile method: {}'.format(method))
    # Get enclosed mass for dependent quantities
    # Loop over particle types
    for iptyp in typlistDEF:
      if iptyp in typlistINC:
        oldcm=inb.get_center(method=centermeth,centergal=galaxy,centertyp=iptyp,move=True)
        idxtot=self.get_index(ptyp=iptyp,pgal=galaxy)
        if   space=='r': rvec=inb.get_var('rxyz',idxtot=idxtot)
        elif space=='R': rvec=inb.get_var('rxy' ,idxtot=idxtot)
        else: raise Exception('Invalid profile space: {}'.format(space))
        if   method=='M'    : 
          mass=inb.get_var('m',idxtot=idxtot)
          iprof,ibins=np.histogram(rvec,bins=Rlist0,weights=mass)
        elif method=='Menc' : 
          mass=inb.get_var('m',idxtot=idxtot)
          mhist,ibins=np.histogram(rvec,bins=Rlist0,weights=mass)
          iprof=np.cumsum(mhist)
        elif method=='spin' :
          mass=inb.get_var('m',idxtot=idxtot)
          Mtot,ibins=np.histogram(rvec,bins=Rlist0,weights=mass) # Selected component
          Menc=np.cumsum(Mtot)
          Jxyz=inb.get_var('Lxyz',idxtot=idxtot)
          Jtot,ibins=mmlmath.myhist(rvec,bins=Rlist0,weights=Jxyz)
          Jenc_x=np.cumsum(Jtot[:,0])
          Jenc_y=np.cumsum(Jtot[:,1])
          Jenc_z=np.cumsum(Jtot[:,2])
          Jenc=np.sqrt(Jenc_x**2+Jenc_y**2+Jenc_z**2)
          G=inb.get_gconst()
          Venc=np.sqrt(G*Menc/profdict[space])
          iprof=Jenc/(np.sqrt(2.)*Menc*Venc*profdict[space])
        elif method=='omega':
          W=inb.get_var('vomega')
          Wtot,ibins=np.histogram(rvec,bins=Rlist0,weights=W)
          ntot,ibins=np.histogram(rvec,bins=Rlist0)
          iprof=Wtot/ntot
        elif method=='Vcirc':
          mass=inb.get_var('m',idxtot=idxtot)
          Mtot,ibins=np.histogram(rvec,bins=Rlist0,weights=mass) # Selected component
          Menc=np.cumsum(Mtot)
          G=self.get_gconst()
          iprof=np.sqrt(G*Menc/profdict[space])
        else: raise Exception('Invalid profile method: {}'.format(method))
        iprof=list(iprof)
      else: iprof=[0.]*len(profdict[space])
      # Add profile to dictionary
      profdict[iptyp]=iprof
    # Add keylist
    profdict['keylist']=[space]+typlistDEF
    # Return dictionary
    return profdict
    
  ##################################################################################################################################
  # METHOD TO RETURN BAR PROPERTIES USING DIFFERENT METHODS
  def get_barsumm(self,**hist_kw):
    barstat={}
    barstat_dens =self.get_barprop(method='dens' ,**hist_kw)
    barstat_scfm2=self.get_barprop(method='scfm2',**hist_kw)
    for ikey in barstat_dens.keys() : barstat['dens_' +ikey]=barstat_dens[ikey]
    for ikey in barstat_scfm2.keys(): barstat['scfm2_'+ikey]=barstat_scfm2[ikey]
    return barstat

  ##################################################################################################################################
  # METHOD TO RETURN BAR STATISTICS
  def get_barstat(self,diskangle=None,method=None,includehist=None,**hist_kw):
    # Set constants
    pgalDEF=1
    ptypDEF='disk'
    nbinsDEF=100
    scaleDEF='linear'
    # Pars input
    diskangle=mmlpars.mml_pars(diskangle,default=(0.,0.),type=tuple,nelements=2)
    method=mmlpars.mml_pars(method,default='dens',type=str)
    includehist=mmlpars.mml_pars(includehist,default=False,type=bool)
    if 'pgal'   not in hist_kw: hist_kw['pgal'  ]=pgalDEF
    if 'ptyp'   not in hist_kw: hist_kw['ptyp'  ]=ptypDEF
    if 'xnbins' not in hist_kw: hist_kw['xnbins']=nbinsDEF
    if 'ynbins' not in hist_kw: hist_kw['ynbins']=nbinsDEF
    if 'xscale' not in hist_kw: hist_kw['xscale']=scaleDEF
    if 'yscale' not in hist_kw: hist_kw['yscale']=scaleDEF
    # Select particles
    idx_part=self.get_index(ptyp=hist_kw['ptyp'],pgal=hist_kw['pgal'])
    inb_part=self
    # Set histogram mode
    if   method=='ellip': histmode='m'
    elif method=='dens' : histmode='m'
    elif method=='scfm2': histmode='scfpotm2rel'
    else: raise Exception('Invalid method: {}'.format(method))
    # Create hist_kw dicts
    histkw_1D=copy.deepcopy(hist_kw)
    histkw_1D['idxtot']=idx_part
    histkw_1D['ymode']=histmode
    if 'xmode'  in histkw_1D: del histkw_1D['xmode']
    if 'xlim'   in histkw_1D: del histkw_1D['xlim']
    histkw_2D=copy.deepcopy(hist_kw)
    histkw_2D['idxtot']=idx_part
    histkw_2D['zmode']=histmode
    if 'xymode' in histkw_2D: del histkw_2D['xymode']
    if 'xlim'   in histkw_2D: del histkw_2D['xlim']
    if 'ylim'   in histkw_2D: del histkw_2D['ylim']
    if 'view'   in histkw_2D: del histkw_2D['view']
    # Correct for disk orientation
    if diskangle[0]!=0. or diskangle[1]!=0.:
      inb_part.rotate(-diskangle[0],axis=[0, 0,1])
      inb_part.rotate(-diskangle[1],axis=[0,-1,0])
    # Proceed based on method
    if method=='ellip':
      barstat=inb_part.get_barstat_ellip(includehist=includehist,**histkw_2D)
    else:
      barstat=inb_part.get_barstat_ampli(includehist=includehist,**histkw_1D)
    # Correct for bar orientation
    inb_part.rotate(-barstat['pa'],axis=[0,0,1])
    # Get bar index
    idx_barx=np.logical_and(idx_part,(abs(self.pos[:,0]) < barstat['x']))
    idx_bary=np.logical_and(idx_part,(abs(self.pos[:,1]) < barstat['y']))
    idx_barz=np.logical_and(idx_part,(abs(self.pos[:,2]) < barstat['z']))
    idx_barxy=np.logical_and(idx_barx,idx_bary)
    idx_bar=np.logical_and(idx_barxy,idx_barz)
    # Handle absence of a bar
    if not any(idx_bar):
      barstat['omega']=0.
      if includehist:
        barstat['vbins']=np.array([0.])
        barstat['vhist']=np.array([0.])
        barstat['whist']=np.array([0.])
      barstat['A']=0.
      barstat['Aalt']=0.
#      raise Exception('There are not any particles in the bar!')
    else:
      # Get bar pattern speed
      rmax=barstat['r']
      dimlim=(-rmax,rmax)
      hist_v,bins_v=self.hist1D(idxtot=idx_bar,ymode='vy',xmode='x',xlim=dimlim,**hist_kw)
      bins_vavg=bins_v[:-1]+(bins_v[1:]-bins_v[:-1])/2.
      hist_w=hist_v/(1000.*bins_vavg) # ans in radians/Myr 
      barstat['omega']=float(hist_w.max())
      if includehist:
        barstat['vbins']=bins_v
        barstat['vhist']=hist_v.data ; barstat['vmask']=hist_v.mask
        barstat['whist']=hist_w.data ; barstat['wmask']=hist_w.mask
      # Get bar amplitude
      if method=='dens':
        if self.verbose: print '[get_barstat] from dens max/min ratio: A={}'.format(barstat['A'])
      else:
        m0=inb_part.scfpot[idx_bar].sum()
        m2=inb_part.scfpotm2[idx_bar].sum()
        if m0==0.0:
          barstat['A']=0.
        else:
          barstat['A']=np.mean(m2/m0)
        if self.verbose: print '[get_barstat] from SCF m2/m0 ratio: A={}'.format(barstat['A'])
      barstat['Aalt']=0.0
    # Return disk to original orientation
    inb_part.rotate(barstat['pa'],axis=[0,0,1])
    if diskangle[0]!=0. or diskangle[1]!=0.:
      inb_part.rotate(diskangle[1],axis=[0,-1,0])
      inb_part.rotate(diskangle[0],axis=[0,0,1])
    # Add bar stat to pNbody object
    self.barstat=barstat
    # Return dictionary
    return barstat
    
  ##################################################################################################################################
  # METHOD TO RETURN BAR STATISTICS BASED ON ELLIPSE FITTING
  def get_barstat_ellip(self,includehist=False,idxtot=None,**extra_kw):
    # Get stat dictionary
    if idxtot==None:
      barstat = mmlbar.get_barstat(self.pos,includehist=includehist)
    else:
      barstat = mmlbar.get_barstat(self.pos[idxtot,:],includehist=includehist)
    # Return dictionary
    return barstat
    

  ##################################################################################################################################
  # METHOD TO RETURN BAR STATISTICS BASED ON AMPLITUDE
  def get_barstat_ampli(self,includehist=False,**hist_kw):
    # Set constants
    edgefact_large=4.0
    edgefact_small=2.0
    nbins=hist_kw['xnbins']
    # Initialize structure
    barstat={}
    # Estimate bar size
    hist_rad,bins_rad=self.hist1D(xmode='rxyz',**hist_kw)
    histidx_barr=np.where(abs(hist_rad) > abs(hist_rad).max()/edgefact_large)
    barstat['r']=float(bins_rad[histidx_barr].max())
    # Get bar angle
    hist_phi,bins_phi=self.hist1D(xmode='phi' ,**hist_kw)
    phirange=bins_phi.max()-bins_phi.min()
    idxmax=hist_phi.argmax()
    barstat['pa']=bins_phi[idxmax]
    if phirange <= math.pi: 
      if idxmax >= hist_kw['xnbins']/2: idxmin=idxmax-(hist_kw['xnbins']/2)
      else:                             idxmin=idxmax+(hist_kw['xnbins']/2)
    else:
      if idxmax >= hist_kw['xnbins']/4: idxmin=idxmax-(hist_kw['xnbins']/4)
      else:                             idxmin=idxmax+(hist_kw['xnbins']/4)
    hist_phi_max=hist_phi[idxmax]
    hist_phi_min=hist_phi[idxmin]
    if hist_phi_max is np.ma.masked: raise Exception('Max of histogram in phi is masked.')
    if hist_phi_min is np.ma.masked: hist_phi_min=0
    # Correct for bar orientation
    self.rotate(-barstat['pa'],axis=[0,0,1])
    # Get bar dimensions (distance of half-max)
    rmax=barstat['r']
    dimlim=(-rmax,rmax)
    hist_x,bins_x=self.hist1D(xmode='x',xlim=dimlim,**hist_kw)
    hist_y,bins_y=self.hist1D(xmode='y',xlim=dimlim,**hist_kw)
    hist_z,bins_z=self.hist1D(xmode='z',xlim=dimlim,**hist_kw)
    halfidx=np.array(range(hist_kw['xnbins']/2,hist_kw['xnbins']-1))
    histidx_barx=np.where(abs(hist_x[halfidx]) > abs(hist_x[halfidx]).max()/edgefact_small)
    histidx_bary=np.where(abs(hist_y[halfidx]) > abs(hist_y[halfidx]).max()/edgefact_small)
    histidx_barz=np.where(abs(hist_z[halfidx]) > abs(hist_z[halfidx]).max()/edgefact_small)
    barstat['x']=float(bins_x[halfidx][histidx_barx].max())
    barstat['y']=float(bins_y[halfidx][histidx_bary].max())
    barstat['z']=float(bins_z[halfidx][histidx_barz].max())
    if includehist:
      barstat['xbins']=bins_x ; barstat['xhist']=hist_x.data ; barstat['xmask']=hist_x.mask
      barstat['ybins']=bins_y ; barstat['yhist']=hist_y.data ; barstat['ymask']=hist_y.mask
      barstat['zbins']=bins_z ; barstat['zhist']=hist_z.data ; barstat['zmask']=hist_z.mask
    # Get bar amplitude
    if hist_phi_max==0.0:
      barstat['A']=0.
    else:
      barstat['A']=(hist_phi_max-hist_phi_min)/hist_phi_max
    # Return disk to original orientation
    self.rotate(barstat['pa'],axis=[0,0,1])
    # Return dictionary
    return barstat

  ##################################################################################################################################
  # METHOD FOR RETURNING INFORMATION ABOUT THE DM HALO PROFILE
  def get_haloprop(self,profile='nfw',**hist_kw):
    hist_kw['ptyp']='halo'
    return self.get_modelprop(profile,**hist_kw)

  ##################################################################################################################################
  # METHOD FOR RETURNING INFORMATION ABOUT THE DISK PROFILE
  def get_diskprop(self,profile='sechdisk',**hist_kw):
    hist_kw['ptyp']='disk'
    return self.get_modelprop(profile,**hist_kw)

  ##################################################################################################################################
  # METHOD FOR RETURNING INFORMATION ABOUT THE BULGE PROFILE
  def get_bulgeprop(self,profile='hernquist',**hist_kw):
    hist_kw['ptyp']='bulge'
    return self.get_modelprop(profile,**hist_kw)

  ##################################################################################################################################
  # METHOD FOR FITTING COMPONENTS WITH DIFFERENT PROFILES
  def get_modelprop(self,profile,**hist_kw):
    # Set constats
    sphmodlist=['nfw','hernquist','jaffe']
    cylmodlist=['sechdisk','expdisk']
    # Pars input
    profile=mmlpars.mml_pars(profile,type=str)
    profile=profile.lower()
    # Select the independent parameters and generate data
    if profile in sphmodlist:
      hist_kw['xmode']='rxyz'
      hist_kw['ymode']='rho'
      hist_kw['xscale']='log'
      y,x=self.hist1D(**hist_kw)
      x=x[:-1]+(x[1:]-x[:-1])/2.
      space='sph'
    elif profile in cylmodlist:
      hist_kw['xymode']='Rz'
      hist_kw['zmode' ]='rho'
      hist_kw['xscale']='log'
      hist_kw['yscale']='linear'
      z,x,y=self.hist2D(**hist_kw)
      x=x[:-1]+(x[1:]-x[:-1])/2.
      y=y[:-1]+(y[1:]-y[:-1])/2.
      X,Y=np.meshgrid(x,y)
      x=(X,Y)
      y=z
      space='cyl'
    else:
      raise Exception('Specified profile not currently supported: {}'.format(profile))
    # Fit model to data
    if profile == 'nfw':
      from astropysics.models import NFWModel as Model
    elif profile == 'hernquist':
      from astropysics.models import HernquistModel as Model
    elif profile == 'expdisk':
      from astropysics.models import ExponentialDiskModel as Model
    elif profile == 'sechdisk':
      from astropysics.models import ExponentialSechSqDiskModel as Model
    else:
      raise Exception('{} profile not currently supported.'.format(profile))
    model=Model()
    modparam = model.fitData(x,y)
    if self.verbose: print '[pNbody.get_modelprop] Check if model parameters make sense.'
    if self.verbose: print profile
    if self.verbose: print model.pardict
#    model.plot()
#    mplib.pyplot.show()
#    raise Exception('pause')
    return model.pardict

  ##################################################################################################################################
  # METHOD TO RETURN 1D HISTOGRAM (mainly parses modes and spaces)
  def hist1D(self,outdict=None,idxtot=None,setxlim=None,setylim=None,limdict=None,**hist_kw):
    '''
    Returns 1D histograms of data
    '''
    # Pars input
    opt=simplot.parsopt_hist1D(hist_kw)
    ptyp=opt['ptyp'] ; pgal=opt['pgal']
    outdict=mmlpars.mml_pars(outdict,default=False,type=bool)
    setxlim=mmlpars.mml_pars(setxlim,default=False,type=bool)
    setylim=mmlpars.mml_pars(setylim,default=False,type=bool)
    # Get data indices selected by type and galaxy
    if idxtot==None: idxtot=self.get_index(ptyp,pgal)
    # Pars xmode
    x=self.get_var(opt['xmode'],idxtot=idxtot)
    # Pars ymode
    if   opt['ymode'] == 'm':   w=self.mass[idxtot]
    elif opt['ymode'] == 'rho': w=self.mass[idxtot]
    elif opt['ymode'] == 'n':   w=np.ones(x.shape)
    else:                       w=self.get_var(opt['ymode'],idxtot=idxtot)
    # Set limits
    if opt['xlim'][0]==opt['xlim'][1]:
      opt['xlim']=self.get_limits(opt['xmode'],opt['xlimpad'],opt['xlimmir'],opt['xscale'],idxtot=idxtot,lim=limdict)
      setxlim=True
    if opt['ylim'][0]==opt['ylim'][1]:
      opt['ylim']=self.get_limits(opt['ymode'],opt['ylimpad'],opt['ylimmir'],opt['yscale'],idxtot=idxtot,lim=limdict)
      setylim=True
    # Create bins
    if   opt['xscale'] == 'linear': xbins=np.linspace(opt['xlim'][0],opt['xlim'][1],opt['xnbins'])
    elif opt['xscale'] == 'log':    xbins=mmlmath.logspace(opt['xlim'][0],opt['xlim'][1],opt['xnbins'],addzero=True)
    else: raise Exception('Unsupported xscale: {}'.format(opt['xscale']))
    # Create histogram
    hist,xbins=np.histogram(x,bins=xbins,weights=w)
    # Mask empty bins
    histn,binsn=np.histogram(x,bins=xbins)
#    xbins=np.ma.masked_where(histn==0,xbins)
    hist =np.ma.masked_where(histn==0,hist ) 
    histn=np.ma.masked_where(histn==0,histn) 
    # Adjust histogram based on method
    if   opt['method']=='hist':
      pass
    elif opt['method']=='mean':
      hist/=histn
    elif opt['method']=='dens':
      xbin1=xbins[:-1] ; xbin2=xbins[1:]
      if   opt['xmode']=='rxyz': # Spherical shells
        vols=(4.0/3.0)*math.pi*(np.power(xbin2,3)-np.power(xbin1,3))
      elif opt['xmode']=='rxy':    # Circular rings
        vols=math.pi*(np.power(xbin2,2)-np.power(xbin1,2))
      else:                      # Linear density
        vols=xbin2-xbin1
      hist/=vols
    else:
      raise Exception('Unrecognized method option: {}'.format(opt['method']))
    # Set limits if singular & unmask if masked
    if opt['xlim'][0]==opt['xlim'][1]: opt['xlim']=(xbins.min(),xbins.max())
    if opt['ylim'][0]==opt['ylim'][1]: opt['ylim']=(hist.data[hist.argmin()],hist.data[hist.argmax()])
    # Set means
    if setxlim: opt['xmean']=xbins.mean()
    if setylim:
      if np.all(hist.mask): opt['ymean']=hist.data.mean()
      else                : opt['ymean']=hist.mean()
    # Return output
    opt['x']=xbins
#    opt['x']=xbins.data
#    opt['xmask']=xbins.mask
    opt['y']=hist.data
    opt['ymask']=hist.mask
    if outdict: return opt
    else:       return hist,xbins
  
  ##################################################################################################################################
  # METHOD TO RETURN 2D HISTOGRAM (mainly parses modes and spaces)
  def hist2D(self,outdict=None,idxtot=None,setxlim=None,setylim=None,setzlim=None,limdict=None,**hist_kw):
    """
    Returns 2D histograms of data
    """
      # Set constants
    viewDict={'x':0,'y':1,'z':2}
    diffDict={'EkLz'    :('Ek'  ,'Lz'  ),
              'EpLz'    :('Ep'  ,'Lz'  ),
              'EtotLz'  :('Etot','Lz'  ),
              'EtotLtot':('Etot','Ltot')}
    # Pars input
    opt=simplot.parsopt_hist2D(hist_kw)
    outdict=mmlpars.mml_pars(outdict,default=False,type=bool)
    setxlim=mmlpars.mml_pars(setxlim,default=False,type=bool)
    setylim=mmlpars.mml_pars(setylim,default=False,type=bool)
    setzlim=mmlpars.mml_pars(setzlim,default=False,type=bool)
    # Pars viewray and rotate
    viewstr=opt['view']
    for iview in viewDict.keys():
      if iview not in viewstr: perpstr=iview
#    if not isinstance(view,str):
#      view=mmlpars.mml_pars(view,type=tuple,nelements=2)
#      self.rotate(angle=-view[0],axis=[0,0,1])
#      self.rotate(angle=-view[1],axis=[1,0,0])
#      viewstr='xy'
#    else:
#      viewstr=view
    # Get index
    if idxtot==None: idxtot=self.get_index(opt['ptyp'],opt['pgal'])
    # Get mode data
    if opt['xymode'] in diffDict: xmode,ymode=diffDict[opt['xymode']]
    else                        : xmode,ymode=opt['xymode'],opt['xymode']
    if opt['package']=='numpy':
      # Pars xymode
      if   opt['xymode']=='pos':
        x=self.get_var(viewstr[0],idxtot=idxtot)
        y=self.get_var(viewstr[1],idxtot=idxtot)
        z=self.get_var(perpstr   ,idxtot=idxtot)
      elif opt['xymode']=='scfpos':
        x=self.scfpos[idxtot,viewDict[viewstr[0]]]
        y=self.scfpos[idxtot,viewDict[viewstr[1]]]
        z=self.scfpos[idxtot,viewDict[perpstr   ]]
      elif opt['xymode']=='vel':
        x=self.get_var('v'+viewstr[0],idxtot=idxtot)
        y=self.get_var('v'+viewstr[1],idxtot=idxtot)
        z=self.get_var('v'+perpstr   ,idxtot=idxtot)
      elif opt['xymode']=='scfvel':
        x=self.scfvel[idxtot,viewDict[viewstr[0]]]
        y=self.scfvel[idxtot,viewDict[viewstr[1]]]
        z=self.scfvel[idxtot,viewDict[perpstr   ]]
      elif opt['xymode']=='Rz':
        x=self.get_var('rxy',idxtot=idxtot)
        y=self.get_var('z',idxtot=idxtot)
        z=None
      elif opt['xymode']=='EkLz':
        x=self.get_var('Ek',idxtot=idxtot)
        y=self.get_var('Lz',idxtot=idxtot)
        z=None
      elif opt['xymode']=='EpLz':
        x=self.get_var('Ep',idxtot=idxtot)
        y=self.get_var('Lz',idxtot=idxtot)
        z=None
      elif opt['xymode']=='EtotLz':
        x=self.get_var('Etot',idxtot=idxtot)
        y=self.get_var('Lz',idxtot=idxtot)
        z=None
      else:
        raise Exception('Unrecognized xymode option: {}'.format(opt['xymode']))
      # Pars zmode
      if   opt['zmode']=='n':   w=np.ones(x.shape)
      elif opt['zmode']=='m':   w=self.mass[idxtot]
      elif opt['zmode']=='rho': w=self.mass[idxtot]
      elif opt['zmode']=='vr':  w=self.vel[idxtot,viewDict[perpstr]]
      elif opt['zmode']=='vt':  w=self.get_var('v'+viewstr[0],idxtot=idxtot)
      else:                     w=self.get_var(opt['zmode'],idxtot=idxtot)
    # Set x, y & z limits
    if opt['xlim'][0]==opt['xlim'][1]:
      opt['xlim']=self.get_limits(xmode,opt['xlimpad'],opt['xlimmir'],opt['xscale'],idxtot=idxtot,lim=limdict)
      setxlim=True
    if opt['ylim'][0]==opt['ylim'][1]:
      opt['ylim']=self.get_limits(ymode,opt['ylimpad'],opt['ylimmir'],opt['yscale'],idxtot=idxtot,lim=limdict)
      setylim=True
    if opt['zlim'][0]==opt['zlim'][1]:
      opt['zlim']=self.get_limits(opt['zmode'],opt['zlimpad'],opt['zlimmir'],opt['zscale'],idxtot=idxtot,lim=limdict)
      setzlim=True
    # Create x bins
    if   opt['xscale']=='linear': xbins=np.linspace(opt['xlim'][0],opt['xlim'][1],opt['xnbins'])
    elif opt['xscale']=='log':    xbins=mmlmath.logspace(opt['xlim'][0],opt['xlim'][1],opt['xnbins'],addzero=True)
    else: raise Exception('Unsupported xscale: {}'.format(opt['xscale']))
    # Create y bins
    if   opt['yscale']=='linear': ybins=np.linspace(opt['ylim'][0],opt['ylim'][1],opt['ynbins'])
    elif opt['yscale']=='log':    ybins=mmlmath.logspace(opt['ylim'][0],opt['ylim'][1],opt['ynbins'],addzero=True)
    else: raise Exception('Unsupported yscale: {}'.format(opt['yscale']))

    # Histogram using pNbody routines
    if opt['package']=='pNbody':
      parkw={'space':opt['xymode'],'mode':opt['zmode'],'view':opt['view'],
             'size':(opt['xlim'][1]-opt['xlim'][0],opt['ylim'][1]-opt['ylim'][0]),
             'shape':(opt['xnbins'],opt['ynbins'])}
      # Create histogram
      param=libutil.extract_parameters([],parkw,self.defaultparameters)
      hist=self.CombiMap(param)
      if opt['zmode']=='n' or opt['zmode']=='m':
        nhist=copy.deepcopy(hist)
      else:
        param['mode']='m'
        nhist=self.CombiMap(param)
      
    # Using numpy histogram
    elif opt['package']=='numpy':
      # Clip to box
      if z!=None:
        boxdepth=np.mean([(opt['xlim'][1]-opt['xlim'][0])/2,(opt['ylim'][1]-opt['ylim'][0])/2])
        boxidx=np.where(abs(z)<=boxdepth)
        x=x[boxidx] ; y=y[boxidx] ; w=w[boxidx]
      # Create histogram
      hist,xbins,ybins=np.histogram2d(x,y,weights=w,bins=[xbins,ybins])
      if opt['zmode']=='n' or opt['zmode']=='m':
        nhist=copy.deepcopy(hist)
      else:
        nhist,nxbins,nybins=np.histogram2d(x,y,bins=[xbins,ybins])
        
    # Mask empty bins
    hist =np.ma.masked_where(nhist==0,hist ) 
    nhist=np.ma.masked_where(nhist==0,nhist) 
    # Adjust histogram based on method
    if   opt['method']=='hist':
      pass
    elif opt['method']=='mean':
      hist/=nhist
    elif opt['method']=='dens':
      if opt['xymode'] == 'Rz':
        X,Y=np.meshgrid(xbins[1:]**2-xbins[:-1]**2,ybins[1:]-ybins[:-1])
        hist/=math.pi*X*Y
      else:
        X,Y=np.meshgrid(xbins[1:]-xbins[:-1],ybins[1:]-ybins[:-1])
        hist/=X*Y
    else:
      raise Exception('Unrecognized method option: {}'.format(opt['method']))
    # Set limits if singular & unmask if masked  
    if opt['xlim'][0]==opt['xlim'][1]: opt['xlim']=(xbins.min(),xbins.max())
    if opt['ylim'][0]==opt['ylim'][1]: opt['ylim']=(ybins.min(),ybins.max())
    if opt['zlim'][0]==opt['zlim'][1]: opt['zlim']=(hist.data.flatten()[hist.argmin()],hist.data.flatten()[hist.argmax()])
    # Set means
    if setxlim: opt['xmean']=xbins.mean()
    if setylim: opt['ymean']=ybins.mean()
    if setzlim:
      if np.all(hist.mask): opt['zmean']=hist.data.mean()
      else                : opt['zmean']=hist.mean()
    # Transpose hist
#    hist=np.transpose(hist)
    # Return output
    opt['x'    ]=xbins
    opt['y'    ]=ybins
#    opt['x'    ]=xbins.data ; opt['xmask']=xbins.mask
#    opt['y'    ]=ybins.data ; opt['ymask']=ybins.mask
    opt['z'    ]=hist.data  ; opt['zmask']=hist.mask
#    mmlio.yorn(str((opt['zmode'],hist.min(),hist.max())))
    if outdict:
      return opt
    else:
      return hist,xbins,ybins

  
  ##################################################################################################################################
  # METHOD TO RETURN 3D HISTOGRAM (mainly parses modes and spaces)
  def hist3D(self,outdict=None,idxtot=None,setxlim=None,setylim=None,setzlim=None,setmlim=None,
             limdict=None,**hist_kw):
    """
    Returns 2D histograms of data
    """
    # Pars input
    opt=simplot.parsopt_hist3D(hist_kw)
    outdict=mmlpars.mml_pars(outdict,default=False,type=bool)
    setxlim=mmlpars.mml_pars(setxlim,default=False,type=bool)
    setylim=mmlpars.mml_pars(setylim,default=False,type=bool)
    setzlim=mmlpars.mml_pars(setzlim,default=False,type=bool)
    setmlim=mmlpars.mml_pars(setmlim,default=False,type=bool)
    # Get index
    if idxtot==None: idxtot=self.get_index(opt['ptyp'],opt['pgal'])
    # Pars space
    if   opt['space']=='pos':
      x=getval(self,'x')
      y=getval(self,'y')
      z=getval(self,'z')
    elif opt['space']=='scfpos':
      x=self.scfpos[:,0]
      y=self.scfpos[:,1]
      z=self.scfpos[:,2]
    elif opt['space']=='vel':
      x=getval(self,'vx')
      y=getval(self,'vy')
      z=getval(self,'vz')
    elif opt['space']=='scfvel':
      x=self.scfvel[:,0]
      y=self.scfvel[:,1]
      z=self.scfvel[:,2]
    elif opt['space']=='sphpos':
      x=np.sqrt(self.x()**2 + self.y()**2 + self.z()**2) # r
      y=np.arctan(self.y()/self.x()) # phi
      z=np.arctan(np.sqrt(self.x()**2 + self.y()**2)/self.z()) # theta
    elif opt['space']=='cylpos':
      x=np.sqrt(self.x()**2 + self.y()**2) # R
      y=np.arctan(self.y()/self.x()) # phi
      z=self.z() # z
    else:
      raise Exception('Unrecognized xymode option: {}'.format(opt['space']))
    # Pars wmode
    if   opt['mode']=='n':   w=np.ones(x.shape)
    elif opt['mode']=='m':   w=self.mass
    elif opt['mode']=='rho': w=self.mass
    else:                    w=self.get_var(opt['mode'])
    # Select data by type
    x=x[idxtot] ; y=y[idxtot] ; z=z[idxtot] ; w=w[idxtot]
    # Set x, y, z & m limits
    if opt['xlim'][0]==opt['xlim'][1]:
      opt['xlim']=self.get_limits(opt['space'],opt['xlimpad'],opt['xlimmir'],opt['xscale'],idxtot=idxtot,lim=limdict)
      setxlim=True
    if opt['ylim'][0]==opt['ylim'][1]:
      opt['ylim']=self.get_limits(opt['space'],opt['ylimpad'],opt['ylimmir'],opt['yscale'],idxtot=idxtot,lim=limdict)
      setylim=True
    if opt['zlim'][0]==opt['zlim'][1]:
      opt['zlim']=self.get_limits(opt['space'],opt['zlimpad'],opt['zlimmir'],opt['zscale'],idxtot=idxtot,lim=limdict)
      setzlim=True
    if opt['mlim'][0]==opt['mlim'][1]:
      opt['mlim']=self.get_limits(opt['mode'],opt['mlimpad'],opt['mlimmir'],opt['mscale'],idxtot=idxtot,lim=limdict)
    # Create x bins
    if   opt['xscale']=='linear': xbins=np.linspace(opt['xlim'][0],opt['xlim'][1],opt['xnbins'])
    elif opt['xscale']=='log':    xbins=mmlmath.logspace(opt['xlim'][0],opt['xlim'][1],opt['xnbins'],addzero=True)
    else: raise Exception('Unsupported xscale: {}'.format(opt['xscale']))
    # Create y bins
    if   opt['yscale']=='linear': ybins=np.linspace(opt['ylim'][0],opt['ylim'][1],opt['ynbins'])
    elif opt['yscale']=='log':    ybins=mmlmath.logspace(opt['ylim'][0],opt['ylim'][1],opt['ynbins'],addzero=True)
    else: raise Exception('Unsupported yscale: {}'.format(opt['yscale']))
    # Create z bins
    if   opt['zscale']=='linear': zbins=np.linspace(opt['zlim'][0],opt['zlim'][1],opt['znbins'])
    elif opt['zscale']=='log'   : zbins=mmlmath.logspace(opt['zlim'][1],opt['zlim'][1],opt['znbins'],addzero=True)
    else: raise Exception('Unsupported zscale: {}'.format(opt['zscale']))
    # Create histogram
    data=np.vstack((x,y,z)).T
    bins=(xbins,ybins,zbins)
    hist,bins=np.histogramdd(data,weights=w,bins=bins)
    if opt['mode']=='n' or opt['mode']=='m':
      nhist=copy.deepcopy(hist)
    else:
      nhist,nbins=np.histogramdd(data,bins=bins)
    # Mask empty bins
    hist =np.ma.masked_where(nhist==0,hist ) 
    nhist=np.ma.masked_where(nhist==0,nhist) 
    # Adjust histogram based on method
    if   opt['method']=='hist':
      pass
    elif opt['method']=='mean':
      hist/=nhist
    elif opt['method']=='dens':
      if opt['xymode'] == 'cylpos':
        X,Y=np.meshgrid(xbins[1:]**2-xbins[:-1]**2,ybins[1:]-ybins[:-1])
        hist/=math.pi*X*Y
      else:
        X,Y=np.meshgrid(xbins[1:]-xbins[:-1],ybins[1:]-ybins[:-1])
        hist/=X*Y
    else:
      raise Exception('Unrecognized method option: {}'.format(opt['method']))
    # Set limits if singular & unmask if masked  
    if opt['xlim'][0]==opt['xlim'][1]: opt['xlim']=(xbins.min(),xbins.max())
    if opt['ylim'][0]==opt['ylim'][1]: opt['ylim']=(ybins.min(),ybins.max())
    if opt['zlim'][0]==opt['zlim'][1]: opt['zlim']=(hist.data.flatten()[hist.argmin()],hist.data.flatten()[hist.argmax()])
    # Set means
    if setxlim: opt['xmean']=xbins.mean()
    if setylim: opt['ymean']=ybins.mean()
    if setzlim:
      if np.all(hist.mask): opt['zmean']=hist.data.mean()
      else                : opt['zmean']=hist.mean()
    # Transpose hist
#    hist=np.transpose(hist)
    # Return output
    opt['x'    ]=xbins
    opt['y'    ]=ybins
#    opt['x'    ]=xbins.data ; opt['xmask']=xbins.mask
#    opt['y'    ]=ybins.data ; opt['ymask']=ybins.mask
    opt['z'    ]=hist.data  ; opt['zmask']=hist.mask
    if outdict:
      return opt
    else:
      return hist,xbins,ybins

  ##################################################################################################################################
  def plotscat(self,plotfile=None,overwrite=None,
               ptyp=None,pgal=None,rmax=None,zdir=None, 
               colormeth=None,sizemeth=None,alphameth=None,
               pntsiz=None,pntmrk=None,
               askuser=None,verbose=None,**input_kw):
    from mpl_toolkits.mplot3d import Axes3D
    # Set constants
    pntsizDEF=1.
    pntmrkDEF='.'
    cmethLIST=['none','galaxy','type'] ; cmethDEF='type'
    smethLIST=['none','mass'] ; smethDEF='none'
    amethLIST=['none','type'] ; amethDEF='type'
    galLIST=self.get_gallist().remove(0)
    typLIST=self.get_typlist()
    # Pars input
    plotfile=mmlpars.mml_pars(plotfile,default='scatter.png',type=str)
    overwrite=mmlpars.mml_pars(overwrite,default=False,type=bool)
    colormeth=mmlpars.mml_pars(colormeth,default=cmethDEF,list=cmethLIST)
    sizemeth=mmlpars.mml_pars(sizemeth,default=smethDEF,list=smethLIST)
    alphameth=mmlpars.mml_pars(alphameth,default=amethDEF,list=amethLIST)
    pntsiz=mmlpars.mml_pars(pntsiz,default=pntsizDEF,type=float,min=0.)
    pntmrk=mmlpars.mml_pars(pntmrk,default=pntmrkDEF,type=str)
    askuser=mmlpars.mml_pars(askuser,default=False,type=bool)
    verbose=mmlpars.mml_pars(verbose,default=False,type=bool)
    plot_kw,input_kw=mmlplot.parsopt(input_kw,outextra=True)
    if not plot_kw['outplot']:
      if not mmlfiles.prep_write(plotfile,overwrite=overwrite,askuser=askuser): return None
    # Select data
    idx=self.get_index(ptyp=ptyp,pgal=pgal,rmax=rmax)
    pos=self.pos[idx,:]
    typ=self.tpe[idx]
    gal=self.gal[idx]
    N=len(gal)
    # Set sizes
    siz=np.array([pntsiz]*N)
    if   sizemeth=='none': pass
    else: raise Exception('Invalid sizemeth: {}'.format(sizemeth))
    # Set colors
#    clr=np.array((plot_kw['fgclr'],)*N)
#    if   colormeth=='none'  : pass
#    elif colormeth=='galaxy':
#      clrdict=mmlplot.get_galcolor()
#      for iclrkey in clrdict.keys():
#        irgb=mmlplot.get_mmlcolor(clrdict[iclrkey])
#        clr[gal==iclrkey,:]=irgb
#    elif colormeth=='type'  :
#      clrdict=mmlplot.get_typcolor()
#      for iclrkey in clrdict.keys():
#        irgb=mmlplot.get_mmlcolor(clrdict[iclrkey])
#        clr[typ==iclrkey,:]=irgb
    if   colormeth=='none'  : pass
    elif colormeth=='galaxy': clrdict=mmlplot.get_galcolor()
    elif colormeth=='type'  : clrdict=mmlplot.get_typcolor(outstr=True)
    else: raise Exception('Invalid colormeth: {}'.format(colormeth))
    # Set alphas
    if   alphameth=='none'  : alpdict={'all':1.0}
    elif alphameth=='galaxy': alpdict=mmlplot.get_galalpha()
    elif alphameth=='type'  : alpdict=mmlplot.get_typalpha(outstr=True)
    else: raise Exception('Invalid alphameth: {}'.format(alphameth))
    # Plot
    fig = mplib.pyplot.figure()
    ax = Axes3D(fig)
    for igal in galLIST:
      for ityp in typLIST:
        # Color
        if   colormeth=='none'  : iclr=plot_kw['fgclr']
        elif colormeth=='galaxy': iclr=mmlplot.get_mmlcolor(clrdict[igal])
        elif colormeth=='type'  : iclr=mmlplot.get_mmlcolor(clrdict[ityp])
        # Alpha
        if   alphameth=='none'  : ialp=1.0
        elif alphameth=='galaxy': ialp=alpdict[igal]
        elif alphameth=='type'  : ialp=alpdict[ityp]
        # Index
        iidx=np.logical_and(gal==igal,typ==ityp)
        if np.any(iidx):
          ax.scatter(pos[iidx,0],pos[iidx,1],pos[iidx,2],
                     zdir=zdir,c=iclr,s=siz,alpha=ialp,
                     marker=pntmrk)#,**plot_kw)
#    ax.scatter(pos[:,0],pos[:,1],pos[:,2],zdir=zdir,c=clr,s=siz,marker=pntmrk)#,**plot_kw)
    # Display image
    if plot_kw['showplot']: fig.show()
    # Return output
    if plot_kw['outplot']: return fig
    else                 : 
      fig.savefig('{}.{}'.format(plotfile,plot_kw['plotext']))
      if verbose:
        mmlio.verbose('Created plot:')
        print '    {}'.format(plotfile)
    return None

  ##################################################################################################################################
  # METHOD TO CREATE SIMPLOT
  def plotsing(self,**input_kw):
    """
    Plots a simple projection of the particles
    """
    out=simplot.plot_sing(self,**input_kw)
    if out != None: return out
    else: return

  ##################################################################################################################################
  # METHOD TO PLOT HISTOGRAMS
  def plot_panels(self,fname=None,overwrite=False,calcdict=None,calcfile=None,unitDict=None,
                  h1DList=None,h2DList=None,typList=None,galList=None,rayList=None,
                  adjustylim=False,adjustxlim=False,**plot_kw):
    """
    Plots panels of simulation info
    """
    # Get calcdict
    calcdict=self.calc(calcfile,overwrite=False,calcdict=calcdict,unitDict=unitDict,
                       h1DList=h1DList,h2DList=h2DList,typList=typList,galList=galList,rayList=rayList,
                       adjustylim=adjustylim,adjustxlim=adjustxlim)
    # Plot calcdict
    plotdict=simplot.plot_mult(fname,overwrite=overwrite,calcdict=calcdict,
                               h1DList=h1DList,h2DList=h2DList,typList=typList,galList=galList,rayList=rayList,
                               **plot_kw)
    # Return output
    return plotdict

  ##################################################################################################################################
  # METHOD TO SAVE HISTOGRAMS TO FILE
  def calc(self,fname=None,overwrite=False,calcdict=None,unitDict=None,
           h1DList=None,h2DList=None,typList=None,galList=None,rayList=None,
           adjustylim=False,adjustxlim=False,skipsave=None,histpackage=None,
           center=None,centergal=None,centertyp=None,centermeth=None,
           limdict=None,**calc_kw):
    '''
    Method to save histograms for a particle set to file
    '''
    # Set constants
    fnameTAG='snapcalc'
    curdir=os.getcwd()
    unitDictDEF={'UnitLength':pNbody.units.Unit_kpc,'UnitMass':pNbody.units.Unit_Ms,
                 'UnitVelocity':pNbody.units.Unit_kms,'UnitTime':pNbody.units.Unit_Myr*1000.}
    unitDictDEF['UnitTime'].symbol='Gyr'
    h1DListDEF=['rxyz_rho','z_n','vz_n','phi_n']
    h2DListDEF=['pos_m','pos_vr','pos_gal','pos_scfpot','pos_scfpotm2']#,'scfpos_scfpotm2']
    typListDEF=self.get_typlist()
    galListDEF=self.get_gallist()
    rayListDEF=['xy','yz','xz']
    hist1D_kw=simplot.parsopt_hist1D() #; hist1D_kw['xlimpad']=0.1 ; hist1D_kw['ylimpad']=0.1
    hist2D_kw=simplot.parsopt_hist2D() ; hist2D_kw['zlimpad']=0.1
    # Pars file input
    skipsave=mmlpars.mml_pars(skipsave,default=False,type=bool)
    if not skipsave:
      fname=mmlpars.mml_pars(fname,default=os.path.join(curdir,fnameTAG),type=str)
      overwrite=mmlpars.mml_pars(overwrite,type=bool)
      if not mmlfiles.prep_write(fname,overwrite): return simcalc.loadcalc(fname)
    # Pars other input
    histpackage=mmlpars.mml_pars(histpackage,default='numpy',list=['numpy','pNbody'])
    unitDict=mmlpars.mml_pars(unitDict,default=unitDictDEF,type=dict)
    h1DList=mmlpars.mml_pars(h1DList,default=h1DListDEF,type=list) ; nh1D=len(h1DList)
    h2DList=mmlpars.mml_pars(h2DList,default=h2DListDEF,type=list) ; nh2D=len(h2DList)
    typList=mmlpars.mml_pars(typList,default=typListDEF,type=list) ; ntyp=len(typList)
    galList=mmlpars.mml_pars(galList,default=galListDEF,type=list) ; ngal=len(galList)
    rayList=mmlpars.mml_pars(rayList,default=rayListDEF,type=list) ; nray=len(rayList)
    # Get info about the galaxy
    statDict=self.get_statdict('general',centergal=1,centertyp='disk',centermeth=centermeth)
    # Center on primary's center of mass by type
    oldcm=self.get_center(method=centermeth,centergal=centergal,centertyp=centertyp,move=center)
    # Convert units
    unitDict=self.convert_units(unitDict)
    # Preallocate
    calcdict=simcalc.initcalc(calcdict,snapdict=self,
                              h1DList=h1DList,h2DList=h2DList,typList=typList,rayList=rayList,galList=galList)
    calcdict['statDict']=statDict
    calcdict['unitDict']=unitDict
    calcdict['ngal']=max(galListDEF)
    for igal in galListDEF: calcdict['com{}'.format(igal)]=self.get_center(method='mass',centergal=igal,move=False)
    # Loop over types
    for ityp in typList:
      # Loo[ over galaxies
      for igal in galList:
        # Loop over 1D histograms
        for ih1D in h1DList:
#          print '[pNbody.calc] check if the object has memory of where it came from [it does]'
#          print calcdict['hist1D'][ih1D][ityp][igal]['xmode']
          iopt=calcdict['hist1D'][ih1D][ityp][igal]
          iopt['ptyp']=ityp
          iopt['pgal']=igal
          iopt['xmode'],iopt['ymode']=ih1D.split('_')
          xmethod,iopt['xscale'],iopt['xlimmir']=simplot.get_histdef(iopt['xmode'])
          ymethod,iopt['yscale'],iopt['ylimmir']=simplot.get_histdef(iopt['ymode'])
          iopt['method']=ymethod
          setxlim=False ; setylim=False
          if adjustxlim or iopt['xlim'][0]==iopt['xlim'][1]:
            iopt['xlim']=self.get_limits(iopt['xmode'],iopt['xlimpad'],iopt['xlimmir'],iopt['xscale'],ptyp='all',pgal=igal,lim=limdict)
            setxlim=True
          if adjustylim or iopt['ylim'][0]==iopt['ylim'][1]:
            iopt['ylim']=self.get_limits(iopt['ymode'],iopt['ylimpad'],iopt['ylimmir'],iopt['yscale'],ptyp='all',pgal=igal,lim=limdict)
            setylim=True
          iopt=self.hist1D(outdict=True,setxlim=setxlim,setylim=setylim,limdict=limdict,**iopt)
          calcdict['hist1D'][ih1D][ityp][igal]=iopt
          del iopt
        # Loop over 2D histograms
        for ih2D in h2DList:
          # Loop over rays
          for iray in rayList:
            # Check for missing keys
            iopt=calcdict['hist2D'][ih2D][ityp][igal][iray]
#            if iopt['zmode'] in ['scfpotm2']:
#              print 'before: {} {} {}: {}'.format(iopt['zmode'],iopt['zscale'],iopt['zlimmir'],iopt['zlim'])
            iopt['package']=histpackage
            iopt['center']=oldcm
            iopt['ptyp']=ityp
            iopt['pgal']=igal
            iopt['view']=iray
            iopt['xymode'],iopt['zmode']=ih2D.split('_')
            xmethod,iopt['xscale'],iopt['xlimmir']=simplot.get_histdef(iopt['xymode'])
            ymethod,iopt['yscale'],iopt['ylimmir']=simplot.get_histdef(iopt['xymode'])
            zmethod,iopt['zscale'],iopt['zlimmir']=simplot.get_histdef(iopt['zmode'])
            iopt['method']=zmethod
            if iopt['xymode']=='pos':
              if igal==0: ipgal=1
              else:       ipgal=0
            else: ipgal=igal
            setxlim=False ; setylim=False ; setzlim=False
            if iopt['xymode']=='pos':
              rlim=self.get_limits(iopt['xymode'],iopt['xlimpad'],iopt['xlimmir'],iopt['xscale'],ptyp=ityp,pgal=ipgal,lim=limdict)
              if adjustxlim or iopt['xlim'][0]==iopt['xlim'][1]:
                iopt['xlim']=rlim
                setxlim=True
              if adjustxlim or iopt['ylim'][0]==iopt['ylim'][1]:
                iopt['ylim']=rlim
                setylim=True
              if adjustylim or iopt['zlim'][0]==iopt['zlim'][1]:
                iopt['zlim']=self.get_limits(iopt['zmode'],iopt['zlimpad'],iopt['zlimmir'],iopt['zscale'],ptyp=ityp,pgal=igal,lim=limdict)
                setzlim=True
            iopt=self.hist2D(outdict=True,setxlim=setxlim,setylim=setylim,setzlim=setzlim,limdict=limdict,**iopt)
#            if iopt['zmode'] in ['scfpotm2']:
#              print 'after : {} {} {}: {}'.format(iopt['zmode'],iopt['zscale'],iopt['zlimmir'],iopt['zlim'])
            calcdict['hist2D'][ih2D][ityp][igal][iray]=iopt
            del iopt
    # Save files
    if not skipsave: simcalc.savecalc(fname,calcdict)
    # Return dictionary
    return calcdict



    
