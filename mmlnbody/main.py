#!/usr/bin/python
import sys,os,shutil,glob,copy,pprint,scipy,math
#import astropysics.models as apymods
import numpy as np
from mmlutils import *
from mmlastro import mmlprofiles
import simlist,simlyze,simstat,simcalc,simhist,simplot,simfile,mmlbuildgal,mmlgadget,mmlscf,mmlgalfit,icgen,compare,mml_pNbody
import pNbody

####################################################################################################################################
# METHOD TO GET MMLSIM OBJECT FROM RUNTAG
def runtag2simobj(runtag=None):
    """
    Returns the mmlsim object for a given runtag
    """
    runlist=simfile.fpar_taglist('icgen','mmlsim')
    if runtag in runlist:
        simdict=simfile.loadfpar('icgen','mmlsim',runtag)
        ifodict=simfile.loadfpar('icgen',simdict['runtyp'],runtag)
        simobj=mmlsim(infodict=ifodict,**simdict)
    else:
        simobj=asksimobj()
    return simobj
def asksimobj(infodict=None,runtyp=None,**simdict):
    """
    Returns an mmlsim object provided user input
    """
    if runtyp not in simlist.LIST_RUNTYPS:
        runtyp=mmlio.askselect('Select a valid run type:',simlist.LIST_RUNTYPS)
    ifodict=simfile.parsfpar('icgen',runtyp,fpar=infodict,askuser=True)
    simdict['runtyp']=runtyp
    simdict['runtag']=simfile.fpar2tag('icgen',runtyp,ifodict)
    simdict=simfile.parsfpar('icgen','mmlsim',fpar=simdict,askuser=True)
    simobj=mmlsim(infodict=ifodict,**simdict)
    return simobj

####################################################################################################################################
####################################################################################################################################
# MMLSIM CLASS
class mmlsim(dict):
    """
    A dictionary subclass with simulation keys. (See icgen.list_fpar for 'mmlsim')
    """

    def __init__(self,infodict={},**inkw):
        """
        Initializes mmlsim objects
        """
        # Get mmlsim dict
        self.update(simfile.parsfpar('icgen','mmlsim',fpar=inkw))
        # Get infodict
        infodict=simfile.parsfpar('icgen',inkw['runtyp'],fpar=infodict)
        setattr(self,'infodict',infodict)
        # Add file dictionary
        if self['runtag'] not in simlist.LIST_RUNTYPS:
            fdict=self.mkfiledict()
            setattr(self,'fdict',fdict)
            ntot=self.get_ntot()
            setattr(self,'ntot',ntot)

    def run(self,mtype,method,**method_kw):
        """
        Handles running different types of methods for runs
        """
        # Pars input
        mtype=mmlpars.mml_pars(mtype,list=simlist.LIST_METHTYP)
        method=mmlpars.mml_pars(method,list=simlist.DICT_METHODS[mtype])
        # Initialize output
        out=None
        # Find correcty method subclass
        if   mtype=='general' : out=general(self,method,**method_kw)
        elif mtype=='compare' : out=compare.run(self,method,**method_kw)
        elif mtype=='buildgal': out=mmlbuildgal.run(self,method,**method_kw)
        elif mtype=='icgen'   : out=icgen.run(self,method,**method_kw)
        elif mtype=='gadget'  : out=mmlgadget.run(self,method,**method_kw)
        elif mtype=='scf'     : out=mmlscf.run(self,method,**method_kw)
        elif mtype=='galfit'  : out=mmlgalfit.run(self,method,**method_kw)
        elif mtype=='simlyze' : out=simlyze.run(self,method,**method_kw)
        elif mtype=='stat'    : out=simstat.run(self,method,**method_kw)
        elif mtype=='calc'    : out=simcalc.run(self,method,**method_kw)
        elif mtype=='hist'    : out=simhist.run(self,method,**method_kw)
        elif mtype=='plots'   : out=simplot.run(self,method,**method_kw)
        elif mtype=='file'    : out=simfile.run(self,method,**method_kw)
        else: raise Exception('Invalid method type: {}'.format(mtype))
        # Return output
        return out

    def mkfiledict(self,**input_kw):
        """
        Returns a dictionary of simulation files and directories
        """
        return simfile.files(simstr=self,**input_kw)

    def get_filelist(self,method,**input_kw):
        """
        Returns a list of relevant simulation files and directories
        """
        return simfile.get_filelist(method,simstr=self,**input_kw)

    def subpbs(self,method,ptype=None,pbsfile=None,schfile=None):
        """
        Submits pbs scripts for different methods
        """
        # Set constants
        methLIST=['buildgal','gadget','scf','galfit']
        typeLIST=simlist.LIST_PTYPS
        # Pars input
        method=mmlpars.mml_pars(method,type=str,list=methLIST)
        if method=='galfit': typeLIST.append('all')
        ptype=mmlpars.mml_pars(ptype,default='disk',type=str,list=typeLIST)
        # Determine files based on method
        files=self.mkfiledict(memcomp='accre',checkaccess=True)[method]
        if method in ['scf','galfit']: files=files[ptype]
        pbsfile=mmlpars.mml_pars(pbsfile,default=files['input']['pbs'],type=str)
        schfile=mmlpars.mml_pars(schfile,default=files['input']['sched'],type=str)
        # Submit run
        mmlio.pbs_submit(pbsfile,schfile=schfile)
        return

    def get_Mvir(self,galaxy=None):
        """
        Returns the virial mass for a given simulation
        """
        if self['runtyp']=='galsim': 
            Mvir=self.infodict['mvir']
        else:
            galaxy=mmlpars.mml_pars(galaxy,default=1,list=[1,2])
            Mvir=simfile.loadfpar('icgen','galsim',self.infodict['gal{}'.format(galaxy)])['mvir']
        return Mvir

    def get_Rvir(self,galaxy=None,**exkw):
        """
        Returns the virial radius for a given simulations (in pc)
        """
        Mvir=self.get_Mvir(galaxy=galaxy)
        exkw['units']='galaxy'
        Rvir=mmlprofiles.mvir2rvir(Mvir,**exkw)
#        Rvir=apymods.NFWModel().Mvir_to_Rvir(Mvir)
        return Rvir

    def get_profpar(self,overwrite=None,**exkw):
        overwrite=mmlpars.mml_pars(overwrite,default=False,type=bool)
        if os.path.isfile(self.fdict['profpar']) and not overwrite:
            out=mmlio.rwdict('R',self.fdict['profpar'])
        else:
            pnobj=self.get_pnbody(ftype='gadget',fext='ic')
            typList=pnobj.get_typlist()
            galList=pnobj.get_gallist() ; galList.remove(0)
            out={}
            for ityp in typList:
                for igal in galList:
                    if ityp=='disk': profid='FLATDISK'
                    else           : profid='HERNQUIST'
                    out[ityp]=pnobj.fitprofile(ptyp=ityp,pgal=igal,profid=profid,**exkw)
            mmlio.rwdict('W',self.fdict['profpar'],out,overwrite=overwrite)
        return out

    def get_flagic(self,fext):
        """
        Determines if an extension is for IC file
        """
        fext=mmlpars.mml_pars(fext,type=str)
        if fext in ['ic','_ic','_000ic']:
            flag_ic=True
        else:
            flag_ic=False
        return flag_ic

    def get_lastsnap(self,ftype=None,ptype=None,outfext=None):
        """
        Determines the full path to the last snapshot for a simulation
        """
        # Pars input
        ftype=mmlpars.mml_pars(ftype,default='gadget',type=str,list=['gadget','scf'])
        ptype=mmlpars.mml_pars(ptype,default='disk',type=str)
        outfext=mmlpars.mml_pars(outfext,default=False,type=bool)
        # Get files
        files=self.fdict
        # Identify snapshot base
        if ftype=='scf':
            snapbase=files['scf'][ptype]['output']['snapbase']
        elif ftype=='gadget':
            snapbase=files['gadget']['output']['snapbase']
        # Search for snapshots
        snaplist=glob.glob(snapbase+'*')
        # Default to IC if no snapshots
        if len(snaplist)==0:
            fext='_000ic'
            lastsnap=self.get_snapname(fext=fext,ftype=ftype,ptype=ptype)
        # Otherwise take last snapshot
        else:
            lastsnap=sorted(snaplist)[-1]
            fext=lastsnap.split(snapbase)[-1]
        # Return output
        if outfext:
            return fext
        else:
            return lastsnap
        
    def get_snapname(self,fext=None,ftype=None,ptype=None,**extra_kw):
        """
        Determines the full path to a simulation snapshot with extension fext
        """
        # Set constants
        typList=simlist.LIST_PTYPS
        fextDEF='_000ic'
        # Pars input
        ftype=mmlpars.mml_pars(ftype,default='gadget',type=str,list=['gadget','scf','buildgal'])
        ptype=mmlpars.mml_pars(ptype,default='disk',type=str)
        fext=mmlpars.mml_pars(fext,type=str)
        # Get files
        files=self.fdict
        # Determine if extension is icfile
        flag_ic=self.get_flagic(fext)
        # Get info on file
        if flag_ic:
            if   ftype=='scf'     : fname=files['scf'][ptype]['output']['snapbase']+fext
            elif ftype=='gadget'  : fname=files['gadget']['input']['ic']
            elif ftype=='buildgal': fname=files['icgen']['bgic']
        else:
            if ftype=='scf'       : fname=files['scf'][ptype]['output']['snapbase']+fext
            elif ftype=='gadget'  : fname=files['gadget']['output']['snapbase']+fext
            elif ftype=='buildgal': fname=files['icgen']['bgic']
        # Return file name and IC flag
        return fname,flag_ic

    def get_loadkeys(self,ftype=None,ptype=None,fname=None,fext=None,idgal2=None,
                     snapdict=None,askuser=None,noidgal2=None,outextra=None,**extra_kw):
        """
        Returns a dictionary of keys for loading snapshot data
        """
        # Set constants
        fextDEF='_000ic'
        ftypList=['gadget','scf','buildgal']
        # Pars input
        askuser=mmlpars.mml_pars(askuser,default=False,type=bool)
        noidgal2=mmlpars.mml_pars(noidgal2,default=False,type=bool)
        outextra=mmlpars.mml_pars(outextra,default=False,type=bool)
        # Get info from pNbody object
        if snapdict!=None:
            fname=snapdict.p_name[0]
            ftype=snapdict.ftype.split('Nbody_')[-1]
            if ftype=='scf': ptype=snapdict.ptype
        # Get info on file type
        if ftype not in ftypList and askuser: ftype=mmlio.askselect('What type of snapshot should be loaded?',ftypList)
        else: ftype=mmlpars.mml_pars(ftype,default='gadget',list=ftypList)
        if ftype=='scf':
            if askuser: ptype=mmlio.askselect('What particle type should be loaded?',simlist.LIST_PTYPS)
            else      : ptype=mmlpars.mml_pars(ptype,default='disk',list=simlist.LIST_PTYPS)
        else: ptype=None
        # Handle absence of fext or fname
        if not isinstance(fext,str) and not isinstance(fname,str) and not askuser:
            raise Exception('Must provide either fext or fname')
        # Get file extension
        if isinstance(fname,str):
            fext=self.get_fext(fname,ftype=ftype,ptype=ptype)
        else:
            fname=None
        if not isinstance(fext,str) and askuser:
            fext=mmlio.askquest('What snapshot extension should be loaded?',default=fextDEF,dtype='str')
        else:
            fext=mmlpars.mml_pars(fext,default=fextDEF,type=str)
        # Get file name
        if not isinstance(fname,str):
            fname,flag_ic=self.get_snapname(fext,ftype=ftype,ptype=ptype)
        # Get idgal2
        if idgal2==None and not noidgal2: idgal2=self.get_idgal2(ftype=ftype)
        # Get ic flag
        flag_ic=self.get_flagic(fext)
        if flag_ic: fext='_000ic'
        # Create loadkeys dictionary
        loadkeys=dict(ftype=ftype,ptype=ptype,fname=fname,fext=fext,idgal2=idgal2,flag_ic=flag_ic)
        # Return loadkeys
        if outextra:
            return loadkeys,extra_kw
        else:
            return loadkeys

    def get_snapdict(self,**input_kw): return self.get_pnbody(**input_kw)

    def get_calcdict(self,**input_kw):
        """
        Loads and returns a calcdict object
        """
        return simcalc.get_calcsnap(self,**input_kw)

    def get_softenings(self,parmfile=None):
        """
        Loads softening lengths for a simulation
        """
        fdict=self.fdict
        parmfile=mmlpars.mml_pars(parmfile,default=fdict['gadget']['input']['param'],type=str)
        params=mmlio.rwdict('R',parmfile)
        softkeylist=['SofteningGas','SofteningHalo','SofteningDisk',
                     'SofteningBulge','SofteningStars','SofteningBndry']
        softdict={}
        for ikey in softkeylist:
            ikey_short=ikey.split('Softening')[-1].lower()
            if ikey in params:
                softdict[ikey_short]=params[ikey]
            else:
                softdict[ikey_short]=0.
        return softdict

    def get_pnbody(self,**input_kw):
        """
        Loads and returns a pNbody object from a simulation snapshot
        """
        # Currently only loads one type at a time for SCF files
        # Set constants
        typList=simlist.LIST_PTYPS
        # Pars input
        loadkeys=self.get_loadkeys(**input_kw)
        ftype  =loadkeys['ftype'  ]
        ptype  =loadkeys['ptype'  ]
        fname  =loadkeys['fname'  ]
        fext   =loadkeys['fext'   ]
        idgal2 =loadkeys['idgal2' ]
        flag_ic=loadkeys['flag_ic']
        # Get files
        files=self.fdict
        # Get supporting files
        if ftype=='scf':
            parmfile=files[ftype]['unitfile']
            makefile=None
            scffiledict={}
        elif ftype=='gadget':
            parmfile=files[ftype]['input']['param']
            makefile=files[ftype]['input']['makefile']
            scffiledict={}
            for ikey in typList:
                scffiledict[ikey]=files['scf'][ikey]['output']['snapbase']+fext
        elif ftype=='buildgal':
            parmfile=files[ftype]['input']['units']
            makefile=None
            scffiledict={}
        else: raise Exception('Invalid file type: {}'.format(ftype))
        # Get snapfile used to set limits
        snapfile0=self.get_snapfile0(ftype=ftype)
        # Return pNbody object
        return pNbody.Nbody(p_name=fname,ftype=ftype,
                            unitsfile=parmfile,makefile=makefile,scffiledict=scffiledict,snapfile0=snapfile0,
                            idgal2=idgal2,flag_ic=flag_ic)

    def get_primary(self,striptag=False):
        """
        Returns the mmlsim object for the primary galaxy
        """
        # Select correct sim dict for isolated primary
        if   self['runtyp']=='galsim': runtag=self['runtag']
        elif self['runtyp']=='intsim': runtag=self.infodict['gal1']
        # Get info from runtag
        galobj=simfile.loadfpar('icgen','galsim',runtag)
        addruntag=galobj['addtag']
        # Strip additional tags as needed to get to base run
        if striptag:
            # Runs without additional tags
            if len(addruntag) == 0:
                runtag0=runtag
            # Runs with additional tags
            else:
                runtag0=runtag.split(addruntag)[0]
        else:
            runtag0=runtag
        # Return simulation dictionary for primary
        simstr0=runtag2simobj(runtag0)
        return simstr0

    def get_snapfile0(self,ftype=None,verbose=None,**extra_kw):
        """
        Returns the comparison snapshot for a simulation
        """
        # Pars input
        ftype=mmlpars.mml_pars(ftype,default='gadget',type=str)
        verbose=mmlpars.mml_pars(verbose,default=False,type=bool)
        # Select primary galaxy isolated run
        simstr0=self.get_primary(striptag=True)
        runtag0=simstr0['runtag']
        # Get the associated file names
        files0=simstr0.fdict
        # Select static IC file
        if   ftype=='buildgal': snap0key='bgic'
        elif ftype=='gadget'  : snap0key='gdic'
        else: raise Exception('Invalid file types: {}'.format(ftype))
        snapfile0=files0['icgen'][snap0key]
        # Return the associated base file
        if verbose: print '{} snapfile0: {} ({})'.format(self['runtag'],snapfile0,runtag0)
        return snapfile0

    def get_idgal2(self,ftype=None,ptype=None):
        """
        Determines the particle ID for a second galaxy
        """
        # Pars input
        ftype=mmlpars.mml_pars(ftype,default='gadget',type=str)
        ptype=mmlpars.mml_pars(ptype,default='disk',type=str)
        ftype=ftype.lower() ; ptype=ptype.lower()
        # Select primary galaxy isolated run
        simstr0=self.get_primary(striptag=False)
        # Get the associated file names
        files0=simstr0.fdict
        # Determine what files to use for each type
        if ftype=='gadget':
            fname=files0['gadget']['input']['ic']
            funit=files0['gadget']['input']['param']
        elif ftype=='scf':
            fname=sorted(glob.glob(files0['scf'][ptype]['output']['snapbase']+'*'))[0]
            funit=files0['scf']['unitfile']
        elif ftype=='buildgal':
            fname=files0['icgen']['bgic']
            funit=files0['buildgal']['input']['units']
        else: raise Exception('Invalid file type: {}'.format(ftype))
        # Load pNbody object for IC file
        nb0=pNbody.Nbody(p_name=fname,ftype=ftype,unitsfile=funit,flag_ic=True)
        # Get max ID for first galaxy
        idgal2=nb0.get_idgal2()
        # Clean up and return ID
        del nb0
        return idgal2

    def get_fext(self,fname=None,ftype=None,ptype=None,**extra_kw):
        """
        Returns the extension for a simulation snapshot
        """
        # Set constants
        fextDEF='_000ic'
        # Pars input
        fname=mmlpars.mml_pars(fname,type=str)
        ftype=mmlpars.mml_pars(ftype,type=str,default='gadget')
        ptype=mmlpars.mml_pars(ptype,type=str,default='disk')
        # Get file dictionary
        filedict=self.fdict
        # Determine base for this type of file
        ftype=ftype.lower()
        if   ftype=='gadget':
            fbase=filedict['gadget']['output']['snapbase']
        elif ftype=='scf':
            fbase=filedict['scf'][ptype]['output']['snapbase']
        elif ftype=='buildgal': 
            fbase=filedict['icgen']['bgic']
        else: raise Exception('Invalid file type: {}'.format(ftype))
        # Determine exception
        if ftype=='buildgal': fext=fextDEF
        else:
            if fbase in fname:
                fext=fname.split(fbase)[-1]
            else:
                fext=fextDEF
            if len(fext) > 10: fext=fextDEF
        # Return file extension
        return fext

    def get_ntot(self):
        """
        Returns the total number of particles in the simulation
        """
        if   self['runtyp']=='galsim': ntot=self.infodict['ntot']
        elif self['runtyp']=='intsim':
            ntot1=simfile.loadfpar('icgen','galsim',self.infodict['gal1'])['ntot']
            ntot2=simfile.loadfpar('icgen','galsim',self.infodict['gal2'])['ntot']
            ntot=ntot1+ntot2
        else: raise Exception('Invalid run type: {}'.format(self['runtyp']))
        return ntot

    def get_rundefaults(self):
        """
        Returns default information for simulation runs
        """
        sclfact=6.
        galfact=100.
        modList=['pos']
        typList=simlist.LIST_PTYPS
        galList=simlist.LIST_PGALS
        defdict={}
        # Get base run                                                                                                                                                          
        simstr0=self.get_primary(striptag=True)
        # Get profile parameters
        profpar=simstr0.get_profpar()
        # Plotting windows
        winddict={}
        for ityp in typList:
            if ityp in profpar: idef=2*(mmlmath.oom(sclfact*profpar[ityp]['rs']),)
            else              : idef=(0.,0.)
            winddict[ityp]=idef
#            if ityp in profpar: winddict[ityp]=mmlio.askquest('Enter valid window for {} {}.'.format(simstr0['runtag'],ityp),default=idef,dtype='tuple')
#            else              : winddict[ityp]=idef
        winddict['visi']=winddict['disk']
        winddict['star']=winddict['disk']
        winddict['gas' ]=winddict['disk']
        defdict['winddict']=winddict
        # Plotting limits
        limdict={}
        for imod in modList:
            limdict[imod]={}
            for ityp in typList:
                limdict[imod][ityp]={}
                for igal in galList:
                    if igal==0 and self['runtyp']=='intsim':
                        xwind=galfact*winddict[ityp][0]
                    else:
                        xwind=winddict[ityp][0]
                    ilim=(-xwind/2.,xwind/2.)
                    limdict[imod][ityp][igal]=ilim
        defdict['limdict']=limdict
        return defdict

    def get_rundefaultsOLD(self):
        """
        Returns default information for simulation runs
        """
        # Get base run                                                                                                                                                          
        simstr0=self.get_primary(striptag=True)
        defdict={}
        # Plotting windows
        halowindowdict={'MW_7M12LH11FG0':(150.,150.),'MW_6p52M11p52LH11FG0':(100.,100.),'MW_6M11LH11FG0':(50.,50.),'MW0B_7M12LH11FG0':(150.,150.),'OTHER':(0.,0.)}
        diskwindowdict={'MW_7M12LH11FG0':( 30., 30.),'MW_6p52M11p52LH11FG0':( 20., 20.),'MW_6M11LH11FG0':(10.,10.),'MW0B_7M12LH11FG0':( 30., 30.),'OTHER':(0.,0.)}
        bulgwindowdict={'MW_7M12LH11FG0':(  6.,  6.),'MW_6p52M11p52LH11FG0':(  4.,  4.),'MW_6M11LH11FG0':( 2., 2.),'MW0B_7M12LH11FG0':(  6.,  6.),'OTHER':(0.,0.)}
        visiwindowdict=diskwindowdict
        starwindowdict=diskwindowdict
        gasdwindowdict=diskwindowdict
        winddict={'gas'    :gasdwindowdict[simstr0['runtag']],
                  'halo'   :halowindowdict[simstr0['runtag']],
                  'disk'   :diskwindowdict[simstr0['runtag']],
                  'bulge'  :bulgwindowdict[simstr0['runtag']],
                  'visi'   :visiwindowdict[simstr0['runtag']],
                  'star'   :starwindowdict[simstr0['runtag']],
                  'all'    :halowindowdict[simstr0['runtag']]}
        defdict['winddict']=winddict
        modList=['pos']
        typList=winddict.keys()
        galList=[0,1,2]
        galfact=100.
        limdict={}
        for imod in modList:
            limdict[imod]={}
            for ityp in typList:
                limdict[imod][ityp]={}
                for igal in galList:
                    if igal==0 and self['runtyp']=='intsim':
                        xwind=galfact*winddict[ityp][0]
                    else:
                        xwind=winddict[ityp][0]
                    ilim=(-xwind/2.,xwind/2.)
                    limdict[imod][ityp][igal]=ilim
        defdict['limdict']=limdict
        # Return default dictionary
        return defdict

####################################################################################################################################
####################################################################################################################################
# WALK THROUGH METHOD
def walk(mtype=None,method=None,runtag=None,verbose=None,**extra_kw):
    """
    Method to walk users through performing different run methods
    """
    # Set constants
    simtypLIST=simlist.LIST_RUNTYPS
    # Get method type if not provided
    mtypeLIST=simlist.LIST_METHTYP
    if mtype not in mtypeLIST: mtype=mmlio.askselect('Select a type of method:',mtypeLIST)
    # Get method if not provided
    methodLIST=simlist.DICT_METHODS[mtype]
    if method not in methodLIST: method=mmlio.askselect('Select a {} method:'.format(mtype),methodLIST)
    method=method.lower()
    # Proceed based on method
    if mtype == 'icgen':
        icgen.run(method,runtag=runtag,**extra_kw)
    else:
        # Get run info
        mmlio.verbose('Getting info on the run...')
        simobj=runtag2simobj(runtag)
        # Call run
        mmlio.verbose('Beginning {} for {} simulation...'.format(method.upper(),simobj['runtyp'].upper()),border=True)
        output=simobj.run(mtype,method,verbose=verbose,**extra_kw)
        print [method,simobj['runtyp'],simobj['runtag']]
        if output==None: return
        else           : return output

####################################################################################################################################
# GENERAL METHODS
def general(simstr,method,**method_kw):
    """
    Handles running different types of methods for runs
    """
    # Set constants
    typListDEF=simlist.LIST_PTYPBASE
    gallist=['both','primary','secondary']
    typListVisible=simlist.LIST_PTYPVISI
    # Pars input
    method=mmlpars.mml_pars(method,type=str)
    method=method.lower()
    # Initialize output
    out=None
    # Find correct method subclass
    # Return mmlsim object
    if   method == 'getsimobj': return simstr
    # Return pNbody object
    elif method == 'getpnbody':
        method_kw['askuser']=True
        out=simstr.get_pnbody(**method_kw)
    # Print info on simulation
    elif method == 'printinfo':
        M1_vir=simstr.get_Mvir(galaxy=1)
        R1_vir=simstr.get_Rvir(galaxy=1)
        print 'SIMSTR:'
        print '    M1_vir={} Msol'.format(M1_vir)
        print '    R1_vir={} kpc'.format(R1_vir)
        if simstr['runtyp']=='intsim':
            M2_vir=simstr.get_Mvir(galaxy=2)
            R2_vir=simstr.get_Rvir(galaxy=2)
            rperi=(R1_vir+R2_vir)*simstr.infodict['rperi']
            print '    M2_vir={} Msol'.format(M2_vir)
            print '    R2_vir={} kpc'.format(R2_vir)
            print '    Rperi={} kpc'.format(rperi)
        fext=mmlio.askquest('What snapshot extension should be loaded?',default='ic',dtype='str')
        pnObj=simstr.get_pnbody(ftype='gadget',fext=fext)
        pnObj.printinfo()
        out=pnObj
    # Check if particles are repeated
    elif method == 'checkrepeat':
        pnobj=general(simstr,'getpnbody')
        out=pnobj.countrepeats()
        print 'N repeats = {}'.format(out)
    # Get profile info
    elif method=='fitprofiles':
        pnobj=general(simstr,'getpnbody')
        out={}
        for ityp in pnobj.get_typlist():
            out[ityp]=pnobj.fitprofile(ptyp=ityp,**method_kw)
            print '{} prof: {}'.format(ityp,out[ityp])
    # Plot pNbody object
    elif method == 'displayplot':
        fext=mmlio.askquest('What snapshot extension should be loaded?',default='ic',dtype='str')
        pnObj=simstr.get_pnbody(ftype='gadget',fext=fext)
        # Galaxy info
        if simstr['simtyp']=='intsim':
            galstr=mmlio.askselect('Which galaxy should be plotted?',gallist)
            galaxy=gallist.index(galstr)
        else:
            galaxy=0
        # Type info
        ptype=mmlio.askselect('What type of particle should be plotted?',['all','visible']+typListDEF)
        if   ptype=='all'    : typList=typListDEF
        elif ptype=='visible': typList=typListVisible
        else                 : typList=[ptype]
        # View window
        view_plane=mmlio.askselect('What plane should be plotted?',['xy','xz'])
        view_window=mmlio.askquest('What should the dimensions of the view be?',default=(0.,0.),dtype='tuple')
        # Center stuff
        center=mmlio.yorn('Should the particles be recentered?')
        centermeth=None
        centertyp=None
        centergal=None
        if center:
            centerptype=mmlio.askselect('What type of particles should the image be centered on?',['all','visible']+typListDEF)
            if   centerptype=='all'     : centertyp=typListDEF
            elif centerptype=='visible' : centertyp=typListVisible
            else                        : centertyp=[centerptype]
            if simstr['simtyp']=='intsim':
                centergalstr=mmlio.askselect('Which galaxy should the image be centered on?',gallist)
                centergal=gallist.index(centergalstr)
            else:
                centergal=0
            centermeth=mmlio.askselect('How should the particles be centered?',mml_pNbody.list_centermethods())
        # File info
        overwrite=None
        fname=None
        if mmlio.yorn('Should the plot be saved?'):
            fname=mmlio.askquest('Where should the plot be saved?',default='{}.displayplot_{}'.format(simstr['runtag'],ptype),dtype='str')
            if os.path.isfile(fname): overwrite=mmlio.yorn('That file already exists. Overwrite?')
            outimage=False
        else:
            outimage=True
        # Plot
        out=simplot.plot_sing(pnObj,fname=fname,overwrite=overwrite,outimage=outimage,
                              typList=typList,galaxy=galaxy,
                              view_plane=view_plane,view_window=view_window,view_palette='light',
                              center=center,centermeth=centermeth,centertyp=centertyp,centergal=centergal)
        # Return output
        if outimage:
            out.show()
            return out
        else:
            print fname
            return pnObj
    # Plot scatter of pNbody object
    elif method == 'scatterplot':
        fext=mmlio.askquest('What snapshot extension should be loaded?',default='ic',dtype='str')
        pnObj=simstr.get_pnbody(ftype='gadget',fext=fext)
        # Galaxy info
        if simstr['runtyp']=='intsim':
            galstr=mmlio.askselect('Which galaxy should be plotted?',gallist)
            pgal=gallist.index(galstr)
        else:
            pgal=0
        # Type info
        ptyp=mmlio.askselect('What type of particle should be plotted?',['all','visi']+typListDEF)
        # Color
        cmethLIST=['none','type','galaxy']
        cmeth=mmlio.askselect('What method should be used to color the scatter points?',cmethLIST)
        # View window
        rmax=mmlio.askquest('How big should the plotted region be?',default=0.,dtype='float')
        # View orientation
        zdir=mmlio.askselect('What cartesian coordinate should be up?',['z','y','x'])
        # Center stuff
        center=mmlio.yorn('Should the particles be recentered?')
        centermeth=None
        centertyp=None
        centergal=None
        if center:
            centertyp=mmlio.askselect('What type of particles should the image be centered on?',['all','visi']+typListDEF)
            if simstr['runtyp']=='intsim':
                centergalstr=mmlio.askselect('Which galaxy should the image be centered on?',gallist)
                centergal=gallist.index(centergalstr)
            else:
                centergal=0
            centermeth=mmlio.askselect('How should the particles be centered?',mml_pNbody.list_centermethods())
            oldcm=pnObj.get_center(centergal=centergal,centertyp=centertyp,method=centermeth,move=center)
        # File info
        overwrite=None
        fname=None
        if mmlio.yorn('Should the plot be saved?'):
            rmaxstr=mmlstring.dec2str(rmax)
            fname=mmlio.askquest('Where should the plot be saved?',default='{}.scatterplot_{}_{}_{}'.format(simstr['runtag'],ptyp,pgal,rmaxstr),dtype='str')
            if os.path.isfile(fname): overwrite=mmlio.yorn('That file already exists. Overwrite?')
            outplot=False
        else:
            outplot=True
        # Plot
        out=pnObj.plotscat(plotfile=fname,overwrite=overwrite,outplot=outplot,zdir=zdir,
                           ptyp=ptyp,pgal=pgal,rmax=rmax,colormeth=cmeth,askuser=True)
        # Return output
        if outplot:
            out.show()
            return out
        else:
            print fname
            return pnObj
    else: raise Exception('Invalid general method: {}'.format(method))
    # Return output
    return out

if __name__ == '__main__':
    walk()

        
