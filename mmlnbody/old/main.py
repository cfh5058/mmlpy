#!/usr/bin/python
import sys,os,shutil,glob,copy,pprint,scipy,math
import astropysics.models as apymods
import numpy as np
from mmlutils import *
import simlist,simlyze,simstat,simcalc,simhist,simplot,simfile,mmlbuildgal,mmlgadget,mmlscf,mmlgalfit,icgen,runsheet,compare,mml_pNbody
import pNbody

####################################################################################################################################
####################################################################################################################################
# MMLSIM CLASS
class mmlsim(mmlclass.mydict):
    """
    A dictionary subclass with keys:
        runtag:      String uniquely identifying the simulation.
        runtyp:      String specifying what type of simulation it is.
        addtag:      String containing additional tag to append to file names.
        memcomp:     String specifying what vpac machine to use for storage on vpac network.
        nprocbg:     Int specifying how many processors BUILDGAL is run on.
        nprocgd:     Int specifying how many processors GADGET is run on.
        nprocscf:    Int specifying how many processors SCF is run on.
        inclbgtag:   Bool specifying if # of BUILDGAL processors tag should be added to files.
        inclgdtag:   Bool specifying if # of GADGET processors tag should be added to files.
        inclscftag:  Bool specifying if # of SCF processors tag should be added to files.
        newstyletag: Bool specifying if tags should be generated using the new style
        ntot:        Total number of particles in the simulation
    See mmlsim.get_form() dictionary for additional details.
    """

    def __init__(self,runtag=None,runtyp=None,addtag=None,
                 memcomp=None,nprocbg=None,nprocgd=None,nprocscf=None,
                 inclbgtag=None,inclgdtag=None,inclscftag=None,newstyletag=None,
                 infodict=None,ntot=None):
        """
        Initializes mmlsim objects
        """
        mmlclass.mydict.__init__(self,runtag=runtag,runtyp=runtyp,addtag=addtag,
                                 memcomp=memcomp,nprocbg=nprocbg,nprocgd=nprocgd,nprocscf=nprocscf,
                                 inclbgtag=inclbgtag,inclgdtag=inclgdtag,inclscftag=inclscftag,
                                 newstyletag=newstyletag,infodict=infodict,ntot=ntot)
        self.pars()
        if self['runtag'] not in simlist.LIST_RUNTYPS:
            fdict=self.mkfiledict()
            setattr(self,'fdict',fdict)

    def get_form(self):
        """
        Returns a dictionary that defines mmlsim key/value pairs.
        """
        if isinstance(self['runtyp'],str): self['runtyp']=self['runtyp'].lower()
        self['runtyp']=mmlpars.mml_pars(self['runtyp'],default='galaxy',type=str,list=simlist.LIST_RUNTYPS)
        self['runtag']=mmlpars.mml_pars(self['runtag'],default='galaxy',type=str)
        # Default info dict stuff
        if   self['runtyp'] == 'galaxy'  : infodictDEF=galsim()
        elif self['runtyp'] == 'interact': infodictDEF=intsim()
        infodictTYP=type(infodictDEF)
        # Create form
        form={
            'runtag'     : mmlpars.parsdict(default='galaxy'   ,type=str                             ),
            'runtyp'     : mmlpars.parsdict(default='galaxy'   ,type=str  ,list=simlist.LIST_RUNTYPS ),
            'ntot'       : mmlpars.parsdict(default=1L         ,type=long ,min=1L                    ),
            'addtag'     : mmlpars.parsdict(default=''         ,type=str                             ),
            'memcomp'    : mmlpars.parsdict(default='vpac38'   ,type=str                             ),
            'nprocbg'    : mmlpars.parsdict(default=1          ,type=int  ,min=1                     ),
            'nprocgd'    : mmlpars.parsdict(default=32         ,type=int  ,min=1                     ),
            'nprocscf'   : mmlpars.parsdict(default=50         ,type=int  ,min=1                     ),
            'inclbgtag'  : mmlpars.parsdict(default=False      ,type=bool                            ),
            'inclgdtag'  : mmlpars.parsdict(default=False      ,type=bool                            ),
            'inclscftag' : mmlpars.parsdict(default=False      ,type=bool                            ),
            'newstyletag': mmlpars.parsdict(default=False      ,type=bool                            ),
            'infodict'   : mmlpars.parsdict(default=infodictDEF,type=infodictTYP                     )
            }
        return form

    def pars(self):
        """
        Pars mmlsim objects
        """
        self=mmlpars.mml_formpars(self,self.get_form())
        return

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
        if self['runtyp']=='galaxy': 
            Mvir=self['infodict']['mvir']
        else:
            galaxy=mmlpars.mml_pars(galaxy,default=1,list=[1,2])
            Mvir=self['infodict']['simstr{}'.format(galaxy)]['infodict']['mvir']
        return Mvir

    def get_Rvir(self,galaxy=None):
        """
        Returns the virial radius for a given simulations
        """
        Mvir=self.get_Mvir(galaxy=galaxy)
        Rvir=apymods.NFWModel().Mvir_to_Rvir(Mvir)
        return Rvir

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
        ftype=mmlpars.mml_pars(ftype,default='gadget',type=str,list=['gadget','scf'])
        ptype=mmlpars.mml_pars(ptype,default='disk',type=str)
        fext=mmlpars.mml_pars(fext,type=str)
        # Get files
        files=self.fdict
        # Determine if extension is icfile
        flag_ic=self.get_flagic(fext)
        # Get info on file
        if flag_ic:
            if ftype=='scf':
                fname=files['scf'][ptype]['output']['snapbase']+fext
            elif ftype=='gadget':
                fname=files['gadget']['input']['ic']
        else:
            if ftype=='scf':
                fname=files['scf'][ptype]['output']['snapbase']+fext
            elif ftype=='gadget':
                fname=files['gadget']['output']['snapbase']+fext
        # Return file name and IC flag
        return fname,flag_ic

    def get_loadkeys(self,ftype=None,ptype=None,fname=None,fext=None,idgal2=None,
                     snapdict=None,askuser=None,noidgal2=None,outextra=None,**extra_kw):
        """
        Returns a dictionary of keys for loading snapshot data
        """
        # Set constants
        fextDEF='_000ic'
        # Pars input
        ftype=mmlpars.mml_pars(ftype,default='gadget',type=str,list=['gadget','scf'])
        ptype=mmlpars.mml_pars(ptype,default='disk',type=str)
        askuser=mmlpars.mml_pars(askuser,default=False,type=bool)
        noidgal2=mmlpars.mml_pars(noidgal2,default=False,type=bool)
        outextra=mmlpars.mml_pars(outextra,default=False,type=bool)
        # Get info from pNbody object
        if snapdict!=None:
            fname=snapdict.p_name[0]
            ftype=snapdict.ftype.split('Nbody_')[-1]
            if ftype=='scf': ptype=snapdict.ptype
        # Handle absence of fext or fname
        if not isinstance(fext,str) and not isinstance(fname,str):
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

    def get_softenings(self):
        """
        Loads softening lengths for a simulation
        """
        fdict=self.fdict
        parmfile=fdict['gadget']['input']['param']
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
        if   self['runtyp'].lower()=='galaxy'  : runtag=self['runtag']
        elif self['runtyp'].lower()=='interact': runtag=self['infodict']['simstr1']['runtag']
        # Get runsheet info from runtag
        runsht=runsheet.get_runsheet(runtag=runtag)
        addruntag=runsht['addruntag']
        # Strip additional tags as needed to get to base run
        if striptag:
            # Runs without additional tags
            if len(addruntag) == 0:
                runtag0=runtag
                runsht0=runsht
            # Runs with additional tags
            else:
                runtag0=runtag.split(addruntag)[0]
                runsht0=runsheet.get_runsheet(runtag=runtag0)
        else:
            runtag0=runtag
            runsht0=runsht
        # Return simulation dictionary for primary
        simstr0=runsheet.runsheet2mmlsim(**runsht0)
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
        else: raise Exception('Invalid file type: {}'.format(ftype))
        # Determine exception
        if fbase in fname:
            fext=fname.split(fbase)[-1]
        else:
            fext=fextDEF
        if len(fext) > 10: fext=fextDEF
        # Return file extension
        return fext

    def get_rundefaults(self):
        """
        Returns default information for simulation runs
        """
        # Get base runsheet
        simstr0=self.get_primary(striptag=True)
        shorttag=runsheet.runtag2runsheet(simstr0['runtag'])
##         # Full mass MW
##         if   'MW_7M12LH11FG0'       in simstr0['runtag']: shorttag='MW1'
##         # Third mass MW
##         elif 'MW_6p52M11p52LH11FG0' in simstr0['runtag']: shorttag='MW3'
##         # Tenth mass MW
##         elif 'MW_6M11LH11FG0'       in simstr0['runtag']: shorttag='MW10'
##         # Other
##         else:
##             raise Exception('Option needs to be added to mmlsim.get_rundefaults for a runtag of {}'.format(simstr0['runtag']))
#            defdict={'addtag':''    ,'nprocgd':128,'nprocbg':100,'memcomp':'vpac38'}
#            if simstr0['runtag']==self['runtag']: defdict['nprocgd']=256
#            defdict={'addtag':'_BG1','nprocgd':64 ,'nprocbg':1  ,'memcomp':'vpac38'}
#            defdict={'addtag':'_BG1','nprocgd':32 ,'nprocbg':1  ,'memcomp':'vpac38'}
##             shorttag='OTHER'
##             defdict={'addtag':''    ,'nprocgd':32 ,'nprocbg':1  ,'memcomp':'vpac38'}
        defdict={}
        # Plotting windows
        halowindowdict={'MW1':(150.,150.),'MW3':(100.,100.),'MW10':(50.,50.),'OTHER':(0.,0.)}
        diskwindowdict={'MW1':( 30., 30.),'MW3':( 20., 20.),'MW10':(10.,10.),'OTHER':(0.,0.)}
        bulgwindowdict={'MW1':(  6.,  6.),'MW3':(  4.,  4.),'MW10':( 2., 2.),'OTHER':(0.,0.)}
        visiwindowdict=diskwindowdict
        starwindowdict=diskwindowdict
        gasdwindowdict=diskwindowdict
        winddict={'gas'    :gasdwindowdict[shorttag],
                  'halo'   :halowindowdict[shorttag],
                  'disk'   :diskwindowdict[shorttag],
                  'bulge'  :bulgwindowdict[shorttag],
                  'visi'   :visiwindowdict[shorttag],
                  'star'   :starwindowdict[shorttag],
                  'all'    :halowindowdict[shorttag]}
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
                    if igal==0 and self['runtyp']=='interact':
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
# GALSIM CLASS
class galsim(mmlclass.mydict):
    """
    Dictionary for passing around info on a galaxy simulation with keys:
        model:    String specifying what model the galaxy is based on.
        haloprof: String specifying what halo profile the galaxy has.
        ntot:     Long sepcifying total number of particles in the sim.    
        mvir:     Float specifying the virial mass of the galaxy's DM halo.
        fgas:     Fraction of gas in the galaxy's disk by mass.
        concen:   Concentration of the galaxy's DM halo.
    """

    def __init__(self,model=None,haloprof=None,ntot=None,mvir=None,fgas=None,concen=None):
        """
        Initializes a new instance of a galsim object
        """
        mmlclass.mydict.__init__(self,model=model,haloprof=haloprof,ntot=ntot,mvir=mvir,fgas=fgas,concen=concen)
        self.pars()

    def get_form(self):
        """
        Returns a dictionary that defines galsim key/value pairs
        """
        modelLIST=mmlbuildgal.LIST_GALMODELS
        hprofLIST=mmlbuildgal.LIST_HALOPROFS
        form={
            'model'    : mmlpars.parsdict(default='MW',type=str  ,list=modelLIST       ),
            'haloprof' : mmlpars.parsdict(default='LH',type=str  ,list=hprofLIST       ),
            'ntot'     : mmlpars.parsdict(default=1L  ,type=long ,min=1L               ),
            'mvir'     : mmlpars.parsdict(default=0.0 ,type=float,min=0.0              ),
            'fgas'     : mmlpars.parsdict(default=0.0 ,type=float,range=[0.0,1.0]      ),
            'concen'   : mmlpars.parsdict(default=11.0,type=float,min=0.0              ),
            'addruntag': mmlpars.parsdict(default=''  ,type=str                        )
            }
        return form
    
    def pars(self):
        """
        Parses galsim objects using get_form info
        """
        self=mmlpars.mml_formpars(self,self.get_form())

    def mkruntag(self,newstyletag=None):
        """
        Creates a runtag for a galaxy simulation from a galsim object
        """
        # Get runtag & return it
        runtag=simfile.get_runtag('galaxy',newstyletag=newstyletag,**self)
        return runtag

    def mkmmlsim(self,addtag=None,memcomp=None,nprocbg=None,nprocgd=None,nprocscf=None,
                 inclbgtag=None,inclgdtag=None,inclscftag=None,newstyletag=None):
        """
        Creates an mmlsim object from a galsim object
        """
        runtag=self.mkruntag(newstyletag=newstyletag)
        runtyp='galaxy'
        return mmlsim(runtag=runtag,runtyp=runtyp,ntot=self['ntot'],
                      addtag=addtag,memcomp=memcomp,nprocbg=nprocbg,nprocgd=nprocgd,nprocscf=nprocscf,
                      inclbgtag=inclbgtag,inclgdtag=inclgdtag,inclscftag=inclscftag,newstyletag=newstyletag,
                      infodict=self)

    

####################################################################################################################################
####################################################################################################################################
# INTSIM CLASS
class intsim(mmlclass.mydict):
    """
    Dictionary for passing around info on an nbody interation simulation with keys:
        runtag1: String uniquely identifying primary galaxy sim.
        runtag2: String uniquely identifying secondary galaxy sim.
        ntot:    Long total number of particles in the simulation.
        rperi:   Float value of pericenter separation in terms of rvir1+rvir2.
        ecc:     Float value of interaction's orbit eccentricity.
        incl:    Float value of interaction's orbit inclination in radians.
    """

    def __init__(self,galsim1=None,galsim2=None,runtag1=None,runtag2=None,ntot=None,rperi=None,ecc=None,incl=None):
        """
        Initializes a new instance of an intsim object
        """
        # Handle galsim input first
        if galsim1!=None and galsim2!=None:
            runtag1=galsim1.mkruntag()
            runtag2=galsim2.mkruntag()
            simstr1=galsim1.mkmmlsim()
            simstr2=galsim2.mkmmlsim()
            ntot=galsim1['ntot']+galsim2['ntot']
        elif isinstance(runtag1,str) and isinstance(runtag2,str):
            runsheet1=runsheet.runtag2runsheet(runtag1)
            runsheet2=runsheet.runtag2runsheet(runtag2)
            simstr1=runsheet.runsheet2mmlsim(runsheet=runsheet1)
            simstr2=runsheet.runsheet2mmlsim(runsheet=runsheet2)
        else:
            mmlio.verbose('WARNING: the returned intsim does not have any info on the interacting galaxies.')
            simstr1=None
            simstr2=None
        # Pars
        mmlclass.mydict.__init__(self,simstr1=simstr1,simstr2=simstr2,runtag1=runtag1,runtag2=runtag2,
                                 ntot=ntot,rperi=rperi,ecc=ecc,incl=incl)
        self.pars()

    def get_form(self):
        """
        Returns a dictionary defining intsim key/value pairs
        """
        simstrDEF=mmlsim() ; simstrTYP=type(simstrDEF)
        form={
            'simstr1': mmlpars.parsdict(default=simstrDEF,type=simstrTYP                   ),
            'simstr2': mmlpars.parsdict(default=simstrDEF,type=simstrTYP                   ),
            'runtag1': mmlpars.parsdict(default='galaxy1',type=str                         ),
            'runtag2': mmlpars.parsdict(default='galaxy2',type=str                         ),
            'ntot'   : mmlpars.parsdict(default=1L       ,type=long ,min=1L                ),
            'rperi'  : mmlpars.parsdict(default=1.5      ,type=float,min=0.0               ),
            'ecc'    : mmlpars.parsdict(default=1.5      ,type=float,min=0.0               ),
            'incl'   : mmlpars.parsdict(default=0.0      ,type=float,range=[0.0,2*np.pi]   )
            }
        return form

    def pars(self):
        """
        Parses intsim objects using get_form info
        """
        self=mmlpars.mml_formpars(self,self.get_form())

##     def getintsuite(self,type='flyby'):
##         '''
##         NAME:
##             mmlnbody.main.intsim.getinsuite
##         PURPOSE:
##             To return a dictionary defining interaction simulation suites
##         CALLING:
##             suite=intsim.getinsuite(type="flyby")
##         KEYWORDS:
##             TYPE:  String identifying what type of interaction suite to return.
##                    ["flyby" is currently the only supported type]
##         OUTPUT:
##             SUITE: Dictionary defining different suites and their corresponding intsim objects.
##         '''
##         if type is 'flyby':
##             suite={}
##             # EXTREME
##             extr_elist=[1.0,1.2,1.4,1.5]
##             suite['extreme']=[]
##             for ie in extr_elist:
##                 iselfcopy=copy.deepcopy(self)
##                 iselfcopy['rperi']=0.1
##                 iselfcopy['ecc']=ie
##                 suite['extreme'].append(iselfcopy)
##             # PARABOLIC
##             para_rlist=[0.5,1.0,1.25]
##             suite['parabolic']=[]
##             for ir in para_rlist:
##                 iselfcopy=copy.deepcopy(self)
##                 iselfcopy['rperi']=ir
##                 iselfcopy['ecc']=1.0
##                 suite['parabolic'].append(iselfcopy)
##             # VARY RPERI
##             varr_rlist=[0.25,0.5,0.75,1.0,1.25]
##             suite['varyrperi']=[]
##             for ir in varr_rlist:
##                 iselfcopy=copy.deepcopy(self)
##                 iselfcopy['rperi']=ir
##                 suite['varyrperi'].append(iselfcopy)
##             # VARY ECC
##             vare_elist=[1.5,2.0,3.0,5.0,7.0]
##             suite['varyecc']=[]
##             for ie in vare_elist:
##                 iselfcopy=copy.deepcopy(self)
##                 iselfcopy['ecc']=ie
##                 suite['varyecc'].append(iselfcopy)
##             # RETROGRADE (VARY RPERI & ECC)
##             suite['retrograde']=[]
##             for irself in suite['varyrperi']:
##                 iselfcopy=copy.deepcopy(irself)
##                 iselfcopy['incl']=180.0
##                 suite['retrograde'].append(iselfcopy)
##             for ieself in suite['varyecc']:
##                 iselfcopy=copy.deepcopy(ieself)
##                 iselfcopy['incl']=180.0
##                 suite['retrograde'].append(iselfcopy)
##             # VARRY INCL
##             vari_ilist=[45.0,90.0]
##             suite['varyincl']=[]
##             for ii in vari_ilist:
##                 iselfcopy=copy.deepcopy(self)
##                 iselfcopy['incl']=ii
##                 suite['varyincl'].append(iselfcopy)
##         else:
##             raise Exception('"flyby" is currently the only supported value for type keyword')

##         return suite

    def mkruntag(self,newstyletag=None):
        """
        Creates a runtag for an interaction simulation from an intsim object
        """
        # Recover runsheets
        self['run1']=runsheet.runtag2runsheet(self['runtag1'])
        self['run2']=runsheet.runtag2runsheet(self['runtag2'])
        # Create and return runtag
        runtag=simfile.get_runtag('interact',newstyletag=newstyletag,**self)
        return runtag

    def mkmmlsim(self,addtag=None,memcomp=None,nprocbg=None,nprocgd=None,nprocscf=None,
                 inclbgtag=None,inclgdtag=None,inclscftag=None,newstyletag=None):
        """
        Generates a mmlsim object from an intsim object
        """
        runtag=self.mkruntag()
        runtyp='interact'
        return mmlsim(runtag=runtag,runtyp=runtyp,ntot=self['ntot'],
                      addtag=addtag,memcomp=memcomp,nprocbg=nprocbg,nprocgd=nprocgd,nprocscf=nprocscf,
                      inclbgtag=inclbgtag,inclgdtag=inclgdtag,inclscftag=inclscftag,newstyletag=newstyletag,
                      infodict=self)

####################################################################################################################################
# WALK THROUGH METHOD
def walk(mtype=None,method=None,simkey=None,runtag=None,verbose=None,**extra_kw):
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
        icgen.run(method,**extra_kw)
    else:
        # Get run info
        mmlio.verbose('Getting info on the run...')
        rapsheet=runsheet.get_runsheet(runsheet=simkey,runtag=runtag)
        # Create mmlsim object from runsheet
        simobj=runsheet.runsheet2mmlsim(**rapsheet)
        # Call run
        mmlio.verbose('Beginning {} for {} simulation...'.format(method.upper(),rapsheet['simtyp'].upper()),border=True)
        output=simobj.run(mtype,method,verbose=verbose,**extra_kw)
        print [method,rapsheet['simtyp'],rapsheet['runsheet'],simobj['runtag']]
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
    if   method == 'getsimstr': return simstr
    # Return pNbody object
    elif method == 'getpnbody':
        ftype=mmlio.askselect('What type of snapshot should be loaded?',['gadget','scf'])
        ptype=None
        if ftype=='scf': ptype=mmlio.askselect('What particle type should be loaded?',simlist.LIST_PTYPS)
        fext=mmlio.askquest('What snapshot extension should be loaded?',default='ic',dtype='str')
        out=simstr.get_pnbody(ftype=ftype,fext=fext,ptype=ptype)
    # Print info on simulation
    elif method == 'printinfo':
        M1_vir=simstr.get_Mvir(galaxy=1)
        R1_vir=simstr.get_Rvir(galaxy=1)
        print 'SIMSTR:'
        print '    M1_vir={} Msol'.format(M1_vir)
        print '    R1_vir={} kpc'.format(R1_vir)
        if simstr['runtyp']=='interact':
            M2_vir=simstr.get_Mvir(galaxy=2)
            R2_vir=simstr.get_Rvir(galaxy=2)
            rperi=(R1_vir+R2_vir)*simstr['infodict']['rperi']
            print '    M2_vir={} Msol'.format(M2_vir)
            print '    R2_vir={} kpc'.format(R2_vir)
            print '    Rperi={} kpc'.format(rperi)
        fext=mmlio.askquest('What snapshot extension should be loaded?',default='ic',dtype='str')
        pnObj=simstr.get_pnbody(ftype='gadget',fext=fext)
        pnObj.printinfo()
        out=pnObj
    # Plot pNbody object
    elif method == 'displayplot':
        fext=mmlio.askquest('What snapshot extension should be loaded?',default='ic',dtype='str')
        pnObj=simstr.get_pnbody(ftype='gadget',fext=fext)
        # Galaxy info
        if simstr['simtyp']=='interact':
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
        winddict=simstr.get_rundefaults()['winddict']
        view_plane=mmlio.askselect('What plane should be plotted?',['xy','xz'])
        view_window=mmlio.askquest('What should the dimensions of the view be?',default=winddict[ptype],dtype='tuple')
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
            if simstr['simtyp']=='interact':
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
        if simstr['runtyp']=='interact':
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
        winddict=simstr.get_rundefaults()['winddict']
        rmax=mmlio.askquest('How big should the plotted region be?',default=winddict[ptyp][0],dtype='float')
        # View orientation
        zdir=mmlio.askselect('What cartesian coordinate should be up?',['z','y','x'])
        # Center stuff
        center=mmlio.yorn('Should the particles be recentered?')
        centermeth=None
        centertyp=None
        centergal=None
        if center:
            centertyp=mmlio.askselect('What type of particles should the image be centered on?',['all','visi']+typListDEF)
            if simstr['runtyp']=='interact':
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

        
