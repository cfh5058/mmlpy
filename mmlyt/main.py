#!/usr/bin/python
import sys,os,shutil,glob,copy,pprint,scipy,math
import numpy as np
from mmlutils import *
from mmlastro import mmlprofiles
import simlist,simfile,mmlgadget

####################################################################################################################################
# METHOD TO GET MMLSIM OBJECT FROM RUNTAG
def runtag2simobj(runtag=None):
    """
    Returns the mmlsim object for a given runtag
    """
    runlist=mmlparam.par_taglist('mmlyt.simlist','mmlsim')
    if runtag in runlist:
        simdict=mmlparam.loadpar('mmlyt.simlist','mmlsim',runtag)
        ifodict=mmlparam.loadpar('mmlyt.simlist',simdict['runtyp'],runtag)
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
    ifodict=mmlparam.parspar('mmlyt.simlist',runtyp,inpar=infodict,askuser=True)
    simdict['runtyp']=runtyp
    simdict['runtag']=mmlparam.par2tag('mmlyt.simlist',runtyp,ifodict)
    simdict=mmlparam.parspar('mmlyt.simlist','mmlsim',inpar=simdict,askuser=True)
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
        self.update(mmlparam.parspar('mmlyt.simlist','mmlsim',inpar=inkw))
        # Get infodict
        infodict=mmlparam.parspar('mmlyt.simlist',inkw['runtyp'],inpar=infodict)
        setattr(self,'infodict',infodict)
        # Add file dictionary
        if self['runtag'] not in simlist.LIST_RUNTYPS:
            fdict,finfo,ftab=self.mkfiledict(retall=True)
            setattr(self,'fdict',fdict)
            setattr(self,'finfo',finfo)
            setattr(self,'ftab' ,ftab )
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
        elif mtype=='gadget'  : out=mmlgadget.run(self,method,**method_kw)
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

    def get_Mvir(self,galaxy=None):
        """
        Returns the virial mass for a given simulation
        """
        if self['runtyp']=='galsim': 
            Mvir=self.infodict['mvir']
        else:
            galaxy=mmlpars.mml_pars(galaxy,default=1,list=[1,2])
            Mvir=mmlparam.loadpar('mmlyt.simlist','galsim',self.infodict['gal{}'.format(galaxy)])['mvir']
        return Mvir

    def get_Rvir(self,galaxy=None,**exkw):
        """
        Returns the virial radius for a given simulations (in pc)
        """
        Mvir=self.get_Mvir(galaxy=galaxy)
        exkw['units']='galaxy'
        Rvir=mmlprofiles.mvir2rvir(Mvir,**exkw)
        return Rvir

    def get_flagic(self,fext):
        """
        Determines if an extension is for IC file
        """
        fext=mmlpars.mml_pars(fext,type=str)
        if fext in ['ic','_ic','_000ic']: flag_ic=True
        else                            : flag_ic=False
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
        noidgal2=mmlpars.mml_pars(noidgal2,default=False,type=bool) ; noidgal2=True
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

    def get_primary(self,striptag=False):
        """
        Returns the mmlsim object for the primary galaxy
        """
        # Select correct sim dict for isolated primary
        if   self['runtyp']=='galsim': runtag=self['runtag']
        elif self['runtyp']=='intsim': runtag=self.infodict['gal1']
        # Get info from runtag
        galobj=mmlparam.loadpar('mmlyt.simlist','galsim',runtag)
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
#        nb0=pNbody.Nbody(p_name=fname,ftype=ftype,unitsfile=funit,flag_ic=True)
#        # Get max ID for first galaxy
#        idgal2=nb0.get_idgal2()
#        # Clean up and return ID
#        del nb0
#        return idgal2

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
            ntot1=mmlparam.loadpar('mmlyt.simlist','galsim',self.infodict['gal1'])['ntot']
            ntot2=mmlparam.loadpar('mmlyt.simlist','galsim',self.infodict['gal2'])['ntot']
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
    mtypeLIST=simlist.LIST_METHTYPS
    if mtype not in mtypeLIST: mtype=mmlio.askselect('Select a type of method:',mtypeLIST)
    # Get method if not provided
    methodLIST=simlist.DICT_METHODS[mtype]
    if method not in methodLIST: method=mmlio.askselect('Select a {} method:'.format(mtype),methodLIST)
    method=method.lower()
    # Proceed based on method
    if   mtype=='buildgal' and method=='mk_ic': mmlbuildgal.mk_ic(runtag=runtag,**extra_kw)
    elif mtype=='gadget'   and method=='mk_ic': mmlgadget.mk_ic(runtag=runtag,**extra_kw)
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
#        out=simstr.get_pnbody(**method_kw)
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
#        pnObj=simstr.get_pnbody(ftype='gadget',fext=fext)
#        pnObj.printinfo()
#        out=pnObj
    # Check if particles are repeated
    elif method == 'checkrepeat':
#        pnobj=general(simstr,'getpnbody')
#        out=pnobj.countrepeats()
        print 'N repeats = {}'.format(out)
    # Get profile info
    elif method=='fitprofiles':
        pnobj=general(simstr,'getpnbody')
        out={}
        for ityp in pnobj.get_typlist():
            out[ityp]=pnobj.fitprofile(ptyp=ityp,**method_kw)
            print '{} prof: {}'.format(ityp,out[ityp])
    else: raise Exception('Invalid general method: {}'.format(method))
    # Return output
    return out

if __name__ == '__main__':
    walk()

        
