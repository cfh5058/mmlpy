#!/usr/bin/python
import sys,os,shutil,glob,copy,pprint,scipy,math,re,pickle
import matplotlib.pyplot as plt
import numpy as np
from mmlutils import *
import simlist,files

LIST_BARKEYS=['fext','time','rmin','rmax','phi','amp']
LIST_BARPROPMETH={'Am2'    :LIST_BARKEYS,
                  'derAm2' :LIST_BARKEYS,
                  'ellipse':LIST_BARKEYS+['ramp','Am4','emax','barFlag'],
                  'scfAm2' :LIST_BARKEYS}
LIST_SINGPLOTMETH=['plot2d','hist2d','pilimg','map2d','ellipse']+['barprop_'+ip for ip in LIST_BARPROPMETH]
LIST_MULTPLOTMETH=['multipanel','gfpanel']
LIST_PLOTMETH=LIST_SINGPLOTMETH+LIST_MULTPLOTMETH

####################################################################################################################################
# METHOD TO GET MMLSIM OBJECT FROM RUNTAG
def allsim(method,runtyp='all',methkw={}):
    """
    Perform some method for all tags
    """
    # Get taglist
    if   runtyp=='all': taglist=mmlparam.par_taglist('pysim.simlist','mmlsim')
    elif runtyp=='gal': taglist=mmlparam.par_taglist('pysim.simlist','galsim')
    elif runtyp=='int': taglist=mmlparam.par_taglist('pysim.simlist','intsim')
    else: raise Exception('Invalid runtyp: {}'.format(runtyp))
    # Loop over runtags
    for tag in taglist:
        sim=loadsimobj(tag)
        sim.getattr(method)(**methkw)
    # Return
    return
    
def loadsimobj(runtag=None,**kwargs):
    """
    Returns the mmlsim object for a given runtag
    """
    runlist=mmlparam.par_taglist('pysim.simlist','mmlsim')
    if runtag in runlist:
        simdict=mmlparam.loadpar('pysim.simlist','mmlsim',runtag,init=True)
        ifodict=mmlparam.loadpar('pysim.simlist',simdict['runtyp'],runtag)
        simobj=mmlsim(infodict=ifodict,**simdict)
    else:
        simobj=asksimobj(runtag=runtag,**kwargs)
    return simobj
def asksimobj(infodict=None,runtyp=None,**simdict):
    """
    Returns an mmlsim object provided user input
    """
    if runtyp not in simlist.LIST_RUNTYPS:
        runtyp=mmlio.askselect('Select a valid run type:',simlist.LIST_RUNTYPS)
    tagstr=simdict.get('runtag',None)
    if infodict!=None: tagstr=None
    ifodict=mmlparam.parspar('pysim.simlist',runtyp,inpar=infodict,tagstr=tagstr,askuser=True)
    simdict['runtyp']=runtyp
    if tagstr: simdict['runtag']=tagstr
    else     : simdict['runtag']=mmlparam.par2tag('pysim.simlist',runtyp,ifodict)
    simdict=mmlparam.parspar('pysim.simlist','mmlsim',inpar=simdict,askuser=True)
    simobj=mmlsim(infodict=ifodict,**simdict)
    return simobj
def get_ntot(runtag=None,**kwargs):
    """
    Returns the number of particles in a simulation
    """
    simobj=loadsimobj(runtag=runtag,**kwargs)
    return simobj.get_ntot()
def ntot2galid(ntot):
    """
    Defines convention for maximum ID given a # of particles
    """
    galid=long(9*mmlmath.oom(float(ntot),nsig=1,method='CEIL'))
    return galid

####################################################################################################################################
####################################################################################################################################
# MMLSIM CLASS
class mmlsim(dict):
    """
    A dictionary subclass with simulation keys.
    """

    def __init__(self,infodict={},**inkw):
        """
        Initializes mmlsim objects
        """
        # Get mmlsim dict
        self.update(mmlparam.parspar('pysim.simlist','mmlsim',inpar=inkw))
        # Get infodict
        infodict=mmlparam.parspar('pysim.simlist',inkw['runtyp'],inpar=infodict)
        setattr(self,'infodict',infodict)
        # Add file dictionary
        if self['runtag'] not in simlist.LIST_RUNTYPS:
            # File lists
            flist_mem={} ; flist_run={}
            for ftype in files.LIST_FILETYPES:
                flist_mem[ftype]=files.get_flist(ftype,compid=self['memcomp'],**self)
                flist_run[ftype]=files.get_flist(ftype,compid=self['runcomp'],**self)
            # Storage files
            fdict_mem,finfo_mem,ftab_mem=files.get_fdict(retall=True,compid=self['memcomp'],**self)
            setattr(self,'fdict',fdict_mem)
            setattr(self,'finfo',finfo_mem)
            setattr(self,'ftab' ,ftab_mem )
            setattr(self,'flist',flist_mem)
            # Runtime files
            fdict_run,finfo_run,ftab_run=files.get_fdict(retall=True,compid=self['runcomp'],**self)
            setattr(self,'fdict_run',fdict_run)
            setattr(self,'finfo_run',finfo_run)
            setattr(self,'ftab_run' ,ftab_run )
            setattr(self,'flist_run',flist_run)
        
    ################################################################################################################################
    ################################################################################################################################
    # CLASS PROPERTIES

    # Number of galaxies
    @property
    def ngal(self):
        if not hasattr(self,'_ngal'):
            if   self['runtyp']=='galsim': self._ngal=1
            elif self['runtyp']=='intsim': self._ngal=2
            else: raise Exception('Unsure how many galaxies are in {}'.format(self['runtyp']))
        return self._ngal
    # Number of particles
    @property
    def ntot(self):
        if not hasattr(self,'_ntot'): self._ntot=self.get_ntot()
        return self._ntot
    # @ntot.setter
    # def ntot(self,val):
    #     if isinstance(val,(long,int,float)):
    #         self._ntot=long(val)
    #     else:
    #         raise Exception('Invalid ntot value')
    # Galaxy ids
    @property
    def galids(self):
        if not hasattr(self,'_galids'): self._galids=self.get_galids()
        return self._galids
    # Limits
    @property
    def limits(self):
        if not hasattr(self,'_limits'): self._limits=self.get_limits()
        return self._limits
    # Families
    @property
    def families(self):
        if not hasattr(self,'_families'):
            from files import buildgal
            # Get list of galaxy info
            if   self['runtyp']=='galsim': ifolist=[self.infodict]
            elif self['runtyp']=='intsim':
                ifolist=[]
                for gal in range(1,self.ngal+1):
                    ifolist.append(mmlparam.loadpar('pysim.simlist','galsim',self.infodict['gal{}'.format(gal)]))
            else: raise Exception('Unsupported runtyp: {}'.format(self['runtyp']))
            # Get list of famlies
            famlist=set(['star'])
            for ifo in ifolist:
                bgpar=buildgal.model2param(ifo)
                for typ in buildgal.LIST_PTYPS:
                    if bgpar[typ+'_mtot']!=0: famlist.add(typ)
            self._families=list(famlist)
        return self._families
    # Galaxies
    @property
    def galaxies(self):
        if not hasattr(self,'_galaxies'): self._galaxies=['galaxy{}'.format(gal) for gal in range(1,self.ngal+1)]
        return self._galaxies
    # Baserun
    @property
    def baserun(self):
        if not hasattr(self,'_baserun'):
            mkdict=files.rw_makelog('R',self.fdict['gadget']['icgen']['makelog'])
            self._baserun=mkdict['baserun']
        return self._baserun

    # INTSIM SPECIFIC PROPERTIES
    # Tperi
    @property
    def Tperi(self):
        if not hasattr(self,'_Tperi'):
            self._Tperi=self.get_Tperi(cmeth='pot',ftype='gadget',verbose=False,askuser=False)
        return self._Tperi
    # Rperi
    @property
    def Rperi(self): 
        if not hasattr(self,'_Rperi'):
            self._Rperi=self.get_Rperi(cmeth='pot',ftype='gadget',verbose=False,askuser=False)
        return self._Rperi

    ################################################################################################################################
    ################################################################################################################################
    # INFO PRINTING METHODS
    def show_galmod(self,galList=None,typList=None,parList=None):
        """Print info on galaxy models"""
        from pysim.files import buildgal
        # Pars input
        if galList is None:
            if   self['runtyp']=='galsim': galList=self['runtag']
            elif self['runtyp']=='intsim': galList=[self.infodict['gal1'],self.infodict['gal2']]
            else: raise Exception('Invalid method for runtyp: {}'.format(self['runtyp']))
        if typList is None: typList=simlist.LIST_PTYPS
        if parList is None: parList=['n','mp','mtot','rscl','zscl']
        # Loop over galaxies
        for g in galList:
            print g
            ix=loadsimobj(g)
            bgpar=buildgal.model2param(ix.infodict)
            # Loop over parameters
            for p in parList:
                print '    '+p
                # Loop over particles types
                for t in typList:
                    if t+'_n' not in bgpar: continue
                    if bgpar[t+'_n']==0: continue
                    print '        {:5s} = {}'.format(t,bgpar[t+'_'+p])
            # Print halo to disk ratio
            print '    Mhalo/Mdisk = {}'.format(bgpar['halo_mtot']/bgpar['disk_mtot'])
        # Print runtag
        print self['runtag']

    ################################################################################################################################
    ################################################################################################################################
    # SNAPSHOT INTERFACE METHOD
    def itersnap(self,method,ftype='gadget',fext=None,nproc=1,**kwargs):
        """
        Perform a method for all snapshots
        """
        # Get module
        module=getattr(self,method)
        # Setup
        args=module(mode='init',loadkw={'ftype':ftype},**kwargs)
        # Locate snapshots
        snapkey,snaptemp,snaplist=files.listsnap(self,ftype,fext=fext)
        # Create function for processing each snapshot
        loadkw=args.pop('loadkw')
        def fiter(ifname):
            # Get & test fext
            ifext=files.get_ext(ifname,snaptemp)
            if not args['testfext'](ifext): return None
            # Add loadkw
            iloadkw=copy.deepcopy(loadkw)
            iloadkw.update(fname=ifname,fext=ifext)
            # Run
            print ifext
            iargs=copy.deepcopy(args)
            iout=module(mode='run',loadkw=iloadkw,**iargs)
            return iout
        # Loop over snapshots
        out=mmlparallel.parmap(fiter,snaplist,nproc=nproc)
        # Post process
        fout=module(mode='end',loadkw=loadkw,**args)
        # Return
        return fout
            
    def allsnaps(self,method=None,ftype='gadget',allfam=None,allgal=None,fext=None,**kwargs):
        """
        Perform a method for all snapshots
        """
        # Pars input
        plotmethLIST=['plot','singplot','multplot']
        methLIST=plotmethLIST
        if method not in methLIST:
            method=mmlio.askselect('Pick a process to run:',methLIST)
        if ftype not in files.LIST_SNAPFORM:
            ftype=mmlio.askselect('Pick a file type:',files.LIST_SNAPFORM)
        if kwargs.get('pmeth',None) in LIST_MULTPLOTMETH: allfam=False ; allgal=False
        if not isinstance(allfam,bool): allfam=mmlio.yorn('Should all families be plotted?')
        if not isinstance(allgal,bool): allgal=mmlio.yorn('Should all galaxies be plotted?')
        # Locate snapshots
        snapkey,snaptemp,snaplist=files.listsnap(self,ftype,fext=fext)
        # Loop over files
        loadkw={}
        loadkw['ftype']=ftype
        kwargs['ftype']=ftype
        idxsortsnap=sorted(range(len(snaplist)),key=lambda x:snaplist[x])
        firstsnap=True
        for i in idxsortsnap:
            loadkw['fname']=snaplist[i]
            # loadkw['ftemp']=snaptemp
            # loadkw['fkey' ]=snapkey
            # Get file extension
            fext=files.get_ext(loadkw['fname'],snaptemp)
            if not fext: 
                if not mmlio.yorn(str(fext)): continue
                fext=snapkey
            loadkw['fext']=fext
            # 335,365
            # try:
            #     if int(float(fext))<335: continue
            # except:
            #     continue
            # print fext
            # Handle first snapshot (ask for user input)
            if firstsnap:
                firstsnap=False
                if   method==    'plot': outargs=self.plotsnap_old(firstsnap=True,**kwargs)
                elif method=='singplot': outargs=self.singplotsnap(firstsnap=True,**kwargs)
                elif method=='multplot': outargs=self.multplotsnap(firstsnap=True,**kwargs)
                else                   : outargs={}
                outargs.setdefault('overwrite',False)
                outargs.setdefault('askuser'  ,False)
                loadkw['pNbody']=outargs.get('pNbody',False)
                if 'snapObj' not in outargs: outargs['snapObj']=self.loadsnap(**loadkw)
                snapObj0=outargs['snapObj']
                # Select findutial galaxy & family
                if method in plotmethLIST and outargs['pmeth'] in LIST_SINGPLOTMETH:
                    fam0=outargs['inpar']['family']
                    gal0=outargs['inpar']['galaxy']
                else:
                    fam0=''
                    gal0=''
                # Get list of families & galaxies
                if allfam:
                    if pNbody: famlist=['']+snapObj0.families()
                    else     : famlist=['']+[f.name for f in snapObj0.families()]
                else         : famlist=[fam0]
                if allgal:
                    if pNbody: gallist=['']+snapObj0.galaxies()
                    else     : gallist=['']+[g.name for g in snapObj0.galaxies()]
                else         : gallist=[gal0]
            # Handle other snapshots
            else: outargs['snapObj']=None
            # Loop over particle families & galaxies
            for fam in famlist:
                for gal in gallist:
                    if method in plotmethLIST and outargs['pmeth'] in LIST_SINGPLOTMETH:
                        outargs['inpar']['family']=fam
                        outargs['inpar']['galaxy']=gal
                    # Get & check file name
                    if method in plotmethLIST:
                        ifile=mmlparam.par2file(self,snapObj0,outargs['pmeth'],plot=True,inpar=outargs['inpar'],fext=fext)['file']
                    if os.path.isfile(ifile) and not outargs['overwrite']: continue
                    # Get snapshot object
                    if not outargs['snapObj']: outargs['snapObj']=self.loadsnap(**loadkw)
                    # Plot
                    try:
                        if method in plotmethLIST:
                            if   outargs['pmeth'] in LIST_SINGPLOTMETH: ifile=self.singplotsnap(fname=ifile,**outargs)
                            elif outargs['pmeth'] in LIST_MULTPLOTMETH: ifile=self.multplotsnap(fname=ifile,**outargs)
                            else: raise Exception('Invalid plot meth: {}'.format(outargs['pmeth']))
                        else: raise Exception('Invalid method: {}'.format(method))
                    except:
                        print outargs
                        print ifile
                        raise
                    print '    '+ifile
        # Animate
        if method in plotmethLIST:
            owanim=mmlio.yorn('Overwrite existing animations?')
            for fam in famlist:
                for gal in gallist:
                    outargs['inpar']['family']=fam
                    outargs['inpar']['galaxy']=gal
                    self.animsnap(snapObj0,outargs['pmeth'],inpar=outargs['inpar'],overwrite=owanim,
                                  fext=files.ftype2snapext(ftype))
                                    

    def loadsnap(self,fext=0,ftype=None,ptype='disk',pgal=1,fname=None,pNbody=False,**kwargs):
        """
        Returns an nbody object for the specified snapshot and type
        """
        # Get file info
        if 'ftemp' in kwargs: raise Exception('ftemp keyword no longer supported')
        if 'fkey'  in kwargs: raise Exception('fkey keyword no longer supported')
        fopt=files.get_snapopt(self,ftype=ftype,fname=fname,fext=fext,ptype=ptype,pgal=pgal)
        # Return if no filename found
        if not fopt['fname']:
            mmlio.verbose('Could not find a file matching your description')
            return None
        # Load pNbody object
        if pNbody:
            import pNbody
            flag_ic=(fopt['ftype']=='gadget' and fopt['fext']=='ic')
            kwargs.update(p_name=fopt['fname'],ftype=fopt['ftype'],flag_ic=flag_ic)
            kwargs.setdefault('idgal2',self.galids[1])
            kwargs.setdefault('unitsfile',files.unitsfile(self.fdict,fopt['ftype']))
            snapObj = pNbody.Nbody(**kwargs)
        # Load modified pynbody object
        else:
            import pysim
            kwargs['filetype']=files.ftype2snapobj(fopt['ftype'])
            kwargs.setdefault('galids',self.galids)
            snapObj = pysim.loadsnap(fopt['fname'],**kwargs)
        #     print snapObj.header.npart
        #     print snapObj.galaxies()
        #     print snapObj.families()
        # mmlio.yorn('?')
        setattr(snapObj,'fext',fopt['fext'])
        return snapObj
            
    ################################################################################################################################
    # ORBIT METHOD
    def orbit_param(self,**kwargs):
        """
        Determines orbit parameters
        """
        from scipy import interpolate,optimize
        from mmlastro import mmlprofiles,mmlconst
        gals=self.galaxies
        fams=['all']#+self.families
        gconst=mmlconst.main('phys',subtag='galaxy')['G']
        # Get orbit from definition
        parPer={}
        M1_sol=loadsimobj(self.infodict['gal1']).infodict['mvir']
        M2_sol=loadsimobj(self.infodict['gal1']).infodict['mvir']
        Mtot_sol=M1_sol+M2_sol
        R1_pc=mmlprofiles.mvir2rvir(M1_sol,units='galaxy')
        R2_pc=mmlprofiles.mvir2rvir(M2_sol,units='galaxy')
        ecc=self.infodict['ecc']
        rtot_pc=R1_pc+R2_pc
        rinit_pc=rtot_pc
        rperi_pc=rtot_pc*self.infodict['rperi']
        minit_sol=Mtot_sol
        mperi_sol=Mtot_sol
        parPer['rinit']=(rinit_pc/1000.)
        parPer['rperi']=(rperi_pc/1000.)
        parPer['minit']=minit_sol
        parPer['mperi']=mperi_sol
        parPer['vinit']=np.sqrt(gconst*((ecc*mperi_sol/rperi_pc)+(minit_sol/rinit_pc)))
        parPer['vperi']=np.sqrt(1.+ecc)*np.sqrt(gconst*mperi_sol/rperi_pc)
        parPer['ecc']=ecc
        # Get orbit from track
        track=self.orbit_track(**kwargs)
        parObs={'tfinal':track['time'][-1],
                'fext_final':track['fext'][-1]}
        tperi={}
        mvar='1menc'
        for f in fams:
            parObs[f+'rinit']=track[f+'rpos'][0]
            parObs[f+'vinit']=track[f+'rvel'][0]
            parObs[f+'minit']=track[f+mvar][0]
            parObs[f+'rfinal']=track[f+'rpos'][-1]
            parObs[f+'vfinal']=track[f+'rvel'][-1]
            parObs[f+'mfinal']=track[f+mvar][-1]
            fofr=interpolate.interp1d(track['time'],track[f+'rpos'],kind='cubic')
            fofv=interpolate.interp1d(track['time'],track[f+'rvel'],kind='cubic')
            fofm=interpolate.interp1d(track['time'],track[f+mvar],kind='cubic')
            tperi[f]=optimize.leastsq(fofr,0.)[0][0]
            parObs[f+'rperi']=float(fofr(tperi[f]))
            parObs[f+'vperi']=float(fofv(tperi[f]))
            parObs[f+'mperi']=float(fofm(tperi[f]))
            vinit_c=np.sqrt(gconst*parObs[f+'minit']/(1000.*parObs[f+'rinit']))
            vperi_c=np.sqrt(gconst*parObs[f+'mperi']/(1000.*parObs[f+'rperi']))
            vfinal_c=np.sqrt(gconst*parObs[f+'mfinal']/(1000.*parObs[f+'rfinal']))
            parObs[f+'einit']=(parObs[f+'vinit']/vinit_c)**2 - 1.
            parObs[f+'eperi']=(parObs[f+'vperi']/vperi_c)**2 - 1.
            parObs[f+'efinal']=(parObs[f+'vfinal']/vfinal_c)**2 - 1.
        parObs['tperi']=tperi['all']
        parObs['fext_peri']=track['fext'][track['time']<=tperi['all']][-1]

        if kwargs.get('verbose',False):
            print self['runtag']
            print 'THEORY:'
            pprint.pprint(parPer)
            print ''
            print 'OBSERVATION:'
            pprint.pprint(parObs)
            print ''
            print 'TIMINGS:'
            pprint.pprint(tperi)

        return parPer,parObs,tperi
        

    def orbit_track(self,overwrite=False,owplot=False,inpar=None,**kwargs):
        """
        Plots the orbit of the secondary
        """
        from pysim.nbody import analysis,array,filt
        from pysim.nbody.family import get_family
        from pysim.nbody.galaxy import get_galaxy
        plotsize=6.
        timunit='Gyr'
        posunit='kpc'
        velunit='km s**-1'
        masunit='Msol'
        tagfmt='{f}{g[6]}{c}{v}'
        kwargs.setdefault('askuser',True)
                                      
        # Get parameter keys
        xyz=['x','y','z']
        famclr={'disk':'c','star':'b','bulge':'r','halo':'g','all':'k'}
        pkeys=['fext','time']
        gals=self.galaxies
        fams=['all']+self.families
        for f in fams:
            pkeys.append(f+'menc')
            pkeys+=[f+c+'pos' for c in ['r']+xyz]
            pkeys+=[f+c+'vel' for c in ['r']+xyz]
            for g in gals:
                pkeys.append(f+g[6]+'menc')
                for v in ['pos','vel']:
                    for c in xyz:
                        pkeys.append(tagfmt.format(f=f,g=g,c=c,v=v))
        # Pars input
        inpar=mmlparam.parspar('pysim.simlist','orbit',inpar=inpar,**kwargs)
        # Load existing prop
        ftable=mmlparam.par2file(self,'pysim.simlist','orbit',inpar=inpar,sing=True)['file']
        if os.path.isfile(ftable) and not overwrite:
            prop=mmlio.rwtable('R',ftable)
            for k in pkeys: 
                if k not in prop['keylist']: overwrite=True
        # Overwrite/intialize prop
        if overwrite or not os.path.isfile(ftable): 
            prop={k:[] for k in pkeys}
            prop['keylist']=pkeys
        # Get snapshots
        snapkey,snaptemp,snaplist=files.listsnap(self,inpar['ftype'])
        # Loop over snapshots
        loadkw={'ftype':inpar['ftype']}
        added=False
        for snap in snaplist:
            fext=files.get_ext(snap,snaptemp)
            if fext in prop['fext'] and not overwrite: continue
            added=True
            loadkw['fname']=snap
            loadkw['fext' ]=fext
            # Load snapshot
            snapObj=self.loadsnap(**loadkw)
            npart=len(snapObj)
            # Get data
            iprop={}
            for f in fams:
                fidx=np.zeros(npart,dtype=bool)
                if f=='all': fidx[:]=True
                else       : fidx[snapObj._get_family_slice(get_family(f))]=True
                # Absolute postions and velocities
                for g in gals:
                    gidx=np.zeros(npart,dtype=bool)
                    gidx[snapObj._get_galaxy_slice(get_galaxy(g))]=True
                    totidx=np.logical_and(gidx,fidx)

                    p,v=analysis.halo.center(snapObj[totidx],mode=inpar['cmeth'],retcen=True,vel=True)
                    p=p.in_units(posunit)
                    v=v.in_units(velunit)

                    for c in xyz:
                        iprop[tagfmt.format(f=f,g=g,c=c,v='pos')]=float(p[xyz.index(c)])
                        iprop[tagfmt.format(f=f,g=g,c=c,v='vel')]=float(v[xyz.index(c)])
                # Relative positions and velocities
                if len(gals)==2: 
                    for c in xyz:
                        iprop[f+c+'pos']=iprop[f+'2'+c+'pos']+iprop[f+'1'+c+'pos']
                        iprop[f+c+'vel']=iprop[f+'2'+c+'vel']+iprop[f+'1'+c+'vel']
                else:
                    for c in xyz:
                        iprop[f+c+'pos']=iprop[f+'1'+c+'pos']
                        iprop[f+c+'vel']=iprop[f+'1'+c+'vel']
                iprop[f+'rpos']=np.sqrt(iprop[f+'xpos']**2+iprop[f+'ypos']**2+iprop[f+'zpos']**2)
                iprop[f+'rvel']=np.sqrt(iprop[f+'xvel']**2+iprop[f+'yvel']**2+iprop[f+'zvel']**2)
                # Enclosed masses
                iprop[f+'menc']=0.
                for g in gals:
                    rcen='{} {}'.format(iprop[f+'rpos'],posunit)
                    ccen=array.SimArray([iprop[f+g[6]+c+'pos'] for c in xyz],posunit)
                    gidx=np.zeros(npart,dtype=bool)
                    gidx[snapObj._get_galaxy_slice(get_galaxy(g))]=True
                    totidx=np.logical_and(gidx,fidx)
                    iprop[f+g[6]+'menc']=float(snapObj[totidx][filt.Sphere(rcen,ccen)]['mass'].sum().in_units(masunit))
                    iprop[f+'menc']+=iprop[f+g[6]+'menc']
            # Add time and fext
            iprop['fext']=fext
            iprop['time']=float(snapObj.properties['time'].in_units(timunit))
            # Add properties to list
            for k in pkeys: prop[k].append(iprop[k])
            # Save properties
            mmlio.rwtable('W',ftable,prop,overwrite=True)
            #print fext
        # Sort arrays
        idxsort=sorted(range(len(prop['time'])),key=lambda k: prop['time'][k])
        prop={k:np.array(prop[k])[idxsort] for k in pkeys}
        # Calculate things
        # for f in fams:
        #     if self['runtyp']=='galsim':
        #         for c in xyz:
        #             prop[f+c+'pos']=prop[f+'1'+c+'pos']
        #             prop[f+c+'vel']=prop[f+'1'+c+'vel']
        #     else:
        #         for c in xyz:
        #             prop[f+c+'pos']=prop[f+'2'+c+'pos']-prop[f+'1'+c+'pos']
        #             prop[f+c+'vel']=prop[f+'2'+c+'vel']-prop[f+'1'+c+'vel']
        #     prop[f+'rpos']=np.sqrt(prop[f+'xpos']**2+prop[f+'ypos']**2+prop[f+'zpos']**2)
        #     prop[f+'rvel']=np.sqrt(prop[f+'xvel']**2+prop[f+'yvel']**2+prop[f+'zvel']**2)
        # Plot properties
        if (added or owplot) and len(prop['time'])>1:
            # Get filename
            plotfile=mmlparam.par2file(self,'pysim.simlist','orbit',inpar=inpar,plot=True,sing=True)['file']
            mmlfiles.mkdirs(os.path.dirname(plotfile))
            # Initialize figure and axes
            plt.clf() ; ncol=3 ; nrow=3
            fig=plt.figure(figsize=(ncol*plotsize,nrow*plotsize))
            axXY=plt.subplot(nrow,ncol,1,aspect='equal')
            axXZ=plt.subplot(nrow,ncol,2,aspect='equal')
            axYZ=plt.subplot(nrow,ncol,3,aspect='equal')
            axTR=plt.subplot(nrow,1,2)
            axTV=plt.subplot(nrow,1,3)
            def plotstart(inax,xarr,yarr,fam,inclstart=True):
                inax.plot(xarr,yarr,'-'+famclr[fam],label=fam)
                if inclstart:
                    inax.plot(xarr[0],yarr[0],'*'+famclr[fam])
            def setaspect(inax,aspect=1.):
                xlim=inax.get_xlim()
                ylim=inax.get_ylim()
                xaspect=xlim[1]-xlim[0]
                yaspect=ylim[1]-ylim[0]
                inax.set_aspect(aspect*xaspect/yaspect)
            # Loop over families
            tlim=[min(prop['time']),max(prop['time'])]
            rmax=0. ; vmax=0.
            for f in fams:
                rmax=max(rmax,prop[f+'rpos'].max())
                vmax=max(vmax,prop[f+'rvel'].max())
                plotstart(axXY,prop[f+'xpos'],prop[f+'ypos'],f)
                plotstart(axXZ,prop[f+'xpos'],prop[f+'zpos'],f)
                plotstart(axYZ,prop[f+'ypos'],prop[f+'zpos'],f)
                plotstart(axTR,prop['time'],prop[f+'rpos'],f,inclstart=False)
                plotstart(axTV,prop['time'],prop[f+'rvel'],f,inclstart=False)

            # Set limits & aspect ratios
            xyzlim=(-rmax,rmax)
            for iax in [axXY,axXZ,axYZ]:
                iax.set_xlim(xyzlim)
                iax.set_ylim(xyzlim)
                setaspect(iax,1.)
            axTR.set_xlim(tlim) ; axTR.set_ylim((0.,rmax)) ; setaspect(axTR,1./3.)
            axTV.set_xlim(tlim) ; axTV.set_ylim((0.,vmax)) ; setaspect(axTV,1./3.)
            # Labels
            tlab='Time [{}]'.format(timunit)
            vlab='Velocity [{}]'.format(velunit)
            rlab='Separation [{}]'.format(posunit)
            xlab='X [{}]'.format(posunit)
            ylab='Y [{}]'.format(posunit)
            zlab='Z [{}]'.format(posunit)
            axs=[axXY,axXZ,axYZ,axTR,axTV]
            xls=[xlab,xlab,ylab,tlab,tlab]
            yls=[ylab,zlab,zlab,rlab,vlab]
            for iax,ix,iy in zip(axs,xls,yls):
                iax.set_xlabel(ix)
                iax.set_ylabel(iy)
            plt.tight_layout()
            plt.savefig(plotfile,bbox_inches='tight')
            mmlio.verbose('Plot:')
            print '    '+plotfile
        # Return
        return prop
        
    def barprop(self,method,overwrite=False,plotsnap=False,owplot=False,owanim=False,fext0=None,**kwargs):
        """
        Plots bar properties
        """
        import matplotlib.pyplot as plt
        from pysim.nbody.analysis import bar
        pkeys=['fext','time','rmin','rmax','phi','amp']
        if   method=='Am2'    : pkeys+=[]
        elif method=='derAm2' : pkeys+=[]
        elif method=='ellipse': pkeys+=['ramp','Am4']
        else: raise Exception('Invalid barprop method: {}'.format(method))
        # Get file type
        if method in ['Am2','derAm2']: ftype='scf'
        else                         : ftype='gadget'
        # Get properties
        if ftype=='scf': kwargs['cmeth']=''
        inpar=mmlparam.parspar('pysim.simlist','barprop_'+method,**kwargs)
        barpar={k:inpar[k] for k in inpar if k not in ['tagstr']}
        # Load existing prop
        ftable=mmlparam.par2file(self,'pysim.simlist','barprop_'+method,inpar=inpar,sing=True)['file']
        if overwrite or not os.path.isfile(ftable): 
            prop={k:[] for k in pkeys}
            prop['keylist']=pkeys
        else:
            prop=mmlio.rwtable('R',ftable)
            pkeys=prop['keylist']
        # Get snapshots
        if not fext0:
            if ftype=='scf': fext0='*_'+inpar['family']+inpar['galaxy'][-1]
            else           : fext0=None
        snapkey,snaptemp,snaplist=files.listsnap(self,ftype,fext=fext0)
        # Loop over snapshots
        loadkw={'ftype':ftype}
        added=False
        for snap in snaplist:
            fext=files.get_ext(snap,snaptemp)
            if fext in prop['fext'] and not overwrite: continue
            added=True
            loadkw['fname']=snap
            loadkw['fext' ]=fext
            # Load snapshot
            snapObj=self.loadsnap(**loadkw)
            # Get snapplot
            snapplot=mmlparam.par2file(self,'pysim.simlist','barprop_'+method,inpar=inpar,plot=True,fext=fext)['file']
            mmlfiles.mkdirs(os.path.dirname(snapplot))
            # Get data
            if   method=='Am2'    : iprop=bar.Am2(snapObj,plot=plotsnap,plotfile=snapplot,**barpar)
            elif method=='derAm2' : iprop=bar.derAm2(snapObj,plot=plotsnap,plotfile=snapplot,**barpar) 
            elif method=='ellipse': iprop=bar.ellipse(snapObj,plot=plotsnap,plotfile=snapplot,**barpar)
            else: raise Exception('Unsupported method: {}'.format(method))
            # Add time and fext
            iprop['fext']=fext
            iprop['time']=float(snapObj.properties['time'].in_units('Gyr'))
            # Add properties to list
            for k in pkeys: prop[k].append(iprop[k])
            # Save properties
            mmlio.rwtable('W',ftable,prop,overwrite=True)
            print fext
        # Plot properties
        if (added or owplot) and len(prop['time'])>1:
            plotfile=mmlparam.par2file(self,'pysim.simlist','barprop_'+method,inpar=inpar,plot=True,sing=True)['file']
            mmlfiles.mkdirs(os.path.dirname(plotfile))
            idxsort=sorted(range(len(prop['time'])),key=lambda k: prop['time'][k])
            plotkeys=[k for k in pkeys if k not in ['time','fext']]
            plt.clf()
            ncol=2 ; nrow=int(np.ceil(float(len(plotkeys))/float(ncol)))
            tlim=[min(prop['time']),max(prop['time'])]
            for idx in range(len(plotkeys)):
                k=plotkeys[idx]
                iax=plt.subplot(ncol,nrow,idx+1)
                iax.plot(np.array(prop['time'])[idxsort],np.array(prop[k])[idxsort])
                iax.set_xlim(tlim)
                iax.set_ylim((min(prop[k]),max(prop[k])))
                iax.set_xlabel('Time (Gyr)')
                iax.set_ylabel(k.capitalize())
            plt.tight_layout()
            plt.savefig(plotfile)
            mmlio.verbose('Plot:')
            print '    '+plotfile
        # Animate
        if (plotsnap or owanim) and len(prop['time'])>1:
            mmlio.verbose('Animation:')
            self.animsnap('pysim.simlist','barprop_'+method,inpar=inpar,overwrite=True)
        # Return
        return prop

    ################################################################################################################################
    # PLOTTING METHODS
    def plotsnap_old(self,pmeth=None,**kwargs):
        """
        Plots a simulation snapshot
        """
        if pmeth not in LIST_PLOTMETH: 
            pmeth=mmlio.askselect('What kind of plot should be created?',LIST_PLOTMETH)
        if   pmeth in LIST_SINGPLOTMETH: return self.singplotsnap(pmeth=pmeth,**kwargs)
        elif pmeth in LIST_MULTPLOTMETH: return self.multplotsnap(pmeth=pmeth,**kwargs)
        return
    def plotsnap(self,pmeth=None,**kwargs):
        """
        Plots a simulation snapshot
        """
        if pmeth not in LIST_PLOTMETH: 
            pmeth=mmlio.askselect('What kind of plot should be created?',LIST_PLOTMETH)
        if   pmeth in LIST_SINGPLOTMETH: return self.plotsnap_sing(pmeth=pmeth,**kwargs)
        elif pmeth in LIST_MULTPLOTMETH: return self.plotsnap_mult(pmeth=pmeth,**kwargs)
        return

    def plotsnap_sing(self,snapObj=None,loadkw={},mode='all',make_plot=True,data=None,
                      pmeth=None,inpar={},plotfile=None,histfile=None,datafile=None,
                      overwrite=None,owdata=False,owhist=False,owanim=True,**kwargs):
        """
        Plot a simulation snapshot
        """
        import nbody
        # Setup
        if mode in ['all','init']:
            # Get input
            if pmeth is None: pmeth=mmlio.askselect('What kind of single plot should be created?',LIST_SINGPLOTMETH)
            if pmeth.startswith('barprop'): barmeth=pmeth.split('barprop_')[-1]
            else                          : barmeth=False
            if make_plot: 
                if not isinstance(overwrite,bool): overwrite=mmlio.yorn('Overwrite existing plots?')
                if not isinstance(owanim   ,bool): owanim   =mmlio.yorn('Overwrite existing animations?')
            else        : 
                overwrite=False
                owanim=False
            # Get loadkw
            loadkw.setdefault('ftype','gadget')
            loadkw.setdefault('fext' ,'000'   )
            if barmeth in ['Am2','derAm2','scfAm2']: loadkw['ftype']='scf'
            if loadkw['ftype']=='scf': kwargs['cmeth']=''
            # Assign histogram and plotting parameters/methods
            inpar=mmlparam.parspar('pysim.simlist',pmeth,inpar=inpar,plot=True,save=True,load=True,**kwargs)
            if   pmeth=='hist2d': inpar['histmeth']='numpy'  ; inpar['plotmeth']='contourf'
            elif pmeth=='pilimg': inpar['histmeth']='pNbody' ; inpar['plotmeth']='pilimg'
            elif barmeth        : inpar['histmeth']='numpy'  ; inpar['plotmeth']='contourf'
            loadkw['pNbody']=(inpar['histmeth']=='pNbody')
            # Get data table (if there is one)
            if barmeth:
                if not datafile: datafile=mmlparam.par2file(self,'pysim.simlist',pmeth,inpar=inpar,sing=True)['file']
                if owdata and os.path.isfile(datafile): os.remove(datafile)
                if not data and os.path.isfile(datafile): data=mmlio.rwtable('R',datafile)
                if data and 'keylist' in data:
                    misskeys=set(LIST_BARPROPMETH[barmeth]).difference(data['keylist'])
                    if len(misskeys)>0: owdata=True
                    if owdata:
                        mmlio.verbose('Missing keys: {}'.format(misskeys))
                        mmlio.verbose('Overwriting data')
                if owdata:
                    data={k:[] for k in LIST_BARPROPMETH[barmeth]}
                    data['keylist']=LIST_BARPROPMETH[barmeth]
            else: 
                datafile=None
                data=None
            # Define function for testing fext
            def testfext(fext):
                if plotfile:
                    ftest=plotfile
                else:
                    ftest=mmlparam.par2file(self,'pysim.simlist',pmeth,inpar=inpar,fext=fext,plot=True)['file']
                    ftest=ftest.replace('simlist',files.ftype2snaptyp(loadkw['ftype'],pNbody=loadkw['pNbody']))
                mmlfiles.mkdirs(os.path.dirname(ftest))
                if data: dataflag=(fext in data['fext'])
                else   : dataflag=True
                return overwrite or not os.path.isfile(ftest) or not dataflag
            # Format arguments
            kwargs=dict(pmeth=pmeth,inpar=inpar,loadkw=loadkw,data=data,make_plot=make_plot,testfext=testfext,
                        overwrite=overwrite,owdata=owdata,owhist=owhist,owanim=owanim,
                        plotfile=plotfile,histfile=histfile,datafile=datafile,**kwargs)
            # Return results/formated input
            if   mode=='all':
                if make_plot and not testfext(loadkw['fext']): return
            else: return kwargs
        # Plot
        if mode in ['all','run']:
            # Load snapshot
            if snapObj:
                if loadkw['pNbody']==issubclass(snapObj.__class__,nbody.snapshot.SimSnap):
                    if loadkw['pNbody']: mmlio.verbose('Snapshot should be pNbody. Reloading...')
                    else               : mmlio.verbose('Snapshot should not be pNbody. Reloading...')
                    snapObj=None
            if not snapObj: snapObj=self.loadsnap(**loadkw)
            # Plot & histogram file
            if not plotfile: plotfile=mmlparam.par2file(self,snapObj,pmeth,inpar=inpar,fext=snapObj.fext,plot=True )['file']
            if not histfile: histfile=mmlparam.par2file(self,snapObj,pmeth,inpar=inpar,fext=snapObj.fext,plot=False)['file']
            # Init plot kw
            plotkw=kwargs.pop('plotkw',{})
            plotkw.update(plotfile=plotfile,histfile=histfile,
                          overwrite=overwrite,owhist=owhist,**inpar)
            tagstr=plotkw.pop('tagstr',None)
            # Galaxy and family
            gal=plotkw['galaxy'] if plotkw['galaxy'] else 'all' 
            fam=plotkw['family'] if plotkw['family'] else 'all'
            # Plot
            if   pmeth=='hist1d': raise Exception('work in progress')
            elif pmeth in ['hist2d','pilimg','map2d','ellipse']:
                # Update special keywords
                plotkw.update(make_plot=make_plot)
                # Limits
                for v in ['x','y','w']:
                    if plotkw[v+'lim'][0]==plotkw[v+'lim'][1]:
                        plotkw[v+'lim']=self.limits[fam][gal].get(plotkw[v],(0.,0.))
                # Update keywords
                xvar=plotkw.pop('x')
                yvar=plotkw.pop('y')
                plotkw.update(wvar=plotkw.pop('w'),avar=plotkw.pop('a'))
                plotkw.setdefault('annotate',True)
                plotkw.setdefault('colorbar',True)
                if pmeth=='ellipse':
                    plotkw['plot_ellipse']=True 
                    plotkw['ellpfile']=mmlparam.par2file(self,snapObj,'ellipse',inpar=inpar,fext=snapObj.fext)['file']
                else:
                    plotkw['plot_ellipse']=False
                    plotkw['ellpfile']=None
                # Histogram
                out=nbody.plot.generic.plot2d(snapObj,xvar,yvar,**plotkw)
            # Bar properties
            elif pmeth.startswith('barprop'):
                from pysim.nbody.analysis import bar,langmm
                barmeth=pmeth.split('barprop_')[-1]
                # Update special keywords
                plotkw.update(owplot=make_plot,filename=datafile,
                              verbose=True,warn=False)
                # Limits
                if 'rlim' in plotkw:
                    if plotkw['rlim'][0]==plotkw['rlim'][1]:
                        plotkw['rlim']=self.limits[fam][gal].get('R',(0.,0.))
                # Get data
                out=langmm.bar(snapObj,barmeth,**plotkw)
            # Error
            else: raise Exception('Invalid plotting method: {}'.format(pmeth))
            # Print info
            if make_plot: print '    '+plotfile
            # Return control
            return out
            # if make_plot: return plotfile
            # else        : return out
        # Animate & plot data table
        if mode in ['end']:
            out=self.animsnap('pysim.simlist',pmeth,inpar=inpar,overwrite=owanim,fext=files.ftype2snapext(loadkw['ftype']),
                              replace=['simlist',files.ftype2snaptyp(loadkw['ftype'],pNbody=loadkw['pNbody'])])
            if pmeth.startswith('barprop'): 
                datafile=mmlparam.par2file(self,'pysim.simlist',pmeth,inpar=inpar,sing=True)['file']
                plotfile=mmlparam.par2file(self,'pysim.simlist',pmeth,inpar=inpar,plot=True,sing=True)['file']
                if os.path.isfile(datafile):
                    mmlio.verbose('Datafile:')
                    print '    '+datafile
                    out=mmlio.plottable(plotfile,datafile,xkey='time',xlab='Time (Gyr)',overwrite=owanim)
            return out

    def plotsnap_mult(self,snapObj=None,loadkw={},mode='all',make_plot=True,
                      pmeth=None,inpar={},overwrite=None,plotfile=None,**kwargs):
        """
        Plot multiple panels for simulation snapshot
        """
        import nbody
        import Image
        import matplotlib
        import matplotlib.pyplot as plt
        plotsize=5.
        cbfrac=0.4
        # Setup
        if mode in ['all','init']:
            # Get input
            if pmeth is None: pmeth=mmlio.askselect('What kind of multipanel plot should be created?',LIST_MULTPLOTMETH)
            if make_plot: 
                if not isinstance(overwrite,bool): overwrite=mmlio.yorn('Overwrite existing plots?')
            else        : overwrite=False
            # Get plotting parameters
            if pmeth=='gfpanel':
                inpar.update(galaxy='',family='')
                # Pars parameter for individual panels
                pmeth0=inpar.pop('pmeth' ,None)
                partg0=inpar.pop('tagstr',None)
                if pmeth0 is None: pmeth0=mmlio.askselect('What kind of single plot should be created?',LIST_SINGPLOTMETH)
                inpar0=mmlparam.parspar('pysim.simlist',pmeth0,inpar=inpar,plot=True,tagstr=partg0,**kwargs)
                if   pmeth0=='hist2d': inpar0['histmeth']='numpy'  ; inpar0['plotmeth']='contourf'
                elif pmeth0=='pilimg': inpar0['histmeth']='pNbody' ; inpar0['plotmeth']='pilimg'
                loadkw['pNbody']=(inpar0['histmeth']=='pNbody')
                kwargs['askuser']=False
                # Add tag for each galaxy/family pair
                gallist=['']+self.galaxies
                famlist=['']+self.families
                if 'star' in famlist: famlist.remove('star')
                pmeth='multipanel'
                inpar={'name':'gfpanel_'+mmlparam.par2tag('pysim.simlist',pmeth0,inpar=inpar0),
                       'ncol':len(gallist),'plotlist':[]}
                if 'tagstr' in inpar0: del inpar0['tagstr']
                for f in famlist:
                    for g in gallist:
                        inpar0.update(galaxy=g,family=f)
                        inpar0=mmlparam.parspar('pysim.simlist',pmeth0,inpar=inpar0,plot=True,save=True,verbose=False)
                        inpar['plotlist'].append(pmeth0+'_'+mmlparam.par2tag('pysim.simlist',pmeth0,inpar=inpar0))
            inpar=mmlparam.parspar('pysim.simlist',pmeth,inpar=inpar,plot=True,save=True,**kwargs)
            # Get loadkw
            loadkw.setdefault('ftype' ,'gadget')
            loadkw.setdefault('fext'  ,'000'   )
            if 'pNbody' not in loadkw: loadkw['pNbody']=mmlio.yorn('Use pNbody?')
            # Define function for testing fext
            def testfext(fext):
                if plotfile:
                    ftest=plotfile
                else:
                    ftest=mmlparam.par2file(self,'pysim.simlist',pmeth,inpar=inpar,fext=fext,plot=True)['file']
                    ftest=ftest.replace('simlist',files.ftype2snaptyp(loadkw['ftype'],pNbody=loadkw['pNbody']))
                return overwrite or not os.path.isfile(ftest)
            # Format arguments
            kwargs=dict(pmeth=pmeth,overwrite=overwrite,inpar=inpar,loadkw=loadkw,
                        make_plot=make_plot,testfext=testfext,plotfile=plotfile,**kwargs)
            # Return results/formated input
            if   mode=='all':
                if make_plot and not testfext(loadkw['fext']): return
            else: return kwargs
        # Plot
        if mode in ['all','run']:
            # Load snapshot
            if snapObj:
                if loadkw['pNbody']==issubclass(snapObj.__class__,nbody.snapshot.SimSnap):
                    if loadkw['pNbody']: mmlio.verbose('Snapshot should be pNbody. Reloading...')
                    else               : mmlio.verbose('Snapshot should not be pNbody. Reloading...')
                    snapObj=None
            if not snapObj: snapObj=self.loadsnap(**loadkw)
            # Plot file
            if not plotfile: plotfile=mmlparam.par2file(self,snapObj,pmeth,inpar=inpar,fext=snapObj.fext,plot=True)['file']
            # Plotkw
            plotkw=kwargs.get('plotkw',{})
            # Plot
            if pmeth=='multipanel':
                # Initialize axes
                ncol=inpar['ncol']
                nrow=int(np.ceil(float(len(inpar['plotlist']))/float(ncol)))
                if not loadkw['pNbody']:
                    fig,axs=plt.subplots(nrow,ncol,squeeze=False,figsize=((1.+cbfrac)*plotsize*ncol,plotsize*nrow))
                if loadkw['pNbody']: nbins=512 ; img=Image.new('RGB',(nbins*ncol,nbins*nrow))
                else               : nbins=100
                # Loop over plots
                for idx,iplot in enumerate(inpar['plotlist']):
                    irow=idx/ncol
                    icol=idx%ncol
                    ipmeth,itag=iplot.split('_',1)
                    print ipmeth,itag,irow,icol
                    if loadkw['pNbody']: iaxs=(img,irow*nbins,icol*nbins)
                    else               : iaxs=axs[irow,icol]#plt.subplot(nrow,ncol,idx)
                    # Load parameters
                    iparfile=mmlparam.par2fname('pysim.simlist',ipmeth,tagstr=itag)
                    if os.path.isfile(iparfile): ippar=mmlparam.loadpar('pysim.simlist',ipmeth,iparfile)
                    else: ippar=mmlparam.parspar('pysim.simlist',ipmeth,plot=True,askuser=True,tagstr=itag)
                    # Plot
                    plotkw['axs']=iaxs
                    x=self.plotsnap_sing(snapObj=snapObj,loadkw=loadkw,make_plot=False,
                                         pmeth=ipmeth,inpar=ippar,plotkw=plotkw)
                # Save plot or return it
                if make_plot:
                    mmlfiles.mkdirs(os.path.dirname(plotfile))
                    if loadkw['pNbody']: img.save(plotfile)
                    else               : 
                        plt.tight_layout()
                        plt.savefig(plotfile)
                else:
                    if loadkw['pNbody']: out=img
                    else               : out=plt.gcf()
            else: raise Exception('Invalid plotting method: {}'.format(pmeth))
            # Print info
            if make_plot: print '    '+plotfile
            # Return control
            if make_plot: return plotfile
            else        : return out
        # Animate
        if mode in ['end']:
            owanim=mmlio.yorn('Overwrite existing animation?')
            out=self.animsnap('pysim.simlist',pmeth,inpar=inpar,overwrite=owanim,fext=files.ftype2snapext(loadkw['ftype']),
                              replace=['simlist',files.ftype2snaptyp(loadkw['ftype'],pNbody=loadkw['pNbody'])])
            return out


    def multplotsnap(self,pmeth=None,snapObj=None,overwrite=None,inpar=None,fname=None,
                     make_plot=True,plotkw={},firstsnap=False,pNbody=None,**kwargs):
        """
        Plots multiple panels for a simulation snapshot
        """
        import nbody
        import matplotlib.pyplot as plt
        import Image
        # Get plotting method
        if pmeth not in LIST_MULTPLOTMETH: 
            pmeth=mmlio.askselect('What kind of multipanel plot should be created?',LIST_MULTPLOTMETH)
        if not isinstance(overwrite,bool):
            overwrite=mmlio.yorn('Overwrite existing plots?')
        # Plot families/galaxies
        if pmeth=='gfpanel':
            # Parameters for each panel
            kwargs['galaxy']=''
            kwargs['family']=''
            pmethsing=kwargs.get('pmethsing',None)
            otkw=self.singplotsnap(pmeth=pmethsing,snapObj=snapObj,overwrite=overwrite,inpar=inpar,firstsnap=True,pNbody=pNbody,**kwargs)
            snapObj=otkw['snapObj']
            pNbody=otkw['pNbody']
            # List of galaxies & families
            if pNbody: 
                famlist=snapObj.families()
                gallist=snapObj.galaxies()
            else:
                famlist=[f.name for f in snapObj.families()]
                gallist=[g.name for g in snapObj.galaxies()]
            famlist=['']+famlist
            if len(gallist)>1: gallist=['']+gallist
            # Initialize parameters for multipanel
            outpar={'name':'gfpanel_'+mmlparam.par2tag(snapObj,otkw['pmeth'],inpar=otkw['inpar']),
                    'ncol':len(gallist),
                    'plotlist':[]}
            # Create plotlist
            ppar=otkw['inpar']
            for f in famlist:
                for g in gallist:
                    ppar.update(galaxy=g,family=f)
                    ppar=mmlparam.parspar(snapObj,otkw['pmeth'],inpar=ppar,plot=True,save=True,verbose=False)
                    outpar['plotlist'].append(otkw['pmeth']+'_'+mmlparam.par2tag(snapObj,otkw['pmeth'],inpar=ppar))
            # Create dictionary for quick replot
            outargs=dict(pmeth='multipanel',overwrite=overwrite,pNbody=pNbody,snapObj=snapObj,inpar=outpar,plotkw=plotkw)
            # Plot
            return self.multplotsnap(firstsnap=firstsnap,fname=fname,make_plot=make_plot,**outargs)
        # Multipanel plot
        elif pmeth=='multipanel':
            if not isinstance(pNbody,bool): pNbody=False
            if not snapObj: snapObj=self.loadsnap(pNbody=pNbody,**kwargs)
            # Get plotting parameters
            ppar=mmlparam.parspar(snapObj,pmeth,inpar=inpar,plot=True,**kwargs)
            # Return info
            if firstsnap:
                ppar0=copy.deepcopy(ppar) ; tagstr=ppar0.pop('tagstr',None)
                outargs=dict(pmeth=pmeth,overwrite=overwrite,pNbody=pNbody,snapObj=snapObj,inpar=ppar0,plotkw=plotkw)
                return outargs
            # Get file name
            if not fname:
                files=mmlparam.par2file(self,snapObj,pmeth,plot=True,inpar=ppar,fext=snapObj.fext)
                fname=files['file']
            if not overwrite and os.path.isfile(fname): return fname
            mmlfiles.mkdirs(os.path.dirname(fname))
            # Initialize axes
            ncol=ppar['ncol']
            nrow=np.ceil(float(len(ppar['plotlist']))/float(ncol))
            if pNbody: nbins=512 ; img=Image.new('RGB',(nbins*ncol,nbins*nrow))
            else     : nbins=100
            # Loop over plots
            for idx,iplot in enumerate(ppar['plotlist']):
                irow=idx%ncol
                icol=idx-irow*ncol
                ipmeth,itag=iplot.split('_',1)
                if pNbody: iaxs=(img,irow*nbins,icol*nbins)
                else     : iaxs=plt.subplot(nrow,ncol,idx,aspect='equal')
                print ipmeth,pNbody
                # Load parameters
                iparfile=mmlparam.par2fname(snapObj,ipmeth,tagstr=itag)
                if os.path.isfile(iparfile): ippar=mmlparam.loadpar(snapObj,ipmeth,iparfile)
                else: ippar=mmlparam.parspar(snapObj,ipmeth,plot=True,askuser=True,tagstr=itag)
                # Plot
                x=self.singplotsnap(pmeth=ipmeth,snapObj=snapObj,inpar=ippar,
                                    nbins=nbins,make_plot=False,axs=iaxs)#,load=True,init=True)
            # Save plot or return it
            if make_plot:
                if pNbody: img.save(fname)
                else     : plt.savefig(fname)
                return fname
            else:
                if pNbody: return img
                else     : return plt.gcf()

    def singplotsnap(self,pmeth=None,snapObj=None,overwrite=None,inpar={},fname=None,
                     make_plot=True,plotkw={},firstsnap=False,axs=None,ftype='gadget',**kwargs):
        """
        Plots a simulation snapshot
        """
        import nbody
        # Get plotting method
        if pmeth not in LIST_SINGPLOTMETH: 
            pmeth=mmlio.askselect('What kind of single plot should be created?',LIST_SINGPLOTMETH)
        if not isinstance(overwrite,bool) and make_plot:
            overwrite=mmlio.yorn('Overwrite existing plots?')
        if   pmeth=='pilimg': pNbody=True
        elif pmeth=='map2d' : 
            if   'histmeth' in inpar : histmeth=inpar['histmeth']
            elif 'histmeth' in kwargs: histmeth=kwargs['histmeth']
            else: histmeth=mmlio.askselect('Choose a histogram method:',['numpy','pNbody'])
            kwargs['histmeth']=histmeth
            pNbody=(kwargs['histmeth']=='pNbody')
        else: pNbody=kwargs.get('pNbody',False)
        kwargs['pNbody']=pNbody
        print 'singplot',pmeth,pNbody
        # Get snapshot object
        if not snapObj: snapObj=self.loadsnap(ftype=ftype,**kwargs)
        else:
            if pNbody and issubclass(snapObj.__class__,nbody.snapshot.SimSnap):
                mmlio.verbose('Snapshot should be pNbody. Reloading...')
                snapObj=self.loadsnap(ftype=ftype,pNbody=pNbody,fext=snapObj.fext)
            elif not pNbody and not issubclass(snapObj.__class__,nbody.snapshot.SimSnap):
                mmlio.verbose('Snapshot should not be pNbody. Reloading...')
                snapObj=self.loadsnap(ftype=ftype,pNbody=pNbody,fext=snapObj.fext)
        # Get plotting parameters
        ppar=mmlparam.parspar(snapObj,pmeth,inpar=inpar,plot=True,**kwargs)
        if   pmeth=='hist2d': ppar['histmeth']='numpy'  ; ppar['plotmeth']='contourf'
        elif pmeth=='pilimg': ppar['histmeth']='pNbody' ; ppar['plotmeth']='pilimg'
        # Return info
        if firstsnap:
            ppar0=copy.deepcopy(ppar) ; tagstr=ppar0.pop('tagstr',None)
            outargs=dict(pmeth=pmeth,overwrite=overwrite,pNbody=pNbody,snapObj=snapObj,inpar=ppar0,plotkw=plotkw)
            return outargs
        # Get file name
        if make_plot:
            if not fname:
                files=mmlparam.par2file(self,snapObj,pmeth,plot=True,inpar=ppar,fext=snapObj.fext)
                fname=files['file']
            if not overwrite and os.path.isfile(fname): return fname
            mmlfiles.mkdirs(os.path.dirname(fname))
        # Plot
        if   pmeth=='hist1d': raise Exception('work in progress')
        elif pmeth in ['hist2d','pilimg','map2d','ellipse']:
            # Files
            histfile=mmlparam.par2file(self,snapObj,'map2d'  ,inpar=ppar,fext=snapObj.fext)['file']
            ellpfile=mmlparam.par2file(self,snapObj,'ellipse',inpar=ppar,fext=snapObj.fext)['file']
            mmlfiles.mkdirs([os.path.dirname(ellpfile),os.path.dirname(histfile)])
            # Limits
            lim=self.limits
            gal=ppar['galaxy'] if ppar['galaxy'] else 'all' 
            fam=ppar['family'] if ppar['family'] else 'all'
            ppar['xlim']=ppar['xlim'] if ppar['xlim'][0]!=ppar['xlim'][1] else lim[fam][gal].get(ppar['x'],None)
            ppar['ylim']=ppar['ylim'] if ppar['ylim'][0]!=ppar['ylim'][1] else lim[fam][gal].get(ppar['y'],None)
            ppar['wlim']=ppar['wlim'] if ppar['wlim'][0]!=ppar['wlim'][1] else lim[fam][gal].get(ppar['w'],None)
            # print lim[fam][gal].keys()
            # print lim[fam][gal]['pot']
            # Keywords
            plotkws=dict(wvar=ppar['w'],avar=ppar['a'],colorbar=True,annotate=True,make_plot=make_plot,axs=axs,
                         plotfile=fname,histfile=histfile,overwrite=overwrite,**ppar)
            if pmeth=='ellipse':
                plotkws['plot_ellipse']=True 
                plotkws['ellpfile']=ellpfile
            else:
                plotkws['plot_ellipse']=False
                plotkws['ellpfile']=None
            # Histogram
            plotkws.update(**plotkw)
            out=nbody.plot.generic.plot2d(snapObj,ppar['x'],ppar['y'],**plotkws)
        else: raise Exception('Invalid plotting method: {}'.format(pmeth))
        # Print info & return control
        if make_plot: return fname
        else        : return out

    def animsnap(self,snapObj,pmeth,inpar=None,tagstr=None,fext='*',overwrite=None,replace=None):
        """
        Creates an animation
        """
        # Pars input
        if isinstance(fext,str): fext=[fext]
        if not isinstance(fext,list):
            raise Exception('fext must be a string or a list of strings.')
        if len(fext)==0: fext=['{*}']
        # Get templates
        fext0='{test}'
        fdict0=mmlparam.par2file(self,snapObj,pmeth,plot=True,inpar=inpar,fext=fext0)
        if replace: 
            ftemp0=fdict0['file'].replace(replace[0],replace[1])
            fanim0=fdict0['anim'].replace(replace[0],replace[1])
        else:
            ftemp0=fdict0['file']
            fanim0=fdict0['anim']
        plottemp=[]
        for iext in fext:
            fdict=mmlparam.par2file(self,snapObj,pmeth,plot=True,inpar=inpar,fext=iext)
            if replace: itemp=fdict['file'].replace(replace[0],replace[1])
            else      : itemp=fdict['file']
            plottemp.append(itemp)
        # Get sorted file list
        plotlist0=mmlfiles.ls(plottemp)
        # Fill in missing frames with duplicates
        frameno=lambda x: int(float(max(re.findall(r'\d+',files.get_ext(x,ftemp0)),key=len)))
        plotlist=[]
        for i in range(len(plotlist0)):
            if i==(len(plotlist0)-1): plotlist.append(plotlist0[i])
            else:
                f0=frameno(plotlist0[i])
                f1=frameno(plotlist0[i+1])
                plotlist+=(f1-f0)*[plotlist0[i]]
        # Make movie
        if len(plotlist)==0:
            return None
        else:
            movestr=mmlplot.ffmpeg(plotlist,fanim0,overwrite=overwrite,rmtemp=True,r=100,
                                   tempDir=os.path.join(os.path.dirname(fanim0),'temp_'+fanim0.split('.')[1]))
            mmlio.verbose('Animation:')
            print '    '+fanim0
            return fanim0

    ################################################################################################################################
    ################################################################################################################################
    # FILE METHODS

    def mvrun(self,**options):
        """
        Move/rename files assosiated with a run  
        """
        # Get destination simobj
        mmlio.verbose('Getting info on new file names...')
        isimstr=copy.deepcopy(self)
        fsimstr=asksimboj(infodict=isimstr.infodict,runtyp=isimstr['runtyp'])
        fname=mmlparam.par2fname('pysim.simlist','mmlsim',inpar=dict(fsimstr))
        # Rename old run files
        fname_new=fname+'_new'
        fname_old=fname+'_old'
        if os.path.isfile(fname): os.rename(fname,fname_old)
        if os.path.isfile(fname_new):
            mmlio.verbose('Using info from new file: {}'.format(fname_new))
            shutil.copy2(fname_new,fname)
        shutil.copy2(fname,fname_new)
        # Move files
        try:
            for ftype in files.LIST_SNAPFORM:
                files.mvfiles(isimstr,fsimstr,ftype,groups='all',**options)
            isimstr.cleantree(**options)
            if os.path.isfile(fname_old): os.remove(fname_old)
            os.remove(fname_new)
        except:
            mmlio.verobse("Error moving files:", sys.exc_info()[0])
            mmlio.verbose("Moving run files back...")
            os.rename(fname,fname_new)
            if os.path.isfile(fname_old): os.rename(fname_old,fname)
            raise
        # Return control
        return

    def mvtree(self,comp_end=None):
        """
        Moves an entire run tree without changing file names
        """
        # Pars input
        if comp_end not in mmlinfo.LIST_COMPUTERS:
            comp_end=mmlio.askselect('What computer should the tree be moved to?',mmlinfo.LIST_COMPUTERS)
        # Get location of sim
        comp_beg=self['memcomp']
        # Get filenames
        simstr1=copy.deepcopy(self) ; simstr1.update(memcomp=comp_beg)
        simstr2=copy.deepcopy(self) ; simstr2.update(memcomp=comp_end)
        files_beg=files.get_fdict(**simstr1)
        files_end=files.get_fdict(**simstr2)
        # Get directories
        memlist = mmlinfo.computers()['memlist']
        dir_beg=files_beg['rundir']
        dir_end=files_end['rundir']
        if comp_beg not in memlist: dir_beg='{}:{}'.format(comp_beg,dir_beg)
        if comp_end not in memlist: dir_end='{}:{}'.format(comp_end,dir_end)
        # Copy the tree
        mmlio.verbose('Source      run tree: '+dir_beg)
        mmlio.verbose('Destination run tree: '+dir_end)
        if mmlio.yorn('Continue moving tree?'):
            # Create directory
            if comp_end in memlist: mmlfiles.mkdirs(os.path.dirname(dir_end))
            else                  : subprocess.call(['ssh',comp_end,'mkdir','-p',os.path.dirname(dir_end.split(':')[-1])])
            # Copy tree
            if comp_beg in memlist and comp_end in memlist: shutil.copytree(dir_beg,dir_end)
            else                                          : subprocess.call(['rsync','-avz',dir_beg,dir_end])
        else: return
        # Remove the old tree
        if comp_beg in memlist and mmlio.yorn('Remove old run tree?'):
            shutil.rmtree(dir_beg)
        # Update the run files
        if mmlio.yorn('Update the run files?'):
            self['memcomp']=comp_end
            mmlparam.savepar('pysim.simlist','galsim',dict(self),overwrite=True)
        return

    def ziprun(self):
        """
        Compress the entire directory tree
        """
        zipcmd=['tar','-zcvf',self.fdict['archiv'],self.fdict['rundir'],'--remove-files']
        subprocess.call(zipcmd)
        shutil.rmtree(self.fdict['rundir'])
        return

    def unziprun(self):
        """
        Uncompresses zipped directory tree
        """
        zipcmd=['tar','-zxvf',self.fdict['archiv'],'-C',self.fdict['simdir']]
        subprocess.call(zipcmd)
        os.remove(self.fdict['archiv'])
        return


    ####################################################################################################################################
    ####################################################################################################################################
    # RUN PROPERTY METHODS

    def get_Mvir(self,galaxy=None):
        """
        Returns the virial mass for a given simulation
        """
        if self['runtyp']=='galsim': 
            Mvir=self.infodict['mvir']
        else:
            galaxy=mmlpars.mml_pars(galaxy,default=1,list=[1,2])
            Mvir=mmlparam.loadpar('pysim.simlist','galsim',self.infodict['gal{}'.format(galaxy)])['mvir']
        return Mvir

    def get_Rvir(self,galaxy=None,**exkw):
        """
        Returns the virial radius for a given simulations (in pc)
        """
        from mmlastro import mmlprofiles
        Mvir=self.get_Mvir(galaxy=galaxy)
        exkw['units']='galaxy'
        Rvir=mmlprofiles.mvir2rvir(Mvir,**exkw)
        return Rvir
    
    def get_Tperi(self,**exkw):
        """
        Returns the time that pericenter occurs at
        """
        if self['runtyp']!='intsim': raise Exception('tperi only valid for interactions (runtyp={})'.format(self['runtyp']))
        parPer,parObs,parTim=self.orbit_param(**exkw)
        return parTim['all']  

    def get_Rperi(self,**exkw):
        """
        Returns the separation at pericenter
        """
        if self['runtyp']!='intsim': raise Exception('rperi only valid for interactions (runtyp={})'.format(self['runtyp']))
        parPer,parObs,parTim=self.orbit_param(**exkw)
        return parObs['allrperi']

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

    def get_galbase(self,galaxy=None,striptag=True):
        """
        Returns the mmlsim object for a selected galaxy
        """
        aliases={'gal1':['gal1','galaxy1','primary'  ,'1',1,None],
                 'gal2':['gal2','galaxy2','secondary','2',2]}
        # Select correct sim dict for isolated primary
        if   self['runtyp']=='galsim': runtag=self['runtag']
        elif self['runtyp']=='intsim': 
            if   galaxy in aliases['gal1']: runtag=self.infodict['gal1']
            elif galaxy in aliases['gal2']: runtag=self.infodict['gal2']
        # Get info from runtag
        galobj=mmlparam.loadpar('pysim.simlist','galsim',runtag)
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
        if self['runtag']==runtag0: simstr0=self
        else                      : simstr0=loadsimobj(runtag0)
        return simstr0
        
    def get_primary(self,striptag=False):
        """
        Returns the mmlsim object for the primary galaxy
        """
        # Select correct sim dict for isolated primary
        if   self['runtyp']=='galsim': runtag=self['runtag']
        elif self['runtyp']=='intsim': runtag=self.infodict['gal1']
        # Get info from runtag
        galobj=mmlparam.loadpar('pysim.simlist','galsim',runtag)
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
        simstr0=loadsimobj(runtag0)
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
        from pysim.files.buildgal import model2param
        if 'timetest' in self['runtag']: ntot=long(10.**(mmlstring.str2dec(self['subtag'].split('on')[0])))
        else:
            if   self['runtyp']=='galsim': ntot=model2param(self.infodict)['total_n']
            elif self['runtyp']=='intsim':
                ntot1=get_ntot(runtag=self.infodict['gal1'])
                ntot2=get_ntot(runtag=self.infodict['gal2'])
                ntot=ntot1+ntot2
            else: raise Exception('Invalid run type: {}'.format(self['runtyp']))
        return ntot

    def get_galids(self):
        """
        Returns the bounding IDs for each galaxy
        """
        from pysim.files.buildgal import model2param
        if   self['runtyp']=='galsim': ntot=[self.get_ntot()]
        elif self['runtyp']=='intsim':
            ntot1=get_ntot(runtag=self.infodict['gal1'])
            ntot2=get_ntot(runtag=self.infodict['gal2'])
            ntot=[ntot1,ntot1+ntot2]
        galids=[0]+[ntot2galid(iN) for iN in ntot]
        if 'timetest' in self['runtag']: galids=[0]
        return galids

    def get_profpar(self,overwrite=False,**exkw):
        from mmlastro import mmlprofiles
        from nbody.array import SimArray
        # Load profile parameters if the exist
        if os.path.isfile(self.fdict['profpar']) and not overwrite:
            fd=open(self.fdict['profpar'],'r')
            out=pickle.load(fd)
            fd.close()
        # Otherwise get it from the IC file
        else:
            import nbody
            pn=self.loadsnap(fext='ic',**exkw)
            if not pn: return None
            rmax=SimArray(self.get_Rvir(),'pc').in_units(pn['pos'].units)
            out={}
            # Loop over families and galaxies present
            for f in pn.families():
                out[f.name]={}
                if f.name=='star': continue
                for g in pn.galaxies():
                    out[f.name][g.name]={}
                    # Select and center data
                    ipn=pn[f][g]
                    #nbody.analysis.angmom.faceon(ipn)
                    # Get profile ID and coordinate lists
                    profid='FLATDISK' if f.name=='disk' else 'HERNQUIST'
                    if   profid=='FLATDISK': coordlist=['R']
                    elif profid=='HIYASHI' : coordlist=['R','z']
                    else: coordlist=['r']
                    # Loop over coordinate list
                    for coord in coordlist:
                        out[f.name][g.name][coord]=nbody.analysis.langmm.fitprofile(ipn,profid,coord=coord,rmax=rmax)
                        print f.name,g.name,profid,coord,out[f.name][g.name][coord]
            out['star']=out['disk']
            # Write to file
            fd=open(self.fdict['profpar'],'w')
            pickle.dump(out,fd)
            fd.close()
        # Return output
        return out

    def get_limits(self,overwrite=False,**exkw):
        """
        Get limits
        """
        sclfact=6.
        mvar=['m','mass']
        pvar=['pos','x','y','z','scfpos']
        rvar=['r','R','rxy','rxyz']
        # Load profile parameters if the exist
        if os.path.isfile(self.fdict['limits']) and not overwrite:
            fd=open(self.fdict['limits'],'r')
            out=pickle.load(fd)
            fd.close()
        # Otherwise get it from the IC file
        else:
            from nbody.array import SimArray
            import nbody
            pn=self.loadsnap(fext='ic',**exkw)
            if not pn: return None
            # Get lists of galaxies
            famlist=pn.families()
            gallist=pn.galaxies()
            ngal=float(len(gallist)) 
            if   ngal==1: galfact=1.
            elif ngal==2: galfact=4.
            else: raise Exception('Set galfact for ngal={}'.format(ngal))
            # Loop over families and galaxies present
            out={'all':{g.name:{} for g in gallist}}
            for f in famlist:
                out[f.name]={'all':{}}
                for g in gallist:
                    ilim={}
                    ifmax=out['all'][g.name]
                    igmax=out[f.name]['all']
                    # Get profpar
                    profpar=self.get_galbase(g.name).get_profpar()
                    profpar[f.name][g.name]=profpar[f.name]['galaxy1']
                    # Select and center data
                    ipn=pn[f][g]
                    nbody.analysis.angmom.faceon(ipn)
                    # Limits
                    ilim['n']=SimArray((1.,float(len(ipn['mass']))/(10.**2)))
                    mlim=SimArray((ipn['mass'].min(),ipn['mass'].sum()/(10.**2)),ipn['mass'].units)
                    print f,g,mlim
                    for var in mvar: ilim[var]=mlim
                    # Position limits
                    pscl=0.
                    for scl in profpar[f.name][g.name].keys():
                        if scl=='z': pscl=max(pscl,profpar[f.name][g.name][scl]['zs'])
                        else       : pscl=max(pscl,profpar[f.name][g.name][scl]['rs'])
                    pscl=float(SimArray(pscl,profpar[f.name][g.name][scl]['xunits']).in_units(ipn['eps'].units))
                    plim=mmlmath.oom(sclfact*pscl)/2.
                    rlim=(ipn['eps'].min(),plim)
                    for var in pvar: ilim[var]=SimArray((-plim,plim),ipn['eps'].units)
                    for var in rvar: ilim[var]=SimArray(rlim,ipn['eps'].units)
                    # Density limits
                    marr=copy.deepcopy(ilim['m']) ; marr[0]*=10.
                    varr=(4./3.)*math.pi*(ilim['r']**3)
                    ilim['rho']=marr[::-1]/varr[::-1]
                    # Get maximums
                    for var in ilim.keys():
                        iunit=ilim[var].units
                        iflim=ilim[var]
                        if   var in pvar: iglim=galfact*ilim[var]
                        elif var in rvar: iglim=SimArray((ilim[var][0],galfact*ilim[var][1]),iunit)
                        else            : iglim=ilim[var]
                        if var in igmax: igmax[var]=SimArray((min(igmax[var][0],iglim[0]),max(igmax[var][1],iglim[1])),iunit)
                        else           : igmax[var]=iglim
                        if var in ifmax: ifmax[var]=SimArray((min(ifmax[var][0],iflim[0]),max(ifmax[var][1],iflim[1])),iunit)
                        else           : ifmax[var]=iflim
                    # Add to output
                    #print f.name,g.name,ilim
                    print ilim['mass']
                    out[f.name][g.name]=ilim
                    out['all'][g.name]=ifmax
                    out[f.name]['all']=igmax
            # Add limits for all galaxies & families
            out['all']['all']={}
            for f in famlist:
                ilim=out[f.name]['all']
                imax=out['all']['all']
                for var in ilim:
                    if var in imax: imax[var]=SimArray((min(imax[var][0],ilim[var][0]),max(imax[var][1],ilim[var][1])),ilim[var].units)
                    else          : imax[var]=ilim[var]
                out['all']['all']=imax
            # Ensure aliases
            if 'dm' in out: out['halo']=out['dm']
            # Write to file
            fd=open(self.fdict['limits'],'w')
            pickle.dump(out,fd)
            fd.close()
        # Return output
        return out

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

    ################################################################################################################################
    ################################################################################################################################
    # REGULARLY USED SUB METHODS
    def stdplots(self,method,family='disk',galaxy='galaxy1',fext='*',
                 overwrite=False,owdata=False,owhist=False,owanim=True,itermeth='itersnap'):
        """
        Creates standard plots
        """
        itermeth='itersnap'
        # Set defaults
        plotfext=fext
        pmeth='map2d'
        ftype='gadget'
        plotkw={}
        inpar=dict(histmeth='numpy',plotmeth='contourf',
                   family=family,galaxy=galaxy,
                   x='x'   ,xscl='lin',xlim=(0.,0.),
                   y='y'   ,yscl='lin',ylim=(0.,0.),
                   w='mass',wscl='log',wlim=(0.,0.),
                   a='',cmeth='pot')
        # PIL image of projected mass
        # disk_galaxy1_x_y_logm_potCenter
        if method=='img':
            pmeth='pilimg'
            inpar.update(w='m',histmeth='pNbody',plotmeth='pilimg')
        # Contourf of projected mass
        # numpy_disk_galaxy1_x_y_logmass_potCenter
        elif method=='ctr': pass
        # Ellipses
        # numpy_disk_galaxy1_x_y_logmass_potCenter
        elif method=='ell':
            pmeth='ellipse'
        # SCF Am2 plots
        elif method=='scf':
            ftype='scf'
            inpar.update(w='Am2',wscl='symlog',wlim=(0.0001,1.),a='N')
            plotkw.update(cmap=plt.get_cmap('jet'),fgcolor=[0,0,0],
                          cbtickloc=[-1,-0.1,-0.01,-0.001,0,0.001,0.01,0.1,1],
                          cbtickfmt='% .3f',
                          annotate_arrow=False)
        # Panel plots
        elif method=='pnl':
            pmeth='gfpanel'
            inpar={'pmeth':'hist2d','tagstr':'x_y_logmass_potCenter'}
            itermeth='itersnap'
        # Bars plots
        elif method.startswith('bar'):
            pmeth=method
            if   method.startswith('bar_'    ): barmeth=method.split('bar_')[-1]
            elif method.startswith('barprop_'): barmeth=method.split('barprop_')[-1]
            else: raise Exception('Poor choice of method naming conventions: {}'.format(method))
            rmax=self.limits[family][galaxy]['R'][1]
            inpar=dict(galaxy=galaxy,family=family,nbins=100,rscl='log',cmeth='pot')
            rlimAm2=(2.*rmax/100.,2.*rmax)
            rlimEll=(2.*rmax/inpar['nbins'],rmax)
            # disk_galaxy1_logr0p3to30_Center_phi_phi10_100bins
            if   barmeth=='Am2':
                ftype='scf'
                inpar.update(maxmeth='phi',phierr=10.,rlim=rlimAm2)
            # disk_galaxy1_logr0p3to30_Center_phi10_100bins
            # disk_galaxy1_logr0p3to30_Center_phi0p1_100bins
            elif barmeth=='derAm2':
                ftype='scf'
                inpar.update(phierr=10.,rlim=rlimAm2)
            #
            elif barmeth=='rot': pass
            # disk_galaxy1_logr0p1to30_potCenter_ell0p1_phi10_100bins_10ellp
            # disk_galaxy1_logr0p1to30_potCenter_ell0p1_phi10_100bins_100ellp
            elif barmeth=='ellipse':
                inpar.update(phierr=10.,ellerr=0.1,nellip=20,ellmin=0.25,rlim=rlimEll,cutmeth='err')
            elif barmeth=='scfAm2':
                ftype='scf'
                inpar.update(phierr=10.,ampmin=0.04,rlim=rlimAm2)
            # Error
            else: raise Exception('Invalid barmeth: {}'.format(barmeth))
            # Adjust rlim for small galaxy
            # if self['runtag']=='MW0B7_MW0B6_Q10RP0p1E1I0' and galaxy=='galaxy2':
            #     if barmeth=='ellipse': inpar.update(rlim=(0.1, 5.0))
            #     else                 : inpar.update(rlim=(1.0,10.0))
        # Error
        else: raise Exception('Invalid method: {}'.format(method))
        # Extension
        if ftype=='scf':
            cmeth=''
            try              : plotfext='{:03d}_{}{}'.format(int(float(fext)),family,galaxy[-1])
            except ValueError: 
                if '_' in fext: plotfext=fext
                else          : plotfext='{}_{}{}'.format(fext,family,galaxy[-1])
        # Bar plots
        if method.startswith('bar_'):
            # Arguments
            kwargs=dict(overwrite=overwrite,owplot=overwrite,wanim=True,
                        plotsnap=True,askuser=True,inpar=inpar,fext0=plotfext)
            # Plot
            out=self.barprop(barmeth,**kwargs)
        # Other plots
        else:
            pprint.pprint(inpar)
            # Arguments
            kwargs=dict(pmeth=pmeth,ftype=ftype,askuser=True,
                        overwrite=overwrite,owdata=owdata,owhist=owhist,owanim=owanim,
                        plotkw=plotkw,fext=plotfext,inpar=inpar,allfam=False,allgal=False,
                        histmeth=inpar.get('histmeth',None),plotmeth=inpar.get('plotmeth',None))
            # Plot
            if   itermeth=='allsnaps': out=self.allsnaps(method='plot',**kwargs)
            elif itermeth=='itersnap': out=self.itersnap(method='plotsnap',**kwargs)
            else: raise Exception('Invalid itermeth: {}'.format(itermeth))
        # Return control
        return out

####################################################################################################################################
####################################################################################################################################
# WALK THROUGH METHOD
# def walk(mtype=None,method=None,runtag=None,verbose=None,**extra_kw):
#     """
#     Method to walk users through performing different run methods
#     """
#     # Set constants
#     simtypLIST=simlist.LIST_RUNTYPS
#     # Get method type if not provided
#     mtypeLIST=simlist.LIST_METHTYPS
#     if mtype not in mtypeLIST: mtype=mmlio.askselect('Select a type of method:',mtypeLIST)
#     # Get method if not provided
#     methodLIST=simlist.DICT_METHODS[mtype]
#     if method not in methodLIST: method=mmlio.askselect('Select a {} method:'.format(mtype),methodLIST)
#     method=method.lower()
#     # Proceed based on method
#     if   mtype=='buildgal' and method=='mk_ic': mmlbuildgal.mk_ic(runtag=runtag,**extra_kw)
#     elif mtype=='gadget'   and method=='mk_ic': mmlgadget.mk_ic(runtag=runtag,**extra_kw)
#     else:
#         # Get run info
#         mmlio.verbose('Getting info on the run...')
#         simobj=loadsimobj(runtag)
#         # Call run
#         mmlio.verbose('Beginning {} for {} simulation...'.format(method.upper(),simobj['runtyp'].upper()),border=True)
#         output=simobj.run(mtype,method,verbose=verbose,**extra_kw)
#         print [method,simobj['runtyp'],simobj['runtag']]
#         if output==None: return
#         else           : return output

# ####################################################################################################################################
# # GENERAL METHODS
# def general(simstr,method,**method_kw):
#     """
#     Handles running different types of methods for runs
#     """
#     # Set constants
#     typListDEF=simlist.LIST_PTYPBASE
#     gallist=['both','primary','secondary']
#     typListVisible=simlist.LIST_PTYPVISI
#     # Pars input
#     method=mmlpars.mml_pars(method,type=str)
#     method=method.lower()
#     # Initialize output
#     out=None
#     # Find correct method subclass
#     # Return mmlsim object
#     if   method == 'getsimobj': return simstr
#     # Return pNbody object
#     elif method == 'getpnbody':
#         method_kw['askuser']=True
# #        out=simstr.get_pnbody(**method_kw)
#     # Print info on simulation
#     elif method == 'printinfo':
#         M1_vir=simstr.get_Mvir(galaxy=1)
#         R1_vir=simstr.get_Rvir(galaxy=1)
#         print 'SIMSTR:'
#         print '    M1_vir={} Msol'.format(M1_vir)
#         print '    R1_vir={} kpc'.format(R1_vir)
#         if simstr['runtyp']=='intsim':
#             M2_vir=simstr.get_Mvir(galaxy=2)
#             R2_vir=simstr.get_Rvir(galaxy=2)
#             rperi=(R1_vir+R2_vir)*simstr.infodict['rperi']
#             print '    M2_vir={} Msol'.format(M2_vir)
#             print '    R2_vir={} kpc'.format(R2_vir)
#             print '    Rperi={} kpc'.format(rperi)
#         fext=mmlio.askquest('What snapshot extension should be loaded?',default='ic',dtype='str')
# #        pnObj=simstr.get_pnbody(ftype='gadget',fext=fext)
# #        pnObj.printinfo()
# #        out=pnObj
#     # Check if particles are repeated
#     elif method == 'checkrepeat':
# #        pnobj=general(simstr,'getpnbody')
# #        out=pnobj.countrepeats()
#         print 'N repeats = {}'.format(out)
#     # Get profile info
#     elif method=='fitprofiles':
#         pnobj=general(simstr,'getpnbody')
#         out={}
#         for ityp in pnobj.get_typlist():
#             out[ityp]=pnobj.fitprofile(ptyp=ityp,**method_kw)
#             print '{} prof: {}'.format(ityp,out[ityp])
#     else: raise Exception('Invalid general method: {}'.format(method))
#     # Return output
#     return out

# if __name__ == '__main__':
#     walk()

        
