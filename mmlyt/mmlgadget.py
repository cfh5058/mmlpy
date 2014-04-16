#!/usr/bin/python
####################################################################################################################################
#
# MEAGAN LANG'S GADGET METHODS
#
####################################################################################################################################
import sys,os,shutil,glob,copy,pprint,scipy,math,struct,fnmatch
import numpy as np
import matplotlib as mplib

LIST_PTYPS=['gas','halo','disk','bulge','stars','bndry']
N_PTYP=len(LIST_PTYPS)
LIST_BLOCKS_STD=['HEAD','POS','VEL','ID','MASS']
LIST_BLOCKS_SPH=['U','RHO','HSML']
LIST_BLOCKS_OPT=['POT','ACCE','ENDT','TSTP']
LIST_BLOCKS=LIST_BLOCKS_STD+LIST_BLOCKS_SPH+LIST_BLOCKS_OPT
MAP_BLOCK2FIELD={'HEAD':'header',
                 'POS' :'position',
                 'VEL' :'velocity',
                 'ID'  :'id',
                 'MASS':'mass',
                 'U'   :'entropy',
                 'RHO' :'density',
                 'HSML':'smoothingLength',
                 'POT' :'potential',
                 'ACCE':'acceleration',
                 'ENDT':'entropyChange',
                 'TSTP':'timeStep'}
MAP_FIELD2BLOCK={ival:ikey for ikey,ival in MAP_BLOCK2FIELD.iteritems()}
                 

LIST_METHODS=['mk_ic','mk_files','setup','retrieve','clean']
LIST_METH_MKIC=['buildgal','snapshot','restart','interact']
from mmlutils import *
import main as mmlyt
import snapshot,simfile

FILELIST_SETUP=['exec','input','ic']
FILELIST_OUTPT=['snapshot','restart','output']
FILELIST_CLEAN=['snapshot','restart']

FEXT_MAKEFILE=['Makefile','makefile','gdmk']
FEXT_PARMFILE=['param','gdpar']

####################################################################################################################################
# GADGET SNAPSHOT CLASS
class NbodySnap_gadget(snapshot.NbodySnap):
    @staticmethod 
    def _can_load(f): return test_snap(f)
    def __init__(self,fname=None,parmfile=None,makefile=None,**exkw):
        # Initialize things
        super(NbodySnap_gadget,self).__init__()
        self._filename=fname
        self._files=[]
        # Initialize from a file
        if self.__class__._can_load(fname):
            # Check for multiple files
            try: fd=open(fname)
            except IOError:
                fd=open(fname+'.0')
                fname+='.0'
            fd.close()
            if fname.endswith('.0'): self._filename=fname.split('.0')[0]
            # Read files
            self._files.append(rw_snap('R',fname,makefile=makefile))
            Nfiles=self._files[0]['HEAD']['Nfiles']
            npart=np.array(self._files[0]['HEAD']['npart'])
            for i in range(1,Nfiles):
                ifname=fname[:-1]+str(i)
                self._files.append(rw_snap('R',fname,makefile=makefile))
                npart+=np.array(self._files[i]['HEAD']['npart'])
            # Setup things 
            self._num_particles=npart.sum()
            self.header=copy.deepcopy(self._files[0]['HEAD'])
            self.header['npart']=npart
            # Initialize stuff
            self._ptype_slice={}
            self._loadable_keys = set([])
            self._ptype_keys=set([])
            self._arrays={}
            self.properties={}
            # Type slices
            for i,ityp in enumerate(LIST_PTYPS):
                self._ptype_slice[ityp]=slice(npart[0:i].sum(),npart[0:(i+1)].sum())
            # Loadable ekys
            for f in self._files:
                self._loadable_keys=self._loadable_keys.union(set(f.keys()))
            self._loadable_keys=[MAP_BLOCK2FIELD[ikey] for ikey in self._loadable_keys]
            
                
        # Initialize from arrays
        else: pass
        
    @staticmethod 
    def _load(fname,parmfile=None,makefile=None):
        if fname==None and hasattr(self,'fname'): fname=self.fname
        if parmfile==None and hasattr(self,'parmfile'): parmfile=self.parmfile
        if makefile==None and hasattr(self,'makefile'): makefile=self.makefile
        if fname==None: raise Exception('File name not provided or initialized')
        # Get file names
        if not isinstance(parmfile,str): parmfile=simfile.search_file(self.fname,FEXT_PARMFILE)
        if not isinstance(makefile,str): makefile=simfile.search_file(self.fname,FEXT_MAKEFILE)
        setattr(self,'parmfile',parmfile)
        setattr(self,'makefile',makefile)
        # Load dictionary of blocks
        nb=rw_snap('R',fname,makefile=makefile)
    @staticmethod 
    def _save(fname,overwrite=False):
        if fname==None: raise Exception('File name not provided or initialized')
        

####################################################################################################################################
####################################################################################################################################
# HIGH LEVEL METHODS REQUIRING MMLSIM
def main(): return mmlnbody.walk(mtype='gadget')

####################################################################################################################################
# METHOD FOR RUNNING DIFFERENT GADGET OPERATIONS
def run(simstr,method,**method_kw):
    """
    Provides interface for running different GADGET operations
    """
    # Set constants
    methLIST=LIST_METHODS
    # Pars input
    method=mmlpars.mml_pars(method,list=methLIST)
    # Initialize default output
    out=None
    # Proceed based on method
    if   method=='mk_ic'     : mk_ic(simstr,**method_kw)
    elif method=='mk_files'  : mk_files(simstr,**method_kw)
    else: raise Exception('Invalid method: {}'.format(method))
    # Return output
    return out

####################################################################################################################################
####################################################################################################################################
# METHODS TO CREATE GADGET FILES

####################################################################################################################################
# METHOD TO CREATE LOCAL COPIES OF GADGET FILES
def mk_files(simstr,verbose=False,askuser=False,**extra_kw):
    """
    Creates local copies of GADGET files
    """
    if verbose: mmlio.verbose('Creating files for a GADGET run locally...',border=True)
    # Get file names
    files=simstr.fdict
    mmlfiles.mkdirs(simstr.fdict['gadget']['input']['dir'])
    # Executable file
    if mmlio.yorn('Create new local copy of executable?'): mk_exec(simstr)
    # Parameter files
    if mmlio.yorn('Create GADGET submission files?'):
        owparam=mmlio.yorn('Overwrite existing files?')
        pardict=mk_param(simstr,overwrite=owparam)
        pbsdict=mk_pbs(simstr,param=pardict,overwrite=owparam)
        if pardict['OutputListOn']==1:
            outlist=mk_outlist(files['gadget']['input']['outlist'],pardict,overwrite=owparam,askuser=askuser)
    # IC file
    if mmlio.yorn('Create new local copy IC file?'      ): mk_ic(simstr,askuser=askuser)
    return
 
####################################################################################################################################
# METHOD TO CREATE GADGET EXECUTABLE FILES
def mk_exec(simstr):
    """
    Creates a local copy of executable files
    """
    # Select default files
    exdeft=simstr.finfo['gadget']['exec']['default']
    mfdeft=simstr.finfo['gadget']['makefile']['default']
    # Check existence
    if not os.path.isfile(exdeft): raise Exception('Executable does not exist: {}'.format(exdeft))
    if not os.path.isfile(mfdeft): raise Exception('Makefile does not exist: {}'.format(mfdeft))
    # Copy files to local directory
    mmlfiles.cpfiles(exdeft,simstr.fdict['gadget']['input']['exec'],overwrite=True)
    mmlfiles.cpfiles(mfdeft,simstr.fdict['gadget']['input']['makefile'],overwrite=True)
    # Return
    return

####################################################################################################################################
# METHOD TO CREATE GADGET PARAMETER FILE
def mk_param(simstr=None,runtag=None,runcomp='accre',overwrite=False,owmake=False,askuser=False):
    """
    Creates a GADGET parameter file and returns content as a dictionary
    """
    if simstr==None: simstr=mmlyt.runtag2simobj(runtag)
    # If the static parameter file dosn't exist, create it
    if overwrite or not os.path.isfile(simstr.fdict['gadget']['icgen']['stat_pm']):
        # Initialize make log
        mklog=mk_makelog(simstr,overwrite=owmake)
        method=mklog['method']
        # Initialize static parameter file
        params=rw_param('R',simstr.finfo['stat_pm']['default'])
        # Get file names
        files_accre=simstr.mkfiledict(memcomp=runcomp,checkaccess=False)
        # Set file names
        params['OutputDir'         ]=files_accre['gadget']['output']['dir']
        params['SnapshotFileBase'  ]=os.path.basename(files_accre['gadget']['output']['snapbase'])
        params['InitCondFile'      ]=files_accre['gadget']['input']['ic']
        params['RestartFile'       ]=os.path.basename(files_accre['gadget']['output']['restbase'])
        params['InfoFile'          ]=os.path.basename(files_accre['gadget']['output']['info'])
        params['TimingsFile'       ]=os.path.basename(files_accre['gadget']['output']['timings'])
        params['CpuFile'           ]=os.path.basename(files_accre['gadget']['output']['cpu'])
        params['EnergyFile'        ]=os.path.basename(files_accre['gadget']['output']['energy'])
        params['ResubmitCommand'   ]=files_accre['gadget']['input']['resub']
        params['OutputListFilename']=files_accre['gadget']['input']['outlist']
        # Set communication buffer
        params['BufferSize'        ]=mmlsim2combuff(simstr)
        # Set softenings
        for ityp in LIST_PTYPS:
            params['Softening'+ityp.capitalize()]=mklog['softs'][ityp]
            params['Softening'+ityp.capitalize()+'MaxPhys']=mklog['softs'][ityp]
        # Save parameters
        rw_param('W',simstr.fdict['gadget']['icgen']['stat_pm'],params,overwrite=overwrite)
        # Check parameter file
        print simstr.fdict['gadget']['icgen']['stat_pm']
        mmlio.yorn('Make changes to the above parameter file now!')
    # Make local copy of parameter file
    mmlfiles.cpfiles(simstr.fdict['gadget']['icgen']['stat_pm'],simstr.fdict['gadget']['input']['param'],overwrite=True)
    # Return output
    return params

####################################################################################################################################
# METHOD TO CREATE PBS SUBMISSION SCRIPT
def mk_pbs(simstr,param=None,overwrite=False):
    """
    Creates a GADGET pbs script and returns content as a dictionary
    """
    # Get file names
    files_accre=simstr.mkfiledict(memcomp='accre',checkaccess=False)
    pbsfile=simstr.fdict['gadget']['input']['pbs']
    if param==None: param=simstr.fdict['gadget']['input']['param']
    # Load existing file if overwrite False
    if not overwrite and os.path.isfile(pbsfile): 
        mmlio.verbose('Loading existing GADGET pbs script:',addval=pbsfile)
        pbsdict=mmlio.rw_pbs('R',pbsfile)
    # Otherwise create file from scratch
    else:
        mmlio.verbose('Creating new GADGET pbs script:',addval=pbsfile)
        # Get parameter file info
        if not isinstance(param,dict):
            param=mmlpars.mml_pars(param,default=simstr.fdict['gadget']['input']['param'],type=str)
            param=rw_param('R',param)
        tcpu=param['TimeLimitCPU']/3600
        mcpu=mmlsim2mcpu(simstr)
        nproc=simstr['nprocev']
        exefile=os.path.basename(files_accre['gadget']['input']['exec'])
        parfile=os.path.basename(files_accre['gadget']['input']['param'])
        # Variables from simstr
        pbsdict={'jobname' :simstr['runtag'],
                 'ppn'     :2,
                 'nodes'   :nproc/ppn,
                 'pmem'    :mcpu,
                 'mem'     :mcpu*nproc,
                 'walltime':tcpu,
                 'outfile' :files_accre['gadget']['output']['runout']}
        pbsdict['cmdlist']=['mpiexec -np {} ./{} {} 0'.format(nproc,exefile,parfile)]
        # Write PBS file
        mmlio.rw_pbs('W',pbsfile,pbsdict,overwrite=overwrite)
    # Return output
    return pbsdict

####################################################################################################################################
# METHODS TO CREATE LIST OF OUTPUT TIMES
def mk_outlist(fname,param,method=None,overwrite=False,askuser=False):
    """
    Creates a file containing a list of output times
    """
    methLIST=['lin','log']
    # Pars input
    fname=mmlpars.mml_pars(fname,type=str)
    param=mmlpars.mml_pars(param,type=dict)
    if method not in methLIST and askuser:
        method=mmlio.askselect('How should the timing of output be scaled?',methLIST)
    else:
        method=mmlpars.mml_pars(method,list=methLIST,default='lin')
    # Recover parameters from param
    nout=int((param['TimeMax']-param['TimeBegin'])/param['TimeBetSnapshot'])
    tbeg=param['TimeBegin']
    tend=param['TimeMax']
    # Create outlist
    if   method=='lin': outlist=list(np.linspace(tbeg,tend,nout))
    elif method=='log': outlist=list(np.logspace(np.log10(tbeg),np.log10(tend),nout))
    else: raise Exception('Invalid method: {}'.format(method))
    # Return output
    return outlist

####################################################################################################################################
# METHOD TO CREATE IC FILES
def mk_ic(simstr=None,runtag=None,askuser=False,overwrite=False,owmake=False):
    """
    Creates IC files
    """
    if simstr==None: simstr=mmlyt.runtag2simobj(runtag)
    # If the static IC file dosn't exist, create it
    if overwrite or not os.path.isfile(simstr.fdict['gadget']['stat_ic']):
        # Initialize make log
        mklog=mk_makelog(simstr,overwrite=owmake)
        method=mklog['method']
        # Create IC file
        if   method=='buildgal': 
            # Load buildgal IC file
            pnb=simstr.loadsnap(simstr.fdict['buildgal']['icgen']['stat_ic'],ftype='buildgal')
            # Convert units
            iunitfile=simstr.fdict['buildgal']['input']['units'  ]
            funitfile=simstr.fdict['gadget'  ]['icgen']['stat_pm']
            pnb.convert_units(oldunit=iunitfile,newunit=funitfile)
            # Reset time
            pnb['time']=0.
            mklog['time0']=pnb['time']
        elif method=='snapshot':
            isimstr=mmlyt.runtag2simobj(mklog['baserun'])
            # Load snapshot file
            pnb=isimstr.loadsnap(fext=mklog['snapext'],ftype='gadget')
            mklog['snaptime']=pnb['time']
            # Reset time
            if mklog['resettime']: pnb['time']=0.
            mklog['time0']=pnb['time']
            # Heat disk
            pnb.heatdisk(mklog['heatfact'])
        elif method=='restart' : raise Exception('Not currently supported.')
        elif method=='interact':
            # Load snapshots
            isimstr1=mmlyt.runtag2simobj(simstr.infodict['gal1'])
            isimstr2=mmlyt.runtag2simobj(simstr.infodict['gal2'])
            pnb1=isimstr1.loadsnap(fext=mklog['snapext1'],ftype='gadget')
            pnb2=isimstr2.loadsnap(fext=mklog['snapext2'],ftype='gadget')
            mklog['snaptime1']=pnb1['time']
            mklog['snaptime2']=pnb2['time']
            # Center galaxies
            pnb1.center(move=True)
            pnb2.center(move=True)
            # Get total mass
            M1=pnb1['mass_tot']*pnb1.convertionFactorTo('Msol')
            M2=pnb2['mass_tot']*pnb2.convertionFactorTo('Msol')
            R1=None ; R2=None
            # Add orbit to second galaxy
            pnb2.add_orbit(M1,M2,R1_pc=R1,R2_pc=R2,**mklog['mkinfo'])
            # Add galaxies
            pnb=pnb1+pnb2
            # Reset time
            pnb['time']=0.0
            mklog['time0']=pnb['time']
        else: raise Exception('Invalid method: {}'.format(method))
        # Save gadget IC file
        pnb.savesnap(simstr.fdict['gadget']['icgen']['stat_ic'],ftype='gadget',flag_ic=True,overwrite=overwrite)
        del pnb
    # Make local copy of IC file
    mmlfiles.cpfiles(simstr.fdict['gadget']['icgen']['stat_ic'],simstr.fdict['gadget']['input']['ic'],overwrite=True)
    # Return
    return

####################################################################################################################################
# METHOD TO CREATE MAKELOG
def mk_makelog(simstr,method=None,overwrite=False):
    """
    Creates a make log
    """
    fmklog=simstr.fdict['gadget']['icgen']['stat_mklog']
    # Read make log
    if os.path.isfile(fmklog) and not overwrite: mklog=simfile.rw_makelog('R',fmklog)
    # Create it
    else:
        # Pars input
        if method not in LIST_METH_MKIC: method=mmlio.askselect('Pick a method for making ICs',LIST_METH_MKIC)
        # Initialize makelog 
        mklog=dict(ftype ='gadget',
                   method=method,
                   runtag=simstr['runtag'],
                   runtyp=simstr['runtyp'],
                   addtag=simstr.infodict['addtag'],
                   ntot  =simstr.get_ntot())
        if   method=='buildgal': mklog.update(baserun=simstr['runtag'],
                                              nproc  =simstr['nprocmk'])
        elif method=='snapshot': mklog.update(baserun=simstr['runtag'].split(mklog['addtag']))
        elif method=='restart' : mklog.update(nproc  =simstr['nprocev'])
        elif method=='interact': mklog.update(baserun=simstr['runtag'])
        else: raise Exception('Invalid method: {}'.format(method))
        simstr0=mmlyt.runtag2simboj(mklog['baserun'])
        mklog['mkinfo']=simstr0.infodict
        # Set defaults
        mklogDEF={}
        if   method=='buildgal': 
            mdict=mmlbuildgal.model2param(simstr0.infodict)
            mklogDEF.update(softs={ityp:mdict.get(ityp+'_eps',0.)/1000. for ityp in LIST_PTYPS})
        elif method=='snapshot': mklog.update(softs=rw_softenings(simstr0['gadget']['icgen']['stat_pm']))
        elif method=='interact':
            mklogDEF.update(softs={})
            simstr1=mmlyt.runtag2simobj(simstr.infodict['gal1'])
            simstr2=mmlyt.runtag2simboj(simstr.infodict['gal2'])
            sdict1=rw_softenings(simstr1['gadget']['icgen']['stat_pm'])
            sdict2=rw_softenings(simstr2['gadget']['icgen']['stat_pm'])
            for ityp in LIST_PTYPS:
                if sdict1[ityp]==sdict2[ityp]: mklogDEF['softs'][ityp]=sdict1[ityp]
                else: raise Exception('The softening lengths are not equal. Decide how to treat this.')
        else: raise Exception('Invalid method: {}'.format(method))
        # Pars makelog
        mklog=mmlparam.parspar('mmlyt.simfile','mklog_{}_{}'.format('gadget',method),inpar=mklog,defpar=mklogDEF,
                               askuser=True,load=False,save=False)
        # Save to file
        simfile.rw_makelog('W',fmklog,mklog,overwrite=overwrite)
    # Return 
    return mklog

####################################################################################################################################
####################################################################################################################################
# METHOD TO TEST FILES
def test_snap(fname):
    """
    Tests validity of a snapshot
    """
    # Try to open file
    try: fd=open(fname)
    except:
        # Try to open first in a set with file extension
        try: fd=open(fname+'.0')
        except: return False
    # Get fileinfo
    out=info_snap(fd)['valid']
    # Close file & return output
    fd.close()
    return out
def info_snap(fd):
    """
    Gets info on the format of the snapshot
    """
    # Initalize output
    out={'valid':False,'version':'None','byteswap':False,'endian':'='}
    # Read first 4 bytes & close file
    (r,) = struct.unpack('=I',fd.read(4))
    # Interpret first integer
    #     Gadget 1: 256 (or byteswapped equivalent)
    #     Gadget 2: 8   (or byteswapped equivalent)
    v={'Gadget1':np.int32(256),'Gadget2':np.int32(8)}
    for iv in v:
        if   r==v[iv]           : out.update(valid=True,version=iv)
        elif r==v[iv].byteswap(): out.update(valid=True,version=iv,byteswap=True)
    # Get endianness of file
    if out['byteswap']:
        if sys.byteorder=='little': out['endian']='>'
        else                      : out['endian']='<'
    # Reset file position & return info
    fd.seek(0,0)
    return out

####################################################################################################################################
# METHOD TO RETURN BLOCK LIST FOR A SNAPSHOT
def get_blocklist(fname,makefile=None):
    """
    Returns a list of blocks for version 1 snapshots
    """
    # Initialize blocklist
    blocklist=copy.deepcopy(LIST_BLOCKS)
    mkopt2block=dict(OUTPUTPOTENTIAL      ='POT' ,
                     OUTPUTACCELERATION   ='ACCE',
                     OUTPUTCHANGEOFENTROPY='ENDT',
                     OUTPUTTIMESTEP       ='TSTP')
    # Check for makefile
    if not isinstance(makefile,str): makefile=simfile.search_file(fname,FEXT_MAKEFILE)
    # Read options from a snapshot
    if isinstance(makefile,str):
        makeopt=rw_makefile('R',makefile)
        for iopt in mkopt2block:
            if not makeopt.has_key(iopt): makeopt[iopt]=False
            if not makeopt[iopt]: blocklist.remove(mkopt2block[iopt])
    # Otherwise remove all optional blocks
    else:
        for iopt in LIST_BLOCKS_OPT: blocklist.remove(iopt)
    return blocklist

####################################################################################################################################
####################################################################################################################################
# METHODS TO READ/WRITE FILES

####################################################################################################################################
# METHOD TO READ/WRITE SNAPSHOTS
def rw_snap(rwid,fname,fdict=None,overwrite=False,makefile=None):
    """
    Reads/writes snapshot files
    """
    if not test_snap(fname): raise Exception('Invalid file type.')
    # Read
    if rwid=='R':
        # Open file
        fd=open(fname)
        # Get info from first 4 bytes
        info=info_snap(fd)
        # Get block list for version 1 snapshots
        if info['version']=='Gadget1': 
            blockidx=0
            blocklist=get_blocklist(fname,makefile=makefile)
        # Continue loading blocks until entire file read
        outdict={}
        while True:
            # Read name & size from block header
            if info['version']=='Gadget2':
                head=fd.read(5*4)
                if len(head)!=5*4: name,blocklen="    ",0
                else:
                    head=struct.unpack(info['endian']+'I4sIII',head)
                    if head[0] != 8 or head[3] != 8 or head[4] != head[2]-8 :
                        raise IOError, "Corrupt header record. Possibly incorrect file format"
                    name,blocklen=head[1],head[2]-8
            # Read size from block header & get name from list
            elif info['version']=='Gadget1':
                head=fd.read(4)
                if len(head)!=4: name,blocklen="    ",0
                else:
                    (blocklen,)=struct.unpack(info['endian']+'I',head)
                    try: name=blocklist[blockidx]
                    except IndexError:
                        warnings.warn("Run out of block names. Using fallbacks: UNK*",RuntimeWarning)
                        name='UNK'+str(blockidx-len(blocklist))
                blockidx+=1
            else: raise IOError,'Unsupprted snapshot version: {}'.format(info['version'])
            tag=name.strip()
            # Break if block is empty
            if blocklen==0: break
            # Header block
            if tag=='HEAD':
                if blocklen != 256: raise IOError, "Mis-sized {} block in {}".format(tag,fname)
                outdict[tag]=rw_snaphead(fd,info['endian'])
                # Update block list based on particles
                if outdict[tag]['npart'][0]==0:
                    for isph in LIST_BLOCKS_SPH: 
                        blocklist.remove(isph)
                if outdict[tag]['nmass'].sum()==0: blocklist.remove('MASS')
            # Body block
            else:
                if   tag=='MASS'           : t_part=outdict['HEAD']['nmass'].sum()
                elif tag in LIST_BLOCKS_SPH: t_part=outdict['HEAD']['npart'][0]
                else                       : t_part=outdict['HEAD']['npart'].sum()
                datatyp=False
                partlen=blocklen/t_part
                if tag=='ID':
                    if   partlen==4: datatyp=np.int32
                    elif partlen==8: datatyp=np.int64
                elif tag in ['POS','VEL','ACCE']:
                    #partlen/=3
                    if   partlen==12: datatyp=np.float32
                    elif partlen==24: datatyp=np.float64
                else:
                    if   partlen==4: datatyp=np.float32
                    elif partlen==8: datatyp=np.float64
                if isinstance(datatyp,bool): raise IOError,"No valid type for {} (partlen={})".format(tag,partlen)
                # Record position & move to end of block
                blockstart=fd.tell()
                extra_len=t_part*partlen
                if extra_len >= 2**32: blocklen=extra_len
                fd.seek(blocklen,1)
                # Add block to output dictionary
                outdict[tag]=dict(start=blockstart,
                                  length=blocklen,
                                  partlen=partlen,
                                  datatyp=datatyp)
            # Read block footer
            foot=fd.read(4)
            if len(foot)!=4: raise IOError, "Could not read block footer"
            (footlen,)=struct.unpack(info['endian']+'I',foot)
            if footlen != blocklen: raise IOError, 'Bad record size for {} in {}'.format(name,fname)
        # Close file
        fd.close()
        # Add mass block
        if 'MASS' not in outdict:
            outdict['MASS']=dict(start=0,length=0,partlen=8,datatyp=np.float64)
        # Return dictionary
        return outdict
    # Write
    elif rwid=='W':
        pass
    # Error
    else: raise Exception('Invalid rwid: {}'.format(rwid))
    return


####################################################################################################################################
# METHOD TO READ/WRITE SNAPSHOT HEADERS
def rw_snaphead(fd,endian='='):
    """
    Reads/writes snapshot header
    """
    # Read data
    if fd.mode=='r':
        data=fd.read(256)
        if len(data) != 256: raise IOError, "Could not read HEAD block in "+fname
        # Construct header structure
        h=dict(npart    = np.zeros(N_PTYP,dtype=np.uint32),
               massarr  = np.zeros(N_PTYP),
               time=0., redshift = 0.,
               flag_sfr=0,flag_feedback=0,
               npartAll = np.zeros(N_PTYP,dtype=np.uint32),
               flag_cooling=0, Nfiles=0,
               BoxSize=0., Omega0=0., OmegaLambda=0., HubbleParam=0.,
               flag_stellarage=0, flag_metals=0,
               NallHW   = np.zeros(N_PTYP,dtype=np.uint32),
               flag_entropy_instead_u=0, flag_doubleprecision=0, 
               flag_ic_info=0, lpt_scalingfactor=0., fill=0)
        # Return if empty
        if data=='': return
        # Construct format string
        fmt=endian+"IIIIIIddddddddiiIIIIIIiiddddiiIIIIIIiiif48s"
        if struct.calcsize(fmt)!=256: raise Exception('Bug: the header format string is not 256 bytes')
        # Unpack data
        (h['npart'][0],h['npart'][1],h['npart'][2],h['npart'][3],h['npart'][4],h['npart'][5],
         h['massarr'][0],h['massarr'][1],h['massarr'][2],h['massarr'][3],h['massarr'][4],h['massarr'][5],
         h['time'],h['redshift'],h['flag_sfr'],h['flag_feedback'],
         h['npartAll'][0],h['npartAll'][1],h['npartAll'][2],h['npartAll'][3],h['npartAll'][4],h['npartAll'][5],
         h['flag_cooling'],h['Nfiles'],h['BoxSize'],h['Omega0'],h['OmegaLambda'],h['HubbleParam'],h['flag_stellarage'],h['flag_metals'],
         h['NallHW'][0],h['NallHW'][1],h['NallHW'][2],h['NallHW'][3],h['NallHW'][4],h['NallHW'][5],
         h['flag_entropy_instead_u'],h['flag_doubleprecision'],h['flag_ic_info'],h['lpt_scalingfactor'],h['fill'])=struct.unpack(fmt,data)
        # Add simply calculated values
        h['nbody']=h['npart'].sum()
        h['nmass']=copy.deepcopy(h['npart'])
        h['nmass'][h['massarr']!=0]=0.
        # Return header dictionary
        return h
    # Write
    elif fd.mode=='w': pass
    # Error
    else: raise Exception('Invalid mode: {}'.format(fd.mode))
    return

####################################################################################################################################
# METHOD TO READ/WRITE PARAMETER FILES
def rw_param(rwid,fname,fdict=None,overwrite=False):
    """
    Reads/writes parameter files
    """
    # Read
    if   rwid=='R': 
        fdict=mmlio.rwdict('R',fname)
        return fdict
    # Write
    elif rwid=='W': 
        mmlio.rwdict('W',fname,fdict,width=29,overwrite=overwrite)
        return
    # Error
    else: raise Exception('Invalid rwid: {}'.format(rwid))
    return

####################################################################################################################################
# METHOD TO READ/WRITE LIST OF OUTPUT TIMES
def rw_outlist(rwid,fname,outlist=None,overwrite=None):
    """
    Reads/writes a list of output times to file
    """
    # Read
    if   rwid=='R': 
        if not os.path.isfile(fname): raise Exception('File name is not valid: {}'.format(fname))
        fid=open(fname,'r')
        outlist=[]
        for line in fid: outlist.append(float(line.strip()))
        fid.close()
    # Write
    elif rwid=='W':
        outlist=mmlpars.mml_pars(outlist,type=list)
        if os.path.isfile(fname) and not overwrite: return outlist
        fid=open(fname,'w')
        for itime in outlist: fid.write(repr(itime)+'\n')
        fid.close()
    # Return output
    return outlist

####################################################################################################################################
# METHOD TO READ SOFTENINGS FROM GADGET PARAMETER FILE
def rw_softenings(fname):
    """
    Reads and returns softening lengths from a GADGET parameter file as a dictionary with type keys
    """
    # Read parameters from parameter file
    params=rw_param('R',fname)
    # Add softenings to dictionary
    soft={ityp:params['Softening'+ityp.capitalize(ityp)] for ityp in LIST_PTYPS}
    # Return dictionary of softening lengths
    return soft

####################################################################################################################################
# METHOD TO READ GADGET MAKE FILE
def rw_makefile(rwid,fname):
    """
    Reads options from a GADGET makefile
    """
    # Read
    if rwid=='R':
        # Check if file is valid
        if not os.path.isfile(fname):
            raise Exception('Makefile name is not valid: {}'.format(fname))
        # Open file
        fid=open(fname,'r')
        # Check each line for options
        makeopt={}
        for line in fid:
            # Handle uncommented options
            if line.startswith('OPT'):
                if line.startswith('OPTIMIZE') or line.startswith('OPTIONS'): continue
                vars=line.split('-D')
                keyval=vars[-1].split('=')
                if len(keyval) == 2:
                    key=keyval[0]
                    val=mmlstring.str2val(keyval[1].split('\n')[0])
                else:
                    key=keyval[0].split('\n')[0]
                    val=True
                makeopt[key]=val
            # Handle commented options
            elif line.startswith('#OPT'):
                if line.startswith('#OPTIMIZE') or line.startswith('#OPTIONS'): continue
                vars=line.split('-D')
                keyval=vars[-1].split('=')
                if len(keyval) == 2:
                    key=keyval[0]
                else:
                    key=keyval[0].split('\n')[0]
                val=False
                makeopt[key]=val
            else:
                continue
        fid.close()
        return makeopt
    # Write
    elif rwid=='W': raise Exception('Why do you want to make a Makefile?')
    # Error
    else: raise Exception('Invalid rwid: {}'.format(rwid))
    return

####################################################################################################################################
# METHOD TO READ/PLOT ENERGY FILES
def rw_energy(fname):
    """
    Reads data from gadget energy files
    """
    if not os.path.isfile(fname):
        raise Exception('Energy file does not exist: {}'.format(fname))
    # Initialize
    colKeys=['t','ei','ep','ek',
             'eiGas','epGas','ekGas',
             'eiHalo','epHalo','ekHalo',
             'eiDisk','epDisk','ekDisk',
             'eiBulge','epBulge','ekBulge',
             'eiStars','epStars','ekStars',
             'eiBndry','epBndry','ekBndry',
             'massGas','massHalo','massDisk',
             'massBulge','massStars','massBndry']
    edict={ikey:[] for ikey in colKeys}
    # Open file
    fid=open(fname,'r')
    # Read input
    for iline in fid:
        colList=iline.split() ; ncols=len(colList)
        if ncols != len(colKeys):
            raise Exception('The number of columns ({}) is not the number expected ({}).'.format(ncols,28))
        for icol in range(ncols):
            edict[colKeys[icol]].append(float(colList[icol]))
    # Return output
    return edict


####################################################################################################################################
####################################################################################################################################
# METHODS FOR PLOTTING GADGET STUFF

####################################################################################################################################
# METHOD FOR PLOTTING GADGET ENERGY STATS
def plot_energy(fname,plotfile=None,overwrite=False,**plot_kw):
    '''
    NAME:
        mmlgadget.plot_energy
    PURPOSE:
        To plot data from a GADGET energy file
    CALLING:
        plot_energy(fname[,plotfile=None,overwrite=False,**plot_kw])
    ARGUMENTS:
        fname:     Str absolute path to GADGET energy file
    KEYWORDS:
        plotfile:  Str absolute path to file where plot should be saved
        overwrite: Bool specifying if existing plot should be overwritten
        Addtional keywords are passed on to the plotting routine (simplot.plot_energy)
    '''
    # Pars input
    if os.path.isfile(plotfile):
        if overwrite:
            os.remove(plotfile)
        else:
            print 'Plot of GADGET energy statistics already exists: {}'.format(plotfile)
            return
    # Read in data
    edict=read_energy(fname)
    # Plot
    t=np.array(edict['t'],float)
    ek=np.array(edict['ek'],float)
    ep=np.array(edict['ep'],float)
    simplot.plot_energy(t,ek,ep,fname=plotfile,overwrite=overwrite,**plot_kw)

####################################################################################################################################
####################################################################################################################################
# METHODS FOR CALCULATING THINGS

####################################################################################################################################
# METHOD TO DETERMINE COMMUNICATION BUFFER
def mmlsim2combuff(simstr):
    """
    Estimates the communication buffer a GADGET simulation will need
    """
    # Get file dictionary
    files=simstr.fdict
    icfile=files['icgen']['gdic']
    # Get buffer size from IC file if it exists
    if os.path.isfile(icfile):
        statinfo=os.stat(icfile)
        combuff=float(statinfo.st_size)/(1.0e6)
    # Otherwise get it by scaling previous results
    else:
        combuff=float(simstr.ntot)*(24.0/560000.0)
    # Return rounded output
    return int(mmlmath.oom(combuff,nsig=2,method='CEIL'))
    
####################################################################################################################################
# METHOD TO DETERMINE MEMORY FROM MMLSIM OBJECT
def mmlsim2mcpu(simstr):
    """
    Estimates the amount of memory a GADGET simulation will need
    """
    # Get file dictionary
    files=simstr.fdict
    parfile=files['icgen']['gdpm']
    # Get parameters from parameter file if it exists
    if os.path.isfile(parfile):
        params=mmlio.rwdict('R',parfile)
        partalloc=float(params['PartAllocFactor'])
        treealloc=float(params['TreeAllocFactor'])
        combuff=int(params['BufferSize'])
        perflag=params['PeriodicBoundariesOn']
        if perflag==1: periodic=True
        else:          periodic=False
    # Otherwise force defaults
    else:
        partalloc=None
        treealloc=None
        combuff=mmlsim2combuff(simstr)
        periodic=None
    # Calculate memory based on GADGET user manual
    mcpu=calcmcpu(simstr.ntot,nproc=simstr['nprocgd'],
                  partalloc=partalloc,treealloc=treealloc,
                  combuff=combuff,periodic=periodic)
    # Force reasonable values as a safety
    if   simstr['runtyp']=='galsim': mcpu=1000. # mb
    elif simstr['runtyp']=='intsim': mcpu=2000. # mb
    else: raise Exception('Invalid run type: {}'.format(simstr['runtyp']))
    # Return rounded output
    return mmlmath.oom(mcpu,nsig=1,method='CEIL')

####################################################################################################################################
# METHOD TO DETERMINE TIME REQUIRED TO RUN GADGET
def mmlsim2tcpu(simstr):
    """
    Estimates the amount of time required to run GADGET on ACCRE
    """
    # Force reasonable values as a safety
    tcpu=mmlio.askquest('What should be the overall time limit?',default=48.,dtype='float')
    # Return output
    return tcpu
    
####################################################################################################################################
# METHOD TO DETERMINE MEMORY REQUIRED TO RUN GADGET
def calcmcpu(ntot,nproc=None,nsph=None,partalloc=None,treealloc=None,combuff=None,
             periodic=None,double=None,treepm=None,pmgrid=None,zoom=None):
    '''
    NAME:
        mmlgadget.calcmcpu
    PURPOSE:
        To determine the amound of memory a GADGET simulation would need.
    CALLING:
        mcpu=calcmcpu(ntot[,nproc=,nsph=,partalloc=,treealloc=,combuff=,
                      periodic=,double=,treepm=,pmgrid=,zoom)
    ARGUMENTS:
        ntot:      Long number of particles
    KEYWORDS:
        nproc:     Int number of processors [DEFAULT=1]
        nsph:      Long number of SPH particles [DEFAULT=0]
        partalloc: Float particle allocation factor [DEFAULT=1.5]
        treealloc: Float tree allocation factor [DEFAULT=0.8]
        combuff:   Int ommunication buffer [DEFAULT=25]
        periodic:  Bool specifying if boundaries are periodic
        double:    Bool specifying if calculations are double floats
        treepm:    Bool specifying if simulation is particle-mesh
        pmgrid:    Int size of particle-mesh grid [DEFAULT=0]
        zoom:      Bool specifying if simulation is a zoom
    OUTPUT:
        mcpu:      Int memory required on each processor (MB)
    '''
    # Pars input
    ntot=mmlpars.mml_pars(ntot,type=long,min=0L)
    nproc=mmlpars.mml_pars(nproc,default=1,type=int,min=0)
    nsph=mmlpars.mml_pars(nsph,default=0L,type=long,min=0L)
    partalloc=mmlpars.mml_pars(partalloc,default=1.5,type=float,min=0)
    treealloc=mmlpars.mml_pars(treealloc,default=0.8,type=float,min=0)
    combuff=mmlpars.mml_pars(combuff,default=25,type=int,min=0)
    periodic=mmlpars.mml_pars(periodic,default=False,type=bool)
    double=mmlpars.mml_pars(double,default=False,type=bool)
    treepm=mmlpars.mml_pars(treepm,default=False,type=bool)
    pmgrid=mmlpars.mml_pars(pmgrid,default=0,type=int)
    zoom=mmlpars.mml_pars(zoom,default=False,type=bool)
    # Additional factors
    if double:
        dbl=2.0
    else:
        dbl=1.0
    # General particle stuff
    M1=dbl*float(ntot)*partalloc*(68.0+64.0*treealloc)
    # SPH particle stuff
    M2=dbl*float(nsph)*partalloc*84.0
    # Communication buffer
    M3=float(nproc)*float(combuff)*(1024.0**2)
    # TreePM stuff
    if treepm:
        if periodic:
            M4=dbl*((float(pmgrid)**3)*12.0-16.0)
            if zoom:
                M5=dbl*((2.0*float(pmgrid))**3)*16.0
            else:
                M5=0.0
        else:
            M4=0.0
            if zoom:
                M5=dbl*((2.0*float(pmgrid))**3)*20.0
            else:
                M5=dbl*((2.0*float(pmgrid))**3)*16.0
    else:
        M4=0.0
        M5=0.0
    # Summ memory & return in MB
    mcpu=(M1+M2+M3+M4+M5)/(1024.0**2)
    # Divide accros processors
    return mcpu/float(nproc)


####################################################################################################################################
####################################################################################################################################
# PROVIDE COMMAND LINE ACCESS
if __name__ == '__main__': main()

def listpar(partyp,**exkw):
    partyp=mmlpars.mml_pars(partyp,type=str)
    if   partyp=='mklog_buildgal': par={'keylist':['nproc','softs','time0'],
                                        'nproc'    :{'form':int  ,'def':0       },
                                        'softs'    :{'form':dict ,'def':{ityp:0. for ityp in LIST_PTYPS}},
                                        'time0'    :{'form':float,'def':0.}}
    elif partyp=='mklog_snapshot': par={'keylist':['snapext','snaptime','resettime','time0','heatfact','softs'],
                                        'snapext'  :{'form':str  ,'def':'_000ic'},
                                        'resettime':{'form':bool ,'def':True    },
                                        'heatfact' :{'form':float,'def':1.0     },
                                        'snaptime' :{'form':float,'def':0.},
                                        'softs'    :{'form':dict ,'def':{ityp:0. for ityp in LIST_PTYPS}},
                                        'time0'    :{'form':float,'def':0.}}
    elif partyp=='mklog_restart' : par={'keylist':['nproc'],
                                        'nproc'    :{'form':int  ,'def':0       }}
    elif partyp=='mklog_interact': par={'keylist':['snapext1','snapext2','snaptime1','snaptime2','time0','softs'],
                                        'snapext1' :{'form':str  ,'def':'_000ic'},
                                        'snapext2' :{'form':str  ,'def':'_000ic'},
                                        'snaptime1':{'form':float,'def':0.},
                                        'snaptime2':{'form':float,'def':0.},
                                        'softs'    :{'form':dict ,'def':{ityp:0. for ityp in LIST_PTYPS}},
                                        'time0'    :{'form':float,'def':0.}}
    elif partyp=='param': par={'OutputDir'                :{'def':'output'          ,'form':str    },
                               'SnapshotFileBase'         :{'def':'snapshot'        ,'form':str    },
                               'SnapFormat'               :{'def':1                 ,'form':[1,2,3]},
                               'NumFilesPerSnapshot'      :{'def':1                 ,'form':int    },
                               'InitCondFile'             :{'def':'galaxy.dat'      ,'form':str    },
                               'ICFormat'                 :{'def':1                 ,'form':[1,2,3]},
                               'RestartFile'              :{'def':'restart'         ,'form':str    },
                               'InfoFile'                 :{'def':'info.txt'        ,'form':str    },
                               'TimingsFile'              :{'def':'timings.txt'     ,'form':str    },
                               'CpuFile'                  :{'def':'cpu.txt'         ,'form':str    },
                               'EnergyFile'               :{'def':'energy.txt'      ,'form':str    },
                               # CPU-time limit and restart options
                               'TimeLimitCPU'             :{'def':172800            ,'form':int    },
                               'ResubmitCommand'          :{'def':resub.sh          ,'form':str    },
                               'ResubmitOn'               :{'def':0                 ,'form':[0,1]  },
                               'CpuTimeBetRestartFile'    :{'def':10800             ,'form':int    },
                               # Simulation specific parameters
                               'TimeBegin'                :{'def':0.0               ,'form':float  },
                               'TimeMax'                  :{'def':10.2269           ,'form':float  },
                               'BoxSize'                  :{'def':1000.0            ,'form':float  },
                               'PeriodicBoundariesOn'     :{'def':0                 ,'form':[0,1]  },
                               'ComovingIntegrationOn'    :{'def':0                 ,'form':[0,1]  },
                               # Cosmological parameters
                               'HubbleParam'              :{'def':0.7               ,'form':float  },
                               'Omega0'                   :{'def':0.3               ,'form':float  },
                               'OmegaLambda'              :{'def':0.7               ,'form':float  },
                               'OmegaBaryon'              :{'def':0.04              ,'form':float  },
                               # Memory allocation
                               'BufferSize'               :{'def':29                ,'form':int    },
                               'PartAllocFactor'          :{'def':3.0               ,'form':float  },
                               'TreeAllocFactor'          :{'def':2.0               ,'form':float  },
                               # Gravitational force accuracy
                               'TypeOfOpeningCriterion'   :{'def':0                 ,'form':[0,1]  },
                               'ErrTolTheta'              :{'def':0.5               ,'form':float  },
                               'ErrTolForceAcc'           :{'def':0.002             ,'form':float  },
                               # Time integration accuracy
                               'MaxSizeTimestep'          :{'def':0.0102269         ,'form':float  },
                               'MinSizeTimestep'          :{'def':0.0               ,'form':float  },
                               'TypeOfTimestepCriterion'  :{'def':0                 ,'form':[0]    },
                               'ErrTolIntAccuracy'        :{'def':0.00125           ,'form':float  },
                               'TreeDomainUpdateFrequency':{'def':0.1               ,'form':float  },
                               'MaxRMSDisplacementFac'    :{'def':0.2               ,'form':float  },
                               # Output of snapshot files
                               'OutputListOn'             :{'def':0                 ,'form':[0,1]  },
                               'OutputListFilename'       :{'def':'output_times.txt','form':str    },
                               'TimeOfFirstSnapshot'      :{'def':0.0               ,'form':float  },
                               'TimeBetSnapshot'          :{'def':0.0102269         ,'form':float  },
                               'TimeBetStatistics'        :{'def':0.00102269        ,'form':float  },
                               'NumFilesWrittenInParallel':{'def':1                 ,'form':int    },
                               # System of units
                               'UnitVelocity_in_cm_per_s' :{'def':1.0e5             ,'form':float  },
                               'UnitLength_in_cm'         :{'def':3.085678e21       ,'form':float  },
                               'UnitMass_in_g'            :{'def':1.989e43          ,'form':float  },
                               'GravityConstantInternal'  :{'def':0.                ,'form':float  },
                               # SPH parameters
                               'DesNumNgb'                :{'def':64                ,'form':int    },
                               'MaxNumNgbDeviation'       :{'def':2                 ,'form':int    },
                               'ArtBulkViscCons'          :{'def':1.0               ,'form':float  },
                               'CourantFac'               :{'def':0.15              ,'form':float  },
                               'InitGasTemp'              :{'def':10000.            ,'form':float  },
                               'MinGasTemp'               :{'def':20.0              ,'form':float  },
                               'MinGasHsmlFractional'     :{'def':0.1               ,'form':float  },
                               # Softenings
                               'SofteningGas'             :{'def':0.0               ,'form':float  },
                               'SofteningHalo'            :{'def':0.1               ,'form':float  },
                               'SofteningDisk'            :{'def':0.05              ,'form':float  },
                               'SofteningBulge'           :{'def':0.05              ,'form':float  },
                               'SofteningStars'           :{'def':0.0               ,'form':float  },
                               'SofteningBndry'           :{'def':0.0               ,'form':float  },
                               'SofteningGasMaxPhys'      :{'def':0.0               ,'form':float  },
                               'SofteningHaloMaxPhys'     :{'def':0.1               ,'form':float  },
                               'SofteningDiskMaxPhys'     :{'def':0.05              ,'form':float  },
                               'SofteningBulgeMaxPhys'    :{'def':0.05              ,'form':float  },
                               'SofteningStarsMaxPhys'    :{'def':0.0               ,'form':float  },
                               'SofteningBndryMaxPhys'    :{'def':0.0               ,'form':float  }}
    else: raise Exception('Invalid parameter type: {}'.format(partyp))
    return par
