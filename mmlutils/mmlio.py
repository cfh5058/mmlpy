####################################################################################################################################
#
# MEAGAN LANG'S I/O METHODS
#
####################################################################################################################################
import mmlclass,mmlstring,mmlpars,mmlinfo,mmlfiles,mmlparam
import os,string,math,copy,re,sys
import numpy as np

DICT_SUBOPT={'pbs'   :['-N {jobname:s}',
                       '-M {email:s}',
                       '-m {emailopt:s}',
                       '-l gpus={gpus:d}:{nodeopt:s}',
                       '-l gpus={gpus:d}',
                       '-l nodes={nodes:d}:ppn={ppn:d}:{nodeopt:s}',
                       '-l nodes={nodes:d}:ppn={ppn:d}',
                       '-l mem={mem:d}mb',
                       '-l pmem={pmem:d}mb',
                       '-l walltime={twall_hrs:d}:{twall_min:d}:{twall_sec:d}',
                       '-o {outfile:s}',
                       '-e {errfile:s}',
                       '-j {logopt:s}'],
             'sbatch':['-J {jobname:s}',
                       '-p {queue:s}',
                       '--mail-user={email:s}',
                       '--mail-type={emailopt:s}',
                       '-n {nproc:d}',
                       '-t {twall_hrs:d}:{twall_min:d}:{twall_sec:d}',
                       '-o {outfile:s}',
                       '-e {errfile:s}',
                       '-A {account:s}']}
DICT_SUBCMD={'pbs'   :['cd $PBS_O_WORKDIR','echo $PBS_O_WORKDIR'],
             'sbatch':['cd $SLURM_SUBMIT_DIR','echo $SLURM_SUBMIT_DIR']}

####################################################################################################################################
# LIST OF PARAMETERS
def listpar(partyp,**exkw):
    """
    Creates a list of parameters
    """
    partyp=mmlpars.mml_pars(partyp,type=str)
    if   partyp=='pbs'   : par={'jobname'  :{'form':str  ,'def':''   },
                                'email'    :{'form':str  ,'def':'meagan.m.lang@vanderbilt.edu'},
                                'emailopt' :{'form':str  ,'def':'bae'},
                                'nproc'    :{'form':int  ,'def':1    },
                                'nodes'    :{'form':int  ,'def':1    },
                                'ppn'      :{'form':int  ,'def':1    },
                                'gpus'     :{'form':int  ,'def':0    },
                                'nodeopt'  :{'form':str  ,'def':''   },
                                'mem'      :{'form':int  ,'def':1000 },
                                'pmem'     :{'form':int  ,'def':1000 },
                                'walltime' :{'form':float,'def':1.0  },
                                'outfile'  :{'form':str  ,'def':'${PBS_JOBNAME}.output'},
                                'errfile'  :{'form':str  ,'def':''   },
                                'logopt'   :{'form':str  ,'def':'oe' },
                                'cmdlist'  :{'form':list ,'def':[]   }}
    elif partyp=='subopt': par={'jobname'  :{'form':str  ,'def':''   },
                                'account'  :{'form':str  ,'def':''   },
                                'queue'    :{'form':str  ,'def':'normal'},
                                'email'    :{'form':str  ,'def':'meagan.m.lang@vanderbilt.edu'},
                                'emailopt' :{'form':str  ,'def':'bae'},
                                'nproc'    :{'form':int  ,'def':1    },
                                'nodes'    :{'form':int  ,'def':1    },
                                'ppn'      :{'form':int  ,'def':1    },
                                'gpus'     :{'form':int  ,'def':0    },
                                'nodeopt'  :{'form':str  ,'def':''   },
                                'mem'      :{'form':int  ,'def':1000 },
                                'pmem'     :{'form':int  ,'def':1000 },
                                'twall'    :{'form':float,'def':1.0  },
                                'twall_hrs':{'form':int  ,'def':1.0  },
                                'twall_min':{'form':int  ,'def':0.0  },
                                'twall_sec':{'form':int  ,'def':0.0  },
                                'outfile'  :{'form':str  ,'def':''   },
                                'errfile'  :{'form':str  ,'def':''   },
                                'logopt'   :{'form':str  ,'def':'oe' },
                                'cmdlist'  :{'form':list ,'def':[]   }}
    else: raise Exception('Invalid parameter type: {}'.format(partyp))
    return par

####################################################################################################################################
# PARS SUBMISSION OPTIONS
def pars_subopt(opt,ftype):
    """
    Parses subopt to allow for multiple submission systems
    """
    opt=mmlpars.mml_pars(opt,type=dict)
    ftype=mmlpars.mml_pars(ftype,list=DICT_SUBOPT.keys())
    # Processors
    opt['ppn']=int(opt['nproc']/opt['nodes']) if 'nproc' in opt and 'nodes' in opt else opt.get('ppn',1)
    if 'nproc' in opt and opt['nproc']<opt['ppn']: opt['ppn']=1
    if 'nproc' in opt: opt['nodes']=opt.get('nodes',int(opt['nproc']/opt['ppn']))
    if 'nodes' in opt: opt['nproc']=opt.get('nproc',int(opt['nodes']*opt['ppn']))
    # Memory
    if   'pmem' not in opt: opt['pmem']=int(opt['mem']/opt['nproc'])
    elif 'mem'  not in opt: opt['mem' ]=int(opt['nproc']*opt['pmem'])
    # Walltime
    if 'twall' in opt:
        opt['twall']=float(opt['twall'])
        hrs=int(np.floor(opt['twall']))
        min=int(np.floor(60.0*(opt['twall'] - hrs)))
        sec=int(np.floor(60.0*(60.0*(opt['twall'] - hrs) - min)))
        opt.update(twall_hrs=hrs,twall_min=min,twall_sec=sec)
    else:
        hrs=float(opt['twall_hrs'])
        min=float(opt['twall_min'])
        sec=float(opt['twall_sec'])
        opt['twall']=hrs+min/60.+sec/(60.*60.)
        for ikey in ['hrs','min','sec']: opt['twall_'+ikey]=int(opt['twall_'+ikey])
    # Pars remaining options
    opt=mmlparam.parspar('mmlutils.mmlio','subopt',inpar=opt,miss=True)
    # Add things you want by default
    needkeys=['email','emailopt','queue','logopt']
    optpar=listpar('subopt')
    for ikey in needkeys: x=opt.setdefault(ikey,optpar[ikey]['def'])
    if ftype=='sbatch':
        if opt['emailopt'].upper() in ['BEGIN','END','FAIL','REQUEUE','ALL']: pass
        elif opt['emailopt']=='bae': opt['emailopt']='all'
        else: 
            verbose('Partial mail options not supported for sbatch. Setting to ALL.')
            opt['emailopt']='all'
    # Return options
    return opt

####################################################################################################################################
# READ/WRITE JOB SUBMISSION FILES
def rw_submit(rwid,fname,ftype,indict=None,overwrite=False,verbflag=False):
    """
    Reads/writes PBS files
    """
    keyreg=r'{(.*?)}' ; valreg='(.*?)'
    # Pars input
    ftype=mmlpars.mml_pars(ftype,list=DICT_SUBOPT.keys())
    # Set dependent variables
    optflag='#'+ftype.upper()
    optdict=copy.deepcopy(DICT_SUBOPT[ftype])
    stdcmds=copy.deepcopy(DICT_SUBCMD[ftype])
    # Read
    if   rwid=='R':
        indict={'cmdlist':[]}
        linelist=rw_bash('R',fname)
        # Separate into option & command lines
        optlines=[] ; cmdlines=[]
        for iline in linelist:
            if iline.startswith(optflag): optlines.append(iline.split(optflag)[-1].strip())
            else                        : cmdlines.append(iline)
        # Pars option lines
        optlines=optlines.join(r' ')
        for optfmt in optdict:
            # Determine if the option is listed
            optreg=re.sub(keyreg,valreg,optfmt)
            optmatch=re.search(optreg,optlines)
            if not optmatch:
                if verbflag: verbose('Option {} not listed'.format(optfmt))
                continue
            # Find keys and fill in dictionary
            optkey=re.findall(keyreg,optfmt)
            optval=list(optmatch.group)[1:]
            for ioptkey,ioptval in zip(optkey,optval):
                ikey,ifmt=ioptkey.split(':')
                ifmt=ifmt.strip()
                if   ifmt=='s': ival=ioptval
                elif ifmt=='d': ival=int(float(ioptval))
                elif ifmt=='f': ival=float(ioptval)
                else: raise Exception('Unknown format {} for key {}'.format(ifmt,ikey))
                indict[ikey]=ival
            # Remove selection from lines
            optlines=re.sub(optreg,'',optlines).strip()
        # Check for skipped options
        if not optlines.isspace() and verbflag:
            verbose('The following lines were not parsed:')
            print optlines
        # Pars command lines
        for iline in linelist:
            # Skip standard commands
            if iline.startswith(tuple(stdcmds)): continue
            # Add non-standard commands
            else: indict['cmdlist'].append(iline)
        # Fill in missing/dependent variables
        outdict=pars_subopt(indict,ftype)
        return outdict
    # Write
    elif rwid=='W':
        outdict=pars_subopt(indict,ftype)
        linelist=[]
        # Loop over expected formats
        for optfmt in optdict:
            try:
                optline=optfmt.format(**outdict)
                optmiss=re.findall(keyreg,optline)
                # Add line to list if all options are filled in
                if len(optmiss)==0:
                    linelist.append('{} {}'.format(optflag,optline))
                else: raise KeyError
            # Skip line if there are missing options
            except KeyError:
                optkey=re.findall(keyreg,optfmt)
                optmiss=[]
                for ikey in optkey:
                    if ikey not in outdict: optmiss.append(ikey)
                if verbflag: verbose('There are keys missing from: {}'.format(optfmt))
                if verbflag: verbose(mmlstring.val2str(optmiss))
        # Add standard lines
        linelist+=stdcmds
        linelist+=['start_time=`date +%s`','']
        linelist+=outdict['cmdlist']
        linelist+=['','wait','echo finished mpiexec',
                   'end_time=`date +%s`',"echo execution time = `expr $end_time - $start_time` s",
                   'exit 0']
        rw_bash('W',fname,linelist,overwrite=overwrite)
        return
    # Error
    else: raise Exception('Invalid rwid: {}'.format(rwid))
    return

####################################################################################################################################
# READ/WRITE PBS FILES
def rw_pbs(rwid,fname,indict=None,overwrite=False):
    """
    Reads/writes PBS files
    """
    pbsstr='#PBS'
    # Read
    if   rwid=='R':
        indict={'cmdlist':[]}
        linelist=rw_bash('R',fname)
        for iline in linelist:
            if   iline.startswith('#PBS'):
                pbsopt=iline.split('#PBS')[-1].lstrip()
                if   pbsopt.startswith('-N'): indict['jobname' ]=pbsopt.split('-N')[-1].strip()
                elif pbsopt.startswith('-M'): indict['email'   ]=pbsopt.split('-N')[-1].strip()
                elif pbsopt.startswith('-m'): indict['emailopt']=pbsopt.split('-m')[-1].strip()
                elif pbsopt.startswith('-l'):
                    if 'walltime' in pbsopt: 
                        wall=pbsopt.split('=')[-1].strip().split(':')
                        indict['walltime']=float(wall[0])+float(wall[1])/60.+float(wall[2])/3600.
                    else:
                        keyvals=pbsopt.split('-l')[-1].strip().split(':')
                        for ikeyval in keyvals:
                            if '=' in ikeyval: ikey,ival=ikeyval.strip('mb').split('=') ; ival=int(float(ival))
                            else             : ikey,ival='nodeopt',ikeyval
                            indict[ikey]=ival
                elif pbsopt.startswith('-j'): indict['logopt' ]=pbsopt.split('-j')[-1].strip()
                elif pbsopt.startswith('-o'): indict['outfile']=pbsopt.split('-o')[-1].strip()
                elif pbsopt.startswith('-e'): indict['errfile']=pbsopt.split('-e')[-1].strip()
            elif iline.startswith('cd $PBS_O_WORKDIR') or iline.startswith('echo $PBS_O_WORKDIR'): continue
            else: indict['cmdlist'].append(iline)
        outdict=mmlparam.parspar('mmlutils.mmlio','pbs',inpar=indict,init=True)
        return outdict
    # Write
    elif rwid=='W':
        outdict=mmlparam.parspar('mmlutils.mmlio','pbs',inpar=indict,init=True)
        # Get hours, minutes, & seconds
        hrs=int(np.floor(outdict['walltime']))
        min=int(np.floor(60.0*(outdict['walltime'] - hrs)))
        sec=int(np.floor(60.0*(60.0*(outdict['walltime'] - hrs) - min)))
        linelist=[]
        if len(outdict['jobname' ])>0: linelist.append(pbsstr+' -N '+outdict['jobname' ])
        if len(outdict['email'   ])>0: linelist.append(pbsstr+' -M '+outdict['email'   ])
        if len(outdict['emailopt'])>0: linelist.append(pbsstr+' -m '+outdict['emailopt'])
        nodeline=pbsstr+' -l nodes={nodes}:ppn={ppn}'.format(**outdict)
        if outdict['gpus']!=0: nodeline+=':gpus={gpus}'.format(**outdict)
        if len(outdict['nodeopt'])>0: nodeline+=':{nodeopt}'.format(**outdict)
        linelist.append(nodeline)
        linelist.append(pbsstr+' -l mem={mem}mb'.format(**outdict))
        linelist.append(pbsstr+' -l pmem={pmem}mb'.format(**outdict))
        linelist.append(pbsstr+' -l walltime={}:{}:{}'.format(hrs,min,sec))
        if len(outdict['outfile'])>0: linelist.append(pbsstr+' -o '+outdict['outfile' ])
        if len(outdict['logopt' ])>0: linelist.append(pbsstr+' -j '+outdict['logopt'  ])
        if len(outdict['errfile'])>0 and outdict['logopt'] not in ['oe','eo']:
            linelist.append(pbsstr+' -e '+outdict['errfile'])
        linelist+=['cd $PBS_O_WORKDIR','echo $PBS_O_WORKDIR']
        linelist+=outdict['cmdlist']
        rw_bash('W',fname,linelist,overwrite=overwrite)
        return
    # Error
    else: raise Exception('Invalid rwid: {}'.format(rwid))
    return


####################################################################################################################################
# SUBMITS PBS FILES TO ACCRE
def pbs_submit(pbsfile,schfile=None):
    """
    Submits .pbs files to ACCRE
    """
    # Check that submit is possible
    if not mmlinfo.hostname2compid()=='bender':
        raise Exception('Python can only submit jobs to ACCRE from BENDER.')
    # Set constants
    curdir=os.getcwd()
    # Pars input
    pbsfile=mmlpars.mml_pars(pbsfile,isfile=True)
    pbsdir=os.path.dirname(pbsfile)
    schfileDEF=os.path.join(pbsdir,'schednum')
    schfile=mmlpars.mml_pars(schfile,ispath=True)
    # Prevent overwrite
    if os.path.isfile(schfile):
        if not yorn('Overwrite existing schedule number file {}?'.format(schfile)): return
    # Change to pbs directory
    os.chdir(pbsdir)
    # Submit run
    os.system('qsub {} >&! {}'.format(os.path.basename(pbsfile),os.path.basename(schfile)))
    # Save scheduler number to file
#    print 'mmlio.pbs_submit: pbsfile={}'.format(pbsfile)
#    print 'mmlio.pbs_submit: schfile={}'.format(schfile)
#    rwfile('W',schfile,[str(schednum)],overwrite=True)
    # Change back to current directory
    os.chdir(curdir)
    return

####################################################################################################################################
# CLASS FOR READING/WRITING PBS SUBMISSION SCRIPTS
####################################################################################################################################
class pbsdict(mmlclass.mydict):
    '''
    NAME:
        mmlio.pbsdict
    PURPOSE:
        To pass around pbs submission script info.
    FORMAT:
        Dictionary with key/value pairs:
            email:     String email address to send notifications to
            mailopt:   String specifying when to send email notification
            ppn:       Int number of processors job will require on each processor 
            nodes:     Int total number of nodes job will require
            pmem:      Int memory required by job on each processor in MB
            mem:       Int total memory required by job in MB
            walltime:  Float length of time the jobs should run for in hours
            simout:    String absolute path to file where output should be saved
            simerr:    String absolute path to file where errors should be saved
            rundir:    String absolute path to directory where input is
            mpipath:   String absolute path to MPI executable
            mpiopt:    String options for MPI exec
            execprep:  String command to execute prior to executable
            execpath:  String absolute path to executable
            execflag:  Int flag for executable
            execinput: String containing input for executable
            bugchk:    String to echo to output
            command:   List of string commands to run.
        See pbsdict.get_form() dictionary for addtional details.
    '''

    def __init__(self,email=None,mailopt=None,ppn=None,nodes=None,pmem=None,mem=None,\
                 walltime=None,simout=None,simerr=None,rundir=None,\
                 mpipath=None,mpiopt=None,command=None,bugchk=None,\
                 execprep=None,execpath=None,execflag=None,execinput=None):
        '''
        NAME:
            mmlio.pbsdict.__init__
        PURPOSE:
            To initialize pbsdict objects.
        CALLING:
            obj=pbsdict(keys)
        KEYWORDS:
            Keywords include all keys listed for pbsdict objects.
        OUTPUT:
            New instance of a pbsdict object.
        '''
        mmlclass.mydict.__init__(self,email=email,mailopt=mailopt,ppn=ppn,nodes=nodes,pmem=pmem,mem=mem,\
                                 walltime=walltime,simout=simout,simerr=simerr,rundir=rundir,\
                                 mpipath=mpipath,mpiopt=mpiopt,command=command,bugchk=bugchk,\
                                 execprep=execprep,execpath=execpath,execflag=execflag,execinput=execinput)
        self.pars()
    
    def get_form(self):
        '''
        NAME:
            mmlio.pbsdict.get_form
        PURPOSE:
            To return a dictionary defining pbsdict key/value pairs.
        CALLING:
            form=pbsdict.get_form()
        OUTPUT:
            Dictionary defining pbsdict key/value pairs.
        '''
        emailDEF='meagan.m.lang@vanderbilt.edu'
        outDEF='myjob.output'
        runDEF='$PBS_O_WORKDIR'
#        runDEF='/scratch/langmm/runs'
        mpiDEF='mpiexec'
#        mpiDEF='/usr/local/mpiexec/latest/x86_64/gcc46/nonet/bin/mpiexec'
        execDEF='Gadget2'
        bugDEF='$PBSNODES'
        form={
            'email'    : mmlpars.parsdict(default=emailDEF ,type=str                 ),
            'mailopt'  : mmlpars.parsdict(default='bae'    ,type=str                 ),
            'ppn'      : mmlpars.parsdict(default=1        ,type=int  ,min=1         ),
            'nodes'    : mmlpars.parsdict(default=1        ,type=int  ,min=1         ),
            'pmem'     : mmlpars.parsdict(default=1000     ,type=int  ,min=0         ),
            'mem'      : mmlpars.parsdict(default=1000     ,type=int  ,min=0         ),
            'walltime' : mmlpars.parsdict(default=1.0      ,type=float,min=1.0/60.0  ),
            'simout'   : mmlpars.parsdict(default=outDEF   ,type=str                 ),
            'simerr'   : mmlpars.parsdict(default='oe'     ,type=str                 ),
            'rundir'   : mmlpars.parsdict(default=runDEF   ,type=str                 ),
            'mpipath'  : mmlpars.parsdict(default=mpiDEF   ,type=str                 ),
            'mpiopt'   : mmlpars.parsdict(default=''       ,type=str                 ),
            'execprep' : mmlpars.parsdict(default=''       ,type=str                 ),
            'execpath' : mmlpars.parsdict(default=execDEF  ,type=str                 ),
            'execflag' : mmlpars.parsdict(default=0        ,type=int                 ),
            'execinput': mmlpars.parsdict(default=''       ,type=str                 ),
            'command'  : mmlpars.parsdict(default=['cmd']  ,type=list                ),
            'bugchk'   : mmlpars.parsdict(default=bugDEF   ,type=str                 )
            }
        return form
    
    def pars(self):
        '''
        NAME:
            mmlio.pbsdict.pars
        PURPOSE:
            To parse pbsdict objects.
        CALLING:
            pbsdict.pars()
        '''
        self=mmlpars.mml_formpars(self,self.get_form())


    def rwpbs(self,rwid,fname,overwrite=False):
        '''
        NAME:
            mmlio.pbsdict.rwpbs
        PURPOSE:
            To read/write pbs submission scripts.
        CALLING:
            pbsdict.rwpbs(rwid,fname[,overwrite=False])
        '''
        # Set constants
        bugchk='$PBSNODES'
        rwidLIST=['R','W']
        # Pars input
        rwid=mmlpars.mml_pars(rwid,type=str,list=rwidLIST)
        # Read
        if rwid == 'R':
            linelist=rw_bash('R',fname)
            self['command']=[]
            # Match lines
            for iline in linelist:
                if iline.startswith('#PBS -M '):
                    self['email']=iline.split('#PBS -M ')[1]
                elif iline.startswith('#PBS -m '):
                    self['mailopt']=iline.split('#PBS -m ')
                elif iline.startswith('#PBS -l nodes='):
                    vars=iline.split('#PBS -l nodes=')[1].split(':ppn=')
                    self['nodes']=int(float(vars[0]))
                    self['ppn']=int(float(vars[1]))
                elif iline.startswith('#PBS -l mem='):
                    self['mem']=int(float(iline.split('#PBS -l mem=')[1].split('mb')[0]))
                elif iline.startswith('#PBS -l mem='):
                    self['pmem']=int(float(iline.split('#PBS -l pmem=')[1].split('mb')[0]))
                elif iline.startswith('#PBS -l walltime='):
                    vars=iline.split('#PBS -l walltime=')[1].split(':')
                    self['walltime']=float(vars[0])+float(vars[1])/60.0+float(vars[2])/3600.0
                elif iline.startswith('#PBS -o '):
                    self['simout']=iline.split('#PBS -o ')[1]
                elif iline.startswith('#PBS -j '):
                    self['simerr']=iline.split('#PBS -j ')[1]
                elif iline.startswith('cd '):
                    self['rundir']=iline.split('cd ')[1]
                elif iline.startswith('echo '):
                    self['bugchk']=iline.split('echo ')[1]
                else:
                    self['command'].append(iline)
        # Write
        elif rwid == 'W':
            # Get hours, minutes, & seconds
            hrs=int(np.floor(self['walltime']))
            min=int(np.floor(60.0*(self['walltime'] - hrs)))
            sec=int(np.floor(60.0*(60.0*(self['walltime'] - hrs) - min)))
            linelist=[
                '#PBS -M {}'.format(self['email'  ]),
                '#PBS -m {}'.format(self['mailopt']),
                '#PBS -l nodes={}:ppn={}'.format(self['nodes'],self['ppn']),
                '#PBS -l mem={}mb'.format(int(self['mem'])),
                '#PBS -l pmem={}mb'.format(int(self['pmem'])),
                '#PBS -l walltime={}:{}:{}'.format(hrs,min,sec),
                '#PBS -o {}'.format(self['simout']),
                '#PBS -j {}'.format(self['simerr']),
                'cd $PBS_O_WORKDIR',
                'echo $PBS_O_WORKDIR',
                'echo {}'.format(self['bugchk'])
                ]
            linelist+=self['command']
            rw_bash('W',fname,linelist,overwrite=overwrite)
            
    def rwmpi(self,rwid,fname,overwrite=False):
        '''
        NAME:
            mmlio.pbsdict.rwmpi
        PURPOSE:
            To read/write MPI pbs submission scripts.
        CALLING:
            pbsdict.rwmpi(rwid,fname[,overwrite=False])
        '''
        # Set constants
        rwidLIST=['R','W']
        # Pars
        rwid=mmlpars.mml_pars(rwid,type=str,list=rwidLIST)
        # Read
        if rwid=='R':
            self.rwpbs('R',fname)
            if len(self['command']) != 1:
                raise Exception('Command is too long to be MPI.')
            # Before and after # processors command
            vars=self['command'][0].split('-np')
            if len(vars) != 2:
                raise Exception('Command does not follow the MPI format'.format(self['command'][0]))
            # MPI commands
            vars_mpi=vars[0].split()
            self['mpipath']=vars_mpi[0]
            if len(vars_mpi) == 1:
                self['mpiopt']=''
            elif len(vars_mpi) == 2:
                self['mpiopt']=vars_mpi[1]
            else:
                raise Exception('Command does not follow the MPI format'.format(self['command'][0]))
            # # processor commands
            vars_npr=vars[1].split('./')
            #self['nproc']=int(float(vars_npr[0].split()))
            if len(vars_npr) !=2 :
                raise Exception('Command does not follow the MPI format'.format(self['command'][0]))
            # Executable commands
            vars_exe=vars_npr[1].split()
            if len(vars_exe) == 2:
                self['execpath']=vars_exe[0]
                self['execinput']=vars_exe[1]
                self['execflag']=0
            elif len(vars_exe) == 3:
                self['execpath']=vars_exe[0]
                self['execinput']=vars_exe[1]
                self['execflag']=int(float(vars_exe[2]))
            else:
                raise Exception('Command does not follow the MPI format'.format(self['command'][0]))
        # Write
        elif rwid=='W':
            nproc=self['ppn']*self['nodes']
            self['command']=['{} {} -np {} ./{} {} {}'.format(self['mpipath'],self['mpiopt'],nproc,
                                                              self['execpath'],self['execinput'],self['execflag'])]
            self.rwpbs('W',fname,overwrite=overwrite)

def rwfile(rwid,fname,linelist=[],overwrite=False):
    '''
    NAME:
        mmlio.rwfile
    PURPOSE:
        To read/write/append lines from/to a file.
    CALLING:
        Read:   linelist=rwfile("R",fname)
        Write:  rwfile("W",fname,linelist[,overwrite=False])
        Append: rwfile("A",fname,linelist)
    ARGUMENTS:
        rwid:      Str specifying what to do (Read-"R", Write-"W", or Append="A")
        fname:     Str absolute path of file to read/write/append from/to
        linelist:  List of lines to write to file (only for write)
    KEYWORDS:
        overwrite: Bool specifying if existing file should be overwritten
    OUTPUT:
        linelist:  List of lines read from file (only for read)
    '''
    # Set constants
    rwidLIST=['R','W','A']
    # Pars input
    rwid=mmlpars.mml_pars(rwid,type=str,list=rwidLIST)
    fname=mmlpars.mml_pars(fname,type=str)
    linelist=mmlpars.mml_pars(linelist,type=list)
    overwrite=mmlpars.mml_pars(overwrite,type=bool)
    # Handle different options
    if rwid=='R':
        if not os.path.isfile(fname):
            raise Exception('File provided for read does not exist: {}'.format(fname))
    elif rwid=='W':
        if os.path.isfile(fname):
            if overwrite: os.remove(fname)
            else: return
        if len(linelist)==0:
            raise Exception('List of lines for write is empty.')
    elif rwid=='A':
        if not os.path.isfile(fname):
            raise Exception('File provided for append does not exist: {}'.format(fname))
        if len(linelist)==0:
            raise Exception('List of lines for append is empty.')
    # Open file
    fid=open(fname,rwid.lower())
    # Read
    if rwid=='R':
        for iline in fid:
            linelist.append(iline)
    # Write/Append
    else:
        for iline in linelist:
            if isinstance(iline,str):
                fid.write(iline+'\n')
            else:
                fid.write(repr(iline)+'\n')
    # Close file
    fid.close()
    # Return output
    if rwid=='R':
        return linelist
    else:
        return

def rw_bash(rwid,fname,linelist=[],overwrite=False):
    '''
    NAME:
        mmlio.rw_bash
    PURPOSE:
        To read/write/append lines from/to a bash script.
    CALLING:
        Read:   linelist=rwfile("R",fname)
        Write:  rwfile("W",fname,linelist[,overwrite=False])
        Append: rwfile("A",fname,linelist)
    ARGUMENTS:
        rwid:      Str specifying what to do (Read-"R", Write-"W", or Append="A")
        fname:     Str absolute path of bash script to read/write/append from/to
        linelist:  List of lines to write to bash script (only for write)
    KEYWORDS:
        overwrite: Bool specifying if existing bash script should be overwritten
    OUTPUT:
        linelist:  List of lines read from bash script (only for read)
    '''
    # Set constants
    rwidLIST=['R','W','A']
    # Pars input
    rwid=mmlpars.mml_pars(rwid,type=str,list=rwidLIST)
    fname=mmlpars.mml_pars(fname,type=str)
    linelist=mmlpars.mml_pars(linelist,type=list)
    overwrite=mmlpars.mml_pars(overwrite,type=bool)
    # Create new line list
    bashhead='#!/bin/sh'
    scriptlines=[bashhead]+linelist
    # Read/write/append
    if rwid=='R':
        scriptlines=rwfile('R',fname)
        if scriptlines[0]==bashhead:
            linelist=scriptlines[1:]
        else:
            linelist=scriptlines
        return linelist
    # Write
    elif rwid=='W':
        scriptlines=[bashhead]+linelist
    # Append
    elif rwid=='A':
        scriptlines=linelist
    # Return
    rwfile('W',fname,scriptlines,overwrite=overwrite)
    return

####################################################################################################################################
# METHOD TO READ/WRITE LISTS TO FILE
def rwlist(rwid,fname,flist=None,overwrite=False,style=None,headline=None,tailline=None):
    """
    Reads/writes a list from/to a file
    """
    # Pars input
    rwid=mmlpars.mml_pars(rwid,type=str,list=['R','W','A'])
    fname=mmlpars.mml_pars(fname,type=str)
    overwrite=mmlpars.mml_pars(overwrite,type=bool)
    if rwid=='W' or rwid=='A':
        if flist==None: raise Exception('Dictionary must be supplied if writing/appending to a file.')
        retflag=mmlfiles.prep_write(fname,overwrite)
        if rwid=='W' and os.path.isfile(fname):
            if not overwrite: return
            else: os.remove(fname)
    if rwid=='R' or rwid=='A':
        if not os.path.isfile(fname):
            raise Exception('File you are trying to read/append does not exist: {}'.format(fname))
    style=mmlpars.mml_pars(style,default='c',type=str,list=['c','fortran'])
    # Set comment character
    if   style == 'c'      : cchar='#'
    elif style == 'fortran': cchar='C'
    # Open file
    fid=open(fname,rwid.lower())
    # Read
    if rwid=='R':
        flist=[]
        for line in fid:
            if line.startswith(cchar):
                pass
            else:
                line=line.split(cchar)[0]
                val=mmlstring.str2val(line)
                flist.append(val)
    # Write/append
    else:
        if isinstance(headline,str): fid.write(cchar+headline+'\n')
        for v in flist:
            if   style == 'c'      : fid.write(mmlstring.val2str(v)+'\n')
            elif style == 'fortran':
                if isinstance(v,bool): fid.write(('.'+mmlstring.val2str(v).upper()+'.')+'\n')
                else                 : fid.write(mmlstring.val2str(v)+'\n')
        if isinstance(tailline,str): fid.write(cchar+tailline+'\n')
    # Close file
    fid.close()
    # Return output
    if rwid=='R': return flist

####################################################################################################################################
# METHOD TO READ/WRITE DICTIONARIES TO FILE
def rwdict(rwid,fname,fdict=None,overwrite=False,width=29,sort=True,style=None,headline=None,tailline=None):
    '''
    NAME:
        mmlio.rwdict
    PURPOSE:
        Reading/writing/appending dictionaries to files.
    CALLING:
        Read:   fdict=rwdict("R",fname)
        Write:  rwdict("W",fname,fdict[,overwrite=False,width=29])
        Append: rwdict("A",fname,fdict[,width=29])
    ARGUMENTS:
        rwid:      Str specifying whether to read ("R"), write ("W"), or append ("A")
        fname:     Str absolute path to file to read/write/append
    KEYWORDS:
        fdict:     Dictionary to write/append to file (required for write & append)
        overwrite: Bool specifying whether to overwrite an existing file
        width:     Int width of key & value columns in characters
        sort:      Bool specifying if output should be sorted by dict keys
        style:     Str specifying what style should be used to read/write dictionary
            "c"      : Comments are # character and keys proceed values
            "fortran": Comments are C character and values proceed keys
        headline:  Str to print as comment at beginning of output
        tailline:  Str to print as comment at end of output
    OUTPUT:
        fdict:     Dictionary read in from file (only output for read)
    '''
    # Pars input
    rwid=mmlpars.mml_pars(rwid,type=str,list=['R','W','A'])
    fname=mmlpars.mml_pars(fname,type=str)
    overwrite=mmlpars.mml_pars(overwrite,type=bool)
    if rwid=='W' or rwid=='A':
        if fdict==None:
            raise Exception('Dictionary must be supplied if writing/appending to a file.')
        retflag=mmlfiles.prep_write(fname,overwrite)
        if rwid=='W' and os.path.isfile(fname):
            if not overwrite: return
            else: os.remove(fname)
    if rwid=='R' or rwid=='A':
        if not os.path.isfile(fname):
            raise Exception('File you are trying to read/append does not exist: {}'.format(fname))
    style=mmlpars.mml_pars(style,default='c',type=str,list=['c','fortran'])
    # Set comment character
    if   style == 'c'      : cchar='#'
    elif style == 'fortran': cchar='C'
    # Open file
    fid=open(fname,rwid.lower())
    # Read
    if rwid=='R':
        fdict={}
        keylist=[]
        for line in fid:
            if line.startswith(cchar):
                pass
            else:
                line=line.split(cchar)[0]
                vars=line.split()
                if len(vars)==2:
                    if style == 'c':
                        key=vars[0]
                        val=mmlstring.str2val(vars[1])
                    elif style == 'fortran':
                        val=mmlstring.str2val(vars[0])
                        key=vars[1]
                else:
                    key=vars[0]
                    val=mmlstring.str2val('')
                fdict[key]=val
                keylist.append(key)
        fdict['keylist']=keylist
    # Write/append
    else:
        if isinstance(headline,str): fid.write(cchar+headline+'\n')
        if fdict.has_key('keylist'):
            klist=fdict['keylist']
        else:
            if sort:
                klist=sorted(fdict.keys(),key=str.lower)
            else:
                klist=fdict.keys()
        for k in klist:
            v=fdict[k]
            if style == 'c':
                fid.write(k.ljust(width)+mmlstring.val2str(v).ljust(width)+'\n')
            elif style == 'fortran':
                if isinstance(v,bool):
                    fid.write(('.'+mmlstring.val2str(v).upper()+'.').ljust(width)+k.ljust(width)+'\n')
                else:
                    fid.write(mmlstring.val2str(v).ljust(width)+k.ljust(width)+'\n')
        if isinstance(tailline,str): fid.write(cchar+tailline+'\n')
    # Close file
    fid.close()
    # Return output
    if rwid=='R': return fdict

####################################################################################################################################
# METHOD TO ASK USER YES OR NO QUESTIONS
####################################################################################################################################
def yorn(queststr,width=None):
    '''
    NAME:
        mmlio.yorn
    PURPOSE:
        Getting answers to yes/no questions from a user.
    CALLING:
        ans=yorn(queststr[,width=50])
    ARGUMENTS:
        queststr: String supplying yes or no question.
    KEYWORDS:
        width:    Width of the string that questions should be in.
    OUTPUT:
        ans:      Boolean ans to yes or not question.
    '''
    # Set constants
    anspos=['Y','y','yes','Yes','YES','1']
    ansneg=['N','n','no','No','NO','0']
    # Pars input
    queststr=mmlpars.mml_pars(queststr,type=str)
    width=mmlpars.mml_pars(width,default=50,type=int)
    # Get answer
    ans=None
    while ans==None:
        ansstr=raw_input(queststr.ljust(width)+' [y/n]: ')
        if ansstr in anspos: ans=True
        if ansstr in ansneg: ans=False
        if ans==None:
            print 'Invalid answer. Try again.'
    return ans

####################################################################################################################################
# METHOD TO ASK QUESTIONS WITH STRING ANSWERS
####################################################################################################################################
def askquest(queststr,width=50,default='',dtype='str'):
    """
    Ask the user to enter information
    """
    ans=None
    while ans==None:
        ansstr=raw_input(queststr.ljust(width)+' ['+repr(default)+']: ')
        if len(ansstr)==0:
            if isinstance(default,str): ansstr=default
            else:                       ansstr=repr(default)
        if dtype=='str':
            ans=string.strip(string.strip(ansstr,'"'),"'")
        else:
            try:
                ans=mmlstring.str2val(ansstr,dtype=dtype)
##                 ansfloat=float(ansstr)
##                 if dtype=='float':
##                     ans=ansfloat
##                 elif dtype=='int':
##                     if int(ansfloat)==ansfloat: ans=int(ansfloat)
##                 elif dtype=='long':
##                     if long(ansfloat)==ansfloat: ans=long(ansfloat)
            except:
                ans=None
        if ans==None:
            print 'Invalid answer. Try again.'
    return ans

####################################################################################################################################
# METHOD TO ASK QUESTION FROM LIST
####################################################################################################################################
def askselect(queststr,anslist,width=50,npercol=25,colpad=5,default=None):
    """
    Method for user to select an answer from a list
    """
    # Pars input
    queststr=mmlpars.mml_pars(queststr,type=str)
    anslist=mmlpars.mml_pars(anslist,type=list)
    if default is None or default not in anslist: default=anslist[0]
    queststr+=' [{}]'.format(default)
    # Print question
    print queststr
    # Get maximum width of ans and number
    maxwid_ans=max(len(str(max(anslist,key=lambda x: len(str(x))))),width/2)
    maxwid_num=len(str(len(anslist)))
    # Count columns and rows
    ncol=len(anslist)/npercol+1
    nrow=min(len(anslist),npercol)
    # Print list of answers
    ipadstr=colpad*' '
    strlist=nrow*['']
    for ians in range(len(anslist)):
        icol=ians/npercol
        irow=ians%npercol
        iansstr=str(anslist[ians]).ljust(maxwid_ans)
        inumstr=str(ians).rjust(maxwid_num)
        strlist[irow]+=str('{} [{}]'.format(iansstr,inumstr)+ipadstr)
    for irow in strlist: print irow
    # Get answer from user
    ans=None
    while ans==None:
        ans=askquest(queststr,width=width,default=anslist.index(default),dtype='int')
        if ans >= len(anslist): ans=None
        if ans==None:
            print 'Invalid answer. Try again.'
    return anslist[ans]
            
####################################################################################################################################
# METHOD TO GET FUNCTION NAME
def get_methname(level=0):
    """
    Method to return a calling function's name
    """
    import inspect
    trace=inspect.stack()[1+level]
    modstr=os.path.basename(trace[1]).split('.py')[0]
    namestr='{}.{}'.format(modstr,trace[3])
    return namestr

####################################################################################################################################
# METHOD TO PRINT INFO
####################################################################################################################################
def verbose(infostr,width=50,addval='',border=False,valpad=50,inclfunc=True,outstring=False):
    '''
    NAME:
        mmlio.verbose
    PURPOSE:
        To display standardized output to the screen.
    CALLING:
        verbose(infostr[,width=50,addval=""])
    ARGUMENTS:
        infostr: Str to print to the screen.
    KEYWORDS:
        width:   Int total width of the printed string with padding
        addval:  Additional variable to print to right of infostr
        border:  Bool specifying if line borders should be placed above/below line
    '''
    import inspect
    # Set constants
    if inclfunc:
        trace=inspect.stack()[1]
#        print inspect.stack()[2]
#        yorn('cotinue?')
        modstr=os.path.basename(trace[1]).split('.py')[0]
        namestr='{}.{}'.format(modstr,trace[3])
        linestr=str(trace[2])
        namewid=max(25,len(namestr))
        linewid=max( 5,len(linestr))
        funcstr='[{}: {}] '.format(namestr.ljust(namewid),linestr.ljust(linewid))
    else:
        funcstr=''
    valpad=max(valpad,len(infostr))
    width=max(width,valpad+len(mmlstring.val2str(addval))+len(funcstr))
    strlist=[]
    if border: strlist.append(mmlstring.linestr(width=width))
    strlist.append(mmlstring.val2str(funcstr+infostr).ljust(valpad)+mmlstring.val2str(addval).rjust(width-valpad))
    if border: strlist.append(mmlstring.linestr(width=width))
    if outstring:
        return strlist
    else:
        for istr in strlist: print istr
        return

####################################################################################################################################
# METHOD TO READ/WRITE TABLES
####################################################################################################################################
def plottable(plotfile,tabdict,xkey=None,ykey=None,xlim=None,xlab=None,ncol=2,overwrite=False,**loadkw):
    """Plot table columns"""
    import matplotlib.pyplot as plt
    validtyps=['float','int','long','complex']
    # Pars input
    loadkw['inarr']=True
    if   isinstance(tabdict,str ): tabdict=rwtable('R',tabdict,**loadkw)
    elif isinstance(tabdict,dict): tabdict.setdefault('keylist',tabdict.keys())
    else: raise Exception('tabdict must be path to file or dictionary table not type: {}'.format(type(tabdict)))
    # Get valid keys
    validkeys=[]
    for k in tabdict['keylist']:
        if not isinstance(tabdict[k],np.ndarray): continue
        for t in validtyps:
            if t in str(type(tabdict[k][0])): 
                validkeys.append(k)
                break
    # Set/parse x variable
    if xkey is None: xkey=validkeys[0]
    if xkey in validkeys: validkeys.remove(xkey)
    else: raise Exception('xkey {} not in list of valid keys: {}'.format(xkey,validkeys))
    # Set/pars y variable
    if ykey is None: ykey=validkeys
    if not isinstance(ykey,list): ykey=[ykey]
    yrmk=[]
    for k in ykey:
        if k not in validkeys:
            yrmk.append(k)
            verbose('ykey {} not in list of valid keys: {}'.format(k,validkeys))
    ykey=set(ykey).difference(yrmk)
    # Sort
    idxsort=sorted(range(len(tabdict[xkey])),key=lambda k: tabdict[xkey][k])
    for k in tabdict['keylist']: tabdict[k]=tabdict[k][idxsort]
    # Plot
    if not mmlfiles.prep_write(plotfile,overwrite): return tabdict
    plt.clf()
    nrow=int(np.ceil(float(len(ykey))/float(ncol)))
    if xlim is None: xlim=[min(tabdict[xkey]),max(tabdict[xkey])]
    if xlab is None: xlab=xkey.capitalize()
    if xlim[0]==xlim[1]: return tabdict
    for idx,y in enumerate(ykey):
        iax=plt.subplot(ncol,nrow,idx+1)
        iax.plot(tabdict[xkey],tabdict[y])
        iax.set_xlim(xlim)
        iax.set_ylim((min(tabdict[y]),max(tabdict[y])))
        iax.set_xlabel(xlab)
        iax.set_ylabel(y.capitalize())
    plt.tight_layout()
    plt.savefig(plotfile)
    verbose('Plot data table:')
    print '    '+plotfile
    # Return
    return tabdict

def rwtable(rwid,fname,tabdict=None,overwrite=False,keylineno=2,typlineno=3,comment='#',sep='   ',keylist=None,
            keystr=None,typstr=None,uniquekey=None,sortkey=None,inarr=False):
    '''
    NAME:
        mmlio.rwtable
    PURPOSE:
        To read/write a table from a file.
    CALLING:
        Read:  tabdict=rwtable("R",fname[,keylineno=,typlineno=,commen=,sep=])
        Write: rwtable("W",fname,tabdict[,overwrite=,keylineno=,typlineno=,comment=,sep=])
    ARGUMENTS:
        rwid:      String specifying if fname should be read ("R") or written ("W")
        fname:     String absolute path of file to read/write
        tabdict:   Dictionary containing table data to write to fname
    KEYWORDS:
        overwrite: Bool specifiying if existing file is/isn"t overwritten
        keylineno: Int number of line with keywords
        typlineno: Int number of line with type info
        comment:   String marking a line as a comment
        sep:       String separating columns in the table
        keylist:   List of keys in order they should be written.
                   [Can also be a key in tabdict]
    OUTPUT:
        tabdict:   Dictionary containing table data read from fname
    '''
    # Pars input
    rwid=mmlpars.mml_pars(rwid,type=str,list=['R','W','A','plot'])
    fname=mmlpars.mml_pars(fname,type=str)
    if   rwid=='R'   : tabdict={}
    elif rwid=='plot': pass
    else             : tabdict=mmlpars.mml_pars(tabdict,type=dict)
    overwrite=mmlpars.mml_pars(overwrite,default=False,type=bool)
    keylineno=mmlpars.mml_pars(keylineno,default=2,type=int)
    typlineno=mmlpars.mml_pars(typlineno,default=3,type=int)
    if typlineno==keylineno:
        raise Exception('Type line number cannot be the same as key line number.')
    comment=mmlpars.mml_pars(comment,default='#',type=str)
    sep=mmlpars.mml_pars(sep,default='   ',type=str)
    if fname not in ['stdout','stderr']:
        mmlfiles.mkdirs(os.path.dirname(fname))
    if rwid=='W': 
        if not mmlfiles.prep_write(fname,overwrite): return
    # Read
    if rwid=='R':
        # Check that file exists
        if not os.path.isfile(fname): raise Exception('File does not exist: {}'.format(fname))
        # Open file
        fid=open(fname,'r')
        linecount=0
        # Loop over lines
        for iline0 in fid:
            iline=iline0.split('\n')[0]
            # Handle comments
            if iline.startswith(comment):
                if linecount==keylineno:
                    keystr=(iline.strip(comment)).split()
                    for ik in keystr: tabdict[ik]=[]
                elif linecount==typlineno:
                    typstr=(iline.strip(comment)).split()
                else: pass
            # Handle body
            else:
                vars=iline.split()
                if keystr==None: keystr=range(len(vars))
                if typstr==None: typstr=len(vars)*[None]
                if len(vars) != len(keystr):
                    raise Exception('Number of columns must match number of keys.')
                if len(vars) != len(typstr):
                    raise Exception('Number of columns must match number of types.')
                for ik in range(len(vars)):
                    # Retype variable
                    ival=mmlstring.str2val(vars[ik],dtype=typstr[ik])
                    if keystr[ik] not in tabdict: tabdict[keystr[ik]]=[]
                    tabdict[keystr[ik]].append(ival)
            # Count
            linecount+=1
        # Close file and return
        fid.close()
        tabdict['keylist']=keystr
        # Sort
        if sortkey in keystr:
            idxsort=sorted(range(len(tabdict[sortkey])),key=lambda k: prop[sortkey][k])
            for k in keystr: tabdict[k]=list(np.array(tabdict[k])[idxsort])
        # Convert to arrays
        if inarr: 
            for k in keystr: tabdict[k]=np.array(tabdict[k])
        return tabdict
    # Write
    elif rwid=='W':
        # Ensure file dosn't get overwritten
        if os.path.isfile(fname) and fname not in ['stdout','stderr']:
            if overwrite: os.remove(fname)
            else        : return
        # Get list of keys
        if isinstance(keylist,list):
            keystr=keylist
        elif 'keylist' in tabdict:
            keystr=tabdict['keylist']
        else:
            keystr=sorted(tabdict.keys())
        # Check lengths of lists
        listbool=np.array([isinstance(tabdict[k],(list,np.ndarray)) for k in keystr])
        if   not np.any(listbool): tabdict={k:[tabdict[k]] for k in keystr}
        lenbody=len(tabdict[keystr[0]])
        for ikey in keystr:
            if isinstance(tabdict[ikey],np.ndarray):
                if len(tabdict[ikey].shape)==1: tabdict[ikey]=list(tabdict[ikey])
                else:
                    raise Exception('Value for key {} is of shape {}, but must be 1D'.format(ikey,tabdict[ikey].shape))
            if not isinstance(tabdict[ikey],list):
                raise Exception('All values must be lists. Value for key {} is {}'.format(ikey,type(tabdict[ikey])))
            if len(tabdict[ikey]) != lenbody:
                raise Exception('All lists must be the same length ({}). List for key {} is {}.'.format(lenbody,ikey,len(tabdict[ikey])))
        # Sort
        if sortkey in keystr:
            idxsort=sorted(range(len(tabdict[sortkey])),key=lambda k: prop[sortkey][k])
            for k in keystr: tabdict[k]=list(np.array(tabdict[k])[idxsort])
        # Determine width of fields and pad things
        labstr=[] ; typstr=[] ; widstr=[]
        for icol in range(len(keystr)):
            ikey=keystr[icol]
            ityplong=str(type(tabdict[ikey][0]))
            if   'float'   in ityplong: ityp='float'
            elif 'int'     in ityplong: ityp='int'
            elif 'long'    in ityplong: ityp='long'
            elif 'complex' in ityplong: ityp='complex'
            elif 'str'     in ityplong: ityp='str'
            elif 'tuple'   in ityplong: ityp='tuple'
            elif 'list'    in ityplong: ityp='list'
            elif 'dict'    in ityplong: ityp='dict'
            elif 'bool'    in ityplong: ityp='bool'
            else: raise Exception('Unidentified type for key {}: {}'.format(ikey,ityplong))
            iwid=len(mmlstring.val2str(max(tabdict[ikey]+[ikey,ityp],key=lambda x: len(mmlstring.val2str(x)))))
            labstr.append(ikey.ljust(iwid))
            typstr.append(ityp.ljust(iwid))
            widstr.append(iwid)
        widhead=np.array(widstr).sum()+len(comment)+len(sep)*len(keystr)
        # Open file
        if   fname=='stdout': fid=sys.stdout
        elif fname=='stderr': fid=sys.stderr
        else                : fid=open(fname,'w')
        # Create header
        lenhead=max(keylineno,typlineno)+1
        for ihead in range(lenhead+1):
            if ihead==keylineno:
                linestr=comment+sep+sep.join(labstr)
            elif ihead==typlineno:
                linestr=comment+sep+sep.join(typstr)
            elif ihead==lenhead:
                linestr=comment*widhead
            else:
                linestr=comment
            fid.write(linestr+'\n')
        # Create body
        for ibody in range(lenbody):
            strlist=[]
            for icol in range(len(keystr)):
                ikey=keystr[icol]
                iwid=widstr[icol]
                ityp=typstr[icol]
                istr=mmlstring.val2str(tabdict[ikey][ibody])
                if '--' in istr.ljust(iwid): 
                    print ikey,ityp,iwid,tabdict[ikey][ibody],type(tabdict[ikey][ibody])
                    raise Exception('Error in creating string.')
                strlist.append(istr.ljust(iwid))
            linestr=' '*len(comment+sep)+sep.join(strlist)
            fid.write(linestr+'\n')
        # Close file and return
        if fname not in ['stdout','stderr']: fid.close()
        return
    # Append
    elif rwid=='A':
        inkw=dict(keylineno=keylineno,typlineno=typlineno,comment=comment,sep=sep,
                  keylist=keylist,keystr=keystr,typstr=typstr)
        # Load existing file if it exists and add data to it
        if os.path.isfile(fname):
            tabdict0=rwtable('R',fname,**inkw)
            # Convert to lists
            if not isinstance(tabdict[uniquekey],list):
                if isinstance(tabdict[uniquekey],np.ndarray):
                    for k in tabdict0['keylist']: tabdict[k]=list(tabdict[k])
                else:
                    for k in tabdict0['keylist']: tabdict[k]=[tabdict[k]]
            if uniquekey not in tabdict0:
                raise Exception('uniquekey ({}) must be a key in the file being appended.'.format(uniquekey))
            if uniquekey not in tabdict:
                raise Exception('uniquekey ({}) must be a key in the entries being added.'.format(uniquekey))
            for irow in range(len(tabdict[uniquekey])):
                if tabdict[uniquekey][irow] in tabdict0[uniquekey]: 
                    if overwrite: 
                        idxmatch=tabdict0[uniquekey].index(tabdict[uniquekey][irow])
                        for ikey in tabdict0['keylist']: del tabdict0[ikey][idxmatch]
                    else: continue
                for ikey in tabdict0['keylist']: tabdict0[ikey].append(tabdict[ikey][irow])
        # Otherwise just create from provided dictionary
        else: tabdict0=tabdict
        # Write updated dictionary to file and return it
        rwtable('W',fname,tabdict0,overwrite=True,**inkw)
        return tabdict0
    # Handle incorrect rwid
    else: raise Exception('Invalid rwid: {}'.format(rwid))
