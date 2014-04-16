###################################################################################################################################
#
# MEAGAN LANG'S GADGET METHODS
#
####################################################################################################################################
import mmlpars,mmlclass,mmlfiles,mmlio,mmlplot,mmlsimplot
import mmlmath as mmath
import matplotlib as mplib
import os,glob,shutil,pNbody,cPickle,copy
import numpy as np
import mmlgadget

####################################################################################################################################
# METHOD TO RETURN CALC FILES
def files(rundir,runid,addid=''):
    '''
    NAME:
        mmlsimcalc.files
    PURPOSE:
        To return files associated with N-body simulation post processesing.
    CALLING:
        fdict=files(rundir,runid,addid="")
    ARGUMENTS:
        rundir: Directory containing all files for an N-body simulation
        runid:  String uniquely identifying a N-body simulation
    KEYWORDS:
        addid:  String to add to the end of runid to create prefix for file names
    OUTPUT:
        fdict:  Dictionary containing the following key/value pairs:
          dir:       Directory containing all simulation calc files
          calclist:  File containing a list of variables calculated from each snapshot
          csnapdir:  Directory containing full diagnostic files for each snapshot
          csnapbase: Base of snapshot diagnostic files
    '''
    pfix_shrt=runid+'.'
    pfix=runid+addid+'.'
    files={'dir':os.path.join(rundir,'calc'+addid),'calcstat':os.path.join(rundir,pfix_shrt+'static_calc')}
#           'clstat':os.path.join(rundir,pfix+'static_calc'),'clevol':os.path.join(rundir,pfix+'evolve_calc')}
    # Add input files
    fkeys=['calclist','genstat','profstat','barstat_dens','barstat_scfm2','barstat_ellip']
    flist=['statlist','genstat','profstat','barstat_dens','barstat_scfm2','barstat_ellip']
    dkeys=['csnapdir']
    dlist=['csnap']
    files=mmlfiles.filedict(files['dir'],dirkeys=dkeys,dirlist=dlist,filekeys=fkeys,filelist=flist,pfix=pfix,fdict=files)
    files['csnapbase']=os.path.join(files['csnapdir'],pfix+'calcsnap')
    # Return file dictionary
    return files

####################################################################################################################################
# METHOD TO PERFORM SIMULATION CALCULATIONS
def analyze(simstr,method='analyze',ftype='gadget',fileopt=None,overwrite=False,shortlist=False,
            calcdict=None,snapfile0=None,idgal2=None,plot_kw={},anim_kw={},
            singlist=['star'],winddict=None,**extra_kw):
    '''
    NAME:
        mmlsimcalc.calc
    PURPOSE:
        To run calculations on a Nbody simulation.
    CALLING:
        calc(simstr[,ftype="gadget"])
    ARGUMENTS:
        simstr:    mmlsim object describing basic parameters of the simulation
    KEYWORDS:
        method:    String specifying a specific part of the analysis to perform.
        ftype:     String specifying what type of simulation it was.
        fileopt:   Dict of file naming options (See icgen.mmlsim.mkfiledict)
        overwrite: Bool specifying whether or not to overwrite existing files
    '''
    # Set constants
    verbwid=150
    singlistVAL=['gas','halo','disk','bulge','star']
    singlistDEF=['star']
    winddictDEF={ising:(6000,6000) for ising in singlistVAL}
    # Pars input
    singlist=mmlpars.mml_pars(singlist,default=singlistDEF,type=list)
    winddict=mmlpars.mml_pars(winddict,default=winddictDEF,type=dict)
    # Set flags
    flag_singonly=False ; flag_multonly=False ; flag_diskonly=False
    if   method == 'calc':
        flag_calc=True  ; flag_plot=False ; flag_anim=False ; flag_stat=False
    elif 'plot' in method:
        flag_calc=False ; flag_plot=True  ; flag_anim=False ; flag_stat=False
        if   method == 'plotgas'  : flag_singonly=True ; singlist=['gas'  ]
        elif method == 'plothalo' : flag_singonly=True ; singlist=['halo' ]
        elif method == 'plotdisk' : flag_singonly=True ; singlist=['disk' ]
        elif method == 'plotbulge': flag_singonly=True ; singlist=['bulge']
        elif method == 'plotstar' : flag_singonly=True ; singlist=['star' ]
        elif method == 'plotsing' : flag_singonly=True
        elif method == 'plotmult' : flag_multonly=True
        elif method == 'plotdisk2': flag_diskonly=True
        else: pass
    elif method == 'anim':
        flag_calc=False ; flag_plot=False ; flag_anim=True  ; flag_stat=False
    elif method == 'stat':
        flag_calc=False ; flag_plot=False ; flag_anim=False ; flag_stat=True
    else:
        flag_calc=True  ; flag_plot=True  ; flag_anim=True  ; flag_stat=True
    # Get file names
    files=simstr.mkfiledict(options=fileopt)
    parmfile=files[ftype]['input']['param']
    makefile=files[ftype]['input']['makefile']
    if ftype=='buildgal':
        raise Exception('Buildgal analysis is a work in progress.')
        snaplist=[]
        icfile=files[ftype]['output']['ic']
    else:
        snaplist=sorted(glob.glob(files[ftype]['output']['snapbase']+'*'),reverse=True)
        icfile=files[ftype]['input']['ic']
    if os.path.isfile(icfile):
        noic=0
        snaplist.insert(0,icfile)
    else:
        noic=1
    # Raise error if no files found
    if len(snaplist)==0:
        raise Exception('There are not any ic/snapshot files for the specified run.')
    # Set scaling for default snapshot
    if   simstr['runtyp'].lower() == 'galaxy'  : adjylim0=False
    elif simstr['runtyp'].lower() == 'interact': adjylim0=True
    # Files that contain default stuff
    snapfile0DEF=simstr.get_snapfile0(ftype=ftype)
    snapfile0=mmlpars.mml_pars(snapfile0,default=snapfile0DEF,type=str)
    calcfile0=files['calc']['calcstat']
    # Make directories
    mmlfiles.mkdirs(os.path.dirname(files['calc']['calclist']))
    mmlfiles.mkdirs(os.path.dirname(files['calc']['csnapbase']))
    mmlfiles.mkdirs(os.path.dirname(files['plots']['energy']))
    mmlfiles.mkdirs(os.path.dirname(files['plots']['statplot']))
    mmlfiles.mkdirs(os.path.dirname(files['plots']['multbase']))
    # Load calc list
    mmlio.verbose('Reading in statistics files',width=verbwid)
    genstatlist =rwstats('R',stattype='gen' ,statfile=files['calc']['genstat' ])
    profstatlist=rwstats('R',stattype='prof',statfile=files['calc']['profstat'])
    barstatlist_dens =rwstats('R',stattype='bar' ,statfile=files['calc']['barstat_dens' ])
    barstatlist_scfm2=rwstats('R',stattype='bar' ,statfile=files['calc']['barstat_scfm2'])
    barstatlist_ellip=rwstats('R',stattype='bar' ,statfile=files['calc']['barstat_ellip'])
#    if os.path.isfile(files['calc']['calclist']):
#        clist=mmlio.rwtable('R',files['calc']['calclist'])
#    else:
#        calckeys=['fname','time','comx','comy','comz','vcomx','vcomy','vcomz',
#                  'clause','sigma',
#                  'disk_A','disk_r','disk_z','halo_A','halo_r','bulge_A','bulge_r','bar_x','bar_y','bar_z','bar_amp','bar_phi','bar_m2']
#        clist={ikey:[] for ikey in calckeys}
    # Ask user what to do
    owcalcstat=False ; owcalc=False
    owgenstat=False ; owprofstat=False ; owbarstat_dens=False ; owbarstat_scfm2=False ; owbarstat_ellip=False
    owmult=False ; owsing=False     ; owdisk=False
    owanim=False ; owstatplot=False ; owenergy=False
    if overwrite:
        if flag_calc:
            owcalcstat=mmlio.yorn('Overwrite existing static calc file for this run?')
            owcalc=True
        if flag_stat:
            owgenstat=True
            owbarstat_dens=True
            owbarstat_scfm2=True
            owbarstat_ellip=True
            owprofstat=True
        if flag_plot:
            owsing=True
            owmult=True
            owdisk=True
            owstatplot=True
            owenergy=True
        if flag_anim:
            owanim=True
    else:
        if flag_calc:
            owcalcstat=mmlio.yorn('Overwrite existing static calc file for this run?')
            owcalc=mmlio.yorn('Overwrite existing calc files for each snapshot?')
        if flag_stat:
            owgenstat=mmlio.yorn('Overwrite existing general statistics entries?')
            owprofstat=mmlio.yorn('Overwrite existing profile statistics entries?')
            owbarstat_dens=mmlio.yorn('Overwrite existing bar statistics entries [density]?')
            owbarstat_scfm2=mmlio.yorn('Overwrite existing bar statistics entries [m=2 SCF]?')
            owbarstat_ellip=mmlio.yorn('Overwrite existing bar statistics entries [ellipticity]?')
        if flag_plot:
            owsing=mmlio.yorn('Overwrite existing sing plots for each snapshot?')
            owmult=mmlio.yorn('Overwrite existing mult plots for each snapshot?')
            owdisk=mmlio.yorn('Overwrite existing disk plots for each snapshot?')
            owstatplot=mmlio.yorn('Create statistics plots?')
            owenergy=mmlio.yorn('Create energy plot?')
        if flag_anim:
            owanim=mmlio.yorn('Overwrite existing animations?')
    # Gadget specific stuff
    if ftype == 'gadget':
        if owenergy:
            energyfile=files[ftype]['output']['energy']
            energyplot=files['plots']['energy']
            mmlgadget.plot_energy(energyfile,plotfile=energyplot,overwrite=owenergy)
            print 'Energy plot:'
            print '    {}'.format(energyplot)
    # Initialize calc dict
    if owcalcstat or not os.path.isfile(calcfile0):
        inb0=simstr.get_pnbody(fname=snapfile0,ftype=ftype,idgal2=idgal2)
#        inb0=pNbody.Nbody(simstr,p_name=snapfile0,ftype=ftype,unitsfile=parmfile,makefile=makefile,idgal2=idgal2)
        calcdict0=inb0.calc(calcfile0,adjustylim=adjylim0,overwrite=True)
        del inb0
    else:
        calcdict0=loadcalc(calcfile0)
    # Loop over snapshots
    if flag_calc or flag_plot:
        for fi in range(len(snaplist)):
            if shortlist and fi >= 2: break
            isnapfile=snaplist[fi]
            mmlio.verbose('Beginning to analyze: ',addval=os.path.basename(isnapfile),
                          border=True,width=verbwid,valpad=50)
            # Make ID string
            if noic==0 and fi==0:
                exts='_000ic'
                icflag=1
            else:
                exts=isnapfile.split(files[ftype]['output']['snapbase'])[1]
                icflag=0
            print exts
            # Create other file names
            icalcfile=files['calc']['csnapbase']+exts
            iscffile=files['scf']['disk']['output']['coefbase']+exts
            imultplot=files['plots']['multbase']+exts+'.'+files['plots']['plotext']
            idiskplot=files['plots']['diskbase']+exts+'.'+files['plots']['plotext']
            isingplot_xy=[]
            isingplot_xz=[]
            for ising in singlist:
                isingplot_xy.append(files['plots'][ising+'base_xy']+exts+'.'+files['plots']['plotext'])
                isingplot_xz.append(files['plots'][ising+'base_xz']+exts+'.'+files['plots']['plotext'])
            # Determine what files will be created
            flag_calcfile=False ; flag_calclist=False
            if flag_calc:
                if owcalc or not os.path.isfile(icalcfile): flag_calcfile=True
                if isnapfile not in  genstatlist['fname'] or owgenstat:  flag_genstat =True
                if isnapfile not in profstatlist['fname'] or owprofstat: flag_profstat=True
                if isnapfile not in  barstatlist_dens[ 'fname'] or owbarstat_dens : flag_barstat_dens =True
                if isnapfile not in  barstatlist_scfm2['fname'] or owbarstat_scfm2: flag_barstat_scfm2=True
                if isnapfile not in  barstatlist_ellip['fname'] or owbarstat_ellip: flag_barstat_ellip=True
#                if isnapfile not in clist['fname']:         flag_calclist=True
            # Determine what plots will be created
            flag_singplot=False ; flag_multplot=False ; flag_diskplot=False
            if flag_plot:
                if owmult or not os.path.isfile(imultplot): flag_multplot=True
                if owdisk or not os.path.isfile(idiskplot): flag_diskplot=True
                if owsing: flag_singplot=True
                else:
                    for ising_xy in isingplot_xy:
                        if not os.path.isfile(ising_xy): flag_singplot=True
                    for ising_xz in isingplot_xz:
                        if not os.path.isfile(ising_xz): flag_singplot=True
                if flag_singonly: flag_multplot=False ; flag_diskplot=False
                if flag_multonly: flag_singplot=False ; flag_diskplot=False
                if flag_diskonly: flag_singplot=False ; flag_multplot=False
            # Check if snapshot or calc file needs to be loaded
            flag_loadsnap=False ; flag_loadcalc=False
            if flag_calcfile or flag_singplot: flag_loadsnap=True
            if flag_genstat or flag_barstat or flag_profstat: flag_loadsnap=True
            if flag_multplot or flag_diskplot:
                flag_loadcalc=True
                if not os.path.isfile(icalcfile):
                    flag_loadsnap=True ; flag_calcfile=True
            # Load snapshot and perform analysis
            if flag_loadsnap:
                inb=simstr.get_pnbody(fname=isnapfile,ftype=ftype,idgal2=idgal2)
#                inb=pNbody.Nbody(simstr,p_name=isnapfile,ftype=ftype,unitsfile=parmfile,makefile=makefile,scffile=iscffile,idgal2=idgal2)
                # Creat calc file
                if flag_calcfile:
                    # Make calcfile
                    calcdict=inb.calc(icalcfile,overwrite=owcalc,calcdict=copy.deepcopy(calcdict0))
                    print 'Created calcfile:'
                    print '    '+icalcfile
                # Plot
                if flag_singplot:
                    for insing in range(len(singlist)):
                        if singlist[insing] == 'star': isingtypList=['disk','bulge','stars']
                        else:                          isingtypList=[singlist[insing]]
                        if owsing or not os.path.isfile(isingplot_xy[insing]):
                            inb.plot(isingplot_xy[insing],overwrite=owsing,typList=isingtypList,galaxy=0,
                                     view_plane='xy',view_window=winddict[singlist[insing]],view_palette='light')
                            print 'Created singular XY {} plot:'.format(singlist[insing])
                            print '    '+isingplot_xy[insing]
                        if owsing or not os.path.isfile(isingplot_xz[insing]):
                            inb.plot(isingplot_xz[insing],overwrite=owsing,typList=isingtypList,galaxy=0,
                                     view_plane='xz',view_window=winddict[singlist[insing]],view_palette='light')
                            print 'Created singular XZ {} plot:'.format(singlist[insing])
                            print '    '+isingplot_xz[insing]
                # Delete pNbody object
                del inb
            else:
                print 'Did not need to load snapshot:'
                print '    '+isnapfile
            # Load calcfile and plot
            if flag_loadcalc:
                calcdict=loadcalc(icalcfile)
##                 # Add things to calclist if flag set
##                 if flag_calclist:
##                     iclist=calcdict['statDict']
##                     for ikey in clist:
##                         clist[ikey].append(iclist[ikey])
##                     mmlio.rwtable('W',files['calc']['calclist'],clist,overwrite=True)
##                     print 'Calc table updated'
##                 # Load info from calclist if flag not set
##                 else:
##                     idxSnapfile=clist['fname'].index(isnapfile)
##                     iclist={}
##                     for ikey in clist:
##                         iclist[ikey]=clist[ikey][idxSnapfile]
##                     print 'Calc table loaded'
                # Make multi-panel plot
                if flag_multplot:
                    mmlsimplot.plot_snapshot(imultplot,calcdict=calcdict,calcdict0=calcdict0,overwrite=owmult,**plot_kw)
                    print 'Created multi-panel plot:'
                    print '    '+imultplot
                # Make disk plot
                if flag_diskplot:
                    mmlsimplot.plot_snaphist2D(idiskplot,typList=['disk'],rayList=['xy','xz'],pgal=0,
                                               calcdict=calcdict,calcdict0=calcdict0,overwrite=owdisk,**plot_kw)
                    print 'Created disk plot:'
                    print '    '+idiskplot
            else:
                print 'Did not need to load calcfile:'
                print '    '+icalcfile

    # Create plots
    if flag_plot:
        # Reiterate energy plot name
        if owenergy:
            mmlio.verbose('Created energy plot',border=True,width=verbwid)
            print '    '+files['plots']['energy']
        # Create stat plots
        if owstatplot:
            mmlio.verbose('Beginning statistics plot',border=True,width=verbwid)
            clist=mmlio.rwtable('R',files['calc']['calclist'])
            mmlsimplot.plot_stats(clist,files['plots']['statplot'],overwrite=owstatplot,**plot_kw)
            print 'Created statistics plot:'
            print '    '+files['plots']['statplot']
    # Create animation
    if flag_anim:
#        if owanim:
        mmlio.verbose('Beginning animations',border=True,width=verbwid)
        anim(files,owanim=owanim,verbwid=verbwid,**anim_kw)
        
    return

####################################################################################################################################
# METHOD TO PERFORM STATISTICS
def rwstats(method,stattype=None,statfile=None,statlist=None,addstatlist=None,owentry=None):
    # Set constants
    stattype_bar=['bar_dens','bar_scfm2','bar_ellip']
    stattypeLIST=['gen','prof']+stattype_bar
    # Pars input
    if isinstance(method,str): method=method.upper()
    method=mmlpars.mml_pars(method,type=str,list=['R','W','A'])
    if isinstance(stattype,str): stattype=stattype.lower()
    stattype=mmlpars.mml_pars(stattype,type=str,list=stattypeLIST)
    owentry=mmlpars.mml_pars(owentry,type=bool,default=False)
    # Get list of keys
    statkeys=['fname','time']
    if stattype=='gen':
        statkeys+=['comx','comy','comz','vcomx','vcomy','vcomz','clause','sigma']
    elif stattype=='prof':
        statkeys+=['disk_A','disk_r','disk_z','halo_A','halo_r','bulge_A','bulge_r']
    elif stattype in stattype_bar:
        statkeys+=['A','Aalt','pa','omega','r','x','y','z']
    # Read
    if method=='R':
        statfile=mmlpars.mml_pars(statfile,type=str)
        if os.path.isfile(statfile):
            statlist=mmlio.rwtable('R',statfile)
        else:
            statlist={ikey:[] for ikey in statkeys}
    # Append
    elif method=='A':
        addstatlist=mmlpars.mml_pars(addstatlist,type=dict)
        if statlist==None: statlist=rwgenstat('R',statfile=statfile)
        if not addstatlist['fname'] in statlist['fname'] or owentry:
            for ikey in statkeys:
                statlist[ikey].append(addstatlist[ikey])
            statlist=rwgenstat('W',statfile=statfile,statlist=statlist)
    # Write
    elif method=='W':
        statfile=mmlpars.mml_pars(statfile,type=str)
        statlist=mmlpars.mml_pars(statlist,type=dict)
        if not os.path.isfile(statfile):
            mmlfiles.mkdirs(os.path.dirname(statfile))
        statlist['keylist']=statkeys
        mmlio.rwtable('W',statfile,statlist,overwrite=True)
    # Return output
    return statlist

####################################################################################################################################
# METHOD FOR ANIMATIONS
def anim(files,owanim=False,verbwid=100,**anim_kw):
    singlistVAL=['gas','halo','disk','bulge','star']
    mmlfiles.mkdirs(files['plots']['animdir'])
    animList=[files['plots']['multbase'],files['plots']['diskbase']]
    for isingVAL in singlistVAL:
        animList.append(files['plots'][isingVAL+'base_xy'])
        animList.append(files['plots'][isingVAL+'base_xz'])
    for ianim in animList:
        # Get file names
        animFname=os.path.join(files['plots']['animdir'],os.path.basename(ianim)+'_anim.'+files['plots']['animext'])
        frameTemp=ianim+'_*.'+files['plots']['plotext']
        # Determine if frame exist
        if len(glob.glob(frameTemp)) == 0: continue
        # Ask if animation should be overwritten
        if not owanim and os.path.isfile(animFname): iowanim=mmlio.yorn('Overwrite {}?'.format(os.path.basename(animFname)),width=verbwid-9)
        else: iowanim=owanim
        # Create animation
        if iowanim or not os.path.isfile(animFname):
            print '    '+animFname
            movestr=mmlplot.ffmpeg(frameTemp,animFname,overwrite=iowanim,rmtemp=True,**anim_kw)

####################################################################################################################################
# METHOD TO LOAD A SIMULATION CALC FILE
def loadcalc(fname):
    fid=open(fname,'r')
    calcdict=cPickle.load(fid)
    return calcdict

####################################################################################################################################
# METHOD TO ANALYZE A SIMULATION
def analyze_classic(simstr,ftype='gadget',fileopt=None,center=None,limopt=None,objdict=None,
                    idgal2=None,overwrite=False,anim_kw={},**plot_kw):
    '''
    NAME:
        mmlsimcalc.analyze
    PURPOSE:
        To run analysis on a Nbody simulation.
    CALLING:
        analyze(simstr[,ftype="gadget"])
    ARGUMENTS:
        simstr:    mmlsim object describing basic parameters of the simulation
    KEYWORDS:
        ftype:     String specifying what type of simulation it was.
        fileopt:   Dict of file naming options (See icgen.mmlsim.mkfiledict)
        overwrite: Bool specifying whether or not to overwrite existing files
    '''
    # Get idgal2
    if simstr['runtyp'].lower() == 'interact':
        if idgal2 == None:
            icgal1=simstr['infodict']['simstr1'].mkfiledict()[ftype]['input']['ic']
            nb1=pNbody.Nbody(icgal1,ftype=ftype)
            idgal2=nb1.get_idgal2()
    # Get file names
    files=simstr.mkfiledict(options=fileopt)
    if ftype=='buildgal':
        raise Exception('Buildgal analysis is a work in progress.')
        snaplist=[]
        icfile=files[ftype]['output']['ic']
    else:
        snaplist=sorted(glob.glob(files[ftype]['output']['snapbase']+'*'),reverse=True)
        icfile=files[ftype]['input']['ic']
    if os.path.isfile(icfile):
        noic=0
        snaplist.insert(0,icfile)
    else:
        noic=1
    # Raise error if no files found
    if len(snaplist)==0:
        raise Exception('There are not any ic/snapshot files for the specified run.')
    # Make directories
    mmlfiles.mkdirs(os.path.dirname(files['calc']['calclist']))
    mmlfiles.mkdirs(os.path.dirname(files['plots']['energy']))
    mmlfiles.mkdirs(os.path.dirname(files['plots']['multbase']))
    mmlfiles.mkdirs(os.path.dirname(files['plots']['animdir']))
    # Load calc list
    if os.path.isfile(files['calc']['calclist']):
        clist=mmlio.rwtable('R',files['calc']['calclist'])
    else:
        calckeys=['fname','time','comx','comy','comz','vcomx','vcomy','vcomz',
                  'clause','sigma',
                  'disk_A','disk_r','disk_z','halo_A','halo_r','bulge_A','bulge_r','bar_x','bar_y','bar_z','bar_amp','bar_phi']
#                  'ekin','epot',
        clist={ikey:[] for ikey in calckeys}
    # Ask if existing snap plots should be overwritten
    if overwrite:
        owmult=True
        owanim=True
    else:
        owmult=mmlio.yorn('Overwrite existing multi-panel plot?')
        owanim=mmlio.yorn('Overwrite existing animations?')
    # Initialize ploting options
#    Im_opt=None ; Cm_opt={} ; Rm_opt=None
#    Iv_opt=None ; Cv_opt={} ; Rv_opt=None
    # Files
    parmfile=files[ftype]['input']['param']
    makefile=files[ftype]['input']['makefile']
    extList=[] ; multplotList=[]
    # Gadget specific stuff
    if ftype == 'gadget':
        owenergy=mmlio.yorn('Create new energy plot?')
        if owenergy:
            energyfile=files[ftype]['output']['energy']
            energyplot=files['plots']['energy']
            mmlgadget.plot_energy(energyfile,plotfile=energyplot,overwrite=True)
            print 'Energy plot:'
            print '    {}'.format(energyplot)
    # Loop over snapshots
    for fi in range(len(snaplist)):
        if fi >= 2: return
        isnapfile=snaplist[fi]
        mmlio.verbose('Beginning to analyze: ',addval=isnapfile,border=True)
        # Make ID string
        if noic==0 and fi==0:
            exts='_ic'
            icflag=1
        else:
            exts=isnapfile.split(files[ftype]['output']['snapbase'])[1]
            icflag=0
        print exts
        extList.append(exts)
        # Create other file names
        imultplot=files['plots']['multbase']+exts+'.'+files['plots']['plotext']
        multplotList.append(imultplot)
        # Check if snapshot needs to be loaded
        loadsnap=False
        if owmult or not os.path.isfile(imultplot):
            flag_multplot=True
        else:
            flag_multplot=False
        if isnapfile not in clist['fname']:
            flag_calclist=True
        else:
            flag_calclist=False
        if flag_multplot or flag_calclist: loadsnap=True
        # Load snapshot and perform analysis
        if loadsnap:
            inb=pNbody.Nbody(p_name=isnapfile,ftype=ftype,unitsfile=parmfile,makefile=makefile,idgal2=idgal2)
            # Add things to calclist if flag set
            if flag_calclist:
                iclist=inb.get_statlist()
                for ikey in clist:
                    clist[ikey].append(iclist[ikey])
                mmlio.rwtable('W',files['calc']['calclist'],clist,overwrite=True)
                print 'Calc table updated'
            # Load info from calclist if flag not set
            else:
                idxSnapfile=clist['fname'].index(isnapfile)
                iclist={}
                for ikey in clist:
                    iclist[ikey]=clist[ikey][idxSnapfile]
                print 'Calc table loaded'
            # Plot positions
            if flag_multplot:
                # Get limits if they dont exists
                if objdict == None and fi !=0:
                    inb_norm=pNbody.Nbody(p_name=snaplist[0],ftype=ftype,unitsfile=parmfile,makefile=makefile,idgal2=idgal2)
                    objdict=inb_norm.plot(multplotList[0],overwrite=True,objdict=objdict,**plot_kw)
                    del inb_norm
                # Plot
                objdict=inb.plot(imultplot,overwrite=owmult,objdict=objdict,**plot_kw)
                print 'Plot created:'
                print '    '+imultplot
            # Delete pNbody object
            del inb
        else:
            print 'Did not need to load snapshot:'
            print '    '+isnapfile
            
    # Create animation
    mmlfiles.mkdirs(files['plots']['animdir'])
    animList=[files['plots']['multbase']]
    print 'Animation:'
    for ianim in animList:
        animFname=os.path.join(files['plots']['animdir'],os.path.basename(ianim)+'_anim.'+files['plots']['animext'])
        frameTemp=ianim+'_*.'+files['plots']['plotext']
        movestr=mmlplot.ffmpeg(frameTemp,animFname,overwrite=owanim,rmtemp=True,**anim_kw)
        print '    '+animFname
    print 'Energy plot:'
    print '    {}'.format(files['plots']['energy'])









#  LocalWords:  contopt
