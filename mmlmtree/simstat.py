#!/usr/bin/python    
###################################################################################################################################
#
# MEAGAN LANG'S SIMSTAT METHODS
#
####################################################################################################################################
import matplotlib
#matplotlib.use('Agg')
import matplotlib.pyplot as plt
import sys,os,shutil,glob,copy,pprint,scipy,math,cPickle
import numpy as np
LIST_METHODS_RUN=['intcut','plotcut','mkvels','halovint','snaptabs','lastints','plotlastints']
LIST_METHODS_INT=['orbit']
LIST_METHODS=LIST_METHODS_RUN+LIST_METHODS_INT
from mmlutils import *
from mmlastro import *
import main as mmlmtree
import simlist,simfile,simutil

def main():
    """
    Provides command line access to STAT methods
    """
    out = mmlnbody.walk(mtype='simstat')
    return out

####################################################################################################################################
####################################################################################################################################
# METHODS FOR GETTING LISTS
def listpar(method,plot=None,**exkw):
    """
    Returns a list of supported parameters
    """
    method=mmlpars.mml_pars(method,type=str)
    plot=mmlpars.mml_pars(plot,default=False,type=bool)
    if   method=='intcut'   : 
        par={'keylist' :['type','mass','q','rsep','z'],
             'type'    :{'def':'ALL'      ,'form':simlist.LIST_INTTYPS},
             'mass'    :{'def':(100.,500.),'form':tuple},
             'q'       :{'def':(1.  ,100.),'form':tuple},
             'rsep'    :{'def':(0.  ,100.),'form':tuple},
             'z'       :{'def':(0.1 ,3.0 ),'form':tuple}}
        if plot: 
            ppar={}
            for ivar in ['x','y','clr','mrk','siz']:
                if ivar in ['x','y']: ivarstr=ivar
                else                : ivarstr=ivar+'var'
                par['keylist']+=[ivarstr,ivar+'lim',ivar+'scl']
                ppar[ivarstr]={'def':'otr','form':['otr']+simlist.LIST_INTPROP}
                ppar[ivar+'lim']={'def':(0.,0.),'form':tuple}
                ppar[ivar+'scl']={'def':'lin','form':mmlplot.LIST_SCALES}
            par.update(ppar)
    elif method=='inttab'   : par={'keylist' :['type','mass','q','rsep','z'],
                                   'type'    :{'def':'ALL'      ,'form':simlist.LIST_INTTYPS},
                                   'mass'    :{'def':(100.,500.),'form':tuple},
                                   'q'       :{'def':(1.  ,100.),'form':tuple},
                                   'rsep'    :{'def':(0.  ,100.),'form':tuple},
                                   'z'       :{'def':(0.1 ,3.0 ),'form':tuple}}
    elif method=='orbit'    : par={'cmeth':{'def':simlist.LIST_COMMETHS[0],'form':simlist.LIST_COMMETHS}}
    elif method=='spin'     : par={'cmeth':{'def':simlist.LIST_COMMETHS[0],'form':simlist.LIST_COMMETHS}}
    elif method=='halovint' : 
        par={'keylist':['cuttag'],
             'cuttag' :{'def':'ALL_mass0to0_q0to0_rsep0to0_z0to0','form':str}}
        if plot: 
            ppar={'x'     :{'def':'halo_tiso','form':['otr']+simlist.LIST_HALOVINT},
                  'y'     :{'def':'halo_spin','form':['otr']+simlist.LIST_HALOVINT},
                  'clrvar':{'def':'int_q'    ,'form':['otr']+simlist.LIST_HALOVINT},
                  'mrkvar':{'def':'int_tag'  ,'form':['otr']+simlist.LIST_HALOVINT},
                  'sizvar':{'def':'halo_mtot','form':['otr']+simlist.LIST_HALOVINT}}
            par['keylist']+=['x','y','clrvar','mrkvar','sizvar']
            par.update(ppar)
    elif method=='lastint'  : 
        par={'keylist':['cuttag','nback'],
             'cuttag':{'def':'ALL_mass0to0_q0to0_rsep0to0_z0to0','form':str},
             'nback':{'def':2,'form':int}}
        if plot: 
            ppar={}
            for ivar in ['x','y','clr','mrk','siz']:
                if ivar in ['x','y']: ivarstr=ivar
                else                : ivarstr=ivar+'var'
                par['keylist']+=[ivarstr,ivar+'lim',ivar+'scl']
                ppar[ivarstr]={'def':'otr','form':['otr']+simlist.LIST_LASTINTKEYS}
                ppar[ivar+'lim']={'def':(0.,0.),'form':tuple}
                ppar[ivar+'scl']={'def':'lin','form':mmlplot.LIST_SCALES}
            par.update(ppar)
    elif method=='snaptab'  : par={'keylist':[]}
    else: raise Exception('Invalid method: {}'.format(method))
    return par
def par2tag(method,inpar,**exkw):
    """
    Returns a file tag given the set of parameters
    """
    import simstat
    method=mmlpars.mml_pars(method,type=str)
    outpar=mmlparam.parspar(simstat,method,inpar)
    if   method=='intcut' :
        outpar['mcut']='{}'.format(mmlstring.lim2str(outpar['mass']))
        outpar['qcut']='{}'.format(mmlstring.lim2str(outpar['q'   ]))
        outpar['rcut']='{}'.format(mmlstring.lim2str(outpar['rsep']))
        outpar['zcut']='{}'.format(mmlstring.lim2str(outpar['z'   ]))
        tagstr='{type}_mass{mcut}_q{qcut}_rsep{rcut}_z{zcut}'.format(**outpar)
    elif method=='orbit'    : tagstr='c{cmeth}'.format(**outpar)
    elif method=='spin'     : tagstr='c{cmeth}'.format(**outpar)
    elif method=='halovint' : tagstr='{cuttag}'.format(**outpar)
    elif method=='lastint'  : tagstr='lastint_n{nback}_{cuttag}'.format(**outpar)
    elif method=='snaptab'  : tagstr='snaptab_{:03d}'
    else: raise Exception('Invalid method: {}'.format(method))
    return tagstr

####################################################################################################################################
####################################################################################################################################
# HIGH LEVEL METHODS REQUIRING MMLSIM

####################################################################################################################################
# METHOD FOR RUNNING DIFFERENT STAT METHODS
def run(simstr,method=None,**mkw):
    """
    Provides interface for running different STAT methods
    """
    if method not in LIST_METHODS: method=mmlio.askselect('Select a valid simstat method:',LIST_METHODS)
    out=None
    # Proceed based on method
    if method in LIST_METHODS_INT:
        inpar=mmlparam.parspar('mmlmtree.simstat',method,inpar=mkw,askuser=True,inclextra=True)
        if not hasattr(simstr,'intlist'): simstr.get_intlist()
        for idx in range(len(simstr.intlist['tag'])):
            iint={ikey:simstr.intlist[ikey][idx] for ikey in intlist['keylist']}
            if   method=='orbit': iout=plot_orbit(simstr,iint,**inpar)
            elif method=='spin' : iout=plot_spin(simstr,iint,**inpar)
            else: raise Exception('Invalid method: {}'.format(method))
    elif method=='halovint': out=halovint(simstr,**mkw)
    elif method=='intcut'  : out=simstr.get_intlist(**mkw)
    elif method=='plotcut' : out=plotcuts(simstr,**mkw)
    elif method=='mkvels'  : mkvels(simstr,**mkw)
    elif method=='snaptabs':
        for s in range(len(simstr.fdict['redlist'])):
            mmlio.verbose('Snapshot {:03d}'.format(s))
            it=snaptab(simstr,s,**mkw)
    elif method in ['lastints','plotlastints']:
        # Pars input
        mkw['par']=mmlparam.parspar('mmlmtree.simstat','lastint',askuser=True,
                                    inpar=mkw.get('par',None),
                                    tagstr=mkw.get('tag',None))
        if 'intlist' not in mkw:
            mkw['intlist']=simstr.get_intlist(cuttag=mkw['par']['cuttag'],mkarray=True)
        # Loop over snapshots
        for s in range(len(simstr.fdict['redlist'])):
            mmlio.verbose('Snapshot {:03d}'.format(s))
            if 'plot' in method: it=plotlastint(simstr,s,**mkw)
            else               : it=lastint(simstr,s,**mkw)
        # Make animation
        if 'plot' in method:
            fdict=snapfiles(simstr,'lastint',mkw['par'],plot=True)
            print fdict['file']
            if len(glob.glob(fdict['file']))>0:
                owanim=mmlio.yorn('Overwrite existing animation?')
                movestr=mmlplot.ffmpeg(fdict['file'],fdict['anim'],overwrite=owanim,rmtemp=True)
                print '    '+fdict['anim']
    else: raise Exception('Invalid method: {}'.format(method))
    # Return output
    return out

####################################################################################################################################
####################################################################################################################################
# FILE CLASSES AND METHODS
                                                                    
####################################################################################################################################
# METHOD TO RETURN DICTIONARY OF SIMCALC FILES FOR CALCULATING SNAPSHOT STATISTICS 
def snapfiles(simstr,method,fpar,plot=None,**exkw):
    """
    Returns a dictionary containing a snapshot hist directory and file base 
    """
    plot=mmlpars.mml_pars(plot,default=False,type=bool)
    if plot:
        exkw['plot']=plot
        if method=='halovint':
            ofpar={}
            for ikey in ['x','y','clrvar','mrkvar','sizvar']: 
                if '_' in fpar[ikey]: ofpar[ikey]=fpar[ikey].split('_')[1]
                else: ofpar[ikey]=fpar[ikey]
            exkw['plotstr0']='plot_{y}vs{x}_{clrvar}CLR{mrkvar}MRK{sizvar}SIZ_'.format(**ofpar)
            exkw['animstr0']='anim_{y}vs{x}_{clrvar}CLR{mrkvar}MRK{sizvar}SIZ_'.format(**ofpar)
        elif method in ['intcut','lastint']:
            addstr=''
            for ivar in ['x','y','clr','mrk','siz']:
                if ivar in ['x','y']: ivarstr=ivar
                else                : ivarstr=ivar+'var'
                if fpar[ivarstr]=='otr': fpar[ivarstr+'str']=''
                else                   : fpar[ivarstr+'str']=fpar[ivarstr]
                if fpar[ivar+'lim'][0]==fpar[ivar+'lim'][1]: fpar[ivar+'limstr']=''
                else: fpar[ivar+'limstr']=mmlstring.lim2str(fpar[ivar+'lim'])
                if fpar[ivar+'scl']=='lin': fpar[ivar+'sclstr']=''
                else: fpar[ivar+'sclstr']=fpar[ivar+'scl']
            if fpar['x']=='otr' and fpar['y']=='otr': raise Exception('Neither x nor y are set')
            xstr='{xsclstr}{xstr}{xlimstr}'.format(**fpar)
            ystr='{ysclstr}{ystr}{ylimstr}'.format(**fpar)
            if len(xstr)!=0 and len(ystr)!=0: addstr+='{}vs{}_'.format(ystr,xstr)
            elif len(xstr)!=0: addstr+='{}_'.format(xstr)
            elif len(ystr)!=0: addstr+='{}_'.format(ystr)
            for ivar in ['clr','mrk','siz']:
                ivstr=str('{'+ivar+'sclstr}{'+ivar+'varstr}{'+ivar+'limstr}').format(**fpar)
                if len(ivstr)!=0: addstr+=ivstr+ivar.upper()
            if not addstr.endswith('_'): addstr+='_'
            exkw['plotstr0']='plot_'+addstr.format(**fpar)
            exkw['animstr0']='anim_'+addstr.format(**fpar)
#            exkw['plotstr0']='plot_{y}vs{x}_{clrvar}CLR{mrkvar}MRK{sizvar}SIZ_'.format(**fpar)
#            exkw['animstr0']='anim_{y}vs{x}_{clrvar}CLR{mrkvar}MRK{sizvar}SIZ_'.format(**fpar)
        else: pass
    if method=='lastint': sing=False
    else                : sing=True
    snfdict=mmlparam.par2file(simstr,'mmlmtree.simstat',method,fpar,sing=sing,**exkw)
    return snfdict

####################################################################################################################################
# METHOD TO RETURN CALC FILES
def files(fdict):
    """
    Returns a dictionary of STAT files & directories
    """
    files={}
    # Return file dictionary
    fdict['stat']=files
    return fdict

####################################################################################################################################
# METHOD TO RETURN LIST OF FILES
def get_filelist(fdict,keylist):
    """
    Returns a list of relavent CALC files
    """
    # Pars input
    fdict=mmlpars.mml_pars(fdict,type=dict)
    keylist=mmlpars.mml_pars(keylist,type=list)
    # Add files
    filelist={}
    # Return output
    return filelist

###################################################################################################################################
# METHOD TO CREATE VELOCITIES
def mkvels(simstr,overwrite=None,**exkw):
    """
    Create groups velocity files
    """
    snapbase=simstr.fdict['snapbase'].split('{')[0]
    snaps=sorted(glob.glob(snapbase+'*'),reverse=True)
    method='sort' # Assumes sequential IDs
#    method=mmlio.askselect('Select a method for getting velocities:',['bigloop','sort'])
    for isnap in snaps:
        isnapnum=int(float(isnap.split(snapbase)[-1]))
        ifids=simstr.fdict['grpids'].format(isnapnum)
        ifvel=simstr.fdict['grpvel'].format(isnapnum)
        if os.path.isfile(ifvel) and not overwrite: continue
        # Get snapshot
        inb=simstr.get_snapshot(snapnum=isnapnum)
        # Get ids
        idids=simfile.readgrpids(ifids)
        # Get index
        if method=='bigloop':
            iidx=np.zeros(idids['ids'].shape,dtype=np.int64)
            for ip in range(len(idids['ids'])): iidx[ip]=inb.getindex(idids['ids'][ip])
            vels=inb.vel[iidx,:]
        elif method=='sort':
            idxnb_for=np.argsort(inb.num)
            print inb.num[idxnb_for[0]],inb.num[idxnb_for[-1]],inb.nbody
            vels=inb.vel[idxnb_for,:][idids['ids']-inb.num.min(),:]
        else: raise Exception('Invalid method: {}'.format(method))
        # Get velocities
        idvel={'Npart':idids['Npart'],'vel':vels}
        # Write velocities to file
        mmlfiles.mkdirs(os.path.dirname(ifvel))
        with open(ifvel,'wb') as file:
            np.array(idvel['Npart'],dtype=np.int32  ).tofile(file)
            np.array(idvel['vel'  ],dtype=np.float32).tofile(file)
        mmlio.verbose('Finished snapshot {}'.format(isnapnum))
    return

###################################################################################################################################
# METHOD TO ADD GROUP INFO TO THE SNAPSHOT TABLES
def snaptab(simstr,snapnum,overwrite=False,verbose=False,**exkw):
    """
    Appends info to the standard snapshot table
    """
    flag_noval=-999.
    # Get keys
    grpkeys=simlist.LIST_TABVARS
    grprads=simlist.LIST_TABRADS
    newkeys=[]
    for r in grprads:
        for k in grpkeys: 
            newkeys.append(k+'_'+r)
    # Get file names
    otabfile=snapfiles(simstr,'snaptab',{},fext=snapnum)['file']
    itabfile=simstr.fdict['tabbase'].format(snapnum)
    z=simstr.fdict['redlist'][snapnum]
    # Read
    if not overwrite and os.path.isfile(otabfile): 
        t=mmlio.rwtable('R',otabfile)
    # Create & write
    else:
        # Read in table & groups
        t=simfile.rwtable('R',itabfile,'snaptab')
        g=simstr.get_groups(snapnum,inclparticles=True)
        # Preallocate for added stuff
        t['keylist']+=newkeys
        t.update(**{k:flag_noval*np.ones(len(t['GroupNum']),dtype=float) for k in newkeys})
        # Loop over groups
        for grpnum in t['GroupNum']:
            # Get physical quatities for particles
            G=g['G'][grpnum]
            mass=np.ones(g['pids'][grpnum].shape,dtype=np.float32)*g['mtot'][grpnum]/g['N'][grpnum]
            ppos=simstr.convert_pos(g['ppos'][grpnum],z)
            pvel=simstr.convert_vel(g['pvel'][grpnum],z)
            # Center
            # cm =self.convert_pos(g['cm' ][grpnum],z)
            # cmv=self.convert_vel(g['cmv'][grpnum],z) 
            cm,cmv=mmlmath.center_mass(mass,ppos,vel=pvel,vflag=True)
            ppos=(ppos-cm ).astype(np.float32)
            pvel=(pvel-cmv).astype(np.float32)
            r=mmlmath.absvect(ppos)
            Rlist=[] ; Nzeroflag=False
            # Loop over radii
            for rkey in grprads:
                R=simstr.convert_pos(t[rkey][grpnum],z)
                Rlist.append(R)
                # Identify particles inside the radius
                ridx=(r<=R)
                N=ridx.sum()
                # Catch empty radii
                if N==0:
                    Nzeroflag=True
                    if verbose:
                        print 'N=0',t['Groupid'][grpnum],len(ridx),rkey,R,(r.min(),r.max())
                    continue
                # Otherwise fill in values
                else:
                    M=mass[ridx].sum()
                    V=np.sqrt(G*M/R)
                    Lxyz=np.sum(mmlmath.angmom(mass[ridx],ppos[ridx],pvel[ridx]),axis=0)
                    L=np.sqrt(np.sum(Lxyz**2))
                    H=L/M
                    spin=mmlmath.spin(H,V,R)
                # Add new values
                t['N_'+rkey][grpnum]=N
                t['m_'+rkey][grpnum]=M
                t['v_'+rkey][grpnum]=V
                t['h_'+rkey][grpnum]=H
                t['spin_'+rkey][grpnum]=spin
            # Plot halos
            ihplot=os.path.join(simstr.fdict['simstat']['dir'],'haloplot','{:03d}'.format(snapnum),'{:07d}.png'.format(t['Groupid'][grpnum]))
            plothalo(ihplot,ppos,Rlist,Rlabs=grprads,overwrite=True)
            if verbose and Nzeroflag: print '    '+ihplot
        # Ensure lists
        for k in t['keylist']: t[k]=list(t[k])
        # Write file
        if len(t['GroupNum'])!=0:
            mmlio.rwtable('W',otabfile,t,overwrite=overwrite)
    # Return table
    return t


###################################################################################################################################
# METHOD TO RETURN INFO ON THE LAST INTERACTION FOR EACH HALO IN A SNAPSHOT
def lastint(simstr,snapnum,overwrite=False,askuser=False,par=None,tag=None,
            intlist=None,**kwargs):
    """
    Gets info on the last interaction for each halo in a snapshot
    """
    flag_noval=-999.
    intkeys=simlist.LIST_INTVARS
    # Get parameters
    par=mmlparam.parspar('mmlmtree.simstat','lastint',inpar=par,tagstr=tag,askuser=askuser)
    # Get filename
    fname=snapfiles(simstr,'lastint',par,fext='_{:03d}'.format(snapnum))['file']
    # Read
    if not overwrite and os.path.isfile(fname): 
        t=mmlio.rwtable('R',fname)
    # Create
    else:
        # Get interaction list if not provided
        if not intlist: intlist=simstr.get_intlist(cuttag=par['cuttag'],mkarray=True)
        isnap=intlist['isnap']
        fsnap=intlist['fsnap']
        # Get snapshot range
        smin=snapnum-par['nback']
        # Get interaction indices
        intidx_past=np.where(np.logical_and(fsnap>=smin,fsnap<snapnum))[0]
        intidx_curr=np.where(np.logical_and(isnap<=snapnum,fsnap>=snapnum))[0]
        #intidx_both=np.logical_or(intidx_past,intidx_curr)
        # Get halo data
        t=snaptab(simstr,snapnum)
        t['keylist']+=intkeys+['LastIntSnap']
        t.update(**{k:flag_noval*np.ones(len(t['GroupNum']),dtype=type(intlist[k][0])) for k in intkeys})
        t['LastIntSnap']=-1*np.ones(len(t['GroupNum']),dtype=int)
        # Loop over halos
        for grpnum in t['GroupNum']:
            # Identify interactions with matching halo IDs
            if   t['Groupid'] in np.hstack(intlist['haloid1'][intidx_curr],intlist['haloid2'][intidx_curr]):
                hidx=intidx_curr[np.where(np.logical_or(intlist['haloid1'][intidx_curr]==t['Groupid'],
                                                        intlist['haloid2'][intidx_curr]==t['Groupid']))]
                hidx=hidx[intlist['q'][hidx]==intlist['q'][hidx].min()]
                t['LastIntSnap'][grpnum]=snapnum
            elif t['Groupid'] in np.hstack(intlist['haloid1'][intidx_past],intlist['haloid2'][intidx_past]):
                hidx=intidx_past[np.where(np.logical_or(intlist['haloid1'][intidx_past]==t['Groupid'],
                                                        intlist['haloid2'][intidx_past]==t['Groupid']))]
                hidx=hidx[intlist['fsnap'][hidx]==intlist['fsnap'][hidx].max()]
                hidx=hidx[intlist['q'][hidx]==intlist['q'][hidx].min()]
                t['LastIntSnap'][grpnum]=intlist['fsnap'][hidx].max()
            else: continue
            # Fill in data
            if len(hidx)==1: hidx=hidx[0]
            else:
                raise Exception('Multiple interactions with the same fsnap and q for grpnum {} in snapshot {}'.format(grpnum,s))
            for k in intkeys: t[grpnum]=intlist[k][hidx]
        # Ensure list
        for k in t['keylist']: t[k]=list(t[k])
        # Write files
        mmlio.rwtable('W',fname,t,overwrite=overwrite)
    # Return info
    return t

###################################################################################################################################
# METHOD FOR PLOTTING LAST INTERACTION STATISTICS
def plotlastint(simstr,snapnum,par={},tag=None,askuser=False,verbose=False,overwrite=False,
                owlastint=None,nbins=100,Nmin=1000,**exkw):
    """
    Plots last interaction info
    """
    # Set constants
    plt.close('all')
    mrkList=['o','s','d','^','v','<','>','*','8','p']
    clrList=['b','g','r','m','k','c','y']
    sizRang=(12.,30.)
    clrMap=plt.get_cmap('jet')
    exart=[]
    # Pars input
    par=mmlparam.parspar('mmlmtree.simstat','lastint',inpar=par,tagstr=tag,askuser=askuser,plot=True)
    if tag==None: tag=mmlparam.par2tag('mmlmtree.simstat','lastint',par)
    # Check for existing file
    plotfile=snapfiles(simstr,'lastint',par,fext='_{:03d}'.format(snapnum),plot=True)['file']
    if not mmlfiles.prep_write(plotfile,overwrite): return par
    # Get data
    grplist=lastint(simstr,snapnum,overwrite=owlastint,askuser=askuser,par=par,tag=tag,**exkw)
    for k in grplist['keylist']: grplist[k]=np.array(grplist[k])
    Ngrp=len(grplist['haloid1'])
    if verbose: mmlio.verbose('Ngrp={}'.format(Ngrp))
    # Get mask
    maskarr=np.array(grplist['N_Rvir'])<Nmin
    if verbose: mmlio.verbose('N<Nmin={}'.format(maskarr.sum()))
    for imkey in mask.keys():
        if imkey not in grplist: continue
        if isinstance(mask[imkey],list): imvallist=mask[imkey]
        else                           : imvallist=[mask[imkey]]
        imarr=np.array(grplist[imkey])
        for imval in imvallist:
            if verbose: mmlio.verbose('N ({}={})={}'.format(imkey,imval,(imarr==imval).sum()))
            maskarr=np.logical_or(maskarr,imarr==imval)
    if verbose: mmlio.verbose('Nmaskarr={}'.format(maskarr.sum()))
    # Mask ranges
    for ivar in ['x','y','clr','siz','mrk']:
        # Select variable
        if ivar in ['x','y']: imkey=par[ivar]
        else                : imkey=par[ivar+'var']
        if imkey not in grplist: continue
        # Mask invalid values
        imaskval=np.logical_not(np.isfinite(np.array(grplist[imkey])))
        maskvar=np.logical_or(maskarr,imaskval)
        # Mask invalid log values
        if   par[ivar+'scl']=='log': imasklog=(np.array(grplist[imkey])<=0)
        elif par[ivar+'scl']=='lin': imasklog=np.array(Ngrp*[False],dtype=bool)
        maskarr=np.logical_or(maskarr,imasklog)
        # Get masked array
        imarr=np.ma.array(grplist[imkey],mask=maskarr)
        if par[ivar+'scl']=='log': imarr=np.ma.log10(imarr)
        # Check limits
        ilim=par[ivar+'lim']
        if par[ivar+'lim'][0]==par[ivar+'lim'][1]: par[ivar+'lim']=(imarr.min(),imarr.max())
        imasklim=np.logical_or(imarr<par[ivar+'lim'][0],imarr>par[ivar+'lim'][1])
        maskarr=np.logical_or(maskarr,imasklog)
        # Print output
        if verbose:
            mmlio.verbose('Limits ({}={})={}'.format(ivar,imkey,par[ivar+'lim']))
            mmlio.verbose('Nmask (val) ({})={}'.format(ivar,imaskval.sum()))
            mmlio.verbose('Nmask (log) ({})={}'.format(ivar,imasklog.sum()))
            mmlio.verbose('Nmask (lim) ({})={}'.format(ivar,imasklim.sum()))
            mmlio.verbose('Nmask (var) ({})={}'.format(ivar,maskarr.sum()))
            mmlio.verbose('Nmask (att) ({})={}'.format(ivar,imarr.mask.sum()))
    if verbose: mmlio.verbose('Nmasklim={}'.format(maskarr.sum()))
    # Count plots
    if maskarr.sum()==Ngrp: raise Exception('All of the data is masked.')
    if verbose: mmlio.verbose('Nplt={}'.format(Ngrp-maskarr.sum()))
    # Get x & y variables
    if   par['x']=='otr' and par['y']=='otr': raise Exception('Neither x nor y variable has been set.')
    elif par['x']=='otr': xlab=par['y'] ; ylab='# of Halos' ; x=np.ma.array(grplist[xlab],mask=maskarr) ; y=None ; xlim=par['ylim'] ; yscl=par['xscl']
    elif par['y']=='otr': xlab=par['x'] ; ylab='# of Halos' ; x=np.ma.array(grplist[xlab],mask=maskarr) ; y=None ; xlim=par['xlim'] ; yscl=par['yscl']
    else:
        xlab=par['x'] ; x=np.ma.array(grplist[xlab],mask=maskarr)
        ylab=par['y'] ; y=np.ma.array(grplist[ylab],mask=maskarr)
    if x!=None and par['xscl']=='log': x=np.ma.log10(x) ; xlab='log '+xlab
    if y!=None and par['yscl']=='log': y=np.ma.log10(y) ; ylab='log '+ylab
    # Print info
    if x!=None and verbose: print xlab,x.min(),x.max()
    if y!=None and verbose: print ylab,y.min(),y.max()
    # Get indices for different splits
    nointIdx=(grplist['LastIntSnap']==-1)
    mergeIdx=np.zeros(len(nointIdx),dtype=bool)
    flybyIdx=np.zeros(len(nointIdx),dtype=bool)
    for tag in [0,1,2,4]: mergeIdx=np.logical_or(mergeIdx,grplist['tag']==tag)
    for tag in [3,5    ]: flybyIdx=np.logical_or(flybyIdx,grplist['tag']==tag)
    # Set up subplots
    plt.close('all')
    f,(nointAx,mergeAx,flybyAx) = plt.subplots(1,3,subplot_kw={'aspect':'equal'})
    # Histogram
    if y==None: 
        nointAx.hist(x[np.logical_and(np.logical_not(maskarr),nointIdx)],bins=nbins,range=xlim,log=(yscl=='log'))
        mergeAx.hist(x[np.logical_and(np.logical_not(maskarr),mergeIdx)],bins=nbins,range=xlim,log=(yscl=='log'))
        flybyAx.hist(x[np.logical_and(np.logical_not(maskarr),flybyIdx)],bins=nbins,range=xlim,log=(yscl=='log'))
    # Scatter plot
    else:
        # Color
        if par['clrvar'] in grplist: clab=par['clrvar'] ; c=np.ma.array(grplist[clab],mask=maskarr)
        else                       : clab='None'        ; c=np.ma.array(len(x)*['b'],mask=maskarr)
        if clab!='None' and par['clrscl']=='log': c=np.log10(c) ; clab='log '+clab
        if not isinstance(c,str): cpar=dict(cmap=clrMap,vmin=par['clrlim'][0],vmax=par['clrlim'][1])
        else                    : cpar={}
        # Size
        if par['sizvar'] in grplist: slab=par['sizvar'] ; s=np.ma.array(grplist[slab],mask=maskarr)
        else                       : slab='None'        ; s=np.ones(x.shape)
        if slab!='None' and par['sizscl']=='log': s=np.log10(s) ; slab='log '+slab
        smin=par['sizlim'][0] ; smax=par['sizlim'][1]
        if smin==smax: s=np.ones(x.shape)*sizRang[0]
        else         : s=(float(sizRang[1]-sizRang[0])/float(smax-smin))*(s-smin)+sizRang[0]
        # Marker
        if par['mrkvar'] in grplist: mlab=par['mrkvar'] ; m=np.ma.array(grplist[mlab],mask=maskarr)
        else                       : mlab='None'        ; m=np.ones(x.shape)
        if mlab!='None' and par['mrkscl']=='log': m=np.log10(m) ; mlab='log '+mlab
        muni=set(m[np.logical_not(maskarr)])
        mbin=float(par['mrklim'][1]-par['mrklim'][0])/float(len(mrkList))
        if len(muni)<=len(mrkList): mwind=muni
        else                      : mwind=[(par['mrklim'][0]+mbin*imrk,par['mrklim'][0]+mbin*(imrk+1)) for imrk in range(len(mrkList))]
        # Plot
        for im,idxm in zip(mwind,range(len(mwind))):
            if isinstance(im,tuple): imidx=np.ma.getdata(np.logical_and(m>=im[0],m<im[1]))
            else                   : imidx=np.ma.getdata(m==im)
            nointMidx=np.logical_and(nointIdx,imidx)
            mergeMidx=np.logical_and(mergeIdx,imidx)
            flybyMidx=np.logical_and(flybyIdx,imidx)
            ip=nointAx.scatter(x[nointMidx],y[nointMidx],s=s[nointMidx],c=c[nointMidx],marker=mrkList[idxm],
                           label='{}'.format(im),alpha=0.7,edgecolors='none',**cpar)
            ip=mergeAx.scatter(x[mergeMidx],y[mergeMidx],s=s[mergeMidx],c=c[mergeMidx],marker=mrkList[idxm],
                           label='{}'.format(im),alpha=0.7,edgecolors='none',**cpar)
            ip=flybyAx.scatter(x[flybyMidx],y[flybyMidx],s=s[flybyMidx],c=c[flybyMidx],marker=mrkList[idxm],
                           label='{}'.format(im),alpha=0.7,edgecolors='none',**cpar)
        # Colorbar
        if clab!='None': cb=plt.colorbar(ip,pad=0.01) ; cb.set_label(clab) #; exart.append(cb)
        # Legend
        if mlab!='None': lgd=plt.legend(title=mlab,loc=3,ncol=3,mode="expand",scatterpoints=1,
                                        borderaxespad=0.,bbox_to_anchor=(0., 1.01, 1., .101)) ; exart.append(lgd)
    plt.xlabel(xlab)
    plt.ylabel(ylab)
    if par['xlim'][0]!=par['xlim'][1]: plt.xlim(par['xlim'])
    if par['ylim'][0]!=par['ylim'][1]: plt.ylim(par['ylim'])
    plt.savefig(plotfile,bbox_extra_artists=tuple(exart),bbox_inches='tight')
    if verbose:
        mmlio.verbose('Created plot:')
        print '     '+plotfile
    # Return the dictionary
    return plotfile

###################################################################################################################################
# METHOD TO RETURN CUT
def mkcuts(simstr,cutpar=None,cuttag=None,overwrite=None,askuser=None):
    """
    Cuts interaction list based on specified params
    """
    # Pars input
    fdict=simstr.fdict
    cutpar=mmlparam.parspar('mmlmtree.simstat','intcut',inpar=cutpar,tagstr=cuttag,askuser=askuser)
    fcut=snapfiles(simstr,'inttab',cutpar)['file']
    if os.path.isfile(fcut) and not overwrite: inttab=simfile.rwtable('R',fcut,'inttab')
    else:
        # Translate cut in redshift to cut in snapshot
        if cutpar['z'][0]!=cutpar['z'][1]:
            zarr0=np.array(fdict['redlist'])
            zidx0=np.logical_and(zarr0>=cutpar['z'][0],zarr0<=cutpar['z'][1])
            if zidx0.sum()==0: raise Exception('No snapshots found in redshift range: {}'.format(cutpar['z']))
            else:
                zint0=np.arange(len(zidx0))[zidx0]
                cutsnap=(zint0.min(),zint0.max())
        else: cutsnap=(0,0)
        # Get z=0 data
        tab_z0=simfile.rwtable('R',fdict['tab0'],'snaptab')
        # Get halos in the correct mass range (@ z=0)
        if cutpar['mass'][0]!=cutpar['mass'][1]:
            marr0=np.array(tab_z0['Mtot'])
            midx0=np.logical_and(marr0>=cutpar['mass'][0],marr0<=cutpar['mass'][1])
            if midx0.sum()==0: raise Exception('No halos were found in the mass range: {}'.format(cutpar['mass']))
            else:
                for ikey in tab_z0['keylist']: tab_z0[ikey]=list(np.array(tab_z0[ikey])[midx0])
        # Get full list of interactions
        inttab=simfile.rwtable('R',fdict['ints'],'inttab')
        # Cut by interaction type
        if cutpar['type']!='ALL':
            tarr=np.array(inttab['Tag'])
            tidx=np.array([False]*len(tarr))
            for ityp in simlist.DICT_INTTYPS[cutpar['type']]: tidx=np.logical_or(tidx,tarr==ityp)
            if tidx.sum()==0: raise Exception('No interactions were found of the type: {}'.format(cutpar['type']))
            else:
                for ikey in inttab['keylist']: inttab[ikey]=list(np.array(inttab[ikey])[tidx])
        # Cut by starting redshift
        if cutsnap[0]!=cutsnap[1]:
            sarr=np.array(inttab['Snapshot'])
            sidx=np.logical_and(sarr>=cutsnap[0],sarr<=cutsnap[1])
            if sidx.sum()==0: raise Exception('No interactions were found in the snapshot range: {}'.format(cutsnap))
            else:
                for ikey in inttab['keylist']: inttab[ikey]=list(np.array(inttab[ikey])[sidx])
        # Cut by ending redshift
        if cutpar['z'][0]!=cutpar['z'][1]:
            zendarr=np.array(inttab['DestructionZ'])
            zendidx=np.logical_and(zendarr>=cutpar['z'][0],zendarr<=cutpar['z'][1])
            if zendidx.sum()==0: raise Exception('No interactions were found ending in the redshift range: {}'.format(cutpar['z']))
            else:
                for ikey in inttab['keylist']: inttab[ikey]=list(np.array(inttab[ikey])[zendidx])
        # Cut by minimum separation
        if cutpar['rsep'][0]!=cutpar['rsep'][1]:
            rarr=np.array(inttab['MinRsep'])
            ridx=np.logical_and(rarr>=cutpar['rsep'][0],rarr<=cutpar['rsep'][1])
            if ridx.sum()==0: raise Exception('No interactions were found with a minimum separation in the range: {}'.format(cutpar['rsep']))
            else:
                for ikey in inttab['keylist']: inttab[ikey]=list(np.array(inttab[ikey])[ridx])
        # Cut by mass ratio
        if cutpar['q'][0]!=cutpar['q'][1]:
            qarr=np.array(inttab['ParentMtot'])/np.array(inttab['SubMtot'])
            qidx=np.logical_and(qarr>=cutpar['q'][0],qarr<=cutpar['q'][1])
            if qidx.sum()==0: raise Exception('No interactions were found with a mass ratio in the range: {}'.format(cutpar['q']))
            else:
                for ikey in inttab['keylist']: inttab[ikey]=list(np.array(inttab[ikey])[qidx])
        # Cut by mass
        if cutpar['mass'][0]!=cutpar['mass'][1]:
            marr=np.array(inttab['ParentHaloid'])
            midx=np.array([False]*len(marr))
            for ihid in tab_z0['FofHaloid']: midx=np.logical_or(midx,marr==ihid)
            if midx.sum()==0: raise Exception('No interactions were found with a mass in the range: {}'.format(cutpar['mass']))
            else:
                for ikey in inttab['keylist']: inttab[ikey]=list(np.array(inttab[ikey])[midx])
        # Write parsed inttab to file
        simfile.rwtable('W',fcut,'inttab',tabdict=inttab,overwrite=overwrite)
        mmlio.verbose('Found {} interactions matching specified parameters.'.format(len(inttab['ParentHaloid'])))
    # Return output
    setattr(simstr,'cutdict',cutpar )
    setattr(simstr,'inttab',inttab)
    return inttab


###################################################################################################################################
# METHOD TO GET HALO DATA FOR INTERACTIONS
def mkcuts_added(simstr,cutpar=None,cuttag=None,askuser=False,verbose=False,overwrite=False,
                 owadd=False,mkarray=False):
    """
    Loads/creates interaction info and adds more info to the table
    """
    flag_noval=-999.
    flag_inval=-666.
    # Pars input
    cutpar=mmlparam.parspar('mmlmtree.simstat','intcut',inpar=cutpar,tagstr=cuttag,askuser=askuser)
    if cuttag==None: cuttag=mmlparam.par2tag('mmlmtree.simstat','intcut',cutpar)
    # Load/create cut data
    fcut=snapfiles(simstr,'intcut',cutpar)['file']
    if verbose: mmlio.verbose(mmlstring.val2str([os.path.isfile(fcut),fcut]))
    if not overwrite and os.path.isfile(fcut): intlist=simfile.rwtable('R',fcut,'intlist')
    else:
        if verbose: mmlio.verbose('Loading table...')
        inttab=mkcuts(simstr,cutpar=cutpar,cuttag=cuttag,askuser=askuser)#,overwrite=overwrite)
        # Transfer variables
        if verbose: mmlio.verbose('Converting data...')
        nint=len(inttab['ParentHaloid'])
        intlist={'keylist':mmlparam.listpar(simlist,'intpar')['keylist']}
        intlist['iz']=list(np.array(simstr.fdict['redlist'])[np.array(inttab['Snapshot'])])
        intlist['fz']=inttab['DestructionZ']
        intlist['isnap']=inttab['Snapshot']
        intlist['haloid1']=inttab['ParentHaloid']
        intlist['haloid2']=inttab['Subhaloid']
        intlist['m1'     ]=inttab['ParentMtot']
        intlist['rvir1'  ]=inttab['ParentRvir']
        intlist['q'      ]=list(np.array(inttab['ParentMtot'])/np.array(inttab['SubMtot']))
        intlist['tag'    ]=inttab['Tag']
        intlist['minrsep']=inttab['MinRsep']
        intlist.update(fsnap=[],rvir1_phys=[])#,m1_z0=[])
        for i in range(nint):
            if intlist['fz'][i] in simstr.fdict['redlist']: intlist['fsnap'].append(simstr.fdict['redlist'].index(intlist['fz'][i]))
            else                                          : intlist['fsnap'].append(simstr.fdict['snap0'])
            # if intlist['haloid1'][i] in tab_z0['Groupid']:
            #     intlist['m1_z0'].append(tab_z0['Mtot'][tab_z0['Groupid'].index(intlist['haloid1'][i])])
            # else: intlist['m1_z0'].append(0.)
            intlist['rvir1_phys'].append(simstr.convert_pos(intlist['rvir1'][i],intlist['iz'][i]))
        # Check for missing keys
        misskeys=[]
        for ikey in intlist['keylist']:
            if ikey not in intlist: misskeys.append(ikey)
        if len(misskeys)>0:
            print misskeys
            mmlio.verbose('The above keys are missing')
    # Check for missing halo fields 
    misspar=False
    for ipar in simlist.LIST_INTHALO:
        if ipar not in intlist: misspar=True
    # Add extra data
    if misspar or owadd:
        # Get relavent data
        isnaparr=np.array(intlist['isnap'])
        fsnaparr=np.array(intlist['fsnap'])
        haloidarr=np.array(intlist['haloid1'])
        # Preallocate for spin data
        for ihpar in simlist.LIST_HALOPROP:
            intlist['i'+ihpar]=flag_noval*np.ones(isnaparr.shape)
            intlist['f'+ihpar]=flag_noval*np.ones(fsnaparr.shape)
        # Loop over snapshots
#        for snapnum in set([31]):
        for snapnum in set(intlist['isnap']+intlist['fsnap']):
            if verbose: mmlio.verbose('Beginning snapshot {}'.format(snapnum))
            # Identify starting and ending snapshots
            idxsnap=np.arange(len(isnaparr))[np.logical_or(isnaparr==snapnum,fsnaparr==snapnum)]
            idx_isnap=(isnaparr[idxsnap]==snapnum) ; idx_iflag=(idx_isnap.sum()!=0)
            idx_fsnap=(fsnaparr[idxsnap]==snapnum) ; idx_fflag=(idx_fsnap.sum()!=0)
            # Continue if there arn't any interactions that begin/end here
            if len(idxsnap)==0: continue
            # Get group numbers
            grpnumarr=np.array(simutil.haloid2groupnum(simstr,snapnum=snapnum,haloid=list(haloidarr[idxsnap])))
            if np.any(grpnumarr==-1): mmlio.yorn('There are missing group numbers.')
            # Load group data
            grps=simstr.get_groups(snapnum=snapnum,inclparticles=True)
            # Add calculated halo data
            iarrdict={'Nvir':[],'mvir':[],'Vvir':[],'hvir':[],'spin_vir':[]}
            for igrpnum,ihaloid,irvir in zip(grpnumarr,list(haloidarr[idxsnap]),np.array(intlist['rvir1_phys'])[idxsnap]):
                iridx=(mmlmath.absvect(grps['ppos_phys'][igrpnum])<=irvir)
                iNvir=iridx.sum()
                iarrdict['Nvir'].append(iNvir)
                ir=mmlmath.absvect(grps['ppos_phys'][igrpnum])
                # Plot halos
                ihplot=os.path.join(simstr.fdict['simstat']['dir'],'haloplot','{:03d}'.format(snapnum),'{:07d}.png'.format(ihaloid))
                plothalo(ihplot,grps['ppos_phys'][igrpnum],irvir,overwrite=True)
                # Handle virial properties
                if iNvir==0:
                    if verbose:
                        print ihaloid,len(iridx),iNvir,iarrdict['mvir'][-1],irvir,ir.min(),ir.max()
                        print '    '+ihplot
                    iarrdict['mvir'].append(flag_noval)
                    iarrdict['Vvir'].append(flag_noval)
                    iarrdict['hvir'].append(flag_noval)
                    iarrdict['spin_vir'].append(flag_noval)
                else:
                    iarrdict['mvir'].append(grps['pmass'][igrpnum][iridx].sum())
                    iarrdict['Vvir'].append(np.sqrt(grps['G'][igrpnum]*iarrdict['mvir'][-1]/irvir))
                    iLxyz=np.sum(mmlmath.angmom(grps['pmass'][igrpnum][iridx],grps['ppos_phys'][igrpnum][iridx],grps['pvel_phys'][igrpnum][iridx]),axis=0)
                    iarrdict['hvir'].append(np.sqrt(np.sum(iLxyz**2))/iarrdict['mvir'][-1])
                    iarrdict['spin_vir'].append(mmlmath.spin(iarrdict['hvir'][-1],iarrdict['Vvir'][-1],irvir))
            # Add read halo data
            for ihpar in simlist.LIST_HALOPROP:
                if ihpar in iarrdict: iarrdict[ihpar]=np.array(iarrdict[ihpar])
                else                : iarrdict[ihpar]=np.array(grps[ihpar])[grpnumarr]
                if idx_iflag: intlist['i'+ihpar][idxsnap[idx_isnap]]=iarrdict[ihpar][idx_isnap]
                if idx_fflag: intlist['f'+ihpar][idxsnap[idx_fsnap]]=iarrdict[ihpar][idx_fsnap]
        # Check and format data
        for ihpar in simlist.LIST_INTHALO:
            # Check that arrays are filled in
            if np.any(intlist[ihpar]==flag_noval): mmlio.verbose('There are missing {} values'.format(ihpar))
            # Add keys and ensure list format
            if ihpar not in intlist['keylist']: intlist['keylist'].append(ihpar)
            intlist[ihpar]=list(intlist[ihpar])
        # Write updated list to file
        fcut=snapfiles(simstr,'intcut',cutpar)['file']
        simfile.rwtable('W',fcut,'intlist',tabdict=intlist,overwrite=True)
    # Add interaction timescale
    intlist['tint']=mmlcosmo.t_lookback(zi=intlist['iz'],zf=intlist['fz'],cosmo='gadget2')
    # Add differences
    for ihpar in simlist.LIST_HALOPROP:
        iarr=np.array(intlist['i'+ihpar]).astype(float)
        farr=np.array(intlist['f'+ihpar]).astype(float)
        diffarr=np.array(len(iarr)*[float('Inf')],dtype=float)
        inz=(iarr!=0)
        diffarr[inz]=(farr[inz]-iarr[inz])/iarr[inz]
        intlist[ihpar+'diff']=list(diffarr)
    # Get angular momentum angle difference
    iLxyzarr=np.vstack((np.array(intlist['iLx']),np.array(intlist['iLy']),np.array(intlist['iLz']))).T
    fLxyzarr=np.vstack((np.array(intlist['fLx']),np.array(intlist['fLy']),np.array(intlist['fLz']))).T
    idfL=np.sum(iLxyzarr*fLxyzarr,axis=1)
    idiL=np.sum(iLxyzarr*iLxyzarr,axis=1)
    fdfL=np.sum(fLxyzarr*fLxyzarr,axis=1)
    intlist['Lphidiff']=list(np.arccos(idfL/np.sqrt(idiL*fdfL)))
    if verbose:
        # print 'nopart',np.sum(np.array(intlist['iNvir'])==0),np.sum(np.array(intlist['fNvir'])==0)
        print 'noval ',np.sum(np.array(intlist['imvir'])==flag_noval),np.sum(np.array(intlist['fmvir'])==flag_noval)
        print 'zero  ',np.sum(np.array(intlist['imvir'])==0),np.sum(np.array(intlist['fmvir'])==0)
    # Convert to arrays
    if mkarray:
        for k in intlist['keylist']: intlist[k]=np.array(intlist[k])
    # Return data
    return intlist

###################################################################################################################################
# METHOD TO PLOT PARTICLE POSITIONS AND RADIUS
def plothalo(fname,pos,R,Rlabs=[],overwrite=False):
    """
    Plots particle positions with an enclosing radius
    """
    plab=['x','y','z']
    plotlist=['xy','xz','yz']
    clrlist=['k','b','r','g','m','o']
    # Check for file
    if os.path.isfile(fname) and not overwrite: return
    mmlfiles.mkdirs(os.path.dirname(fname))
    # Pars input
    if not isinstance(R,list): R=[R]
    if len(Rlabs)!=len(R): Rlabs=[]
    # Initialize stuff
    plt.close('all')
    f,axs=plt.subplots(1,3,subplot_kw={'aspect':'equal'})
    # Get limits
    pmax=max(np.abs(pos).max(),max(R))
    plim=(-pmax,+pmax)
    # Create circle(s)
    phi=np.linspace(0,2.*np.pi,100)
    x_circ=[] ; y_circ=[]
    for iR in R:
        x_circ.append(iR*np.cos(phi))
        y_circ.append(iR*np.sin(phi))
    # Plot
    for idx,plot in enumerate(plotlist):
        xvar=plot[0] ; yvar=plot[1]
        # Get data
        x=pos[:,plab.index(xvar)]
        y=pos[:,plab.index(yvar)]
        # Plot points
        axs[idx].scatter(x,y,c=clrlist[0])
        # Plot circles
        for idxR in range(len(R)):
            iRlab=Rlabs[idxR] if len(R)==len(Rlabs) else None
            axs[idx].plot(x_circ[idxR],y_circ[idxR],clrlist[idxR+1],label=iRlab)
        # Set labels and limits
        axs[idx].set_xlabel(xvar)
        axs[idx].set_ylabel(yvar)
        axs[idx].set_xlim(plim)
        axs[idx].set_ylim(plim)
        # Hide tick labels
        axs[idx].xaxis.set_ticklabels([])
        if idx!=0: axs[idx].yaxis.set_ticklabels([])
    # Add legend to last axs
    if len(Rlabs)!=0: axs[-1].legend()
    # Save figure
    f.savefig(fname)
    # Return control
    return

###################################################################################################################################
# METHOD FOR PLOTTING INTERACTION CUTS
def plotcuts(simstr,cutpar=None,cuttag=None,askuser=None,verbose=False,overwrite=None,
             owcut=None,nbins=None,mask={},Nmin=1000,**exkw):
    """
    Plots interaction info
    """
    # Set constants
    plt.close('all')
    mrkList=['o','s','d','^','v','<','>','*','8','p']
    clrList=['b','g','r','m','k','c','y']
    sizRang=(12.,30.)
    clrMap=plt.get_cmap('jet')
    exart=[]
    # Pars input
    cutpar=mmlpars.mml_pars(cutpar,default={},type=dict)
    askuser=mmlpars.mml_pars(askuser,default=False,type=bool)
    overwrite=mmlpars.mml_pars(overwrite,default=False,type=bool)
    nbins=mmlpars.mml_pars(nbins,default=100,type=int,min=1)
    cutpar=mmlparam.parspar('mmlmtree.simstat','intcut',inpar=cutpar,tagstr=cuttag,askuser=askuser,plot=True)
    if cuttag==None: cuttag=mmlparam.par2tag('mmlmtree.simstat','intcut',cutpar)
    # Check for existing file
    plotfile=snapfiles(simstr,'intcut',cutpar,plot=True)['file']
    if not mmlfiles.prep_write(plotfile,overwrite): return cutpar
    # Get data
    intlist=simstr.get_intlist(cutpar=cutpar,cuttag=cuttag,overwrite=owcut,verbose=verbose,askuser=askuser,**exkw)
    Nint=len(intlist['haloid1'])
    if verbose: mmlio.verbose('Nint={}'.format(Nint))
    # Get mask
    maskarr=np.logical_or(np.array(intlist['iN'])<Nmin,np.array(intlist['fN'])<Nmin)
    if verbose: mmlio.verbose('N<Nmin={}'.format(maskarr.sum()))
    for imkey in mask.keys():
        if imkey not in intlist: continue
        if isinstance(mask[imkey],list): imvallist=mask[imkey]
        else                           : imvallist=[mask[imkey]]
        imarr=np.array(intlist[imkey])
        for imval in imvallist:
            if verbose: mmlio.verbose('N ({}={})={}'.format(imkey,imval,(imarr==imval).sum()))
            maskarr=np.logical_or(maskarr,imarr==imval)
            # maskarr=np.logical_or(maskarr,abs(imarr-imval)<=1.0e-3)
    if verbose: mmlio.verbose('Nmaskarr={}'.format(maskarr.sum()))
    # Mask ranges
    for ivar in ['x','y','clr','siz','mrk']:
        # Select variable
        if ivar in ['x','y']: imkey=cutpar[ivar]
        else                : imkey=cutpar[ivar+'var']
        if imkey not in intlist: continue
        # Mask invalid values
        imaskval=np.logical_not(np.isfinite(np.array(intlist[imkey])))
        maskvar=np.logical_or(maskarr,imaskval)
        # Mask invalid log values
        if   cutpar[ivar+'scl']=='log': imasklog=(np.array(intlist[imkey])<=0)
        elif cutpar[ivar+'scl']=='lin': imasklog=np.array(Nint*[False],dtype=bool)
        maskarr=np.logical_or(maskarr,imasklog)
        # Get masked array
        imarr=np.ma.array(intlist[imkey],mask=maskarr)
        if cutpar[ivar+'scl']=='log': imarr=np.ma.log10(imarr)
        # Check limits
        ilim=cutpar[ivar+'lim']
        if cutpar[ivar+'lim'][0]==cutpar[ivar+'lim'][1]: cutpar[ivar+'lim']=(imarr.min(),imarr.max())
        imasklim=np.logical_or(imarr<cutpar[ivar+'lim'][0],imarr>cutpar[ivar+'lim'][1])
        maskarr=np.logical_or(maskarr,imasklog)
        # Print output
        if verbose:
            mmlio.verbose('Limits ({}={})={}'.format(ivar,imkey,cutpar[ivar+'lim']))
            mmlio.verbose('Nmask (val) ({})={}'.format(ivar,imaskval.sum()))
            mmlio.verbose('Nmask (log) ({})={}'.format(ivar,imasklog.sum()))
            mmlio.verbose('Nmask (lim) ({})={}'.format(ivar,imasklim.sum()))
            mmlio.verbose('Nmask (var) ({})={}'.format(ivar,maskarr.sum()))
            mmlio.verbose('Nmask (att) ({})={}'.format(ivar,imarr.mask.sum()))
    if verbose: mmlio.verbose('Nmasklim={}'.format(maskarr.sum()))
    # Count plots
    if maskarr.sum()==Nint: raise Exception('All of the data is masked.')
    if verbose: mmlio.verbose('Nplt={}'.format(Nint-maskarr.sum()))
    # Get x & y variables
    if   cutpar['x']=='otr' and cutpar['y']=='otr': raise Exception('Neither x nor y variable has been set.')
    elif cutpar['x']=='otr': xlab=cutpar['y'] ; ylab='# of Halos' ; x=np.ma.array(intlist[xlab],mask=maskarr) ; y=None ; xlim=cutpar['ylim'] ; yscl=cutpar['xscl']
    elif cutpar['y']=='otr': xlab=cutpar['x'] ; ylab='# of Halos' ; x=np.ma.array(intlist[xlab],mask=maskarr) ; y=None ; xlim=cutpar['xlim'] ; yscl=cutpar['yscl']
    else:
        xlab=cutpar['x'] ; x=np.ma.array(intlist[xlab],mask=maskarr)
        ylab=cutpar['y'] ; y=np.ma.array(intlist[ylab],mask=maskarr)
    if x!=None and cutpar['xscl']=='log': x=np.ma.log10(x) ; xlab='log '+xlab
    if y!=None and cutpar['yscl']=='log': y=np.ma.log10(y) ; ylab='log '+ylab
    # Print info
    if x!=None and verbose: print xlab,x.min(),x.max()
    if y!=None and verbose: print ylab,y.min(),y.max()
    # Histogram
    if y==None: plt.hist(x[np.logical_not(maskarr)],bins=nbins,range=xlim,log=(yscl=='log'))
    # Scatter plot
    else:
        # Color
        if cutpar['clrvar'] in intlist: clab=cutpar['clrvar'] ; c=np.ma.array(intlist[clab],mask=maskarr)
        else                          : clab='None'           ; c=np.ma.array(len(x)*['b'],mask=maskarr)
        if clab!='None' and cutpar['clrscl']=='log': c=np.log10(c) ; clab='log '+clab
        if not isinstance(c,str): cpar=dict(cmap=clrMap,vmin=cutpar['clrlim'][0],vmax=cutpar['clrlim'][1])
        else                    : cpar={}
        # Size
        if cutpar['sizvar'] in intlist: slab=cutpar['sizvar'] ; s=np.ma.array(intlist[slab],mask=maskarr)
        else                          : slab='None'           ; s=np.ones(x.shape)
        if slab!='None' and cutpar['sizscl']=='log': s=np.log10(s) ; slab='log '+slab
        smin=cutpar['sizlim'][0] ; smax=cutpar['sizlim'][1]
        if smin==smax: s=np.ones(x.shape)*sizRang[0]
        else         : s=(float(sizRang[1]-sizRang[0])/float(smax-smin))*(s-smin)+sizRang[0]
        # Marker
        if cutpar['mrkvar'] in intlist: mlab=cutpar['mrkvar'] ; m=np.ma.array(intlist[mlab],mask=maskarr)
        else                          : mlab='None'           ; m=np.ones(x.shape)
        if mlab!='None' and cutpar['mrkscl']=='log': m=np.log10(m) ; mlab='log '+mlab
        muni=set(m[np.logical_not(maskarr)])
        mbin=float(cutpar['mrklim'][1]-cutpar['mrklim'][0])/float(len(mrkList))
        if len(muni)<=len(mrkList): mwind=muni
        else                      : mwind=[(cutpar['mrklim'][0]+mbin*imrk,cutpar['mrklim'][0]+mbin*(imrk+1)) for imrk in range(len(mrkList))]
#        if len(muni)>len(mrkList): raise Exception('Variable {} has more values ({}) than the marker list ({})'.format(mlab,len(muni),len(mrkList)))
        # Plot
        for im,idxm in zip(mwind,range(len(mwind))):
            if isinstance(im,tuple): imidx=np.ma.getdata(np.logical_and(m>=im[0],m<im[1]))
            else                   : imidx=np.ma.getdata(m==im)
            ip=plt.scatter(x[imidx],y[imidx],s=s[imidx],c=c[imidx],marker=mrkList[idxm],
                           label='{}'.format(im),alpha=0.7,edgecolors='none',**cpar)
        # Colorbar
        if clab!='None': cb=plt.colorbar(ip,pad=0.01) ; cb.set_label(clab) #; exart.append(cb)
        # Legend
        if mlab!='None': lgd=plt.legend(title=mlab,loc=3,ncol=3,mode="expand",scatterpoints=1,
                                        borderaxespad=0.,bbox_to_anchor=(0., 1.01, 1., .101)) ; exart.append(lgd)
    plt.xlabel(xlab)
    plt.ylabel(ylab)
    if cutpar['xlim'][0]!=cutpar['xlim'][1]: plt.xlim(cutpar['xlim'])
    if cutpar['ylim'][0]!=cutpar['ylim'][1]: plt.ylim(cutpar['ylim'])
    plt.savefig(plotfile,bbox_extra_artists=tuple(exart),bbox_inches='tight')
    if verbose:
        mmlio.verbose('Created plot:')
        print '     '+plotfile
    # Return the dictionary
    return plotfile

###################################################################################################################################
# METHOD TO GET INTERACTION AND HALO DATA
def halovint(simstr,snapnum=None,inpar=None,askuser=None,overwrite=None,
             exclflags=None,plot=None,owplot=None,nbins=None,**exkw):
    """
    Loads/creates halo+interaction info
    """
    # Set constants
    mrkList=['.','o','s','d','+','x','*','v','^','<','>']
    clrList=['b','g','r','m','k','c','y']
    sizRang=(1.,12.)
    clrMap=plt.get_cmap('jet')
    # Pars input
    inpar=mmlpars.mml_pars(inpar,default={},type=dict)
    askuser=mmlpars.mml_pars(askuser,default=False,type=bool)
    overwrite=mmlpars.mml_pars(overwrite,default=False,type=bool)
    if not isinstance(plot,bool) and askuser: plot=mmlio.yorn('Plot?')
    else: plot=mmlpars.mml_pars(plot,default=False,type=bool)
    if plot and askuser and not isinstance(owplot,bool): owplot=mmlio.yorn('Overwrite plot?')
    else: owplot=mmlpars.mml_pars(owplot,default=False,type=bool)
    nbins=mmlpars.mml_pars(nbins,default=100,type=int,min=1)
    exclflags=mmlpars.mml_pars(exclflags,default=[-1,-3],type=list)
    if not isinstance(snapnum,int):
        if askuser: snapnum=mmlio.askquest('Enter a valid snapshot number:',default=len(simstr.fdict['redlist'])-2,dtype='int')
        else      : snapnum=len(simstr.fdict['redlist'])-1
    if 'cuttag' not in inpar: inpar['cuttag']=None
    cutpar=mmlparam.parspar('mmlmtree.simstat','intcut',tagstr=inpar['cuttag'],askuser=askuser,load=True,save=True,init=True)
    if inpar['cuttag']==None: inpar['cuttag']=mmlparam.par2tag('mmlmtree.simstat','intcut',cutpar)
    outpar=mmlparam.parspar('mmlmtree.simstat','halovint',inpar=inpar,tagstr=inpar['cuttag'],askuser=askuser,plot=plot,load=True,save=True,init=True)
    # Redshift
    fz=simstr.fdict['redlist'][snapnum]
    # File names
    fname=snapfiles(simstr,'halovint',outpar,fext='_{:03d}'.format(snapnum))['file']
    # Load file if it exists
    if os.path.isfile(fname) and not overwrite: outvar=mmlio.rwtable('R',fname)
    else:
        # Load groups
        grps=simstr.get_groups(snapnum=snapnum,inclparticles=True)
        grpnums=range(len(grps['N']))
        haloids=simutil.haloid2groupnum(simstr,snapnum=snapnum,groupnum=grpnums,askuser=askuser)
        # Load interactions
        ints=simstr.get_intlist(cutpar=cutpar,cuttag=outpar['cuttag'],askuser=askuser)
        # Preallocate output dictionary
        intkeys=copy.deepcopy(simlist.LIST_INTPROP)+['flag']
        halokeys=['id','num']+copy.deepcopy(simlist.LIST_HALOPROP)
        outvar={'keylist':[]}
        for ikey in halokeys:
            ihalokey='halo_'+ikey
            outvar['keylist'].append(ihalokey)
            outvar[ihalokey]=[]
        for ikey in intkeys:
            iintkey='int_'+ikey
            outvar['keylist'].append(iintkey)
            outvar[iintkey]=[]
        # Loop over haloids
        for ihaloid,igrpnum in zip(haloids,grpnums):
            # Get halo properties
            if igrpnum==-1: ihaloprop={ikey:-1 for ikey in halokeys}
            else          : ihaloprop={'id'  :ihaloid,
                                       'num' :igrpnum,
                                       'cm'  :list(grps['cm'  ][igrpnum]),
                                       'cmv' :list(grps['cmv' ][igrpnum]),
                                       'mtot':float(grps['mtot'][igrpnum]),
                                       'mgas':float(grps['mgas'][igrpnum]),
                                       'spin':float(grps['spin'][igrpnum])}
            # Get interaction properties
            iintprop,iintflag=simutil.haloid2lastint(simstr,snapnum=snapnum,haloid=ihaloid,intlist=ints,askuser=askuser,outflags=True)
            iintprop['flag']=iintflag
            # Get time since last interaction
            if iintprop['iz'] in [-1]:
                ihaloprop['tiso']=iintprop['iz']
            else:
                ihaloprop['tiso']=mmlcosmo.t_lookback(zi=iintprop['iz'],zf=fz,cosmo='gadget2')
            # Add properties to dictionary
            for ikey in halokeys: outvar['halo_'+ikey].append(ihaloprop[ikey])
            for ikey in intkeys : outvar['int_' +ikey].append( iintprop[ikey])
        # Save the dictionary
        mmlio.rwtable('W',fname,outvar,overwrite=overwrite)
    # Plot
    if plot:
        # Files
        plotfile=snapfiles(simstr,'halovint',outpar,fext='_{:03d}'.format(snapnum),plot=plot)['file']
        if not mmlfiles.prep_write(plotfile,owplot): return outvar
        # Get mask
        maskval=np.array(outvar['int_flag'])
        print 'Ntot={}'.format(len(maskval))
        maskarr=np.array(len(maskval)*[False])
        for imval in exclflags:
            print 'N ({})={}'.format(imval,(maskval==imval).sum())
            maskarr=np.logical_or(maskarr,maskval==imval)
        if maskarr.sum()==len(maskarr): raise Exception('All of the data is masked.')
        print 'Nplt={}'.format(len(maskarr)-maskarr.sum())
        mmlio.yorn('?')
        # Get x & y variables
        if   outpar['x']=='otr' and outpar['y']=='otr': raise Exception('Neither x nor y variable has been set.')
        elif outpar['x']=='otr': xlab=outpar['y'] ; ylab='# of Halos' ; x=np.ma.array(outvar[xlab],mask=maskarr) ; y=None
        elif outpar['y']=='otr': xlab=outpar['x'] ; ylab='# of Halos' ; x=np.ma.array(outvar[xlab],mask=maskarr) ; y=None
        else:
            xlab=outpar['x'] ; x=np.ma.array(outvar[xlab],mask=maskarr)
            ylab=outpar['y'] ; y=np.ma.array(outvar[ylab],mask=maskarr)
        # Histogram
        if y==None: plt.hist(x[np.logical_not(maskarr)],bins=nbins)
        # Scatter plot
        else:
            # Color
            if outpar['clrvar'] in outvar: clab=outpar['clrvar'] ; c=np.ma.array(outvar[clab],mask=maskarr)
            else                         : clab='None'           ; c=np.ma.array(len(x)*['b'],mask=maskarr)
            if not isinstance(c,str): cpar=dict(cmap=clrMap,vmin=c.min(),vmax=c.max())
            else                    : cpar={}
            # Size
            if outpar['sizvar'] in outvar: slab=outpar['sizvar'] ; s=np.ma.array(outvar[slab],mask=maskarr)
            else                         : slab='None'           ; s=np.ones(x.shape)
            smin=s.min() ; smax=s.max()
            if smin==smax: s=np.ones(x.shape)*sizRang[0]
            else         : s=(float(sizRang[1]-sizRang[0])/float(smax-smin))*(s-smin)+sizRang[0]
            # Marker
            if outpar['mrkvar'] in outvar: mlab=outpar['mrkvar'] ; m=np.ma.array(outvar[mlab],mask=maskarr)
            else                         : mlab='None'           ; m=np.ones(x.shape)
            muni=set(m[np.logical_not(maskarr)])
            print muni
            if len(muni)>len(mrkList): raise Exception('Variable {} has more values ({}) than the marker list ({})'.format(mlab,len(muni),len(mrkList)))
            # Plot
            for im,idxm in zip(muni,range(len(muni))):
                imidx=np.ma.getdata(m==im)
                ip=plt.scatter(x[imidx],y[imidx],s=s[imidx],c=c[imidx],marker=mrkList[idxm],
                               label='{}={}'.format(mlab,im),alpha=0.7,**cpar)
            # Legend
            if mlab!='None': plt.legend(title=mlab)
            # Colorbar
            if clab!='None': cb=plt.colorbar(ip) ; cb.set_label(clab)
        plt.xlabel(xlab)
        plt.ylabel(ylab)
        plt.savefig(plotfile)
        mmlio.verbose('Create plot:')
        print '     '+plotfile
    # Return the dictionary
    return outvar

###################################################################################################################################
###################################################################################################################################
# PROVIDE COMMAND LINE ACCESS
if __name__ == '__main__': main()
