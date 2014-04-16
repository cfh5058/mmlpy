#!/usr/bin/python
import pysim
import os,copy,pprint,glob
import numpy as np
import matplotlib.pyplot as plt
from mmlutils import *
from pysim.files import gadget,scf,buildgal

PERF_DIR='/scratch/langmm/tests/performance'


def load(ftype,ftable=None,**exkw):
    """Load performance stats"""
    # Pars input
    if ftable is None: ftable=os.path.join(PERF_DIR,'timings_{}.txt'.format(ftype))
    exkw.setdefault('compid','stampede')
    # Load data from table
    data=mmlio.rwtable('R',ftable)
    keylist=data['keylist']
    data={k:np.array(data[k]) for k in keylist}
    # Pars data
    exclkeys=[]
    for k in exkw:
        if k not in keylist: continue
        exclkeys.append(exkw[k])
        mmlio.verbose('Selecting only data points with {} = {}'.format(k,exkw[k]))
        kidx=(data[k]==exkw[k]) ; nk=kidx.sum()
        errstr='{} entries found with {} = {}'.format(nk,k,exkw[k])
        if nk==0: raise Exception(errstr)
        else    : mmlio.verbose('    '+errstr)
        data={k:data[k][kidx] for k in keylist}
    exclstr='_'.join(exclkeys)
    N=len(data['npart'])
    # Scale SCF time by number of processors
    # if ftype=='scf':
    #     data['time']*=data['nproc']
    # Add calculated statistics
    addkeys=['Ts','Tp','To','S']
    data['Ts']=np.zeros(N)
    for n in set(data['npart']):
        nidx=(data['npart']==n)
        midx=np.argmin(data['nproc'][nidx])
        data['Ts'][nidx]=data['time'][nidx][midx]/data['nproc'][nidx][midx]
    data['Tp']=data['time']
    data['To']=data['nproc']*data['Tp']-data['Ts']
    data['S' ]=data['Ts']/data['Tp']
    # Return data
    data['keylist']=keylist+addkeys
    return data,exclstr

def predict(ftype,npart=None,nproc=None,compid=None,runtyp=None,
            simstr=None,ftable=None,**exkw):
    """Predict performance based on previous runs"""
    from sklearn import linear_model
    exp_npart=1.0 ; exp_nproc=0.25 
    def sclX(x):
        x=x.reshape((-1,2))
        x[:,0]=np.log10(x[:,0])
        # x[:,0]=(x[:,0]**exp_npart)*np.log10(x[:,0])
        # x[:,1]=np.log10(x[:,1])
        # x[:,1]=(x[:,1]**exp_nproc)*np.log10(x[:,1])
        return x
    def sclY(y,undo=False):
        if undo: y=10.**y
        else   : y=np.log10(y)
        return y
    # Pars input
    if simstr:
        runtypDEF=simstr['runtyp']
        npartDEF=simstr.ntot
        if   ftype=='buildgal':
            nprocDEF=simstr['nprocmk']
            compidDEF=simstr['mkcomp']
        elif ftype=='gadget':
            nprocDEF=simstr['nprocev']
            compidDEF=simstr['runcomp']
        elif ftype=='scf':
            nprocDEF=simstr['nprocscf']
            compidDEF=simstr['scfcomp']
        else: raise Exception('Unsupported file type: {}'.format(ftype))
    else:
        runtypDEF=None
        npartDEF=None
        nprocDEF=None
        compidDEF=None
    if runtyp is None: runtyp=runtypDEF
    if  npart is None:  npart= npartDEF
    if  nproc is None:  nproc= nprocDEF
    if compid is None: compid=compidDEF
    # Get model based on existing data
    data,exclstr=load(ftype,ftable=ftable,compid=compid,runtyp=runtyp)
    print 'nsamples  = {}'.format(len(data['npart']))
    print 'nfeatures = {}'.format(2)
    print 'ntargets  = {}'.format(1)
    XTrain=sclX(np.vstack([data['npart'],data['nproc']]).T)
    yTrain=sclY(data['time'])
    print 'X.shape = {}'.format(XTrain.shape)
    print 'y.shape = {}'.format(yTrain.shape)
    # Fit linear regression
    cfl=linear_model.LinearRegression()
    cfl.fit(XTrain,yTrain,n_jobs=1)
    def pfunc(npart,nproc,test=False):
        X=sclX(np.array([npart,nproc],dtype=float).T)
        y=cfl.predict(X)
        if not test: y=sclY(y,undo=True)
        if len(y)==1: y=y[0]
        return y
    # RMS
    y0=yTrain
    X0=np.vstack([data['npart'],data['nproc']]).T
    rms=np.sqrt(np.mean((y0-pfunc(X0[:,0],X0[:,1],test=True))**2))
    print 'rms = {}'.format(rms)
    # Return prediction
    if npart and nproc: return pfunc(npart,nproc)
    # Return function
    else              : return pfunc

def plot(ftype,ftable=None,fplot=None,fitflag=True,
         exp_npart=1.0,exp_nproc=0.25,**exkw):
    """Plot performance"""
    clrs=['b','r','g','m','c']
    syms=['','x','^','s']
    lins=['-',':',';','.']
    import matplotlib.pyplot as plt
    from matplotlib.font_manager import FontProperties
    from scipy.interpolate import splev as interp1d
    fontP = FontProperties()
    fontP.set_size('small')
    axsiz=4.3
    fig_kw={'figsize':(10,5)}
    subplot_kw=None
    # Pars input
    exkw.setdefault('compid','stampede')
    linkey=exkw.get('linkey','compid')
    symkey=exkw.get('symkey','runtyp')
    sclfact=exkw.get('sclfact',1./3600.)#'npart')
    fitfunc=exkw.get('fitfunc',None)
    # Load data from table
    data,exclstr=load(ftype,ftable=ftable,**exkw)
    keylist=data['keylist']
    N=len(data['time'])
    nlim={'npart':(min(data['npart'])   ,max(data['npart'])   ),
          'nproc':(min(data['nproc'])/1.,max(data['nproc'])*1.)}
    # Get timing units and label
    unitset=set(data['unit'])
    if len(unitset)==1: tunit=list(unitset)[0]
    else: raise Exception('Multiple time units: {}'.format(unitset))
    if   sclfact=='npart': tunit+='/particle'
    elif sclfact==1./3600.:
        if   tunit=='s/Gyr'     : tunit='hrs/Gyr'
        elif tunit=='s/snapshot': tunit='hrs/snapshot'
        else: raise Exception('Invalid tunit: {}'.format(tunit))
    tlab='Time ({})'.format(tunit)
    # Loop over keys
    prop=[]
    linset=set(data[linkey])
    symset=set(data[symkey])
    for ilin,klin in enumerate(linset):
        for isym,ksym in enumerate(symset):
            kidx=np.logical_and(data[linkey]==klin,data[symkey]==ksym)
            klab=''
            if len(linset)!=1: klab+=klin+' '
            if len(symset)!=1: klab+=ksym+' '
            dk={k:data[k][kidx] for k in keylist}
            # Time vs. npart
            ivar=dict(lab=klab,lin=lins[ilin],sym=syms[isym],sclfact=sclfact,tunit=tunit,
                      nlim=nlim)
            prop+=plotprop(0,dk,xvar='npart',yvar='time'    ,**ivar)
            prop+=plotprop(1,dk,xvar='nproc',yvar='time'    ,**ivar)
            #prop+=plotprop(2,dk,xvar='npart',yvar='overhead',**ivar)
            #prop+=plotprop(3,dk,xvar='npart',yvar='speedup' ,**ivar)
            #if fitfunc: axs[iax].plot(np.array(x),fitfunc(np.array(x)),':k')
    # Initialize plotting stuff
    amax=max([p['ax'] for p in prop])+1
    ncol=2
    nrow=int(np.ceil(float(amax)/float(ncol)))
    fig_kw.update(figsize=(ncol*axsiz*1.3,nrow*axsiz))
    fig,axs=plt.subplots(nrow,ncol,subplot_kw=subplot_kw,squeeze=False,**fig_kw)
    axs=[a for alist in axs for a in alist]
    # Plot
    for i,p in enumerate(prop):
        iax=axs[p['ax']]
        # Plot
        iax.plot(p['x'],p['y'],p['style'],color=p['color'],marker=p['marker'],label=p['label'],lw=3)
        # Fit NlogN to scaling with Npart
        if p['xvar']=='npart' and p['yvar']=='time' and fitflag:
            exp=exp_npart
            ffit=lambda N: p['y'][0]*((N**exp)*np.log10(N))/((p['x'][0]**exp)*np.log10(p['x'][0]))
            iax.plot(p['x'],ffit(p['x']),'--k')
        # Fit logN to scaling with Nproc
        if p['xvar']=='nproc' and p['yvar']=='time' and fitflag:
            exp=exp_nproc
            ffit=lambda N: p['y'][0]*((N**exp)*np.log10(N))/((p['x'][0]**exp)*np.log10(p['x'][0]))
            iax.plot(p['x'],ffit(p['x']),'--k')
        # Legend
        if i==0: legkw={'bbox_to_anchor':(0., 1.02, 1., .102),'loc':3,
                        'ncol':3,'mode':'expand','borderaxespad':0.}
        else   : legkw={'bbox_to_anchor':(0., 1.02, 1., .102),'loc':3,
                        'ncol':3,'mode':'expand','borderaxespad':0.}
        iax.legend(prop=fontP,numpoints=1,handletextpad=0.3,handlelength=1.2,**legkw)
        # Labels, limits, and scales
        if p['xlab']!=None: iax.set_xlabel(p['xlab'])
        if p['ylab']!=None: iax.set_ylabel(p['ylab'])
        if p['xlim']!=None: iax.set_xlim(p['xlim'])
        if p['ylim']!=None: iax.set_ylim(p['ylim'])
        if p['xscl']!=None: iax.set_xscale(p['xscl'])
        if p['yscl']!=None: iax.set_yscale(p['yscl'])
    # Plot contour
    # iax=axs[-1]
    # xbin=np.linspace(np.log10(nlim['npart'][0]),
    #                  np.log10(nlim['npart'][1]),11,endpoint=True)
    # ybin=np.linspace(nlim['nproc'][0],nlim['nproc'][1],11,endpoint=True)
    # Z,xbin,ybin=np.histogram2d(np.log10(data['npart']),
    #                            np.array(data['nproc']),
    #                            [xbin,ybin],weights=np.log10(data['time']))
    # Y,X = np.meshgrid((ybin[:-1]+ybin[1:])/2.,
    #                   (xbin[:-1]+xbin[1:])/2.)
    # iax.contourf(X,Y,Z.T)
    # iax.set_xlabel('logNpart')
    # iax.set_ylabel('Nproc')
    # Save
    if fplot is None: fplot=os.path.join(PERF_DIR,'timings_{}_{}'.format(ftype,exclstr))
    plt.tight_layout(rect=[0,0,1,0.90],h_pad=2.0)
    plt.savefig(fplot)
    print '    '+fplot
    # Return
    return

# Plot scaling with problem size or number of processors
def plotprop(ax,d,xvar='npart',yvar='time',sclfact=None,tunit='s',
             nlim={},**exkw):
    """Plot time vs. npart for processor number"""
    from matplotlib import colors,cm
    syms=['.','x','^','s']
    lins=['-',':',';','.']
    if   xvar=='npart': nvar='nproc'
    elif xvar=='nproc': nvar='npart'
    else: raise Exception('Invalid x variable: {}'.format(xvar))
    nlist=set(d[nvar])
    if   nvar=='npart': cmap=cm.get_cmap('cool')
    elif nvar=='nproc': cmap=cm.get_cmap('cool')
    lab=exkw.get('lab','')
    lin=exkw.get('lin',lins[0])
    sym=exkw.get('sym',syms[0])
    p=[]
    for i,n in enumerate(sorted(nlist)):
        # Get x data and label
        idx=(d[nvar]==n)
        x=d[xvar][idx]
        xlab=xvar.title()
        if   nvar=='npart': ilab=lab+'log{} = {:3.2f}'.format(nvar.title(),np.log10(n))
        elif nvar=='nproc': ilab=lab+'{} = {}'.format(nvar.title(),n)
        else: raise Exception('Unsupported nvar: {}'.format(nvar))
        if   xvar=='npart': xscl='log'
        elif xvar=='nproc': xscl='linear'
        else: raise Exception('Unsupported xvar: {}'.format(xvar))
        # Get y data and label
        if   yvar=='time':
            y=d['time' ][idx]
            if sclfact=='npart': 
                y/=d['npart'][idx]
            elif isinstance(sclfact,(float,int,np.ndarray)):
                y*=sclfact
            ylab='Time ({})'.format(tunit)
            yscl='log'
        elif yvar=='overhead':
            y=d['To'   ][idx]
            ylab='Overhead ({})'.format(tunit)
            yscl='log'
        elif yvar=='speedup':
            y=d['S'    ][idx]
            ylab='Speedup'
            yscl=None
        else: raise Exception('Invalid variable: {}'.format(yvar))
        # Sort by increasing x
        idxsort=np.argsort(x)
        x=x[idxsort]
        y=y[idxsort]
        # Add properties to list
        nlim.setdefault(xvar,(min(x),max(x)))
        nlim.setdefault(nvar,(min(nlist),max(nlist)))
        if   nvar=='npart': norm=colors.LogNorm(*nlim[nvar])
        elif nvar=='nproc': norm=colors.Normalize(*nlim[nvar])
        else: raise Exception('Unsupported nvar: {}'.format(nvar))
        clr=cmap(norm(n))
        p.append({'xvar':xvar,'x':x,'xlab':xlab,'xlim':nlim[xvar],'xscl':xscl,
                  'yvar':yvar,'y':y,'ylab':ylab,'ylim':None      ,'yscl':yscl,
                  'ax':ax,'style':lin,'marker':sym,'color':clr,'label':ilab})
    return p

def files(baserun,npart):
    """Returns dictionary of files for test runs"""
    # Get associated base simulation
    sim=pysim.loadsim(baserun)
    # General options
    fopt={}
    fopt['runtyp']=sim['runtyp']
    fopt['runtag']='{}_{}'.format(baserun,mmlstring.dec2str(np.log10(npart)))
    fopt['topdir']=PERF_DIR
    fopt['simdir']=os.path.join(fopt['topdir'],fopt['runtyp'])
    fopt['rundir']=os.path.join(fopt['simdir'],fopt['runtag'])
    # Gadget
    fopt['gadget']={'dir':os.path.join(fopt['rundir'],'gadget')}
    fopt['gadget']['ic']=os.path.join(fopt['gadget']['dir'],'{}.ic'.format(fopt['runtag']))
    fopt['gadget']['pm']=os.path.join(fopt['gadget']['dir'],'{}.pm'.format(fopt['runtag']))
    fopt['gadget']['mk']=os.path.join(fopt['gadget']['dir'],'{}.mk'.format(fopt['runtag']))
    # SCF
    fopt['scf']={'dir':os.path.join(fopt['rundir'],'scf')}
    fopt['scf']['ic']=os.path.join(fopt['scf']['dir'],'{}.000'.format(fopt['runtag']))
    # Return options
    return fopt

def make(ftype,baserun,npart,overwrite=False):
    """Creates/copies over the necessary files for a downsampled run"""
    sim=pysim.loadsim(baserun)
    # Downsample for gadget

        

def mksim(ftype,baserun,npart,nproc,compid):
    """Create test sim object"""
    # Load base run
    sim0=pysim.loadsim(runtag=baserun)
    # Create label for new run
    runtag='timetest'
    subtag='{}on{}'.format(mmlstring.dec2str(np.log10(npart)),nproc)
    #subtag='{}on{}_{}'.format(mmlstring.dec2str(np.log10(npart)),nproc,baserun)
    # Create infodict
    ifo=copy.deepcopy(sim0.infodict) ; ifotag=ifo.pop('tagstr',None)
    ifo=mmlparam.parspar('pysim.simlist',sim0['runtyp'],inpar=ifo,tagstr=runtag,
                         askuser=False,init=True,save=True,overwrite=True,verbose=False)
    # Create sim
    sim={'runtag':runtag,'subtag':subtag,'runtyp':sim0['runtyp']}
    for k in ['nprocmk','nprocev','nprocscf']: sim[k]=nproc
    for k in ['mkcomp' ,'runcomp','scfcomp' ]: sim[k]=compid
    sim=mmlparam.parspar('pysim.simlist','mmlsim',inpar=sim,tagstr=runtag,
                         askuser=False,init=True,save=True,overwrite=True,verbose=False)
    sim['infodict']=ifo
    sim=pysim.loadsim(**sim)
    return sim

# RUN AND RETRIEVE TESTS
def test(method,ftype,npart,nproc,compid,baserun,parmrun=None,**exkw):
    """Do timing tests"""
    # Pars input
    if parmrun is None: parmrun=baserun
    sim=mksim(ftype,baserun,npart,nproc,compid)
    # Ask
    if not mmlio.yorn('Do {} {} for {}?'.format(method,ftype,sim['runtag']+sim['subtag'])): return None
    # Submit run
    if   method=='sub':
        mkw={}
        if   ftype=='gadget'  : mkw={'method':'timetest','baserun':baserun,'parmrun':parmrun,'snapext':'ic','ntot':npart}
        elif ftype=='buildgal': pass
        elif ftype=='scf'     : 
            snstat=sim.fdict['gadget']['output']['snapshot'].format(0)
            if True:#not os.path.isfile(snstat+'*'):
                icstat=sim.fdict['gadget']['input']['ic']
                # if not os.path.isfile(icstat+'*'):
                #     fmv=[sim.fdict['gadget']['input']['ic'],
                #          sim.fdict['gadget']['input']['makefile'],
                #          sim.fdict['gadget']['input']['param']]
                #     for fdst in fmv:
                #         fsrc=fdst.replace('on{}'.format(nproc),'on{}'.format(256))
                #         mmlfiles.mkdirs(os.path.dirname(fdst))
                #         for ifsrc in glob.glob(fsrc+'*'):
                #             print fdst+ifsrc.split(fsrc)[-1]
                #             mmlfiles.cp(ifsrc,fdst+ifsrc.split(fsrc)[-1])
                for fic in glob.glob(icstat+'*'):
                    mmlfiles.cp(fic,snstat+fic.split(icstat)[-1],overwrite=True)
        pysim.files.subrun(sim,ftype,overwrite=True,owmake=True,mkw=mkw)
    # End run
    elif method=='end':
        pysim.files.endrun(sim,ftype)
    elif method=='record': pass
    # Error
    else: raise Exception('Invalid method: {}'.format(method))
    # Get performance specs
    if method in ['end','record']:
        save(ftype,sim,overwrite=True,**exkw)
    # Return sim
    return sim

def calc(ftype,sim,**exkw):
    """Calculate performance statistics"""
    # Use module specfic function if it exists
    fmod=pysim.files.ftype2module(ftype)
    if hasattr(fmod,'performance'): return fmod.performance(sim,**exkw)
    # Handle generic performance measure from output file
    fout=sim.fdict[ftype]['output']['runout']
    print '    '+fout
    # Get default units
    if   ftype=='scf'     : unitDEF='s/snapshot'
    elif ftype=='buildgal': unitDEF='s/ic'
    else: raise Exception('Not default unit for file type: {}'.format(ftype))
    unitDEF=exkw.get('unit',unitDEF)
    # Get user input
    time=mmlio.askquest('Enter {} execution time from above file.'.format(ftype),default=0.,dtype='float')
    unit=mmlio.askquest('Enter {} units for time just entered.'.format(ftype),default=unitDEF,dtype='string')
    # Scale by number of processors
    if ftype=='scf': unit*=float(sim['nprocscf'])
    else: raise Exception('Add option {} here.'.format(ftype))
    # Return time and units
    return time,unit

def save(ftype,sim,ftable=None,overwrite=False,**exkw):
    """Save performance data to file"""
    from . import ftype2nproc,ftype2compid
    # Pars input
    if ftable is None: ftable=os.path.join(PERF_DIR,'timings_{}.txt'.format(ftype))
    # Collect data
    time,unit=calc(ftype,sim,**exkw)
    npart=sim.ntot
    nproc=ftype2nproc(ftype,sim)
    compid=ftype2compid(ftype,sim)
    baserun=sim.baserun
    # Append table
    keylist=['tag','baserun','runtyp','compid','npart','nproc','time','unit']
    entry={'tag':'{}_{}on{}_{}'.format(compid,mmlstring.dec2str(np.log10(npart)),nproc,baserun),
           'baserun':baserun,'runtyp':sim['runtyp'],'compid':compid,'npart':npart,'nproc':nproc,'time':time,'unit':unit}
    mmlio.rwtable('A',ftable,entry,keylist=keylist,uniquekey='tag',overwrite=overwrite)
    # Return control
    return time,unit
