#!/usr/bin/python
import pysim
import os,copy,pprint
import numpy as np
import matplotlib.pyplot as plt
from mmlutils import *
from pysim.files import gadget,scf,buildgal


def load_performance(ftype,ftable=None,**exkw):
    """Load performance stats"""
    # Pars input
    if ftable is None: ftable=os.path.expanduser('~/timings_{}.txt').format(ftype)
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

def pred_performance(ftype,npart=None,nproc=None,compid=None,runtyp=None,
                     simstr=None,ftable=None,**exkw):
    """Predict performance based on previous runs"""
    from sklearn import linear_model
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
    data,exclstr=load_performance(ftype,ftable=ftable,compid=compid,runtyp=runtyp)
    print 'nsamples  = {}'.format(len(data['npart']))
    print 'nfeatures = {}'.format(2)
    print 'ntargets  = {}'.format(1)
    XTrain=np.vstack([data['npart'],data['nproc']]).T
    XTrain[0]*=np.log10(XTrain[0])
    XTrain[1] =np.log10(XTrain[1])
    yTrain=data['time']
    #yTrain=np.log10(data['time'])
    print 'X.shape = {}'.format(XTrain.shape)
    print 'y.shape = {}'.format(yTrain.shape)
    # Log fit
    exp_npart=1.0 ; exp_nproc=0.25
    
    # Fit linear regression
    cfl=linear_model.LinearRegression()
    cfl.fit(XTrain,yTrain,n_jobs=1)
    def pfunc(npart,nproc):
        X=np.array([npart,nproc],dtype=float).T
        X[0]*=np.log10(X[0])
        X[1] =np.log10(X[1])
        y1=cfl.predict(X)
        
        return y1
    # Return prediction
    if npart and nproc: return pfunc(npart,nproc)
    # Return function
    else              : return pfunc

def plot_performance(ftype,ftable=None,fplot=None,**exkw):
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
    data,exclstr=load_performance(ftype,ftable=ftable,**exkw)
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
        if tunit=='s/Gyr': tunit='hrs/Gyr'
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
            prop+=plotprop(2,dk,xvar='npart',yvar='overhead',**ivar)
            prop+=plotprop(3,dk,xvar='npart',yvar='speedup' ,**ivar)
            #if fitfunc: axs[iax].plot(np.array(x),fitfunc(np.array(x)),':k')
    # Initialize plotting stuff
    amax=max([p['ax'] for p in prop])+1
    ncol=2
    nrow=int(np.ceil(float(amax)/float(ncol)))
    fig_kw.update(figsize=(ncol*axsiz*1.3,nrow*axsiz))
    fig,axs=plt.subplots(nrow,ncol,subplot_kw=subplot_kw,**fig_kw)
    axs=[a for alist in axs for a in alist]
    # Plot
    for i,p in enumerate(prop):
        iax=axs[p['ax']]
        # Plot
        iax.plot(p['x'],p['y'],p['style'],color=p['color'],marker=p['marker'],label=p['label'],lw=3)
        # Fit NlogN to scaling with Npart
        if p['xvar']=='npart' and p['yvar']=='time':
            ffit=lambda N: p['y'][0]*(N*np.log10(N))/(p['x'][0]*np.log10(p['x'][0]))
            iax.plot(p['x'],ffit(p['x']),'--k')
        # Fit logN to scaling with Nproc
        if p['xvar']=='nproc' and p['yvar']=='time':
            exp=0.25
            ffit=lambda N: p['y'][0]*((N**exp)*np.log10(N))/((p['x'][0]**exp)*np.log10(p['x'][0]))
            iax.plot(p['x'],ffit(p['x']),'--k')
        # Legend
        if i==0: legkw={'bbox_to_anchor':(0., 1.02, 1., .102),'loc':3,
                        'ncol':3,'mode':'expand','borderaxespad':0.}
        else   : legkw={'bbox_to_anchor':(0., 1.02, 1., .102),'loc':3,
                        'ncol':3,'mode':'expand','borderaxespad':0.}
        iax.legend(prop=fontP,numpoints=1,handletextpad=0.3,handlelength=1.3,**legkw)
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
    if fplot is None: fplot=os.path.expanduser('~/timings_{}_{}').format(ftype,exclstr)
    plt.tight_layout(rect=[0,0,1,0.97],h_pad=2.0)
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

def calc_performance(method,ftype,npart,nproc,compid,baserun,parmrun=None,ftable=None):
    """Do timing tests"""
    if ftable is None: ftable=os.path.expanduser('~/timings_{}.txt').format(ftype)
    # Pars input
    if parmrun is None: parmrun=baserun
    sim0=pysim.loadsim(runtag=baserun)
    # Label run
    runtag='timetest'
    subtag='{}on{}'.format(mmlstring.dec2str(np.log10(npart)),nproc)
    # Ask
    if not mmlio.yorn('Do {} {} for {}?'.format(method,ftype,runtag+subtag)): return sim0
    # Create infodict
    ifo=copy.deepcopy(sim0.infodict) ; ifotag=ifo.pop('tagstr',None)
    ifo=mmlparam.parspar('pysim.simlist',sim0['runtyp'],inpar=ifo,tagstr=runtag,
                         askuser=False,init=True,save=True,overwrite=True,verbose=False)
    # Create sim
    sim={'runtag':runtag,'subtag':subtag,
         'runtyp':sim0['runtyp'],
         'mkcomp' :compid,'runcomp':compid,'scfcomp' :compid,
         'nprocmk':nproc ,'nprocev':nproc ,'nprocscf':nproc }
    sim=mmlparam.parspar('pysim.simlist','mmlsim',inpar=sim,tagstr=runtag,
                         askuser=False,init=True,save=True,overwrite=True,verbose=False)
    sim['infodict']=ifo
    sim=pysim.loadsim(**sim)
    # Submit run
    if   method=='sub':
        mkw={}
        if   ftype=='gadget'  : mkw={'method':'timetest','baserun':baserun,'parmrun':parmrun,'snapext':'ic','ntot':npart}
        elif ftype=='buildgal': pass
        elif ftype=='scf'     : mmlfiles.cp(sim.fdict['gadget']['input']['ic'],sim.fdict['gadget']['output']['snapshot'].format(0))
        pysim.files.subrun(sim,ftype,overwrite=True,owmake=True,mkw=mkw)
    # End run
    elif method=='end':
        # Collect files
        pysim.files.endrun(sim,ftype)
    # Get performance specs
    elif method in ['end','record']:
        # Collect data
        time,unit=pysim.files.performance('record',sim,ftype,overwrite=True)
        # Append table
        keylist=['tag','baserun','runtyp','compid','npart','nproc','time','unit']
        entry={'tag':'{}_{}on{}_{}'.format(compid,mmlstring.dec2str(np.log10(npart)),nproc,baserun),
               'baserun':baserun,'runtyp':sim['runtyp'],'compid':compid,'npart':npart,'nproc':nproc,'time':time,'unit':unit}
        mmlio.rwtable('A',ftable,entry,keylist=keylist,uniquekey='tag',overwrite=True)
    # Error
    else: raise Exception('Invalid method: {}'.format(method))
    # Return sim
    return sim
