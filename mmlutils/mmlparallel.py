#!/usr/bin/python 
####################################################################################################################################
#
# MEAGAN LANG'S PARALLEL METHODS
#
####################################################################################################################################
import numpy as np
import copy
import multiprocessing
from multiprocessing import Process, Pipe
from itertools import izip

####################################################################################################################################
# METHODS FOR MULTIPROCESSING
# def spawn(f):
#     def fun(q_in,q_out):
#         while True:
#             i,x = q_in.get()
#             if i == None:
#                 break
#             q_out.put((i,f(x)))
#     return fun

# def parmap(f, X, nproc = multiprocessing.cpu_count()/3):
#     q_in   = multiprocessing.Queue(1)
#     q_out  = multiprocessing.Queue()

#     proc = [multiprocessing.Process(target=spawn(f),args=(q_in,q_out)) for _ in range(nproc)]
#     for p in proc:
#         p.daemon = True
#         p.start()

#     sent = [q_in.put((i,x)) for i,x in enumerate(X)]
#     [q_in.put((None,None)) for _ in range(nproc)]
#     res = [q_out.get() for _ in range(len(sent))]

#     [p.join() for p in proc]

#     return [x for i,x in sorted(res)]

def spawn(f):
    def fun(ppipe,cpipe,x):
        ppipe.close()
        cpipe.send(f(x))
        cpipe.close()
    return fun

def parmap(f,X,nproc=multiprocessing.cpu_count()/3):
    out=[]
    if nproc==1:
        for x in X: out.append(f(x))
    else:
        nX=len(X)
        nB=int(np.ceil(float(nX)/float(nproc)))
        print 'nX={}'.format(nX)
        print 'nB={}'.format(nB)
        for iB in range(nB):
            beg=nproc*iB
            end=min(nX,nproc*(iB+1))
            print '(iB,beg,end)=({},{},{})'.format(iB,beg,end)
            iX=X[beg:end]
            pipe=[Pipe() for x in iX]
            proc=[Process(target=spawn(f),args=(p,c,x)) for x,(p,c) in izip(iX,pipe)]
            [p.start() for p in proc]
            out+=[p.recv() for (p,c) in pipe]
            [p.join() for p in proc]
    return out

####################################################################################################################################
# METHOD FOR PLOTING PERFORMANCE STATISTICS
def plot_performance(nproc,npart,times,fname='performance.png',timelabel='Time/particle (s)',extra=None,fitfunc=None):
    """
    Plot performance statistics
    """
    clrs=['r','b','g','m','c']
    # Initialize plotting
    import matplotlib.pyplot as plt
    from matplotlib.font_manager import FontProperties
    from scipy.interpolate import splev as interp1d
    fontP = FontProperties()
    fontP.set_size('small')
    fig_kw={'figsize':(10,5)}
    subplot_kw=None
    fig,axs=plt.subplots(1,2,subplot_kw=subplot_kw,**fig_kw)
    #axs=[item for sublist in axs for item in sublist]
    xlab=[None]*len(axs) ; ylab=[None]*len(axs)
    xlim=[None]*len(axs) ; ylim=[None]*len(axs)
    xscl=[None]*len(axs) ; yscl=[None]*len(axs)
    npartlim=(min(npart)/5,max(npart)*5)
    sclarr=np.array(npart)
    sclarr=1.
    # Get missing values
    fun=[None]*len(nproc)
    times_fill=copy.deepcopy(times)
    for i,p in enumerate(nproc):
        x=np.array(npart)
        y=np.array(times[i])/x#sclarr
        # NlogN scaling
        fun[i]=lambda inx: y[y>0].max()*(inx*np.log10(inx))/(x[y>0].max()*np.log10(x[y>0].max()))
        # Linear interp
        # fun[i]=lambda inx: y[y>0].min()+inx*((y[y>0].max()-y[y>0].min())/(x[y>0].max()-x[y>0].min()))
        # f = interp1d(x[y>0],y[y>0])
        # Times
        if np.any(y==0): 
            y[y==0]=fun[i](x[y==0])
            times_fill[i]=list(y*x)#sclarr)
        
    # Plot scaling with problem size
    iax=0
    for i,p in enumerate(nproc):
        x=np.array(npart)
        y=np.array(times[i])/sclarr
        print p,(np.log10(y[2])-np.log10(y[1]))/(np.log10(x[2])-np.log10(x[1]))
        axs[iax].plot(x,y,'-x'+clrs[i],label='Nproc={}'.format(p))
        if np.any(y==0): 
            y_fill=np.array(times_fill[i])/sclarr
            axs[iax].plot(x,y_fill,'--'+clrs[i])#,label='Nproc={}'.format(p))
    if fitfunc: axs[iax].plot(np.array(npart),fitfunc(np.array(npart)),':k')
    if extra!=None:
        x=extra['npart']
        if isinstance(sclarr,float): extscl=1.
        else                       : extscl=extra['npart']
        y=extra['time']/extscl
        axs[iax].plot(x,y,'ok')#,label='Nproc={}'.format(extra['nproc']))
    xlab[iax]='Npart'
    ylab[iax]=timelabel
    xlim[iax]=npartlim
    xscl[iax]='log'
    yscl[iax]='log'
    # Get variables
    Ts=np.array(times[nproc.index(1)])
    Ts_fill=np.array(times_fill[nproc.index(1)])
    # # Plot overhead
    # iax=1
    # for i,p in enumerate(nproc):
    #     if p==1: continue
    #     x=np.array(npart)
    #     Tp=np.array(times[i])
    #     To=p*Tp-Ts
    #     axs[iax].plot(x,To,'x-'+clrs[i],label='Nproc={}'.format(p))
    # xlab[iax]='Npart'
    # ylab[iax]='Overhead (s)'
    # xlim[iax]=npartlim
    # xscl[iax]='log'
    # yscl[iax]='log'
    # Plot speedup
    iax=1
    for i,p in enumerate(nproc):
        if p==1: continue
        x=np.array(npart)
        Tp=np.array(times[i])
        S=Ts/Tp
        axs[iax].plot(x[S!=0],S[S!=0],'x-'+clrs[i],label='Nproc={}'.format(p))
        if np.any(S==0): 
            S_fill=np.array(Ts_fill)/Tp
            axs[iax].plot(x,S_fill,'--'+clrs[i])#,label='Nproc={}'.format(p))
    xlab[iax]='Npart'
    ylab[iax]='Speedup'
    xlim[iax]=npartlim
    xscl[iax]='log'
    # # Plot efficiency
    # iax=3
    # for i,p in enumerate(nproc):
    #     if p==1: continue
    #     x=np.array(npart)
    #     Tp=np.array(times[i])
    #     S=Ts/Tp
    #     E=S/p
    #     axs[iax].plot(x,E,'x-'+clrs[i],label='Nproc={}'.format(p))
    # xlab[iax]='Npart'
    # ylab[iax]='Efficiency'
    # xlim[iax]=npartlim
    # xscl[iax]='log'
    # Legends,labels,limits
    for i,iax in enumerate(axs):
        if i==0: legkw={'bbox_to_anchor':(0., 1.02, 1., .102),'loc':3,
                        'ncol':3,'mode':'expand','borderaxespad':0.}
        else   : legkw={'bbox_to_anchor':(0., 1.02, 1., .102),'loc':3,
                        'ncol':2,'mode':'expand','borderaxespad':0.}
        iax.legend(prop=fontP,numpoints=1,handletextpad=0.3,handlelength=1.3,**legkw)
        if xlab[i]!=None: iax.set_xlabel(xlab[i])
        if ylab[i]!=None: iax.set_ylabel(ylab[i])
        if xlim[i]!=None: iax.set_xlim(xlim[i])
        if ylim[i]!=None: iax.set_ylim(ylim[i])
        if xscl[i]!=None: iax.set_xscale(xscl[i])
        if yscl[i]!=None: iax.set_yscale(yscl[i])
        # x0,x1 = iax.get_xlim()
        # y0,y1 = iax.get_ylim()
        # iax.set_aspect((x1-x0)/(y1-y0))
    # Save
    plt.tight_layout(rect=[0,0,1,0.97],h_pad=2.0)
    plt.savefig(fname)
    # Return
    return
