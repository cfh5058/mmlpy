#!/usr/bin/python
import pysim
import os,copy,pprint
import numpy as np
import matplotlib.pyplot as plt
from mmlutils import *
from pysim.files import performance,gadget,scf,buildgal
LIST_GDPTYPS=gadget.LIST_PTYPS

#method='end'
#method='sub'
#method='gadget4scf'
#method='plot'
method='calc'
#method='performance'
#ftypeList=['gadget']
ftypeList=['scf']

compid='stampede'
#compid='accre'

testrun_ic={'gadget':'MW0B8_MW0B8_Q1RP0p1E1I0_TEST',
            'scf'   :'MW0B_8M12LH11FG0'}
testrun_pm={'gadget':'MW0B8_MW0B8_Q1RP0p1E1I0_TEST',
            'scf'   :'MW0B_8M12LH11FG0'}
testrun_tp={'gadget':'intsim',
            'scf'   :'galsim'}

dt_snap=0.01 # Gyr
tsim_iso=3.0 # Gyr
tsim_int=5.0 # Gyr
Npart_MW1=1.0e8

# Isolated stuff
Npart_iso={'name' :['MW1','MW3','MW10'],
           'fact' :[1.0  ,3.0  ,10.0  ],
           'nsim' :[1    ,1    ,1     ],
           'tsim' :3*[tsim_iso]}
Npart_iso['npart']=[Npart_MW1/f for f in Npart_iso['fact']]
Npart_iso['nsnap']=[t/dt_snap for t in Npart_iso['tsim']]
# Interactions
Npart_int={'name' :['q1' ,'q3' ,'q10' ],
           'fact' :[1.0  ,3.0  ,10.0  ],
           'nsim' :[8    ,5    ,5     ],
           'tsim' :3*[tsim_int]}
Npart_int['npart_gal']=[[Npart_MW1,Npart_MW1/f] for f in Npart_int['fact']]
Npart_int['npart']=[Npart_MW1*(1.+1./f) for f in Npart_int['fact']]
Npart_int['nsnap']=[t/dt_snap for t in Npart_int['tsim']]

Npart_test={'gadget':[long(2e6),long(2e7),long(2e8),
                      long(1e6),long(1e7),long(1e8)],
            'scf'   :[long(1e6),long(1e7),long(1e8)]}

Nproc_test={'gadget':[64,128,256],
            'scf'   :[10, 50,100]}
Nproc_scf=50
Nproc_gd=256

if method=='gadget4scf':
    method='sub'
    compid='accre'
    ftypeList=['gadget']
    testrun_ic['gadget']=testrun_ic['scf']
    testrun_pm['gadget']=testrun_pm['scf']
    Nproc_test['gadget']=Nproc_test['scf']
    Npart_test['gadget']=Npart_test['scf']


# Gadget Done on ACCRE
# 6p3on64,6p3on128,6p3on256
# 7p3on64,7p3on128
# 8p3on64,8p3on128

# Error on 8p3on256 due to insufficient memory

# times={'4on1' :{'scf':5. ,'gadget':18.}, # done
#        '4on10':{'scf':6. ,'gadget':10.}, # done
#        '4on50':{'scf':20.,'gadget':31.}, # done
#        '5on1' :{'scf':6. ,'gadget':1015.}, # done
#        '5on10':{'scf':7. ,'gadget':300.}, # done
#        '5on50':{'scf':13.,'gadget':217.20/0.01}, #139.}, # done
#        '6on1' :{'scf':28.,'gadget':3128.+3136.+73473.+18541.},
#        '6on10':{'scf':22.,'gadget':3082.+3115.+38425.}, # done
#        '6on50':{'scf':27.,'gadget':3097.+15798.}, # done
#        '7on1' :{'scf':67.,'gadget':0.}, # won't run
#        '7on10':{'scf':53.,'gadget':3467.+76313.+76498.},
#        '7on50':{'scf':55.,'gadget':3277.+74024.+73839.},
#        # New runs
#        '8on64' :{'scf':0. ,'gadget':0. },
#        '8on128':{'scf':0. ,'gadget':0. },
#        '8on256':{'scf':0. ,'gadget':0. }}

# Set up plot
if method=='plot': 
    for ftype in ftypeList:
        fname='/home/langmm/code/python/mine/pysim/performance/perf_{}.png'.format(ftype)
        exp_npart=1.
        exp_nproc=0.25
        if ftype=='scf':
            exp_npart=0.25
            exp_nproc=0.75
        performance.plot(ftype,fplot=fname,compid=compid,baserun=testrun_ic[ftype],runtyp=testrun_tp[ftype],
                         exp_npart=exp_npart,exp_nproc=exp_nproc)
# Calculate costs
elif method=='calc':
    sbuf='    '
    s2hr=1.0/3600.0
    # BUILDGAL
    if 'buildgal' in ftypeList:
        print 'BUILDGAL'
        tbg=0.
    # SCF
    if 'scf' in ftypeList:
        print 'SCF'
        fscf_iso=performance.predict('scf',compid=compid,runtyp='galsim')
        fscf_int=fscf_iso
        #fscf_int=performance.predict('scf',compid=compid,runtyp='galsim')
        print 1*sbuf+'Isolation'
        tscf_iso=0.
        for i in range(len(Npart_iso['name'])): 
            Nsim=Npart_iso['nsim'][i]
            tscf=np.ceil(s2hr*Npart_iso['nsnap'][i]*fscf_iso(Npart_iso['npart'][i],Nproc_scf))
            tscf_iso+=tscf*Nsim
            print '{}{:4s} {} ({} x {})'.format(2*sbuf,Npart_iso['name'][i],tscf*Nsim,tscf,Nsim)
        print 1*sbuf+'Interactions'
        tscf_int=0.
        for i in range(len(Npart_int['name'])):
            Nsim=Npart_int['nsim'][i]
            tscf=0.
            for gNpart in Npart_int['npart_gal'][i]:
                tscf+=np.ceil(s2hr*Npart_int['nsnap'][i]*fscf_iso(gNpart,Nproc_scf))
            tscf_int+=tscf*Nsim
            print '{}{:4s} {} ({} x {})'.format(2*sbuf,Npart_int['name'][i],tscf*Nsim,tscf,Nsim)
        print 'Total={} SUs'.format(tscf_iso+tscf_int)
    # GADGET
    if 'gadget' in ftypeList:
        print 'GADGET'
        fgd_iso=performance.predict('gadget',compid=compid,runtyp='intsim')
        fgd_int=fgd_iso
        #fgd_int=performance.predict('gadget',compid=compid,runtyp='galsim')
        print 1*sbuf+'Isolation'
        tgd_iso=0.
        for i in range(len(Npart_iso['name'])):
            Nsim=Npart_iso['nsim'][i]
            tgd=np.ceil(s2hr*Npart_iso['tsim'][i]*fgd_iso(Npart_iso['npart'][i],Nproc_gd))
            tgd_iso+=tgd*Nsim
            print '{}{:4s} {} ({} x {})'.format(2*sbuf,Npart_iso['name'][i],tgd*Nsim,tgd,Nsim)
        print 1*sbuf+'Interactions'
        tgd_int=0.
        for i in range(len(Npart_int['name'])):
            Nsim=Npart_int['nsim'][i]
            tgd=np.ceil(s2hr*Npart_int['tsim'][i]*fgd_int(Npart_int['npart'][i],Nproc_gd))
            tgd_int+=tgd*Nsim
            print '{}{:4s} {} ({} x {})'.format(2*sbuf,Npart_int['name'][i],tgd*Nsim,tgd,Nsim)
        print 'Total={} SUs ({} w/o isolation)'.format(tgd_iso+tgd_int,tgd_int)
    
# Submitting and retrieving files
else: 
    # Loop over file types
    for ftype in ftypeList:
        # Loop over particle number
        for idxpart,iNpart in enumerate(Npart_test[ftype]):
            # Loop over number of processors
            for idxproc,iNproc in enumerate(Nproc_test[ftype]):
                # Get test sim
                sim=performance.test(method,ftype,iNpart,iNproc,compid,testrun_ic[ftype],testrun_pm[ftype])


    # elif method=='oldplot':
    #     if   imeth=='gadget': 
    #         tlab='Wall-clock Time/Gyr (s)'
    #         fitf=lambda N: data[-1][0]*(N*np.log10(N))/(Npart_test[0]*np.log10(Npart_test[0]))
    #     elif imeth=='scf'   : 
    #         tlab='Wall-clock Time/snap (s)'
    #         fitf=lambda N: data[-1][0]*(np.log10(N))/(np.log10(Npart_test[0]))
    #     else: tlab=None
    #     if   imeth=='gadget': 
    #         simext=pysim.loadsim('MW0B7_MW0B7_Q1RP0p1E1I0')
    #         tscl=(3.085678e21/1e5)/(3.15569e16)
    #         cpustats=gadget.rw_cpustat('R',simext.fdict['gadget']['output']['cpu'])
    #         itsimarr=np.array(cpustats['tsim'])
    #         isteparr=np.array(cpustats['step'])
    #         idxend=np.arange(len(isteparr))[isteparr<=128].max()
    #         #idxend=np.arange(len(itsimarr))[itsimarr<=0.0102269].max()-1
    #         exttime=cpustats['tot'][idxend]/((cpustats['tsim'][idxend]/tscl))
    #         extproc=simext['nprocev']
    #         extpart=simext.ntot
    #         ext={'time':exttime,'nproc':extproc,'npart':extpart}
    #         print exttime
    #         print extpart
    #         print extproc
    #     else:
    #         ext=None
            
    #     # print imeth
    #     # pprint.pprint(data)
    #     mmlparallel.plot_performance(Nproc_test,Npart_test,data,fname=fname,
    #                                  timelabel=tlab,extra=ext,fitfunc=fitf)
    #     print '    '+fname
        # for iNpart in Npart_test:
        #     ilogN=int(np.log10(iNpart))
        #     x=np.array(Nproc_test)
        #     y=np.array(plotlist[ilogN],dtype=float)
        #     iaxs.plot(x,y,label='log(Npart)={}'.format(ilogN))
        # # Labels
        # plt.xlabel('Nproc')
        # ylab='Time'
        # if perpart: ylab+='/particle'
        # if imeth=='scf': ylab+='/snap'
        # else:
        #     if perstep: ylab+='/step'
        #     if pertime: ylab+='/Gyr'
        # plt.ylabel(ylab+' (s)')
        # iaxs.set_xscale('log')
        # iaxs.set_yscale('log')
        # # Legend
        # iaxs.legend(loc=1)
        # fname='/home/langmm/code/python/mine/runscripts/timetest_'+imeth+'.png'
        # plt.savefig(fname)
        # plt.cla()
        # print '    '+fname


# Print info
# elif method=='oldcalc':
#     s2hr=1.0/3600.0
#     dt_snap=0.01 # Gyr
#     time_iso=3.0 # Gyr
#     time_int=5.0 # Gyr
#     Nsnap_iso=time_iso/dt_snap
#     Nsnap_int=time_int/dt_snap
#     Nsim_q1 =8
#     Nsim_q3 =5
#     Nsim_q10=5
#     # BUILDGAL
#     print 'BUILDGAL'
#     npart0_scf=1.0e6
#     t0_bgal_accre=(5.*3600.+10.*60.+8.)
#     t0_bgal_stamp=t0_bgal_accre*(2.3/2.7)
#     tscf_scaleN=lambda N: t0_bgal_stamp*((N/npart0_scf)**2)
#     ts_bgal_MW1 =s2hr*tscf_scaleN(Npart_MW1 )
#     ts_bgal_MW3 =s2hr*tscf_scaleN(Npart_MW3 )
#     ts_bgal_MW10=s2hr*tscf_scaleN(Npart_MW10)
#     print '    t (MW1 ) = {} hrs'.format(ts_bgal_MW1 )
#     print '    t (MW3 ) = {} hrs'.format(ts_bgal_MW3 )
#     print '    t (MW10) = {} hrs'.format(ts_bgal_MW10)
#     print 'Total={} hrs'.format(ts_bgal_MW1+ts_bgal_MW3+ts_bgal_MW10)
#     # SCF
#     print 'SCF'
#     Nproc_scf=10
#     npart0_scf=1.0e7
#     t0_scf=times['7on10']['scf']
#     tscf_scaleN=lambda N: 10.**(np.log10(t0_scf)+0.3*(np.log10(N)-np.log10(npart0_scf)))#t0_scf*(np.log10(npart0_scf)/np.log10(N))
#     print 'Isolation'
#     print s2hr*Nproc_scf*tscf_scaleN(Npart_MW1 )
#     tscf_MW1 =np.ceil(Nsnap_iso*s2hr*Nproc_scf*tscf_scaleN(Npart_MW1 ))
#     tscf_MW3 =np.ceil(Nsnap_iso*s2hr*Nproc_scf*tscf_scaleN(Npart_MW3 ))
#     tscf_MW10=np.ceil(Nsnap_iso*s2hr*Nproc_scf*tscf_scaleN(Npart_MW10))
#     print '    MW1   {}'.format(tscf_MW1 )
#     print '    MW3   {}'.format(tscf_MW3 )
#     print '    MW10  {}'.format(tscf_MW10)
#     print 'Interactions'
#     tscf_q1 =np.ceil(Nsnap_int*s2hr*Nproc_scf*tscf_scaleN(Npart_MW1 + Npart_MW1 ))
#     tscf_q3 =np.ceil(Nsnap_int*s2hr*Nproc_scf*tscf_scaleN(Npart_MW1 + Npart_MW3 ))
#     tscf_q10=np.ceil(Nsnap_int*s2hr*Nproc_scf*tscf_scaleN(Npart_MW1 + Npart_MW10))
#     print '    q=1   {} ({} x {})'.format(tscf_q1 *Nsim_q1 ,tscf_q1 ,Nsim_q1 )
#     print '    q=3   {} ({} x {})'.format(tscf_q3 *Nsim_q3 ,tscf_q3 ,Nsim_q3 )
#     print '    q=10  {} ({} x {})'.format(tscf_q10*Nsim_q10,tscf_q10,Nsim_q10)
#     print 'Total={} SUs'.format(tscf_MW1+tscf_MW3+tscf_MW10+tscf_q1 *Nsim_q1+tscf_q3 *Nsim_q3+tscf_q10*Nsim_q10)
#     methLIST=[]
#     # GADGET
#     print 'GADGET'
#     npart0_gd=19800000
#     Nproc_gd=256
#     t0_gd=82472.6262858 # s
#     tgd_scaleN=lambda N: t0_gd*N*np.log10(N)/(npart0_gd*np.log10(npart0_gd))
#     print 'Isolation'
#     print s2hr*tgd_scaleN(Npart_MW1 )*Nproc_gd
#     tgd_MW1 =np.ceil(time_iso*s2hr*tgd_scaleN(Npart_MW1 )*Nproc_gd)
#     tgd_MW3 =np.ceil(time_iso*s2hr*tgd_scaleN(Npart_MW3 )*Nproc_gd)
#     tgd_MW10=np.ceil(time_iso*s2hr*tgd_scaleN(Npart_MW10)*Nproc_gd)
#     print '    MW1   {}'.format(tgd_MW1 )
#     print '    MW3   {}'.format(tgd_MW3 )
#     print '    MW10  {}'.format(tgd_MW10)
#     print 'Interactions'
#     tgd_q1 =np.ceil(time_int*s2hr*Nproc_gd*tgd_scaleN(Npart_MW1 + Npart_MW1 ))
#     tgd_q3 =np.ceil(time_int*s2hr*Nproc_gd*tgd_scaleN(Npart_MW1 + Npart_MW3 ))
#     tgd_q10=np.ceil(time_int*s2hr*Nproc_gd*tgd_scaleN(Npart_MW1 + Npart_MW10))
#     print '    q=1   {} ({} x {})'.format(tgd_q1 *Nsim_q1 ,tgd_q1, Nsim_q1)
#     print '    q=3   {} ({} x {})'.format(tgd_q3 *Nsim_q3 ,tgd_q3, Nsim_q3)
#     print '    q=10  {} ({} x {})'.format(tgd_q10*Nsim_q10,tgd_q10,Nsim_q10)
#     print 'Total={} SUs'.format(tgd_MW1+tgd_MW3+tgd_MW10+tgd_q1 *Nsim_q1+tgd_q3 *Nsim_q3+tgd_q10*Nsim_q10)
    
