#!/usr/bin/python
import matplotlib.pyplot as plt
import os,math,array,sys
from pylab import *
from numpy import *

def problem2():
    print "Problem 2"

    # Load in observation data from file
    ngc=[] ; r_e=[] ; sigma=[] ; I=[]
    f = open('/home/langmm/utils/python/mine/astro353/ps3/fplane.dat','r')
    i=0
    for line in f:
        if not line.startswith('#'):
            vars=line.split()
            ngc.append(float(vars[0]))
            r_e.append(float(vars[1]))
            sigma.append(float(vars[2]))
            I.append(float(vars[3]))
            i+=1
    nobs=i-1
    f.close
    ngc=array(ngc[:])
    r_e=array(r_e[:])
    sigma=array(sigma[:])
    I=array(I[:])

    # Fit
    logreFIT_sigma,sig_rvsigma=model(log10(sigma),log10(r_e)) ; reFIT_sigma=10.**logreFIT_sigma
    logreFIT_I,sig_rvI=model(log10(I),log10(r_e)) ; reFIT_I=10.**logreFIT_I

    # Plot individual variables
    plotname1='/home/langmm/utils/python/mine/astro353/ps3/astro353_ps3_problem2a.png'
    fig=plt.figure()
    ax1=subplot(2,1,1)
    ax1.plot(sigma,r_e,'bo',label='data')
    ax1.plot(sigma,reFIT_sigma,'b--',label='fit (scat={:3.3f})'.format(sig_rvsigma))
    ax1.set_xlabel('$\sigma$ (km/s)')
    ax1.set_ylabel('$r_e$ (kpc)')
    ax1.set_xscale('log') ; ax1.set_yscale('log')
    ax1.legend(loc=2)
    ax2=subplot(2,1,2)
    ax2.plot(I,r_e,'r^',label='data')
    ax2.plot(I,reFIT_I,'r--',label='fit (scat={:3.3f})'.format(sig_rvI))
    ax2.set_xlabel('$I$ (L$_{\odot}$/kpc$^2$)')
    ax2.set_ylabel('$r_e$ (kpc)')
    ax2.set_xscale('log') ; ax2.set_yscale('log')
    ax2.legend(loc=3)
    fig.tight_layout()
    fig.savefig(plotname1)
    plt.close(fig)

    # Get fundamental plane fit
    x=1.24 ; y=-0.82
    fund=(sigma**x)*(I**y)
    logreFIT_fund,sig_fund=model(log10(fund),log10(r_e)) ; reFIT_fund=10.**logreFIT_fund

    # Plot fundamental plane
    plotname2='/home/langmm/utils/python/mine/astro353/ps3/astro353_ps3_problem2b.png'
    fig=plt.figure()
    ax1=subplot(1,1,1)
    ax1.plot(fund,r_e,'gs',label='data')
    ax1.plot(fund,reFIT_fund,'g--',label='fit (scat={:3.3f})'.format(sig_fund))
    ax1.set_xlabel('$\sigma^{1.24}I^{-0.82}$')
    ax1.set_ylabel('$r_e$ (kpc)')
    ax1.set_xscale('log') ; ax1.set_yscale('log')
    ax1.legend(loc=2)
    fig.tight_layout()
    fig.savefig(plotname2)
    plt.close(fig)

def problem3():
    print "Problem 3"

    # Load in observation data from file
    sigma=[] ; rb=[] ; mub=[]
    f = open('/home/langmm/utils/python/mine/astro353/ps3/faber1997.dat','r')
    i=0
    for line in f:
        if not line.startswith('#'):
            vars=line.split()
            sigma.append(10.**float(vars[0]))
            rb.append(10.**float(vars[1]))
            mub.append(float(vars[2]))
            i+=1
    nobs=i-1
    f.close
    sigma=array(sigma[:])
    rb=array(rb[:])
    mub=array(mub[:])

    # Plot
    plotname1='/home/langmm/utils/python/mine/astro353/ps3/astro353_ps3_problem3.png'
    fig=plt.figure()
    ax1=fig.add_subplot(2,2,1)
    logrb_fit1,rb_sig1=model(log10(sigma),log10(rb)) ; rb_fit1=10.**logrb_fit1
    ax1.plot(sigma,rb,'bo',label='data')
    ax1.plot(sigma,rb_fit1,'b--',label='fit (scat={:3.3f})'.format(rb_sig1))
    ax1.set_xlabel('$\sigma_0$')
    ax1.set_ylabel('$r_{b}$')
    ax1.set_xscale('log') ; ax1.set_yscale('log')
    ax1.legend(loc=4)
#    ax1.legend(bbox_to_anchor=(0., 1.02, 1., .102), loc=3,
#               ncol=1, mode="expand", borderaxespad=0.)
    ax2=fig.add_subplot(2,2,2)
    logrb_fit2,rb_sig2=model(mub,log10(rb)) ; rb_fit2=10.**logrb_fit2
    ax2.plot(mub,rb,'r^',label='data')
    ax2.plot(mub,rb_fit2,'r--',label='fig (scat={:3.3f})'.format(rb_sig2))
    ax2.set_xlabel('$\mu_{b}$')
    ax2.set_ylabel('$r_{b}$')
    ax2.set_yscale('log')
    ax2.legend(loc=4)
#    ax2.legend(bbox_to_anchor=(0., 1.02, 1., .102), loc=3,
#               ncol=1, mode="expand", borderaxespad=0.)
    ax3=fig.add_subplot(2,2,3)
    logsig_fit,sig_sig1=model(mub,log10(sigma)) ; sig_fit=10.**logsig_fit
    ax3.plot(mub,sigma,'gs',label='data')
    ax3.plot(mub,sig_fit,'g--',label='fig (scat={:3.3f})'.format(sig_sig1))
    ax3.set_xlabel('$\mu_{b}$')
    ax3.set_ylabel('$\sigma_0$')
    ax3.set_yscale('log')
    ax3.legend(loc=2)
#    ax3.legend(bbox_to_anchor=(0., 1.02, 1., .102), loc=3,
#               ncol=1, mode="expand", borderaxespad=0.)
    ax4=fig.add_subplot(2,2,4)
    fund=10.**(0.766*log10(rb)-0.257*mub+2.97)
    logsig_fit2,sig_sig2=model(log10(fund),log10(sigma)) ; sig_fit2=10.**logsig_fit2
    ax4.plot(fund,sigma,'m+',label='data')
    ax4.plot(fund,sig_fit2,'m--',label='fig (scat={:3.3f})'.format(sig_sig2))
    ax4.set_xlabel('$r_{b}^{0.766}10^{-0.257\mu_{b}+2.97}$')
    ax4.set_ylabel('$\sigma_0$')
    ax4.set_xscale('log') ; ax4.set_yscale('log')
    ax4.legend(loc=2)
    fig.tight_layout()
    fig.savefig(plotname1)
    plt.close(fig)

def model(x,y):
    pfit=polyfit(x,y,1)
    yfit=pfit[1]+pfit[0]*x
    stddev=sqrt(((y-yfit)**2).sum()/float(len(x)-2))
    return yfit,stddev


def main():
#    problem1()
    problem2()
    problem3()
#    problem5()


if __name__ == '__main__':
    main()
