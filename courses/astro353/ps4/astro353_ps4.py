#!/usr/bin/python
import matplotlib.pyplot as plt
import os,math,array,sys
from pylab import *
from numpy import *

def problem1():
    print "Problem 1"
    

def problem2():
    print "Problem 2"

    # Determine mass in white dwarfs at different times
    nbin=100
    mmin_tot=0.1
    mmax_tot=100.
    t=linspace(0.00001,13.7,nbin)
    mmin_wd=(t/10.)**(-1./2.5)
    nwd=(1./1.35)*((mmin_wd**(-1.35)) - (mmax_tot**(-1.35)))
    mwd=0.6*nwd
    mtot=(1./0.35)*((mmin_tot**(-0.35)) - (mmax_tot**(-0.35)))
    mwd_mtot=mwd/mtot

    # Plot individual variables
    plotname1='/home/langmm/utils/python/mine/astro353/ps3/astro353_ps4_problem2c.png'
    fig=plt.figure()
    ax1=subplot(1,1,1)
    ax1.plot(t,mwd_mtot,label='data')
    ax1.set_ylabel('$M_{wd}/M_{tot}$')
    ax1.set_xlabel('t (Gyr)')
    ax1.set_xscale('linear')
    fig.tight_layout()
    fig.savefig(plotname1)
    plt.close(fig)
    print '    {}'.format(plotname1)

def problem3():
    print "Problem 3"

def main():
#    problem1()
    problem2()
#    problem3()


if __name__ == '__main__':
    main()
