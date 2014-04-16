#!/usr/bin/python
import matplotlib.pyplot as plt
import os,math,array,sys
from pylab import *
from numpy import *

def problem1():
    print "Problem 1"
    

def problem3():
    print "Problem 3"

    # Determine mass in white dwarfs at different times
    nbin=100
    rmin=0.001
    rmax=1000.
    r_b=logspace(log10(rmin),log10(rmax),nbin)
    rho_rho0=(r_b**2 + 1)**(-5./2.)
    sig_sig0=(r_b**2 + 1)**(-2.0)
    Rc=sqrt(sqrt(2.)-1)

    # Plot individual variables
    plotname1='/home/langmm/utils/python/mine/astro353/ps5/astro353_ps5_problem3c.png'
    fig=plt.figure()
    ax1=subplot(2,1,1)
    ax1.plot(r_b,rho_rho0,label=r'$\rho/\rho_0$')
    ax1.set_xlabel('$r/b$')
    ax1.set_ylabel(r'$\rho(r)/\rho(0)$')
    ax1.set_xscale('log')
    ax1.set_yscale('log')
    ax1.set_ylim((1.0e-16,10))
    ax1.axvline(Rc,ls='--',label='$R_c$')
    ax1.legend()
    ax2=subplot(2,1,2)
    ax2.plot(r_b,sig_sig0,label='$\Sigma/\Sigma_0$')
    ax2.set_xlabel('$R/b$')
    ax2.set_ylabel('$\Sigma(R)/\Sigma(0)$')
    ax2.set_xscale('log')
    ax2.set_yscale('log')
    ax2.set_ylim((1.0e-13,10))
    ax2.axvline(Rc,ls='--',label='$R_c$')
    ax2.legend()
    fig.tight_layout()
    fig.savefig(plotname1)
    plt.close(fig)
    print '    {}'.format(plotname1)

def problem4():
    print "Problem 4"

    # Determine density & potential at range of radii
    nbin=100
    rmin=0.001
    rmax=10000.
    clrList=['b','r','g']
    etaList=[1.,2.,3.]
    rhoList_bhof=[] ; potList_bhof=[] ; rhoList_bhon=[] ; potList_bhon=[]
    r=logspace(log10(rmin),log10(rmax),nbin)
    for ieta in etaList:
        rhoList_bhof.append(eta_model_rho(r,eta=ieta,bhon=False))
        potList_bhof.append(eta_model_pot(r,eta=ieta,bhon=False))
        rhoList_bhon.append(eta_model_rho(r,eta=ieta,bhon=True))
        potList_bhon.append(eta_model_pot(r,eta=ieta,bhon=True))

    # Plot density & potential of eta model vs. radius
    for ifig in range(2):
        if ifig==0:
            plotname1='/home/langmm/utils/python/mine/astro353/ps5/astro353_ps5_problem4b.png'
            rhoList=rhoList_bhof
            potList=potList_bhof
        else:
            plotname1='/home/langmm/utils/python/mine/astro353/ps5/astro353_ps5_problem4c.png'
            rhoList=rhoList_bhon
            potList=potList_bhon
        fig=plt.figure()
        ax1=subplot(2,1,1)
        for irho in range(len(rhoList)):
            ax1.plot(r,rhoList[irho],label=r'$\eta={}$'.format(etaList[irho]),color=clrList[irho])
        ax1.set_xlabel('$r$')
        ax1.set_ylabel(r'$\rho(r)$')
        ax1.set_xscale('log')
        ax1.set_yscale('log')
        ax1.set_ylim((1.0e-16,1.0e6))
        ax1.legend()
        ax2=subplot(2,1,2)
        for ipot in range(len(potList)):
            ax2.plot(r,abs(potList[ipot]),label=r'$\eta={}$'.format(etaList[ipot]),color=clrList[ipot])
        ax2.set_xlabel('$r$')
        ax2.set_ylabel(r'$|\Phi_{\eta}(r)|$')
        ax2.set_xscale('log')
        ax2.set_yscale('log')
        ax2.set_ylim((1.0e-7,10))
        ax2.legend()
        fig.tight_layout()
        fig.savefig(plotname1)
        plt.close(fig)
        print '    {}'.format(plotname1)

def eta_model_rho(r,eta=1,bhon=False):
    rho=(eta/4.*math.pi)/((r**(3.-eta))*((1.+r)**(1.+eta)))
    return rho
def eta_model_pot(r,eta=1,bhon=False):
    if eta==1:
        pot=-log(1.+1./r)
#        pot=-1./(r+1.)
#    elif eta==2:
#        pot=-r/((r+1.)**2.)
#    elif eta==3:
#        pot=-(r**2)/((r+1.)**3.)
    else:
        pot=(1./(eta-1.))*((r**(eta-1.))/((1.+r)**(eta-1.)) - 1.)
#        raise Exception('Unsupported eta value. Integrate numerically.')
    if bhon: pot+=(-0.02/r)
    return pot

def main():
#    problem1()
#    problem2()
#    problem3()
    problem4()


if __name__ == '__main__':
    main()
