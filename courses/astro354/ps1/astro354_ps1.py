#!/usr/bin/python
import numpy as np
import matplotlib.pyplot as plt
import scipy.integrate as integrate

def main():
    print "Astro 354: PS1"
    # Set constants
    mmin=0.08
    mmax=age2mmax(np.array([0.0]))
    ltot=1.0
    # Determine constant in terms of Ltot
    k=get_imfcoef(ltot,mmin,mmax[0])
    ltot=get_imflumin(k,mmin,mmax) # Check that integration constant sets this to 1
    print '    Ltot={}'.format(ltot)
    print '    C={}'.format(k)

    # PROBLEM 1
    print "Problem 1"
    # Get magnitude differences as a funciton of time
    tlist=np.linspace(0.,10.,100.)
    mmaxlist=age2mmax(tlist)
    llist=get_imflumin(k,mmin,mmaxlist)
    magdiff=-2.5*np.log(llist/ltot)
    # Plot
    plotname1='/home/langmm/utils/python/mine/courses/astro354/ps1/astro354_ps1_problem1.png'
    fig1=plt.figure()
    ax11=plt.subplot(1,1,1)
    ax11.plot(tlist,magdiff)
    ax11.set_xlabel('Time (Gyr)')
    ax11.set_ylabel('M(t)-M(0)')
    fig1.tight_layout()
    fig1.savefig(plotname1)
    plt.close(fig1)
    print '    {}'.format(plotname1)

    # PROBLEM 2
    print "Problem 2"
    # Determine constant in terms of Ltot
    colorlist=integratecolor(mmin,mmaxlist)
    # Plot
    plotname2='/home/langmm/utils/python/mine/courses/astro354/ps1/astro354_ps1_problem2.png'
    fig2=plt.figure()
    ax21=plt.subplot(1,1,1)
    ax21.plot(tlist,colorlist)
    ax21.set_xlabel('Time (Gyr)')
    ax21.set_ylabel('g-r Color')
    fig2.tight_layout()
    fig2.savefig(plotname2)
    plt.close(fig2)
    print '    {}'.format(plotname2)

    # PROBLEM 3
    print 'Problem 3'
    # Get g-r and M-M0 at specified times
    tlist_label=np.array([0.01,0.1,1.,2.,5.,10.])
    labellist=['10 Myr','100 Myr','1 Gyr','2 Gyr','5 Gyr','10 Gyr']
    mmaxlist_label=age2mmax(tlist_label)
    llist_label=get_imflumin(k,mmin,mmaxlist_label)
    magdiff_label=-2.5*np.log(llist_label/ltot)
    colorlist_label=integratecolor(mmin,mmaxlist_label)
    # Plot
    plotname3='/home/langmm/utils/python/mine/courses/astro354/ps1/astro354_ps1_problem3.png'
    fig3=plt.figure()
    ax31=plt.subplot(1,1,1)
    ax31.plot(colorlist_label,magdiff_label,'bo')
    ax31.set_xlabel('g-r Color')
    ax31.set_ylabel('M(t)-M(0)')
    for label, x, y in zip(labellist,list(colorlist_label),list(magdiff_label)):
        plt.annotate(label,
                     xy=(x,y),xytext=(-20,20),
                     textcoords = 'offset points', ha = 'right', va = 'bottom',
                     bbox = dict(boxstyle = 'round,pad=0.5', fc = 'yellow', alpha = 0.5),
                     arrowprops = dict(arrowstyle = '->', connectionstyle = 'arc3,rad=0'))
    fig3.tight_layout()
    fig3.savefig(plotname3)
    plt.close(fig3)
    print '    {}'.format(plotname3)
    
    return

# Get coefficient for Salpeter IMF
def get_imfcoef(ltot,mmin,mmax):
    exp=(3.5-2.35)+1
    k=(mmax**exp - mmin**exp)/exp
    return ltot/k
# Get mass from Salpeter IMF
def get_imfmass(k,mmin,mmax):
    exp=(1.-2.35)+1.
    nele=len(mmax)
    mminlmmax=(mmin< mmax)
    mmingmmax=(mmin>=mmax)
    mtot=k*(mmax**exp - mmin**exp)/exp
    mtot[mmingmmax]=0.
    return mtot
# Get luminosity from Salpeter IMF
def get_imflumin(k,mmin,mmax):
    exp=(3.5-2.35)+1.
    nele=len(mmax)
    ltot=np.zeros(nele)
    minlmax=(mmin>=mmax)
    mingmax=(mmin< mmax)
    ltot=ltot=k*(mmax**exp - mmin**exp)/exp
    ltot[minlmax]=0.
    return ltot
# Get most massive galaxy alive
def age2mmax(age):
    mmax=50.0
    nele=len(age)
    ageg0=(age > 0)
    mass=np.ones(nele)*mmax
    mass[ageg0]=(age[ageg0]/10.0)**(-1./2.5)
    massgmmax=(mass>mmax)
    mass[massgmmax]=mmax
    return mass

# Get color from IMF (involves integration)
def integratecolor(mmin,mmax):
    nele=len(mmax)
    numer=np.zeros(nele)
    for iele in range(nele):
        numer[iele] = np.exp(0.06)*integrate.quad(lambda x: colorintfunc(x),mmin,mmax[iele])[0]
    exp = -2.35 + 1.
    denom = (mmax**exp - mmin**exp)/exp
    Lg_Lr=numer/denom
    gr_color=-2.5*np.log10(Lg_Lr)
    return gr_color

# Function for color integration
def colorintfunc(M):
    f=((M+2.)**(-0.92))*(M**(-1.43))
    return f

if __name__ == '__main__':
    main()
