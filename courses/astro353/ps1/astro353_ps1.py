#!/usr/bin/python
import matplotlib.pyplot as plt
import os,math,array,sys
from pylab import *
from numpy import *

def problem1():
    sig_sb=5.67e-8 # W m^-2 T^-4
    Lsol=3.84e26 # W
    pc2m=3.086e16 # m
    Lmw=(2.0e10)*Lsol # W
    Rpc=5000.0
    Rm=Rpc*pc2m
    Teff=(Lmw/(4*math.pi*(Rm**2)*sig_sb))**0.25

    print "Problem 1"
    print "   Lmw=%(Lmw)f W" % {'Lmw':Lmw}
    print "   R=%(Rm)f m" % {'Rm':Rm}
    print "   Teff=%(Teff)f K" % {'Teff':Teff}
    

def intnmag(aprmags,magbins):
    nmag=len(magbins)
    cumhist=zeros([nmag],dtype=float)
    histmag=histogram(aprmags,bins=magbins)
    for i in range(0,len(histmag)):
        cumhist[i+1]=sum(histmag[0:i])
    return cumhist

    
def problem2():
    print "Problem 2"
    plotnameA="astro353_ps1_2a.png"
    plotnameB="astro353_ps1_2b.png"

    # Set constants
    xsol=0.0
    zsol=0.0
    absmag=5.0
    nstars=30000

    # Load in observation data from file
    amagobs=zeros([nstars],dtype=float)
    glatobs=zeros([nstars],dtype=float)
    glonobs=zeros([nstars],dtype=float)
    f = open('/home/langmm/utils/python/astro353/data/starcounts.dat','r')
    i=0
    for line in f:
        if not line.startswith('#'):
            vars=line.split()
            glonobs[i]=float(vars[0])*(math.pi/180.0)
            glatobs[i]=float(vars[1])*(math.pi/180.0)
            amagobs[i]=float(vars[2])
            i+=1
    nobs=i-1
    f.close

    # Load in real data from file
    ndat=18
    amagdat=zeros([ndat],dtype=float)
    z_histdat=zeros([ndat],dtype=float)
    x_histdat=zeros([ndat],dtype=float)
    f = open('/home/langmm/utils/python/astro353/data/NmAPQ3.txt','r')
    i=0
    for line in f:
        if not line.startswith('#'):
            vars=line.split()
            amagdat[i]=float(vars[0])
            z_histdat[i]=10**float(vars[1])
            x_histdat[i]=10**float(vars[2])
            i+=1
    ndat=i
    f.close

    # Create model in galactic coordinats
    nuni=1e6
    duni=8.0*random.ranf(size=[nuni])
    glonuni=2.0*math.pi*random.ranf(size=[nuni])
    glatuni=arcsin(2.0*random.ranf(size=[nuni])-1.0)
    amaguni=absmag+5.0*log10(100.0*duni)

    # Convert coordinates
    dobs=(10**((amagobs-absmag)/5.0))/100.0
    xobs=dobs*cos(glatobs)*cos(glonobs)-xsol
    yobs=dobs*cos(glatobs)*sin(glonobs)
    zobs=dobs*sin(glatobs)

    xuni=duni*cos(glatuni)*cos(glonuni)-xsol
    yuni=duni*cos(glatuni)*sin(glonuni)
    zuni=duni*sin(glatuni)
    
    # Plot data
    plt.figure(1)
    ax1=plt.subplot(211,aspect='equal')
    ax1.plot(xobs,zobs,'bo',xsol,zsol,'r^')
    ax1.legend(('Stars (from data)','Sun'),'upper left')
    plt.xlabel('x (kpc)')
    plt.ylabel('z (kpc)')
    plt.ylim=(ax1.get_xlim())

    ax2=plt.subplot(212,aspect='equal',sharex=ax1,sharey=ax1)
    ax2.plot(xuni,zuni,'bo',xsol,zsol,'r^')
    ax2.legend(('Stars (uniform)','Sun'),'upper left')
    plt.xlabel('x (kpc)')
    plt.ylabel('z (kpc)')
    plt.savefig(plotnameA)
    plt.close(1)
    print "   %(plotnameA)s" % {'plotnameA':plotnameA}

    # Find correct data along pencil ray
    rpen=1.0 # kpc
    xcutobs = ((xobs>=0)-1)+ \
              ((yobs>=-rpen)-1)+((yobs<=rpen)-1)+ \
              ((zobs>=-rpen)-1)+((zobs<=rpen)-1)
    zcutobs = ((xobs>=-rpen)-1)+((xobs<=rpen)-1)+ \
              ((yobs>=-rpen)-1)+((yobs<=rpen)-1)+ \
              ((zobs>=0)-1)
    xcutuni = ((xuni>=0)-1)+ \
              ((yuni>=-rpen)-1)+((yuni<=rpen)-1)+ \
              ((zuni>=-rpen)-1)+((zuni<=rpen)-1)
    zcutuni = ((xuni>=-rpen)-1)+((xuni<=rpen)-1)+ \
              ((yuni>=-rpen)-1)+((yuni<=rpen)-1)+ \
              ((zuni>=0)-1)

    x_amagobs=amagobs[where(xcutobs==0)]
    z_amagobs=amagobs[where(zcutobs==0)]
    x_amaguni=amaguni[where(xcutuni==0)]
    z_amaguni=amaguni[where(zcutuni==0)]
    x_amagdat=amagdat
    z_amagdat=amagdat
    
    # Create bins for histogram
    nmagbins=100
    magmin=array([x_amagobs.min(),x_amaguni.min(),x_amagdat.min(),
            z_amagobs.min(),z_amaguni.min(),z_amagdat.min()])
    magmax=array([x_amagobs.max(),x_amaguni.max(),x_amagdat.max(),
            z_amagobs.max(),z_amaguni.max(),z_amagdat.max()])
    magrange=[magmin.min(),magmax.max()]
    xnorm_histobs,x_binobs=histogram(x_amagobs,bins=nmagbins,range=magrange,normed=True)
    znorm_histobs,z_binobs=histogram(z_amagobs,bins=nmagbins,range=magrange,normed=True)
    xnorm_histuni,x_binuni=histogram(x_amaguni,bins=nmagbins,range=magrange,normed=True)
    znorm_histuni,z_binuni=histogram(z_amaguni,bins=nmagbins,range=magrange,normed=True)
    x_binobs=x_binobs[1:]
    z_binobs=z_binobs[1:]
    x_binuni=x_binuni[1:]
    z_binuni=z_binuni[1:]
    x_bindat=amagdat
    z_bindat=amagdat
    xnorm_histdat=x_histdat/(x_histdat.sum())
    znorm_histdat=z_histdat/(z_histdat.sum())

    # Plot histogram
    # x real vs. uni
    figure(2)
    ax1=plt.subplot(221)
    ax1.plot(x_bindat,xnorm_histdat,'bo',
             x_binuni,xnorm_histuni,'r^')
    ax1.legend(('Real Counts','Model (uniform)'),'upper left')
    plt.xlabel('apparent magnitude')
    plt.ylabel('N(m)/Ntot (positive x)')
    # x real vs. model
    ax2=plt.subplot(222,sharex=ax1,sharey=ax1)
    ax2.plot(x_bindat,xnorm_histdat,'bo',
              x_binobs,xnorm_histobs,'r^')
    ax2.legend(('Real Counts','Model (from file)'),'upper left')
    ax2.ylim=(ax1.get_ylim())
#    plt.xlabel('apparent magnitude')
#    plt.ylabel('log N(m) (positive x)')
    # z real vs. uni
    ax3=plt.subplot(223,sharex=ax1,sharey=ax1)
    ax3.plot(z_bindat,znorm_histdat,'bo',
             z_binuni,znorm_histuni,'r^')
    ax3.legend(('Real Counts','Model (uniform)'),'upper left')
    plt.ylim=(ax1.get_ylim())
    plt.xlabel('apparent magnitude')
    plt.ylabel('N(m)/Ntot (positive z)')
    # z real vs. model
    ax4=plt.subplot(224,sharex=ax1,sharey=ax1)
    ax4.plot(z_bindat,znorm_histdat,'bo',
             z_binobs,znorm_histobs,'r^')
    ax4.legend(('Real Counts','Model (from file)'),'upper left')
    ax4.ylim=(ax1.get_ylim())
#    plt.xlabel('apparent magnitude')
#    plt.ylabel('log N(m) (positive z)')
    plt.savefig(plotnameB)
    plt.close(2)
    print "   %(plotnameB)s" % {'plotnameB':plotnameB}

def dmag(mag,absmag):
    dpc=10.0*(10**((mag-absmag)/5.0))
#    print "D=%(dpc)f pc" % {'dpc':dpc}
    return dpc

def phi(mmin,mmax,absmag,starden_pc,fov_deg):
    omega_sph=41253.0 # deg^2
    dmin=dmag(mmin,absmag)
    dmax=dmag(mmax,absmag)
    Vshell=(4.0*math.pi/3.0)*(dmax**3 - dmin**3)
    Vfov=Vshell*(fov_deg/omega_sph)
    Nfov=starden_pc*Vfov
    phi_mag=Nfov

#    print "Phi(m)=%(phi_mag)f stars" % {'phi_mag':phi_mag}
    return phi_mag

def problem3():
    print "Problem 3"

    # Set constants
    fov_arcmin=15.0*15.0
    fov_deg=fov_arcmin/(60*60)
    absmag=5.0
    starden_pc=1.0/25.0
    nnum=13

    # a)
    dmax=dmag(20.0,absmag)
    print "a) dmax=%(dmax)f pc" % {'dmax':dmax}

    # b)
    n1920=phi(19,20,absmag,starden_pc,fov_deg)
    print "b) n(19<m<20)=%(n1920)f" % {'n1920':n1920}

    # c)
    n1213=phi(12,13,absmag,starden_pc,fov_deg)
    print "c) n(12<m<13)=%(n1213)f" % {'n1213':n1213}

    # d)
    dmin=starden_pc**(-1.0/3.0)
    mmin=absmag+5*log10(dmin/10.0)
    print "d) mmin=%(mmin)f" % {'mmin':mmin}
    print "   dmin=%(dmin)f pc" % {'dmin':dmin}

    # e)
    fname="astro353_ps1_3e.png"
    # Calculate number of stars
    maglist=linspace(8,20,nnum)
    numlist=zeros((nnum))
    for im in range(0,nnum):
        m1=maglist[im]
        m2=m1+1
        numlist[im]=phi(m1,m2,absmag,starden_pc,fov_deg)
    lognumlist=log10(numlist)
    slope=(lognumlist[3]-lognumlist[2])/(maglist[3]-maglist[2])
    intercept=lognumlist[3]-slope*maglist[3]
    magone=(0.0-intercept)/slope
    print "   mone=%(magone)f" % {'magone':magone}

    # Plot
    plt.figure(2)
    plt.subplot(211)
    plt.plot(maglist,numlist)
    plt.xlabel('m')
    plt.ylabel('Phi(m)')

    plt.subplot(212)
    plt.plot(maglist,lognumlist)
    plt.xlabel('m')
    plt.ylabel('LOG[Phi(m)]')
    plt.savefig(fname)
    plt.close(2)
    print "e) %(file)s" % {'file':fname}

def problem5():
    print "Problem 5"
    fname="astro353_ps1_5.png"
    
    # Load in observation data from file
    nstars_data=3641
    parallax=zeros([nstars_data],dtype=float)
    parallax_err=zeros([nstars_data],dtype=float)
    aprvmag=zeros([nstars_data],dtype=float)
    bvcolor=zeros([nstars_data],dtype=float)
    bvcolor_err=zeros([nstars_data],dtype=float)
    f = open('/home/langmm/utils/python/astro353/data/hipparcos.dat','r')
    i=0
    for line in f:
        if not line.startswith('#'):
            vars=line.split()
            parallax[i]=float(vars[0])/1000
            parallax_err[i]=float(vars[1])/1000
            aprvmag[i]=float(vars[2])
            bvcolor[i]=float(vars[3])
            bvcolor_err[i]=float(vars[4])
            i+=1
    f.close
#    print i

    # Load in theorhetical data from file
    nstars_theo=63
    mass_theo=zeros([nstars_theo],dtype=float)
    absvmag_theo=zeros([nstars_theo],dtype=float)
    bvcolor_theo=zeros([nstars_theo],dtype=float)
    f = open('/home/langmm/utils/python/astro353/data/zams_02.dat','r')
    i=0
    for line in f:
        if not line.startswith('#'):
            vars=line.split()
            mass_theo[i]=float(vars[0])
            absvmag_theo[i]=float(vars[1])
            bvcolor_theo[i]=float(vars[2])
            i+=1
    f.close
#    print i

    # Calculate distance
    d=1.0/parallax
    d_err=sqrt(((-1/parallax)**4)*(parallax_err**2))

    # Calculate absolute magnitudes
    absvmag=aprvmag-5*log(d/10)/log(10)
    absvmag_err=sqrt(((-50/(d*log(10)))**2)*(d_err**2))

    # Plot
    plt.figure(1)
    ax1=plt.subplot(111)
    plt.errorbar(bvcolor,absvmag,xerr=bvcolor_err,yerr=absvmag_err,
                 fmt='o',ecolor='g',label='Data')
    plt.errorbar(bvcolor_theo,absvmag_theo,fmt='r-',label='Theory')
    ax1.invert_yaxis()
    plt.legend()
#    plt.legend(('label1','label2'),'upper right')
    plt.xlabel('B-V')
    plt.ylabel('Mv')
    plt.savefig(fname)
    plt.close(1)
    print "   %(file)s" % {'file':fname}


def problem6():
    print "Problem 6"
    fname="astro353_ps1_6.png"
    Ebv=0.35
    
    # Load in observation data from file
    nstars_data=140
    aprvmag=zeros([nstars_data],dtype=float)
    bvcolor=zeros([nstars_data],dtype=float)
    f = open('/home/langmm/utils/python/astro353/data/myades.dat','r')
    i=0
    for line in f:
        if not line.startswith('#'):
            vars=line.split()
            aprvmag[i]=float(vars[0])
            bvcolor[i]=float(vars[1])
            i+=1
    f.close
#    print i

    # Load in theorhetical data from file
    nstars_theo=63
    mass_theo=zeros([nstars_theo],dtype=float)
    absvmag_theo=zeros([nstars_theo],dtype=float)
    bvcolor_theo=zeros([nstars_theo],dtype=float)
    f = open('/home/langmm/utils/python/astro353/data/zams_02.dat','r')
    i=0
    for line in f:
        if not line.startswith('#'):
            vars=line.split()
            mass_theo[i]=float(vars[0])
            absvmag_theo[i]=float(vars[1])
            bvcolor_theo[i]=float(vars[2])
            i+=1
    f.close
    # print i

    # Find the V-band extinction
    Av=3.1*Ebv
    print "   Av=%(Av)f" % {'Av':Av}

    # Correct the V-band magnitudes
    V0=aprvmag-Av

    # Correct the B-V color
    bv0=bvcolor-Ebv

    # Correct the magnitude for distance
    V_dcorr=4.5
    dcorr=10*(10**(V_dcorr/5))
    print "   d=%(dcorr)f pc" % {'dcorr':dcorr}

    # Plot
    plt.figure(1)
    ax1=plt.subplot(111)
    plt.plot(bv0,V0,'bo',label='Observations')
    plt.plot(bv0,V0-V_dcorr,'g+',label='Corrected Observations')
    plt.plot(bvcolor_theo,absvmag_theo,'r-',label='Theory')
    ax1.invert_yaxis()
    plt.legend()
    plt.xlabel('B-V')
    plt.ylabel('Mv,V')
    plt.savefig(fname)
    plt.close(1)
    print "   %(file)s" % {'file':fname}


def main():
#    problem1()
#    problem2()
#    problem3()
#    problem5()
    problem6()


if __name__ == '__main__':
    main()
