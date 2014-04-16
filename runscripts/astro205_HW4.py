#!/usr/bin/python
import numpy as np

mh=1.67e-27 # kg
Msol=1.99e30 # kg
G=6.67e-11 #kms
k=1.38e-23 #kms

def Mj1(n_m,T_K):
    Mj_kg=4.*((k*T_K/(G*mh))**(3./2.))*((n_m*mh)**(-1./2.))
    return Mj_kg

def Mj2(r_m,T_K):
    Mj_kg=3.*k*T_K*r_m/(G*mh)
    return Mj_kg

def Rj1(n_m,T_K):
    Rj_m=np.sqrt(15.*k*T_K/(8.*np.pi*G*(mh**2)*n_m))
    return Rj_m

def Rj2(M_kg,T_K):
    Rj_m=G*mh*M_kg/(2*k*T_K)
    return Rj_m

def rhoj2(M_kg,T_K):
    rho_kgm=(3./(4.*np.pi*(M_kg**2)))*((3.*k*T_K/(G*mh))**3.)
    return rho_kgm

def tff(rho_kgm):
    tff_s=np.sqrt(3.*np.pi/(32.*G*rho_kgm))
    return tff_s

def printstat(title,R_pc,n,T1,T2):
    T=(T1+T2)/2. # K
    R=3.09e16*R_pc
    rho=mh*n
    M=4.*np.pi*(R**3)*rho/3.
    print ''
    print title
    print 'T = {} K'.format(T)
    print 'M = {} Msol'.format(M/Msol)
    print 'M = {} kg'.format(M)
    print 'Mj1 = {} kg ({})'.format(Mj1(n,T),M>Mj1(n,T))
    print 'Mj2 = {} kg ({})'.format(Mj2(R,T),M>Mj2(R,T))
    print 'R = {} m'.format(R)
    print 'Rj1 = {} m ({})'.format(Rj1(n,T),R<Rj1(n,T))
    print 'Rj2 = {} m ({})'.format(Rj2(M,T),R<Rj2(M,T))
    print 'rho = {} kg/m**3'.format(rho)
    print 'rhoj2 = {} kg/m**3 ({})'.format(rhoj2(M,T),rho>rhoj2(M,T))
    print 'nj2 = {} 1/m**3 ({})'.format(rhoj2(M,T)/mh,n>(rhoj2(M,T)/mh))
    print 'tff = {} s'.format(tff(rho))
    print 'tff = {} yr'.format(tff(rho)/3.16e7)

# Problem 2
print ''
print ''
print 50*'-'
print 'Problem 2'
print 50*'-'
printstat('Part a (Tavg)',1.0,1.0e10, 30.,100.)
# printstat('Part a (Tmin)',1.0,1.0e10, 30., 30.)
# printstat('Part a (Tmax)',1.0,1.0e10,100.,100.)
printstat('Part b (Tavg)',0.1,1.0e12, 30.,200.)
# printstat('Part b (Tmin)',0.1,1.0e12, 30., 30.)
# printstat('Part b (Tmax)',0.1,1.0e12,200.,200.)


print ''
