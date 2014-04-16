#!/usr/bin/python
from mmlutils import *
import numpy as np
import math,os,shutil,copy,cPickle,sys,types,pprint



def list_stabmeths():
    return ['op73','efstat82','toomre81']

####################################################################################################################################
# DISK STABILITY OSTRIKER & PEEBLES (1973)
def stability_op73(Trot,W):
    """
    Determines the stability of a disk based on OSTRIKER & PEEBLES (1973)
    """
    tOP_max=0.14 # +/- 0.03
    tOP=Trot/np.abs(W)
    mmlio.verbose('tOP={}'.format(tOP))
    if tOP > tOP_max:
        crit=False
    else:
        crit=True
    return crit

####################################################################################################################################
# DISK STABILITY Efstathiou et al.(1982)
def stability_efstat82(Vpeak,Md,Rd,G):
    alpha_min=1.1
    alpha=Vpeak/np.sqrt(G*Md/Rd)
    mmlio.verbose('alpha={}'.format(alpha))
    if alpha < alpha_min:
        crit=False
    else:
        crit=True
    return crit
    
####################################################################################################################################
# DISK STABILITY Julian & Toomre (1966), Toomre (1981)
def stability_toomre81(R,sigr,kappa,rhoS,G):
    Q_min=2.0
    Q=R*(kappa**2)/(4.*np.pi*G*rhoS)
    mmlio.verbose('Q={}'.format(Q))
    if Q < Q_min:
        crit=False
    else:
        crit=True
    return crit
