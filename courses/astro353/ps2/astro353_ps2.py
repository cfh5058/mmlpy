#!/usr/bin/python
import matplotlib.pyplot as plt
import os,math,array,sys
from pylab import *
from numpy import *

def problem4(r):
    u=7.66*((r/2.0)**0.25)
    part1=exp(-u)
    tot=0.0
    for n in range(8):
        tot+=part1*(1.0/float(math.factorial(n)))*(u**n)
    print tot

def problem5():
    mu_disk=24.5
    r_disk=60.0
    mu_bulg=22.5
    r_bulg=20.0

    part1=2.0/7.22
    part2=(r_disk/r_bulg)**2.0
    part3=10**(-(mu_disk-mu_bulg)/2.5)
    frac=1.0/(1.0+(part1*part2*part3))
    print frac

def main():
    problem4(11.152)
    problem5()


if __name__ == '__main__':
    main()
