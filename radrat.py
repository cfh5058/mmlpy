#!/usr/bin/python
import math

def main(x,c=10.):
    print x/c
    print mass(x)/mass(c)
    return


def mass(x):
    M=math.log(1.+x)-x/(1.+x)
    return M

if __name__ == '__main__':
    main()
