#!/usr/bin/python

####################################################################################################################################
# METHOD TO CALL GENERIC IMF
####################################################################################################################################
def main(imf,m,integral=False,**extra_kw):
    '''
    NAME:
        mml_imf.main
    PURPOSE:
        To return results for chosen IMF
    CALLING:
        dn_dm=main(imf,m[,integral=,**extra_kw])
    ARGUMENTS:
        imf:      Str specifying what IMF to use
        m:        Float mass to compute IMF at
    KEYWORDS:
        integral: Bool specifying if the integral value should be returned.
        Additional keywords are spaced to the specific IMF method that is called.
    OUTPUT:
        dn_dm:
    '''
    imf=mpart.mml_pars(imf,type=str,list=['salpeter','scalo','kroupa','chabrier'])
    if imf == 'salpeter':
        dn_dm=salpeter(m,integral=integral,**extra_kw)
    elif imf == 'scalo':
        dn_dm=scalo(m,integral=integral,**extra_kw)
    elif imf == 'kroupa':
        dn_dm=kroupa(m,integral=integral,**extra_kw)
    elif imf == 'chabrier':
        dn_dm=chabrier(m,integral=integral,**extra_kw)
    else:
        raise Exception('Unsupported imf value of {}'.format(imf))
    return dn_dm

####################################################################################################################################
# Method to compute Salpeter IMF
def salpeter(m,integral=False,alpha=2.35):
    '''
    the Salpeter 1955 IMF: dn/dm ~ m^-2.35
    (adopted from Adam Ginsburg"s imf.py)
    '''
    if integral: alpha -= 1
    return m**(-alpha)

####################################################################################################################################
# Method to compute Scalo IMF
def scalo(m,integral=False):
    '''
    '''
    exp1=1.83
    exp2=3.27
    exp3=2.45
    if integral:
        exp1 += 1 ; exp2 -= 1 ; exp3 -= 1
    zeta = 0.                                                  * (m < 0.08)
    zeta += (m**+exp1 / 0.2**+exp1 * 0.2**-exp2) * (m >= 0.08) * (m < 0.2 )
    zeta += (m**-exp2)                           * (m >= 0.2 ) * (m < 10.0)
    zeta += (m**-exp3 / 10.**-exp3 * 10.**-exp2) * (m >= 10.0)
    return zeta

####################################################################################################################################
# Method to compute Kroupa IMF
def kroupa(m,integral=False):
    '''
    Kroupa 2001 IMF (http://arxiv.org/abs/astro-ph/0009005, http://adsabs.harvard.edu/abs/2001MNRAS.322..231K)
    (adopted from Adam Ginsburg"s imf.py)
    '''
    exp1 = 0.3
    exp2 = 1.3
    exp3 = 2.3
    if integral:
        exp1 += 1 ; exp2 -= 1 ; exp3 -= 1
    zeta =  0.                                                   * (m < 0.01)
    zeta += (m**+exp1 / 0.08**+exp1 * 0.08**-exp2) * (m >= 0.01) * (m < 0.08)
    zeta += (m**-exp2)                             * (m >= 0.08) * (m < 0.5 )
    zeta += (m**-exp3 / 0.50**-exp3 * 0.50**-exp2) * (m >= 0.5 )
    return zeta

####################################################################################################################################
# Method to compute Chabrier IMF
def chabrier(m,integral=False):
    '''
    Chabrier 2003 IMF
    http://adsabs.harvard.edu/abs/2003PASP..115..763C
    (only valid for m < 1 msun)
    not sure which of these to use...
    integral is NOT IMPLEMENTED
    (adopted from Adam Ginsburg"s imf.py)
    '''
    if integral: raise Exception("Chabrier integral NOT IMPLEMENTED")
    # This system MF can be parameterized by the same type of lognormal form as
    # the single MF (eq. [17]), with the same normalization at 1 Msun, with the
    # coefficients (Chabrier 2003)
    return 0.86 * np.exp(-1*(np.log10(m)-np.log10(0.22))**2/(2*0.57**2))
    # This analytic form for the disk MF for single objects below 1 Msun, within these uncertainties, is given by the following lognormal form (Chabrier 2003):
    return 0.158 * np.exp(-1*(np.log10(m)-np.log10(0.08))**2/(2*0.69**2))


if __name__ == '__main__':
    main()
