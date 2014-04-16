####################################################################################################################################
#
# STRING CLASSES & METHODS
#
####################################################################################################################################
import string,copy
import mmlpars,mmlio

####################################################################################################################################
# METHOD TO RETURN LONG LINE
####################################################################################################################################
def linestr(width=100):
    return width*'-'

####################################################################################################################################
# METHOD TO CONVERT LIMITS TO STRING
def lim2str(lim,joinstr='to',**exkw):
    """
    Converts limits to string
    """
    lim=mmlpars.mml_pars(lim,type=tuple,nelements=2)
    limstr='{}{}{}'.format(dec2str(lim[0],**exkw),joinstr,dec2str(lim[1],**exkw))
    return limstr

####################################################################################################################################
# METHOD TO CREATE STRING
####################################################################################################################################
def val2str(val):
    '''
    NAME:
        mmlstring.val2str
    PURPOSE:
        To create a standardized string for any object.
    CALLING:
        valstr=val2str(val)
    ARGUMENTS:
        val:    Variable of any type
    OUTPUT:
        valstr: String version of val
    '''
    # Return the input if it's a string
    if   isinstance(val,str  ): valstr=val
    # Handle types where spaces are added
    elif isinstance(val,tuple): valstr=repr(val).replace(', ',',')
    elif isinstance(val,list ): valstr=repr(val).replace(', ',',')
    elif isinstance(val,dict ): valstr=repr(val).replace(', ',',').replace(': ',':')
    # Otherwise use repr()
    else: valstr=repr(val)
    # Return output
    return valstr

####################################################################################################################################
# METHOD TO RECOVER VALUE FROM STRING
####################################################################################################################################
def str2val(valstr,dtype=None):
    '''
    NAME:
        mmlstring.str2val
    PURPOSE:
        Recovering variables from a string.
    CALLING:
        val=str2val(valstr)
    ARGUMENTS:
        valstr: String to decipher
    KEYWORDS:
        dtype:  String specifying the data type
    OUTPUT:
        val:    Deciphered string
    '''
    valstr_strip=valstr.strip()
    # Get data type if not provided
    if not isinstance(dtype,str):
        if   valstr_strip.startswith('"') and valstr_strip.endswith('"'):
            valstr_strip=valstr_strip.strip('"')
            dtype='str'
        elif valstr_strip.lower() == 'nan': dtype='float'
        elif valstr_strip.endswith('L') and valstr_strip.rstrip('L').isdigit(): dtype='long'
        elif valstr_strip.isdigit():
            if int(valstr_strip) == long(valstr_strip): dtype='int'
            else                                      : dtype='long'
        elif valstr_strip.startswith("'") and valstr_strip.endswith("'"):
            valstr_strip=valstr_strip.strip("'")
            dtype='str'
        elif valstr_strip.lower() == 'none'          : dtype='none'
        elif valstr_strip.strip('.').lower() in ['true','false']: dtype='bool'
        elif valstr_strip.startswith('(') and valstr_strip.endswith(')'):
            if ',' in valstr: dtype='tuple'
            elif 'j' in valstr: dtype='complex'
            else: dtype='str'
        elif valstr_strip.startswith('[') and valstr_strip.endswith(']'): dtype='list'
        elif valstr_strip.startswith('{') and valstr_strip.endswith('}'): dtype='dict'
    # If you know the data type
    if isinstance(dtype,str):
        dtype=dtype.lower()
        # Floats
        if   dtype=='s' or 'str' in dtype:
            if   valstr_strip.startswith("'") and valstr_strip.endswith("'"): val=valstr_strip.strip("'")
            elif valstr_strip.startswith('"') and valstr_strip.endswith('"'): val=valstr_strip.strip('"')
            else: val=valstr_strip
        elif dtype=='f' or 'float'   in dtype: val=float(valstr_strip)
        elif dtype=='d' or 'double'  in dtype: val=float(valstr_strip)
        elif dtype=='i' or 'int'     in dtype: val=int(float(valstr_strip.replace('L','')))
        elif dtype=='l' or 'long'    in dtype: val=long(float(valstr_strip.replace('L','')))
        elif dtype=='c' or 'complex' in dtype: val=complex(valstr_strip.replace(' ',''))
        elif dtype=='n' or 'none'    in dtype: val=None
        elif dtype=='b' or 'bool'    in dtype:
            if   'true'  in valstr_strip.strip('.').lower(): val=True
            elif 'false' in valstr_strip.strip('.').lower(): val=False
            else: raise Exception('Unrecognised string for bool type: {}'.format(valstr))
        elif 'tuple' in dtype:
            vallist=valstr_strip[1:-1].split(',')
            val=[]
            for ivalstr in vallist:
                val.append(str2val(ivalstr))
            val=tuple(val)
        elif 'list' in dtype:
            vallist=valstr_strip[1:-1].split(',')
            val=[]
            for ivalstr in vallist:
                val.append(str2val(ivalstr))
        elif 'dict' in dtype:
            vallist=valstr_strip[1:-1].split(',')
            val={}
            for ivalstr in vallist:
                ikey=str2val(ivalstr.split(':')[0],dtype='str')
                ival=str2val(ivalstr.split(':')[1])
                val[ikey]=ival
        else: raise Exception('Unsupported type string: {}'.format(dtype))
    # If you don't know the data type try a few things
    else:
        # Try as complex
        try:
            # Try as complex
            if 'j' in valstr_strip:
                val=complex(valstr_strip.replace(' ',''))
            # Try as flaot
            else:
                val=float(valstr_strip)
        # Otherwise return as string
        except ValueError:
            val=valstr_strip
    # Return output
    return val

####################################################################################################################################
def dec2str(numval,formstr=None,decrep=None,negrep=None,nsig=None,trimzero=None):
    """
    Converts decimals into strings
    """
    # Pars input
    decrep=mmlpars.mml_pars(decrep,default='p',type=str)
    negrep=mmlpars.mml_pars(negrep,default='n',type=str)
    nsig=mmlpars.mml_pars(nsig,default=3,type=int,min=1)
    formstr=mmlpars.mml_pars(formstr,default='.{}g'.format(nsig),type=str)
    trimzero=mmlpars.mml_pars(trimzero,default=True,type=bool)
    # Skip if integer/long
    if isinstance(numval,(int,long)): raise Exception('[mmlstring.dec2str] This function is not intended for integers or longs.')
    # Round numval if within precision
    if trimzero:
        precval=10.**(-(nsig-1))
        if numval == 0:
            numval=int(numval)
        else:
            if int(numval) == 0:
                if abs(numval) < precval:               numval=int(numval)
            else:
                if abs(numval % int(numval)) < precval: numval=int(numval)
    # Set default format for integers/longs
    if isinstance(numval,(int,long)): formstr='d'
    # Create string
    decstr=(r'{:'+formstr+'}').format(numval)
    # Trim zeros
    if trimzero and '.' in decstr:
        if 'e' in decstr or 'E' in decstr:
            if 'e' in decstr: decbas,decexp=decstr.split('e')
            if 'E' in decstr: decbas,decexp=decstr.split('E')
            if 'e' and 'E' in decstr:
                raise Exception('[mmlstring.dec2str] Both e and E are in the string: {}'.format(decstr))
        else:
            decbas=copy.deepcopy(decstr)
            decexp=''
        decbas.rstrip('0')
        decbas.rstrip('.')
        decstr=decbas+decexp
    # Create tables to replace .- characters and remove + .
    rtab='+ '
    itab='.-'
    ftab=decrep+negrep
    trantab=string.maketrans(itab,ftab)
    # Return translated string
    return decstr.translate(trantab,rtab)

####################################################################################################################################
def str2dec(decstr,decrep=None,negrep=None):
    """
    Converts decimal strings into floats
    """
    # Pars input
    decrep=mmlpars.mml_pars(decrep,default='p',type=str)
    negrep=mmlpars.mml_pars(negrep,default='n',type=str)
    # Create tables to reinstate .- characters
    itab=decrep+negrep
    ftab='.-'
    trantab=string.maketrans(itab,ftab)
    # Translate string and convert to float
    numval=float(decstr.translate(trantab))
    return numval

####################################################################################################################################
# DECIMAL
####################################################################################################################################
def decimal(numval,formstr=None,decrep='p',negrep='n'):
    '''
    NAME:
        mmlstring.decimal
    PURPOSE:
        To generate a string from a number w/ or w/o decimals.
        Decimals and negative signs are replaced while plus signs and sapces are removed.
    CALLING:
        str=decimal(numval,formstr=None,decrep="p",negrep="n")
    ARGUMENTS:
        numval:  An integer, long, or float.
    KEYWORDS:
        formstr: A string that is used to format string.
                 Default for int/long: "d"
                 Default for float:    ".3g"
        decrep:  String to replace decimal points with.
        negrep:  String to replace negative sign with.
    OUTPUT:
        decstr:  Formated string version of numval.
    '''
    # Round numval if within precision
    if numval == 0:
        numval=int(numval)
    else:
        if int(numval) == 0:
            if numval < 0.01:                 numval=int(numval)
        else:
            if (numval % int(numval)) < 0.01: numval=int(numval)
    
    # Set default format
    if formstr==None:
        if isinstance(numval,(int,long)):
            formstr='d'
        else:
            formstr='.3g'

    # Create string
    decstr=(r'{:'+formstr+'}').format(numval)

    # Create tables to replace .- characters and remove + .
    rtab='+ '
    itab='.-'
    ftab=decrep+negrep
    trantab=string.maketrans(itab,ftab)

    # Return transplate string
    return decstr.translate(trantab,rtab)

    
