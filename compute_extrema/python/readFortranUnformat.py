__author__    = "Marco Mancini (MM)"
__modifiers__ = ''
__date__      = '02/06/2015'
__version__   = '1.0'


"""Module containing a simple function to read unformatted fortran files."""

import numpy as np
import struct

def readFortranUnformat(filename,endian='>',dtype='f',ntype=4):
    """
    Function which permits to read "unformatted" fortran file.
    The file has to contains homogenous arrays (in the sense that only
    a type is permitted in a record, and any record must have the same
    number of components.



    Parameters
    ----------
    filename : str
          name of the file to read.
    endian : character, optional
        Specify the endian-ness of the file.  Possible values are
        '>', '<', '@' and '='.  See the documentation of Python's
        struct module for their meanings.  The deafult is '>'.
    dtype : str
        Specify the type to read. May be a simple character or a string.
        Possible values are : 'f','d','i' and 'l'. The deafult is 'd'.
    ntype : int
        Number of byte of the header. The deafult is 4.

    """
    if endian in '<>@=':
        if endian == '@': endian = '='
    else:
        raise ValueError('Cannot set endian-ness to: '+endian)

    use_struct = len(dtype) > 1

    for ii in dtype:
        if ii not in 'ildf': raise ValueError('Cannot set type to: '+dtype)


    try:
        ffile = open(filename)
    except:
        raise IOError("Cannot open the file: "+filename)


    formato=endian+dtype
    len_data = np.fromstring(ffile.read(ntype),dtype=endian+"i")[0]
    result = []
    if use_struct :
        while True:
            data = ffile.read(len_data)
            if data=='' : break
            result.append(struct.unpack(dtype,data))
            ffile.read(2*ntype)
    else :
        while True:
            data = ffile.read(len_data)
            if data=='' : break
            result.append(np.fromstring(data,dtype=formato))
            ffile.read(2*ntype)

    ffile.close()

    return np.array(result)



if __name__ == '__main__':
    result = readFortranUnformat("bo.dat",endian='=')
