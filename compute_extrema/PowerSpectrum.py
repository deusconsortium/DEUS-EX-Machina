from pylab import *
import numpy as num
from matplotlib.colors import LogNorm
import sys
import os
import python.readFortranUnformat as FF
import python.Tools as tool
import python.Physics as phy


BINARY_FILE = True
DATA_PATH = 'data/'

if len(sys.argv) == 4:
	boxlen = sys.argv[1]
	npart = sys.argv[2]
	cosmo = sys.argv[3]
else:
	boxlen = raw_input("enter simulation boxlen : ")
	npart = raw_input("enter simulation npart : ")
	cosmo = raw_input("enter cosmo name (XXXXXwX)  : ")


Pfile = 'P_'+boxlen+'_n'+npart+'_'+cosmo+'.txt'

simu = phy.Simu(int(boxlen),int(npart),cosmo)
simu.SpectrumCheck(a0 = simu.af,P_file = Pfile)
