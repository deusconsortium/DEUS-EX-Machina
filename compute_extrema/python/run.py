#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
#  sans titre.py
#  
#  Copyright 2015 pdefromont <pdefromont@advil.obspm.fr>
#  
#  This program is free software; you can redistribute it and/or modify
#  it under the terms of the GNU General Public License as published by
#  the Free Software Foundation; either version 2 of the License, or
#  (at your option) any later version.
#  
#  This program is distributed in the hope that it will be useful,
#  but WITHOUT ANY WARRANTY; without even the implied warranty of
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#  GNU General Public License for more details.
#  
#  You should have received a copy of the GNU General Public License
#  along with this program; if not, write to the Free Software
#  Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston,
#  MA 02110-1301, USA.
#  
#  

import sys
import os

#~ 
#~ if len(sys.argv) == 4:
	#~ boxlen = sys.argv[1]
	#~ npart = sys.argv[2]
	#~ cosmo = sys.argv[3]
#~ 
#~ else:
	#~ boxlen = raw_input("enter simulation boxlen : ")
	#~ npart = raw_input("enter simulation npart : ")
	#~ cosmo = raw_input("enter cosmo name (XXXXXwX)  : ")

def char_title(n):
	r = ''
	if n>1000:
		r='0'
	elif n>100:
		r='00'
	elif n>10:
		r='000'
	else:
		r='0000'
	return r+str(n)
	
boxlen = 648
npart = 512
output = 13
cosmo = 'rpcdmw5'

####

raffinement = 2
nfiles = 64
nb_cores = 8
####



def writefile():
	f = open('../namelist/run_extremes','w')
	f.write('cat > parameters_'+char_title(output)+'.nml <<EOF\n')
	f.write('&parameters\n')
	f.write('nfact='+str(raffinement)+'\n')
	f.write('cic=2\n')
	f.write('h0=72.\n')
	f.write('read_alafof=.true.\n')
	f.write('ramses_read_part=2 \n')
	f.write("rep='/data_bingo/Babel/boxlen"+str(boxlen)+"_n"+str(npart)+"_"+str(cosmo)+"/cube_"+char_title(output)+"/'\n")
	f.write("infofile='info_"+char_title(output)+".txt'\n")
	f.write("inputfiles='fof_bboxlen"+str(boxlen)+"_n"+str(npart)+"_"+str(cosmo)+"_cube_'\n")
	f.write("nchar='"+char_title(output)+"'\n")
	f.write('nfiles='+str(nfiles)+'\n')
	f.write('iogroupsize=0\n')
	f.write("outputfile='run_extremes_"+char_title(output)+".txt'\n")
	f.write('/\n')
	f.write('EOF\n')
	f.write('\n#RUN\n')
	f.write('mpirun -np '+str(nb_cores)+' ~/code/powergrid_challenge/powergrid parameters_'+char_title(output)+'.nml\n')	
	f.close()

writefile()
