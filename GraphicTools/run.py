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

from DEUSProfile import*

test = DEUSGraphics(1100.)
r1_tab = [15.,20.,50.,100.]
colors = ['r','b','g','k']

test.setDoPlot(False)

figure(1)
grid(True)

#Mp constant

xlabel('r in [Mpc/h]')
ylabel('$f(r)$')
title("simu1 = 5184_n2048 and simu2 = 2592_n1024") 
for i in range(size(r1_tab)):
	test.Load('work_5184_2048_output14')
	r0,f0,df,v,dv = test.getMeanProfile(r1_tab[i])
	#if f is not None:
		#errorplot(r,f,df,color = colors[i],label = str(solve(r,f,1.0)))
	plot(r0,f0,marker='o',color = colors[i],label = "simu1 "+str(r1_tab[i]))

	test.Load('work_2592_1024_output24')
	r,f,df,v,dv = test.getMeanProfile(r1_tab[i])
	#if f is not None:
		#errorplot(r,f,df,marker = '*',linestyle = '--',color = colors[i],label = str(solve(r,f,1.0)))
	x,y = fraction(r0,f0,r,f)
	plot(r,f,marker = '*',linestyle = '-',color = colors[i],label = "simu2 "+str(r1_tab[i]))

"""
#varier Mp
r1 = 20.
xlabel('r in [Mpc/h]')
ylabel('$\\Delta(r)+1$')
test.Load('work_2592_1024_output24')
r0,f0,df,v,dv = test.getMeanProfile(r1)
#plot(x,y,marker='+',color = 'r',label ="2592_1024")
plot(r0,f0,marker='+',color = 'r',label ="2592_1024")
plot(r0,1.0+0.*r0,'k')

test.Load('work_648_1024_output24')
r,f,df,v,dv = test.getMeanProfile(r1)
fraction
x,y = fraction(r0,f0 - 1.0,r,f - 1.0,20)
#~ plot(r0,f0,marker='+',color = 'b',label ="648_1024")
plot(r,f,marker='+',color = 'b',label ="648_1024")

test.Load('work_2592_2048_output60')
r,f,df,v,dv = test.getMeanProfile(r1)
x,y = fraction(r0,f0 - 1.0,r,f - 1.0,20)
plot(r,f,marker='+',color = 'g',label ="2592_2048")
#plot(r,f,marker='+',color = 'g',label ="2592_2048")
"""

legend()
show()


