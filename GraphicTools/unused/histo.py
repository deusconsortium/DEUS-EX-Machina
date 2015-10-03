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

def Legend(to_add):
	legend_tab.append(to_add)

def Show(legend_pos = 1,title = ''):
	global legend_tab
	if size(legend_tab) > 0:
		legend(legend_tab,legend_pos,title = title)
	show()
	del legend_tab[:]

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

def file_number():
	i=0
	while True:
		fname = DATA_PATH + 'minimum_'+char_title(i)
		if BINARY_FILE:
			fname+= '.dat'
		else:
			fname+= '.txt'
		if not os.path.isfile(fname):
			break
		else:
			i+=1
	return i

def Load(number = -1):
	if number == -1:
		n = file_number()
		print('founded '+str(n)+' files')
		X,Y,Z,D,S,Sm = Load(0)
		for i in range(n-1):
			x,y,z,d,s,sm = Load(i+1)
			X= num.concatenate([X,x])
			Y = num.concatenate([Y,y])
			Z = num.concatenate([Z,z])
			D = num.concatenate([D,d])
			S = num.concatenate([S,s])
			Sm = num.concatenate([Sm,sm])
		print('loading done')
		return X,Y,Z,D,S,Sm
	elif number >=0 and number < file_number():
		name =  DATA_PATH + 'minimum_'+char_title(number) 
		print('loading file '+name)
		if BINARY_FILE:
			result = FF.readFortranUnformat(name+'.dat',dtype='f')
		else:
			result = num.loadtxt(name+'.txt', usecols=(0,1,2,3,4,5), unpack=False)
		return result[:,0],result[:,1],result[:,2],result[:,3],result[:,4],result[:,5]

def Info(data):
	X,Y,Z,D,S,Sm = data
	xmin = min(X)
	xmax = max(Y)
	ymin = min(Y)
	ymax = max(Y)
	zmin = min(Z)
	zmax = max(Z)
	pV = 0
	dm = 0.0
	for i in range(size(D)):
		dm+= D[i]
		if D[i] == 0.0:
			pV+=1
	dm = dm/size(D)
	pV = 100.*pV/size(D)
	print('---info---\n- x in ['+str(xmin)+', '+str(xmax)+']\n- y in ['+str(ymin)+', '+str(ymax)+']\n- z in ['+str(zmin)+', '+str(zmax)+']\n- empty cells : '+str(round(pV,4))+' % of cells\n- mean density 1+d = '+str(round(dm,4))+'\n----------')
	
def mask(data,tab,seuil = 0.,seuil_moyen = 0.,dmax = 1000.):
	print('masking data such as seuil > '+str(seuil)+' and d+1 <= '+str(dmax))
	x,t,z,d,s,sm = data
	mask_seuil = (s>=seuil)
	mask_seuil_m = (sm>=seuil_moyen)
	mask_minimum = (d<=dmax)
	mask1 = num.logical_and(mask_seuil,mask_minimum)
	mask = num.logical_and(mask1,mask_seuil_m)
	tab_masque = tab[mask]	
	print('resulting data is '+str(100.*float(size(tab_masque))/float(size(tab)))+' % of the total data')
	return tab_masque

def histo(data,tab,seuil = [],bins = 100,window = [0.0,1.0],dmax = 1000.,normalized = True,doplot = True,no_mask = False):
	if doplot:
		figure(1)
		grid(True)
	if no_mask == True:
		print('treating data without mask')
		y0,x0 = num.histogram(tab,bins,range = (window[0],window[1]))
		x = num.empty(size(y0))
		y = num.empty(size(y0))
		Dx = (max(x0) - min(x0))/(max(tab) - min(tab))
		for i in range(size(y)):
			x[i] = 0.5*(x0[i+1]+x0[i])
			if normalized:
				y[i] = float(y0[i])/(float(size(tab))*(x0[i+1]-x0[i]))
			else:
				y[i] = float(y0[i])
		if doplot:
			plot(x - 1.,y,linestyle = '-', marker = '')
	elif seuil == []:
		D = mask(data,tab,dmax=dmax)
		y0,x0 = num.histogram(D,bins,range = (window[0],window[1]))
		x = num.empty(size(y0))
		y = num.empty(size(y0))
		Dx = (max(x0) - min(x0))/(max(D) - min(D))
		for i in range(size(y)):
			x[i] = 0.5*(x0[i+1]+x0[i])
			if normalized:
				y[i] = float(y0[i])/(float(size(D))*(x0[i+1]-x0[i]))
			else:
				y[i] = float(y0[i])
		if doplot:
			plot(x - 1.,y,linestyle = '-', marker = '')
	else:
		x = num.empty((size(seuil),bins))
		y = num.empty((size(seuil),bins))
		for s in range(size(seuil)):
			D = mask(data,tab,seuil = seuil[s],dmax=dmax)
			y0,x0 = num.histogram(D,bins,range = (window[0],window[1]))
			Dx = (max(x0) - min(x0))/(max(D) - min(D))
			for i in range(size(y0)):
				x[s,i] = 0.5*(x0[i+1]+x0[i])
				if normalized:
					y[s,i] = float(y0[i])/(float(size(D))*(x0[i+1]-x0[i]))
				else:
					y[s,i]= float(y0[i])
			N = float(size(seuil))
			if doplot:
				if N > 1:
					plot(x[s] - 1.,y[s],linestyle = '-', marker = '',color = [float(s)/(N-1),0.0,1.0-float(s)/(N-1)])
				else:
					plot(x[s] - 1.,y[s],linestyle = '-', marker = '',color = 'b')
				Legend(str(round(seuil[s],2)))
	if doplot:
		Show(title = 'seuil')
	return x-1.0,y

def clear0(tab):
	mini = num.min(tab[num.nonzero(tab)])
	tab[tab == 0]=mini
	return tab

legend_tab = list()
X,Y,Z,D,S,Sm = Load(0)
data = [X,Y,Z,D,S,Sm]
#~ seuil = num.linspace(0.0,0.05,10)
x,y = histo(data,D,doplot = False,normalized = True,window = [0.9,1.1],no_mask = False)
#~ Info(data)

"""
theory
"""

theo = phy.Simu(5184,2048,'lcdmw5')
theo.SpectrumCheck(P_file = 'P_5184_n2048_lcdmw5.txt')
d,N = theo.P_Voids_d_a()
d0,N0 = theo.P_Voids_d0(use_af = True)

figure(1)
grid(True)
title('Minimum density for z = '+str(1./theo.af - 1.))
plot(x,y,'bo',label = 'data')
plot(d,N,'r-',label = 'evolved')
plot(d0,N0,'r--',label = 'linear')
xlabel('$\\delta$')
ylabel('$N_v(\\delta)$')
legend(loc = 2)
show()


"""
======== GRAPHES DE Nv(z) EVOLUTION VS LINEAR ==========
"""
#~ 
#~ s1 = phy.Simu(648,256,'rpcdmw5')
#~ a = num.logspace(-3.,-2.,5)
#~ figure(1)
#~ grid(True)
#~ d0,N0 = s1.P_Voids_d0()
#~ plot(d0,N0,'bo')
#~ d1,N1 = s1.P_Voids_d_a(a)
#~ for i in range(size(a)):
	#~ s2 = phy.Simu(648,256,'rpcdmw5',a0 = a[i])
	#~ d2,N2 = s2.P_Voids_d0()
	#~ plot(d2,N2,'--',color = tool.RBColor(i,size(a)))
	#~ plot(d1[i],N1[i],'-',color = tool.RBColor(i,size(a)),label = str(1./a[i]-1.))
	#~ legend(loc = 2,title = 'z')
#~ show()


#~ simu = phy.Simu(648,256,'rpcdmw5',a0 = 0.0121622)
#~ simu = phy.Simu(648,256,'rpcdmw5')
#~ simu.SpectrumCheck()
#~ s0 = simu.s0
#~ simu0 = phy.Simu(648,256,'rpcdmw5',a0 = 0.0121622)
#~ d0,N0 = simu0.P_Voids_d0(nmax = 200)
#~ d,N = simu.P_Voids_d_a(af=0.0121622)
#~ figure(1)
#~ grid(True)
#~ plot(x,y,'bo')
#~ Legend('data z = 81')
#~ plot(d,N,'r-')
#~ plot(d0,N0,'r--')
#~ plot(x,exp(-0.5*(x/s0)**2.)/sqrt(2.*pi*s0**2.),'r-')
#~ Show()



#~ simu = phy.Simu(648,1024,'rpcdmw5')
#~ d,N = simu.P_Voids_d_a(0.301,nmax = 100)
#~ d,N = simu.P_Voids_d0(nmax = 100)
#~ Legend('theory')
#~ 
#~ xlim(-1.0,0.0)
#~ plot(d,N/max(N),'go')
#~ plot(x,y/max(y),'b-')
#~ show()


"""
======== GRAPHES DE D(a) ==========
"""

#~ D_a_file = "/efiler2/bingo_save/Babel/data/mpgrafic_input_lcdmw5.dat"
#~ a,dD,D = num.loadtxt(D_a_file, usecols=(0,3,2), unpack=True)
#~ 
#~ for i in range(size(a)):
	#~ if a[i] > 0.01:
		#~ n = i
		#~ break
#~ 
#~ a2 = num.empty(size(a) - n)
#~ for i in range(size(a)-n):
	#~ a2[i] = a[i+n]
#~ 
#~ w =-1.
#~ Wm0 = 0.2573
#~ 
#~ d0 = -0.0001
#~ phi = phy.Evolve_single_psi(1.+d0,a2,w,Wm0)
#~ Dt = 1.0 +3.*(1.-phi)/d0
#~ figure(1)
#~ grid(True)
#~ xscale('log')
#~ plot(a,D/max(D),'b--')
#~ 
#~ plot(a2,Dt/max(Dt),'r-')
#~ show()
