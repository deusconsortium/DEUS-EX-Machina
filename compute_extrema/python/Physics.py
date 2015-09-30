#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
#  Physics.py
#  
#  Copyright 2015 pdefromont <pdefromont@oin>
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

from Tools import*

A_out = 1./(2.*pi**2.)
A_in = 1./(2.*pi)
ZCMB = 1100.
ACMB = 1./(ZCMB + 1.)

def Eta(a,w,Wm0):
	alpha = a*(ZCMB + 1.)
	
	if Wm0 == 1.:
		return -1./3. + alpha/5. + 2./(15.*alpha**(3./2.))
	
	else:
		e0 = Epsilon0(w,Wm0)
		
		eta1 = hyp2f1(1.-1./(3.*w),1.+1./(2.*w),(9.+1./w)/6.,-1./e0)
		eta2_a = hyp2f1(1./2.+1./(3.*w),(w-1.)/(2.*w),3./2.-1./(6.*w),-alpha**(3.*w)/e0)
		eta2_1 = hyp2f1(1./2.+1./(3.*w),(w-1.)/(2.*w),3./2.-1./(6.*w),-1./e0)
		eta3 = hyp2f1(3./2.-1./(2.*w),3./2.+1./(3.*w),5./2.-1./(6.*w),-1./e0)
		eta4_a = hyp2f1(-1./(3.*w),1./(2.*w),(3.+1./w)/6.,-alpha**(3.*w)/e0)
		eta4_1 = hyp2f1(-1./(3.*w),1./(2.*w),(3.+1./w)/6.,-1./e0)
		
		zeta_a = 6.*(1.-9.*w)*alpha**((3.*w-1.)/2.)*eta1*eta2_a + (1.+3.*w)*((6.+3.*w-9.*w**2.)*eta3 + (1.-3.*w)*(1.-9.*w)*e0*eta2_1)*eta4_a
		zeta_1 = 6.*(1.-9.*w)*eta1*eta2_1 + (1.+3.*w)*((6.+3.*w-9.*w**2.)*eta3 + (1.-3.*w)*(1.-9.*w)*e0*eta2_1)*eta4_1
		
		return (zeta_a/zeta_1-1.)/3.
		
def d0_of_d_linear(d,af,w,Wm0):
	f = 1. + d
	A = Eta(af,w,Wm0)
	B = (-27.*A**5.*f**2. - 27.*A**6.*f**2. + sqrt(
  108.*A**9.*f**3. + (-27.*A**5.*f**2. - 27.*A**6.*f**2.)**2.))**(1./3.)
	f0 = 1. + 1./A - 2**(1./3.)/B + B/(3.*2.**(1./3.)*A**3.*f)
	return f0 - 1.

def Epsilon0(w,Wm0):
	if Wm0 == 1.0:
		return 0.0
	return (1.-Wm0)/Wm0*(ZCMB+1.)**(3.*w)


def sqrt_Omega_m_y(w,e0,y):
	if e0 == 0.:
		return 1.
	else:
		return tanh(3.*w/sqrt(2.)*y + arcsinh(1./sqrt(e0)))

def dynamic_y(t,Y,param):
	w,e0,f0 = param
	return [Y[1], -Y[1]/(sqrt(2.)*sqrt_Omega_m_y(w,e0,t)) + Y[0] - f0/Y[0]**2.]

def jac_y(t,Y,param):
	w,e0,f0 = param
	return [[0.,1.],[1.+2.*f0/Y[0]**3.,-1./(sqrt(2.)*sqrt_Omega_m_y(w,e0,t))]]

def y_to_alpha(w,e0,y):
	if e0 == 0.:
		return exp(sqrt(2.)*y)
	else:
		return e0**(1./(3.*w))*(sinh(arcsinh(1./sqrt(e0)) + 3.*w*y/sqrt(2.)))**(2./(3.*w))

def alpha_to_y(w,e0,alpha):
	if e0 == 0.:
		return log(alpha)/sqrt(2.)
	else:
		return sqrt(2.)/(3.*w)*(arcsinh(alpha**(3.*w/2.)/sqrt(e0)) - arcsinh(1./sqrt(e0)))

def dynamic_y_II(v,t,param):
	w,e0,f0 = param
	return [v[1], -v[1]/(sqrt(2.)*sqrt_Omega_m_y(w,e0,t)) + v[0] - f0/v[0]**2.]


def Evolve_single_psi_old(f0,a_tab,w,Wm0):
	if f0 == 1.:
		if size(a_tab) > 1:
			psi = num.empty(size(a_tab))
			for i in range(size(psi)):
				psi[i] = 1.
			return psi
		else:
			return 1.0
	
	e0 = Epsilon0(w,Wm0)
	y0 = [1.,0.]
	t0 = 0.
	
	print('evolving psi until a = '+str(max(a_tab))+' with : \n\te0 = '+str(e0)+'\n\tw = '+str(w)+'\n\tWm0 = '+str(Wm0))
	
	if size(a_tab) > 1:
		psi = num.zeros(size(a_tab))
		psi[0] = 1.0
		k = 0
		af = max(a_tab)
		tf = alpha_to_y(w,e0,af*(ZCMB + 1.))
		
		r = integrate.ode(dynamic_y,jac_y).set_integrator('lsoda',method ='bdf',with_jacobian = True)
		param = [w,e0,f0]
		r.set_initial_value(y0,t0).set_f_params(param).set_jac_params(param)
		
		dy = tf/10000.
		ym = 0.0
		psim = 1.0
		while r.successful() and k < size(a_tab):
			yz = alpha_to_y(w,e0,a_tab[k]*(ZCMB + 1.))
			yp = r.t
			psip = r.y[[0]]
			if yp >= yz:
				psi[k] = psim + (yz - ym)*(psip - psim)/(yp - ym)
				k = k+1
			ym = r.t
			psim = r.y[[0]]
			r.integrate(r.t + dy)
		return psi
	else:
		af = a_tab
		tf = alpha_to_y(w,e0,af*(ZCMB + 1.))
		
		r = integrate.ode(dynamic_y,jac_y).set_integrator('lsoda',method ='bdf',with_jacobian = True)
		param = [w,e0,f0]
		r.set_initial_value(y0,t0).set_f_params(param).set_jac_params(param)
		
		dy = tf/10000.
		ym = 0.0
		psim = 1.0
		psi = 1.0
		done = False
		while r.successful() and done == False:
			yp = r.t
			psip = r.y[[0]]
			if yp >= tf:
				done = True
				psi = psim + (tf - ym)*(psip - psim)/(yp - ym)
			ym = yp
			psim = psip
			r.integrate(r.t + dy)
		
		return psi

def Evolve_single_psi(f0,a_tab,w,Wm0,coeff = 1.0):
	if f0 == 1.:
		if size(a_tab) > 1:
			psi = num.empty(size(a_tab))
			for i in range(size(psi)):
				psi[i] = 1.
			return psi
		else:
			return 1.0
	
	e0 = Epsilon0(w,Wm0)
	y0 = [1.,coeff*(1.-f0)/3.]
	
	alpha_tab = num.zeros(size(a_tab)+1)
	alpha_tab[0] = 1.0
	if size(a_tab)==1:
		alpha_tab[1] = a_tab*(ZCMB + 1.0)
	else:
		for i in range(size(a_tab)):
			alpha_tab[i + 1] = a_tab[i]*(ZCMB + 1.)
	
	#~ print('evolving psi from alpha = '+str(alpha_tab[0])+' to alpha = '+str(alpha_tab[size(alpha_tab) - 1])+' with : \n\te0 = '+str(e0)+'\n\tw = '+str(w)+'\n\tWm0 = '+str(Wm0))
	
	params = [w,e0,f0]
	t = alpha_to_y(w,e0,alpha_tab)
	sol = integrate.odeint(dynamic_y_II,y0,t,args=(params,))
	retour = sol[:,0]
	retour = num.delete(retour,0)
	return retour



def Nu_of_r1(r1_tab,P0,k,factor = 0.5):
	print('\ncomputing <v|r1> ...')
	nu = num.zeros(size(r1_tab))
	for i in range(size(r1_tab)):
		r1 = r1_tab[i]
		if r1> 0.0:
			P = P0*W(factor*k*r1)**2.
			b2 = Beta_r1(r1,P,k)
			g = gamma(P,k)
			s0 = sqrt(sigma2_0(0,P,k))
			x = sqrt((b2**2. + g**2. - 2.*b2*g**2.)/(b2 - g**2.)**2.)
			nu[i] = g*(1.-b2)*K(x)/(g**2.- b2)
		else:
			nu[i] = 0.0
	return nu,nu*s0

def d_of_d0(d0,af,w,Wm0,linear = False):
	if size(d0)==1:
		if linear:
			n = Eta(ad,w,Wm0)
			return (1.+d0)/(1.-d0*n)**3. -1.
		else:
			psi = Evolve_single_psi(1. + d0,af,w,Wm0)
			return (1.+d0)/psi**3. -1.
	else:
		d = num.zeros(size(d0))
		for i in range(size(d0)):
			d[i] = d_of_d0(d0[i],af,w,Wm0,linear)
		return d

def sigma2B(r,n,P,k,R=0.0):
	if size(r)>1:
		sigma2 = num.zeros(r.size)
		for i in range(r.size):
			temp = num.zeros(k.size)
			for l in range(k.size):
				x = A_in*k[l]*r[i]
				if R == 0.:
					temp[l] = k[l]**(2. + 2.*n)*P[l]*W(x)
				else:
					y = A_in*k[l]*R
					temp[l] = k[l]**(2. + 2.*n)*P[l]*W(x)*sinc(y)
				
			sigma2[i] = A_out*Integrate(temp,k)
	else:
		temp = num.zeros(k.size)
		for l in range(k.size):
			x = A_in*k[l]*r
			if R == 0.:
				temp[l] = k[l]**(2. + 2.*n)*P[l]*W(x)
			else:
				y = A_in*k[l]*R
				temp[l] = k[l]**(2. + 2.*n)*P[l]*W(x)*sinc(y)
		sigma2 = A_out*Integrate(temp,k)
	return sigma2

def sigma2_0(n,P,k):
	temp = num.empty(k.size)
	for l in range(k.size):
		temp[l] = P[l]*k[l]**(2.+2.*n)
	
	sigma2 = A_out*Integrate(temp,k,from_zeros = True)
		
	return sigma2
	
def s20_tab(P,k,nmax=2,to_print = True):
	s2_tab = num.zeros(nmax + 1)
	for i in range(nmax + 1):
		s2_tab[i] = sigma2_0(i,P,k)
		if to_print:
			print('- sigma_'+str(i)+'^2(0) \t= '+str(s2_tab[i]))
	return s2_tab

def Beta_r1(r1,P,k):
	if size(r1) > 1:
		beta = num.empty(size(r1))
		for i in range(size(r1)):
			beta[i] = sigma2B(r1[i],0,P,k)*sigma2_0(1,P,k)/(sigma2B(r1[i],1,P,k)*sigma2_0(0,P,k))
	else:
		beta = sigma2B(r1,0,P,k)*sigma2_0(1,P,k)/(sigma2B(r1,1,P,k)*sigma2_0(0,P,k))
	return beta

def gamma(P,k):
	return sigma2_0(1,P,k)/sqrt(sigma2_0(0,P,k)*sigma2_0(2,P,k))

def N_voids_cum_infinity(r1,P0,k,factor = 0.5,delta_c = 0.0,nu_c = 0.):
	N = 10
	power0 = log10(r1)
	power = num.linspace(power0,2.7,N)
	r1_tab = num.zeros(N)
	for i in range(N):
		r1_tab[i] = 10.**power[i]
	nv,ncum,Ninf = N_Voids_r1(r1_tab,P0,k,factor = factor,delta_c = delta_c,nu_c = nu_c,end = True)
	return ncum[size(ncum)-1]

	
def N_Voids_v(v_tab,g,nmax=20):
	print('computing Nv_0(v) with g = '+str(g)+' and '+str(nmax)+' terms in the Cn series')
	Lv = L(v_tab,g,nmax)
	return exp(-v_tab**2./(2.*(1-g**2.)))*Lv/((2.*pi)**2.*sqrt(2.*(1.-g**2.)))

def N_Voids_d_a(af,w,Wm0,g,s0,vwindow=[-10.,0.],Npoints = 100,nmax=200,coeff=1.0):
	print('\n***\ncomputing the evolved Nv(d) profile ...')
	v = num.linspace(vwindow[0],vwindow[1],Npoints)
	d0  = s0*v
	Nv0 = N_Voids_v(v,g,nmax)
	print('evolving profile form z '+str(ZCMB)+' to z = '+str(1./af - 1.))
	if size(af)==1:
		psi = num.zeros(size(d0))
		for i in range(size(psi)):
			psi[i] = Evolve_single_psi(1. + d0[i],af,w,Wm0,coeff)
		d = (1. + d0)/psi**3. -1.
		A = 1./(1. - 3.*D(psi,d0)*(1.0+d0)/psi)	
		return d,A*(1.+d0)/((1.+d)*s0)*Nv0
		#~ return d,Nv0/(D(d,d0)*s0)
	else:
		psi = num.zeros((size(af),size(d0)))
		A = num.zeros((size(af),size(d0)))
		Nv = num.zeros((size(af),size(d0)))
		d = num.zeros((size(af),size(d0)))
		
		p = num.zeros(size(af))
		for i in range(size(d0)):
			p = Evolve_single_psi(1. + d0[i],af,w,Wm0,coeff)
			for a in range(size(af)):
				psi[a,i] = p[a] 
				d[a,i] = (1. + d0[i])/(psi[a,i])**3. -1.
		for a in range(size(af)):			
			A[a] = 1./(1. - 3.*D(psi[a,:],d0)*(1.0 + d0)/psi[a,:])
			Nv[a] = A[a]*(1.+d0)/((1.+d[a,:])*s0)*Nv0	
		return d,Nv

def N_Voids_r1(r1_tab,P0,k,Reff = -1.,factor = 0.5,delta_c = 0.0,nu_c = 0.,end = False):
	nv = num.empty(size(r1_tab))
	
	print('\ncomputing Nv(r1) ...')
	Ncum = num.zeros(size(r1_tab))
	if Reff == -1.:
		for i in range(size(r1_tab)):
			r1 = r1_tab[i]
			if r1 == 0.0:
				nv[i] = 0.
			else:
				P = P0*W(factor*k*r1)**2.
				b2 = Beta_r1(r1,P,k)
				g = gamma(P,k)
				x = sqrt((b2**2. + g**2. - 2.*b2*g**2.)/(b2 - g**2.)**2.)
				y = (sigma2(r1,0,P,k)/sigma2B(r1,0,P,k) - sigma2(r1,1,P,k)/sigma2B(r1,1,P,k))
				s0 = sqrt(sigma2_0(0,P,k))
				z = (g**2. - b2)/((1. - b2)*g*s0)
				if delta_c == 0. and nu_c == 0.:
					Ix = I(x)
					nv[i] = 3.*b2/(2.*r1*pi**2.*sqrt(2.))*g*sqrt(1. - g**2.)/(g**2. - b2)**2.*y*Ix
					#print(x,Ix)
					#print('Nv('+str(r1)+') = '+str(nv[i]))
				else:
					if nu_c > 0.:
						Jx = J(x,z*s0*nu_c)
					else:
						Jx = J(x,delta_c*z)
						
					nv[i] = 3.*b2/(2.*r1*pi**2.*sqrt(2.))*g*sqrt(1. - g**2.)/(g**2. - b2)**2.*y*Jx
					#print('Nv('+str(r1)+') = '+str(nv[i]))
				if nv[i] < 0.0:
					nv[i] = 0.0
				if end == False:
					Ncum = Primitive(nv,r1_tab,from_zeros = True)
				else:
					Ncum = Primitive(nv,r1_tab)
		
		if end == False:
			Ninf = N_voids_cum_infinity(r1_tab[size(r1_tab)-1],P0,k,factor,delta_c,nu_c)
			print('\n'+str(Ninf)+' to add to '+str(Ncum[size(Ncum)-1]))
			Ninf = Ninf + Ncum[size(Ncum)-1]
			
		else: 
			Ninf = 1.0
	else:
		P = P0*W(k*Reff)**2.
		for i in range(size(r1_tab)):
			r1 = r1_tab[i]
			b2 = Beta_r1(r1,P,k)
			g = gamma(P,k)
			x = sqrt(b2**2. + g**2. - 2.*b2*g**2.)/(b2 - g**2.)
			nv[i] = g*(g**2. - 1.)*3.*b2/((2.*pi)**2.*sqrt(2.*(1.-g**2.))*(g**2. - b2)**2.*r1)*(sigma2(r1,0,P,k)/sigma2B(r1,0,P,k) - sigma2(r1,1,P,k)/sigma2B(r1,1,P,k))*I(x)
	print('done')

	return nv,Ncum,Ninf

def sigma2(r,n,P,k,R = 0.0):
	if size(r)>1:
		sigma2 = num.zeros(r.size)
		for i in range(r.size):
			temp = num.zeros(k.size)
			for l in range(k.size):
				x = A_in*k[l]*r[i]
				if R== 0.:
					temp[l] = k[l]**(2. + 2.*n)*P[l]*sinc(x)
				else:
					y = A_in*k[l]*R
					temp[l] = k[l]**(2. + 2.*n)*P[l]*sinc(x)*sinc(y)
				
			sigma2[i] = A_out*Integrate(temp,k)
	else:
		temp = num.zeros(k.size)
		for l in range(k.size):
			x = A_in*k[l]*r
			if R == 0.:
				temp[l] = k[l]**(2. + 2.*n)*P[l]*sinc(x)
			else:
				y = A_in*k[l]*R
				temp[l] = k[l]**(2. + 2.*n)*P[l]*sinc(x)*sinc(y)
			
		sigma2 = A_out*Integrate(temp,k)
	return sigma2


class Simu :
	def __init__(self,boxlen=648,npart=256,cosmo='rpcdmw5',a0=ACMB,LoadFromFile = True):
		self.load(boxlen,npart,cosmo,a0,LoadFromFile)
		
	def load(self,boxlen,npart,cosmo,a0=ACMB,LoadFromFile = True):
		self.boxlen = boxlen
		self.npart = npart
		self.cosmo = cosmo
		self.a0 = a0
		self.af = a0
		if self.cosmo == 'lcdmw5' or self.cosmo == 'lcdmw7':
			self.w = -1.
			self.Wm0 = 0.2573
		if self.cosmo == 'rpcdmw5':
			self.w = -.87
			self.Wm0 = 0.23
		if self.cosmo == 'wcdmw5':
			self.w = -1.2
			self.Wm0 = 0.275
		
		if(LoadFromFile == True):
			data = num.loadtxt('data/info_simu.txt')
			self.af = data[0]
			self.Wm0 = data[1]
			self.npart = int(data[2])
		
		self.loadPK()
		print('simu well loaded.\n\tboxlen = '+str(self.boxlen)+'\n\tnpart = '+str(self.npart)+'\n\tcosmo = '+str(self.cosmo)+'\n\ta0 = '+str(self.a0)+'\n\taf = '+str(self.af)+'\n\tWm0 = '+str(self.Wm0)+'\n\tw = '+str(self.w)+'\n\tRf = '+str(self.Rf) +' [Mpc/h]'+'\n\ts0 = '+str(self.s0)+'\n\tg = '+str(self.g)+'\n\n')
			
	def loadPK(self):
		print('\n***\nloading spectrum at z = '+str(1./self.a0 -1.))
		
		P_file = "/efiler2/bingo_save/Babel/data/pk_"+self.cosmo+'.dat'
		D_a_file = "/efiler2/bingo_save/Babel/data/mpgrafic_input_"+self.cosmo+".dat"
			
		try:
			k,P = num.loadtxt(P_file, usecols=(0, 1), unpack=True)
			a,D,dD = num.loadtxt(D_a_file, usecols=(0, 2,3), unpack=True)
		except IOError:
			print("!!! ERROR : no file " + P_file +" or "+D_a_file+" !!!")
		
		self.D0 = Solve(D,a,self.a0)
		self.Df = Solve(D,a,self.af)
		self.D1 = Solve(D,a,1.0)
		self.Coeff = Solve(dD,a,self.a0)

		P = P*(self.D0/self.D1)**2. 
		print('P(k,z) well loaded')
		
		#linearisation de P,k
		print('linearising P(k)')
		kmin = 2.*pi/(float(self.boxlen))
		kmax = float(self.npart)*kmin/2.
		ku = num.linspace(kmin,kmax, self.npart)
		Pu = num.zeros(size(ku))
		for i in range(size(ku)):
			Pu[i] = Solve(P,k,ku[i])		
		
		Rf = num.loadtxt('data/info.txt')

		Rbox = float(self.boxlen)
		Rf = Rf*Rbox
		self.Rf = Rf
		Rcell = Rbox/float(self.npart)
		self.k = ku
		self.cic_correction = 0.86
		self.P = Pu*W(self.cic_correction*self.k*Rcell)**2.*exp(-(self.k*Rf)**2.)

		self.s0 = sqrt(sigma2_0(0,self.P,self.k))
		self.g = gamma(self.P,self.k)
		print('***\n')
	
	def SpectrumCheck(self,a0 = ACMB,P_file = 'powergridcicfof_00001.txt',cic_correction = -1.):
		print('\n***\nspectrum cheking at z = '+str(1./a0 - 1.)+'\n')
			
		if(cic_correction != -1.):
			self.P = self.P/W(self.cic_correction*self.k*Rcell)**2.
			self.cic_correction == cic_correction
			self.P = self.P*W(self.cic_correction*self.k*Rcell)**2.
			
		D_a_file = "/efiler2/bingo_save/Babel/data/mpgrafic_input_"+self.cosmo+".dat"
			
		try:
			k,P = num.loadtxt(P_file, usecols=(0, 1), unpack=True)
			a,D,dD = num.loadtxt(D_a_file, usecols=(0, 2,3), unpack=True)
		except IOError:
			print('!!ERROR !! No non-linear spectrum founded')
		
		D1 = Solve(D,a,a0)
		
		figure(1)
		grid(True)
		xscale('log')
		xlabel('k')
		ylabel('$k^6P(k)$')
		yscale('log')
		title('$k^6P_{lin}(k)\\times W(\\alpha kR_{cell})$ with $\\alpha = $'+str(self.cic_correction))
		plot(k,k**6.*P,'b+',label = 'non-linear')
		plot(self.k,self.k**6.*self.P*(D1/self.D0)**2.*exp((self.k*self.Rf)**2.),'r-',label='linear')
		g0 = gamma(P,k)
		s0 = sqrt(sigma2_0(0,P,k))
		print('spectrum comparison : gDATA = '+str(g0)+' vs gTHEO = '+str(self.g))
		print('spectrum comparison : s0DATA = '+str(s0)+' vs s0THEO = '+str(self.Df/self.D0*self.s0)+'\n***\n')
		legend(loc = 2,title = 'z = '+str(1./self.af- 1.))
		show()
	
	def P_Voids_d0(self,use_af = False,vwindow=[-5,0.],Npoints = 100,nmax=200):
		v = num.linspace(vwindow[0],vwindow[1],Npoints)
		if use_af == True:
			s0 = self.Df/self.D0*self.s0
		else:
			s0 = self.s0
		d0  = s0*v
		Nv0 = N_Voids_v(v,self.g,nmax)
		return d0,(40.*sqrt(5.)*pi**(3./2.))/(29. - 6.*sqrt(6.))*Nv0/s0

	def P_Voids_d_a(self,af = -1.,vwindow=[-5.,0.],Npoints = 50,nmax=200):
		if af == -1.:
			af = self.af
		d,N = N_Voids_d_a(af,self.w,self.Wm0,self.g,self.s0,vwindow,Npoints,nmax,self.Coeff)
		return d,(40.*sqrt(5.)*pi**(3./2.))/(29. - 6.*sqrt(6.))*N
