from pylab import *
import numpy as num
import sys
import os

from scipy.special import erfc
from scipy.special import erf
from scipy.special import hyp2f1

def f(x):
	return (x**3.-3.*x)*(erf(sqrt(5./2.)*x)+erf(sqrt(5./2.)*x/2.))/2. + sqrt(2./(5.*pi))*((31.*x**2./4.+8./5.)*exp(-5.*x**2./8.)+(x**2./2.-8./5.)*exp(-5.*x**2./2.))

def derivative(x,y):
	dy = num.zeros(size(y))
	
	for i in range(size(x) - 2):
		dy[i + 1] = (y[i + 2] - y[i])/(x[i + 2] - x[i])
	
	dy[0] = (y[1] - y[0])/(x[1]-x[0])
	dy[size(dy)-1] = (y[size(y) - 1] -  y[size(y) - 2])/(x[size(y) - 1] - x[size(y) - 2])
	return dy

def solve(x,y,y0):
		if size(x) == size(y):
			if y[0] > y0:
				for i in range(size(y) - 1):
					if y[i] >= y0 and y[i+1] <= y0:
						return x[i + 1] + (y[i+1] - y0)*(x[i] - x[i+1])/(y[i+1] - y[i])
			else:
				for i in range(size(y) - 1):
					if y[i] <= y0 and y[i+1] >= y0:
						return x[i + 1] + (y[i+1] - y0)*(x[i] - x[i+1])/(y[i+1] - y[i])
		return None


def Wth(x):
	return 3.*(num.sin(x) - x*num.cos(x))/(math.pow(x,3.0))

def integrate(x,y,imax = -1,from_zeros = True):
	Int = 0.0
	if imax == -1:
		imax = size(x)

	if from_zeros == True and x[0] > 0.0:
		Int += 0.5*y[0]*x[0]
		
	for i in range(imax - 1):
		Int += 0.5*(y[i+1] + y[i])*(x[i+1] - x[i])
	return Int

#physical functions

def Epsilon0(w,Wm0,zcmb):
	if Wm0 == 1.0:
		return 0.0
	return (1.-Wm0)/Wm0*(zcmb+1.)**(3.*w)

def Eta(a,w,Wm0,d_logD_d_loga,zcmb):
	s0 = d_logD_d_loga
	alpha = a*(zcmb + 1.)
	
	if Wm0 == 1.:
		return -1./3. + alpha*(1.+2.*s0/3.)/5. + 2.*(1.-s0)/(15.*alpha**(3./2.))
	
	else:
		e0 = Epsilon0(w,Wm0,zcmb)
		
		eta1 = hyp2f1(1.-1./(3.*w),1.+1./(2.*w),(9.+1./w)/6.,-1./e0)
		eta2_a = hyp2f1(1./2.+1./(3.*w),(w-1.)/(2.*w),3./2.-1./(6.*w),-math.pow(alpha,3.*w)/e0)
		eta2_1 = hyp2f1(1./2.+1./(3.*w),(w-1.)/(2.*w),3./2.-1./(6.*w),-1./e0)
		eta3 = hyp2f1(3./2.-1./(2.*w),3./2.+1./(3.*w),5./2.-1./(6.*w),-1./e0)
		eta4_a = hyp2f1(-1./(3.*w),1./(2.*w),(3.+1./w)/6.,-math.pow(alpha,3.*w)/e0)
		eta4_1 = hyp2f1(-1./(3.*w),1./(2.*w),(3.+1./w)/6.,-1./e0)
		
		zeta_a = math.pow(alpha,(3.*w-1.)/2.)*eta2_a*(9.*w - 1.)*(6.*eta1 - 2.*e0*s0*eta4_1*(3.*w+1.)) + 3.*(1.+3.*w)*(w-1.)*(2.+3.*w)*eta3*eta4_a + e0*(1.+3.*w)*(1.+2.*s0-3.*w)*(9.*w-1.)*eta2_1*eta4_a
		zeta_1 = eta2_1*(9.*w - 1.)*(6.*eta1 - 2.*e0*s0*eta4_1*(3.*w+1.)) + 3.*(1.+3.*w)*(w-1.)*(2.+3.*w)*eta3*eta4_1 + e0*(1.+3.*w)*(1.+2.*s0-3.*w)*(9.*w-1.)*eta2_1*eta4_1
		
		return (zeta_a/zeta_1-1.)/3.

class SplineFit:
	def __init__(self,x,y):
		if size(x) is not size(y):
			print 'splineFit : x and y have not the same size !'
		else:
			self._x = x
			
			p = num.zeros(size(x))
			p[0] = (y[1]-y[0])/(x[1]-x[0])
			p[size(p) - 1] = (y[size(y)-1] - y[size(y)-2])/(x[size(x)-1]-x[size(x)-2])
			for i in range(size(p) - 2):
				p[i+1]=(y[i+2]-y[i])/(x[i+2]-x[i])
			self._c = num.zeros(((size(x) - 1),4))
			for i in range(size(x) - 1):
				self._c[i][0] = y[i]
				self._c[i][1] = p[i]
				self._c[i][2] = 3.0*(y[i+1]-y[i]) - (p[i+1]+2.*p[i])*(x[i+1]-x[i])
				self._c[i][3] = -2.0*(y[i+1]-y[i]) + (p[i] + p[i+1])*(x[i+1]-x[i])
	def __call__(self, x):
		if size(x) > 1:
			y = num.zeros(size(x))
			for i in range(size(x)):
				y[i] = self.__call__(x[i])
			return y
		else:
			index = -1
			for i in range(size(self._x) - 1):
				if self._x[i] <= x and self._x[i+1] >= x:
					index = i
					break
			if index >= 0:
				X = (x-self._x[index])/(self._x[index+1] - self._x[index])
				return self._c[index][0] + self._c[index][1]*(x - self._x[index]) + self._c[index][2]*math.pow(X,2) + self._c[index][3]*math.pow(X,3)
			else:
				return None
			
#integration scheme

def alpha_to_y(w,e0,alpha):
	if e0 == 0.:
		return num.log(alpha)/num.sqrt(2.)
	else:
		return num.sqrt(2.)/(3.*w)*(num.arcsinh(math.pow(alpha,3.*w/2.)/num.sqrt(e0))) - num.arcsinh(1./num.sqrt(e0))

def y_to_alpha(w,e0,y):
	if e0 == 0.:
		return num.exp(num.sqrt(2.)*y)
	else:
		return math.pow(e0,1./(3.*w))*math.pow(num.sinh(num.arcsinh(1./num.sqrt(e0)) + 3.*w*y/num.sqrt(2.)),2./(3.*w))

def sqrt_Omega_m_y(w,e0,y):
	if e0 == 0.:
		return 1.
	else:
		return num.tanh(3.*w/num.sqrt(2.)*y + num.arcsinh(1./num.sqrt(e0)))

def dynamic_y(v,t,param):
	w,e0,f0 = param
	return [v[1], -v[1]/(num.sqrt(2.)*sqrt_Omega_m_y(w,e0,t)) + v[0] - f0*math.pow(v[0],-2)]

def Evolve_single_psi(r0,f0,a_tab,w,Wm0,zcmb,s0):
	if f0 == 1.:
		if size(a_tab) > 1:
			psi = num.empty(size(a_tab))
			for i in range(size(psi)):
				psi[i] = 1.
			return psi
		else:
			return 1.0
	
	e0 = Epsilon0(w,Wm0,zcmb)
	y0 = [1.,s0*(1. - f0)/3.]
	alpha = a_tab*(zcmb + 1.)
	alpha_tab = num.zeros(size(a_tab)+1)
	alpha_tab[0] = 1.0
	if size(a_tab)==1:
		alpha_tab[1] = a_tab*(zcmb + 1.0)
	else:
		for i in range(size(a_tab)):
			alpha_tab[i + 1] = a_tab[i]*(zcmb + 1.)
		
	params = [w,e0,f0]
	t = alpha_to_y(w,e0,alpha_tab)
	sol = integrate.odeint(dynamic_y,y0,t,args=(params,))
	retour,v = sol.T
	retour = num.delete(retour,0)
	v = num.delete(v,0)
	
	a = alpha_tab[size(alpha_tab)-1]
	return retour,r0*num.sqrt(Wm0/2.)*(zcmb + 1.0)**(1./2.)/(100000.*math.pow(a,3.*w/2.-1.)*C0)*v

def evolveProfile(r0,f0,af,w,Wm0,zcmb,dlogD_dloga_cmb):
	s0 = dlogD_dloga_cmb
	if size(af) == 1:
		r = num.zeros(size(f0))
		f = num.zeros(size(f0))
		v = num.zeros(size(f0))
		for i in range(size(f0)):
			psi,v[i] = Evolve_single_psi(r0[i],f0[i],af,w,Wm0,zcmb,s0)
			r[i] = r0[i]*psi
			if psi > 0. and psi < 1.0e+5:
				f[i] = f0[i]/(psi**3.)
			else:
				f[i] = None
	else:
		r = num.zeros((size(af),size(r0)))
		f = num.zeros((size(af),size(r0)))
		v = num.zeros((size(af),size(r0)))
		for i in range(size(r0)):
			psi,v[i] = Evolve_single_psi(r0[i],f0[i],af,w,Wm0,zcmb,s0)
			for j in range(size(af)):
				r[j][i] = r0[i]*psi[j]
				if psi[j] > 0. and psi[j] < 1.0e+5:
					f[j][i] = f0[i]/(psi[j]**3.)
				else:
					f[j][i] = None
	return r,f,v
