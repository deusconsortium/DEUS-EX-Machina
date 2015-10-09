from pylab import *
import numpy as num
import sys
import os

from scipy import integrate as Nintegrate
from scipy.special import erfc
from scipy.special import erf
from scipy.special import hyp2f1


#tedious functions

def f(x):
	return (x**3.-3.*x)*(erf(sqrt(5./2.)*x)+erf(sqrt(5./2.)*x/2.))/2. + sqrt(2./(5.*pi))*((31.*x**2./4.+8./5.)*exp(-5.*x**2./8.)+(x**2./2.-8./5.)*exp(-5.*x**2./2.))

def I_dn_dr1(x):
	return -(-(8.*x**4.+2.*x**2.+5.0)/(x**4.*(5.+4*x**2.)**2.)*9.*sqrt(5./(2.*pi)) + 3.*(x**2. - 1.)/(x**5.*sqrt(2.*pi))*(arctan(sqrt(5.)/x) + arctan(sqrt(5.)/(2.*x))))

def Wg(a,b):
	return exp(-(a+b)**2./2.)-exp(-(a-b)**2./2.) + b/2.*sqrt(2.*num.pi)*(erf((a-b)/sqrt(2.)) + erf((a+b)/sqrt(2.)))

#useful functions

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
	W =  3.*(num.sin(x) - x*num.cos(x))/(num.power(x,3))
	nan_pos = num.isnan(W)
	W[nan_pos] = 1.0 
	return W

def integrate(x,y,imax = -1,from_zeros = True,method = 'quadric'):
	Int = 0.0
	if imax == -1:
		imax = size(x)

	if from_zeros == True and x[0] > 0.0:
		Int += 0.5*y[0]*x[0]
	
	if method == 'spline':
		fit = SplineFit(x,y)
		Int += fit.Integrate(x[imax - 1])
	else:
		for i in range(imax - 1):
			Int += 0.5*(y[i+1] + y[i])*(x[i+1] - x[i])
	
	return Int
	
	
def getPowerLawFitCoefficients(x,y,order = 2):
	S = num.empty([order + 1,order + 1])
	for i in range(order + 1):
		for j in range(order + 1):
			if i == 0 and j == 0:
				S[0,0] = size(x)
			else:
				S[i,j] = num.sum(x**(i+j)) 	
	P = num.empty(order + 1)
	for i in range(order + 1):
		P[i] = num.sum(y*(x**i))
	A = dot(inv(S),P)
	return A

def gaussianDensitySmoothing(r,f0,R,end_point_to_remove = 5,r1_factor = 2):
	if R == 0.0 :
		return f0
	
	print 'GAUSSIAN smoothing on '+str(R)+' [Mpc/h]'
	#print 'extending f0 by asymptotic behavior'
	
	#continuing artificially f0 by asymptotic fit
	r1 = solve(r,f0,1.0)
	index = 0
	size_0 = size(r)
	for i in range(size(r)):
		if r[i] >= r1:
			index = i
			break
	quality = size(r) - 1 - end_point_to_remove - r1_factor*index
	if quality < 3:
		print 'warning: not enough point to get AsymptoticBehavior'
	else:
		xu = num.zeros(quality)
		yu = num.zeros(quality)
		for i in range(quality):
			xu[i] = r[r1_factor*index + i]
			yu[i] = f0[r1_factor*index + i]
		A = getPowerLawFitCoefficients(log(xu), log(num.fabs(1. - yu)),1)
		r_big = num.linspace(r[0],2*(r[size(r)-1]-r[0]),2*size(r) - 1)
		f_big = num.zeros(size(r_big))
		for i in range(size(r_big)):
			if i < size(r) - end_point_to_remove:
				f_big[i] = f0[i]
			else:
				if yu[0] > 1.:
					f_big[i] = 1.0 + exp(A[0] + A[1]*log(r_big[i]))
				else:
					f_big[i] = 1.0 - exp(A[0] + A[1]*log(r_big[i]))
		#print 'asymptotic behaviour founded such as |1-f| ~ r^'+str(A[1])
		f0 = f_big
		r = r_big
	
	x = r
	fg = num.zeros(size(x))
	d0 = f0 -r/3.*derivative(r,f0)
	
	for i in range(size(x)):
		fg[i] = 3.*R/(sqrt(2*num.pi)*r[i]**3.)*integrate(x,x*d0*Wg(r[i]/R,x/R))
	
	fr = num.zeros(size_0)
	for i in range(size_0):
		fr[i] = fg[i]
	return fr

#physical functions

def getDensity(r,f):
	d = num.zeros(size(r))
	n = size(r)-1
	d[n]=(f[n]*r[n]**3. - f[n-1]*r[n-1]**3)/(r[n]**3. - r[n-1]**3)
	d[n-1]=d[n]
	for i in range(n-1):
		x = r[n-1-i]/r[n-2-i]
		d[n-2-i] = -d[n-1-i]*(1.+x*(2.+3.*x))/(3.+x*(2.+x)) + 4.*(x**3.*f[n-1-i]-f[n-2-i])/(x-1.)
	return d

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
		if size(x) == size(y):
			self._x = x
	
			p = num.zeros(size(x))
			p[0] = (y[1]-y[0])/(x[1]-x[0])
			p[size(p) - 1] = (y[size(y)-1] - y[size(y)-2])/(x[size(x)-1]-x[size(x)-2])
			for i in range(size(p) - 2):
				p[i+1]=(y[i+2]-y[i])/(x[i+2]-x[i])
			self._c = num.zeros(((size(x) - 1),4))
			for i in range(size(x) - 1):
				self._c[i][0] = y[i]
				self._c[i][1] = p[i]*(x[i+1] - x[i])
				self._c[i][2] = 3.0*(y[i+1]-y[i]) - (p[i+1]+2.*p[i])*(x[i+1]-x[i])
				self._c[i][3] = -2.0*(y[i+1]-y[i]) + (p[i] + p[i+1])*(x[i+1]-x[i])
		else:
			print 'splineFit : x and y have not the same size ! : '+str(size(x))+' vs '+str(size(y))
			
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
				return self._c[index][0] + self._c[index][1]*X + self._c[index][2]*X**2. + self._c[index][3]*X**3.
			else:
				return None
			
	def Integrate(self,x = None):
		if x is None:
			I = 0.0
			for i in range(size(self._x)):
				D = self._x[i+1] - self._x[i]
				I += D*(self._c[i][0] + self._c[i][1]/2. + self._c[i][2]/3. + self._c[i][3]/4.)
			return I
		else:
			if size(x) > 1:
				I = num.zeros(size(x))
				for i in range(size(x)):
					I[i] = self.Integrate(x[i])
				return I
			else:
				index = -1
				for i in range(size(self._x) - 1):
					if self._x[i] <= x and self._x[i+1] >= x:
						index = i
						break
				if index >= 0:
					I = 0.0
					for i in range(index):
						D = self._x[i+1]-self._x[i]
						I += D*(self._c[i][0] + self._c[i][1]/2. + self._c[i][2]/3. + self._c[i][3]/4.)
					X = (x - self._x[index])/(self._x[index+1]-self._x[index])
					D = self._x[index+1]-self._x[index]
					I += D*(self._c[index][0]*X + self._c[index][1]*X**2./2. + self._c[index][2]*X**3./3. + self._c[index][3]*X**4./4.)
					return I
				else:
					return self.Integrate(self)
#plotting function

def plotProfiles(x,various_y,labels):
	figure(1)
	
	grid(True)	
	if various_f.ndim == 1:
		plot(r,various_f,linestyle = '-', color = 'b',label = labels)
	elif various_f.ndim == 2:
		N = num.size(various_y,0)
		for i in range(N):
			plot(r,various_f[i],linestyle = '-', color = [float(i)/(N-1),0.0,1.0-float(i)/(N-1)],label = labels[i])
	legend()
	show()				


#integration scheme

def alpha_to_y(w,e0,alpha):
	if e0 == 0.:
		return num.log(alpha)/num.sqrt(2.)
	else:
		return num.sqrt(2.)/(3.*w)*(num.arcsinh(num.power(alpha,3.*w/2.)/num.sqrt(e0))) - num.arcsinh(1./num.sqrt(e0))

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

def evolvePsi(r0,f0,a_tab,w,Wm0,zcmb,s0):
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
	sol = Nintegrate.odeint(dynamic_y,y0,t,args=(params,))
	retour,v = sol.T
	retour = num.delete(retour,0)
	v = num.delete(v,0)
	
	a = alpha_tab[size(alpha_tab)-1]
	return retour,r0*num.sqrt(Wm0/2.)*(zcmb + 1.0)**(1./2.)/(299792458.0*20000.*num.power(a,3.*w/2.-1.))*v

def evolveProfile(r0,f0,af,w,Wm0,zcmb,dlogD_dloga_cmb):
	print 'evolving parameters : \n\ta = '+str(af)+'\n\tw = '+str(w)+'\n\tWm0 = '+str(Wm0)+'\n\tzcmb = '+str(zcmb)+'\n\tdlogD_dloga_cmb = '+str(dlogD_dloga_cmb)
	s0 = dlogD_dloga_cmb
	if size(af) == 1:		
		r = num.zeros(size(r0))
		f = num.zeros(size(r0))
		v = num.zeros(size(r0))
		for j in range(size(r0)):
			psi, v[j] = evolvePsi(r0[j],f0[j],af,w,Wm0,zcmb,s0)
			r[j] = r0[j]*psi
			if psi > 1.0e-2 and psi < 1.0e+5:
				f[j] = f0[j]/(psi**3.)
			else:
				f[j] = None
	else:
		r = num.zeros((size(af),size(r0)))
		f = num.zeros((size(af),size(r0)))
		v = num.zeros((size(af),size(r0)))
		for i in range(size(r0)):
			psi,v[i] = evolvePsi(r0[i],f0[i],af,w,Wm0,zcmb,s0)
			for j in range(size(af)):
				r[j][i] = r0[i]*psi[j]
				if psi[j] > 0. and psi[j] < 1.0e+5:
					f[j][i] = f0[i]/(psi[j]**3.)
				else:
					f[j][i] = None
	return r,f,v
