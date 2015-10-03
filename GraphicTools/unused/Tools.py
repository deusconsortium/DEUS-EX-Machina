from pylab import *
import numpy as num
from scipy.optimize import curve_fit
from scipy.special import hyp2f1
from scipy.special import gamma
from scipy.special import erfc
from scipy.special import erf
from scipy import integrate
import sys
import os

def P(i,x):
	if i == 1:
		return 2500. + 3125.*x**2. -750.*x**4.+2240.*x**6.+256.*x**8.
	elif i == 2:
		return 100. -65.*x**2.-8.*x**4.
	elif i == 3:
		return 3.*x*sqrt(5.)*(5.+2.*x**2.+8.*x**4.)

def K(x):
	Up = P(1,x)*sqrt(5.+x**2.)+P(2,x)*sqrt(5.+4.*x**2.)*(5.+4.*x**2.)**2.
	Down = 15.*x*sqrt((5.+x**2.)*(5.+4.*x**2.))*((x**2.-1.)*(5.+4.*x**2.)**2.*(arctan(sqrt(5.)/x)+arctan(sqrt(5.)/(2.*x)))-P(3,x))
	
	return -sqrt(2.*pi/5.)*Up/Down

def f(x):
	return (x**3.-3.*x)*(erf(sqrt(5./2.)*x)+erf(sqrt(5./2.)*x/2.))/2. + sqrt(2./(5.*pi))*((31.*x**2./4.+8./5.)*exp(-5.*x**2./8.)+(x**2./2.-8./5.)*exp(-5.*x**2./2.))

def I(x):
	return -(-(8.*x**4.+2.*x**2.+5.0)/(x**4.*(5.+4*x**2.)**2.)*9.*sqrt(5./(2.*pi)) + 3.*(x**2. - 1.)/(x**5.*sqrt(2.*pi))*(arctan(sqrt(5.)/x) + arctan(sqrt(5.)/(2.*x))))

def Jlim(a,b):
	return 3./a**5.*(a*b*(1.+a**2.*(b**2./3.-1.))*exp(-(a*b)**2./2.)+(1-a**2.)*sqrt(pi/2.)*erfc(a*b/sqrt(2.)))

def J(A,x0):
	J0 = 0.0
	if x0 <= 3.5:
		#zone entre x0 et 3.5
		print('x0 = '+str(x0)+', computing the integral numerically ...')
		N = 30
		x = num.linspace(x0,3.5,N)
		y = f(x)*x*exp(-(x*A)**2./2.)
		J0 = Integrate(y,x)
		#zone entre 3.5 et l'infini
		J0 = J0 + Jlim(A,3.5)
	else:
		J0 = Jlim(A,x0)
	return J0
	
def Cn(n,g):
	A = 1.-g**2.
	return g**n/(gamma(n+1.)*(1.-g**2.)**n)*1./(5.*sqrt(5.*pi))*(-1.)**n*2.**(-2. + n/2.)*gamma((1. + n)/2.)*(2.*sqrt(A)**(1.+n)*(-16.*(1. + 5.*A)**(-(1./2.) - n/2.) + 2.**(5. + n)*(4. + 5.*A)**(-(1./2.) - n/2.) + 5.*A*(31.*2.**(2. + n)*(1. + 5.*A)**((3. + n)/2.) + (4. + 5.*A)**((3. + n)/2.))*(4. + 25.*(A + A**2.))**(1./2.*(-3. - n))*(1. + n)) + 25.*sqrt(A)**(3. + n)*(1. + n)*(-6.*hyp2f1(1./2., (3. + n)/2., 3./2., -5.*A) - 3.*hyp2f1(1./2., (3. + n)/2., 3./2., -((5.*A)/4.)) + A*(3. + n)*(2.*hyp2f1(1./2., (5. + n)/2., 3./2., -5.*A) + hyp2f1(1./2., (5. + n)/2., 3./2., -((5.*A)/4.)))))
	
def L(v,g,n):
	if size(v)==1:
		Li = 0.
		for i in range(n+1):
			if i ==0:
				Li+= Cn(0,g)
			else:
				Li+= v**i*Cn(i,g)
		return Li
	else:
		Lt = num.zeros(size(v))
		for j in range(size(v)):
			Lt[j] = L(v[j],g,n)
		return Lt



def Sum(X):
	s = 0.0
	for k in range(size(X)):
		s += X[k]
	return s

def PowerLaw_fit(x,y,order = 2):
	S = num.empty([order + 1,order + 1])
	for i in range(order + 1):
		for j in range(order + 1):
			if i == 0 and j == 0:
				S[0,0] = size(x)
			else:
				S[i,j] = Sum(x**(i+j)) 	
	P = num.empty(order + 1)
	for i in range(order + 1):
		P[i] = Sum(y*(x**i))
	A = dot(inv(S),P)
	return A

def AB(x,y):
	a = num.empty(size(x)-1)
	b = num.empty(size(x)-1)
	s1 = (y[2]-y[0])/(x[2]-x[0])
	a[0] = -s1 + 2.*(y[1]-y[0])/(x[1]-x[0])
	b[0] = s1/(x[1]-x[0])+(y[0]-y[1])/(x[1]-x[0])**2.
	for i in range(size(x)-1):
		if i > 0:
			a[i] = a[i-1] + 2.*b[i-1]*(x[i]-x[i-1])
			b[i] = (y[i+1] - y[i])/(x[i+1]-x[i])**2. - a[i]/(x[i+1]-x[i])
	return a,b

def Continuity_coeffs(x,y,n):
	if n > 0 and n < size(x)-2:
		A = num.empty(4)
		Delta = x[n+1] - x[n]
		s0 = (y[n+1]-y[n-1])/(x[n+1]-x[n-1])
		s1 = (y[n+2]-y[n])/(x[n+2]-x[n])
		A[0] = y[n]
		A[1] = s0
		A[2] = 3.*(y[n+1]-y[n])/Delta**2. - (s1 + 2.*s0)/Delta
		A[3] = 2.*(y[n] - y[n+1])/Delta**3. + (s1 + s0)/Delta**2.
		
		return A
	if n == 0:
		s1 = (y[2]-y[0])/(x[2]-x[0])
		d = x[1] - x[0]
		A = num.empty(4)
		A[0] = y[0]
		A[1] = -s1 + 2.*(y[1]-y[0])/d
		A[2] = s1/d - (y[1] - y[0])/d**2.
		A[3] = 0.0
		return A
	if n == size(x) -2:
		s1 = (y[n+1]-y[n-1])/(x[n+1]-x[n-1])
		d = x[n+1] - x[n]
		A = num.empty(4)
		A[0] = y[n]
		A[1] = s1
		A[2] = (y[n+1] - y[n])/d**2. -s1/d
		A[3] = 0.0
		return A
	print('\n\nERROR !! n ='+str(n)+'\n\n')
	return -1

def D(y,x):
	d = num.empty(size(x))
	for i in range(size(x)):
		if i==0:
			d1 = (y[2] - y[0])/(x[2] - x[0])
			d[0] = 2.*(y[1] - y[0])/(x[1] - x[0]) - d1
		if i== size(x)-1:
			d[i] = (y[i]-y[i-1])/(x[i]-x[i-1])
		if (i > 0 and i< size(x)- 1):
			d[i] = (y[i+1]-y[i-1])/(x[i+1]-x[i-1])
	return d

def Minimum(t):
	if size(t) > 1:
		mini = t[0]
		for i in range(size(t)):
			if t[i]<mini:
				mini = t[i]
	else:
		mini = t
	return mini

def Maximum(t):
	if t.size > 1:
		maxi = t[0]
		for i in range(t.size):
			if t[i]>maxi:
				maxi = t[i]
	else:
		maxi = t
	return maxi


def sinc(x):
	if size(x)>1:
		y = num.empty(size(x))
		for i in range(size(x)):
			y[i] = sinc(x[i])
		return y
	else:
		if x == 0.0:
			return 1.0
		return sin(x)/x
	
def W2(k,Rcell,algo = 'CIC'):
	if algo == 'CIC':
		return sinc(k*Rcell/(2.*sqrt(2.)))**6.
	elif algo == 'TSC':
		return sinc(k*Rcell/(2.*sqrt(3.)))**9.
	else:
		return 0.0
	

def W(x):
	if x.size == 1:
		if x == 0.0:
			return 1.0
		return 3.*(sin(x) - x*cos(x))/x**3.
	else:
		w = num.empty(x.size)
		for i in range(x.size):
			if x[i] == 0.:
				w[i] = 1.
			else:
				w[i] = 3.*(sin(x[i]) - x[i]*cos(x[i]))/x[i]**3.
		return w

def Solve(x,F,a,options = 'none'):
	#print('calling Solve function with x,F two arrays of dimensions : ('+str(x.size)+', '+str(F.size)+') and a = '+str(a))
	k = 0
	if F[0] > a:
		for i in range(x.size):
			if F[i] < a and options != 'up':
				k = i
				#print('Solve done, returning x0 = '+str(((x[k-1]-x[k])*a + F[k-1]*x[k] - F[k]*x[k-1])/(F[k-1] - F[k])))
				return ((x[k-1]-x[k])*a + F[k-1]*x[k] - F[k]*x[k-1])/(F[k-1] - F[k])
			if F[i] == a:
				#print('Solve done, returning x0 = '+str(x[i]))
				return x[i]
		#print('Solve done but wihtout solution !! Returning -1')
		return -1.0
	else:
		for i in range(x.size):
			if F[i] >= a:
				k = i
				#print('Solve done, returning x0 = '+str(((x[k-1]-x[k])*a + F[k-1]*x[k] - F[k]*x[k-1])/(F[k-1] - F[k])))
				return ((x[k-1]-x[k])*a + F[k-1]*x[k] - F[k]*x[k-1])/(F[k-1] - F[k])
			if F[i] == a and options != 'down':
				#print('Solve done, returning x0 = '+str(x[i]))
				return x[i]
		#print('Solve done but wihtout solution !! Returning -1')
		return -1.0


def Get_d1(d,r,r1,method = 0):
	i0 = -1
	for i in range(size(r)):
		if r[i] >= r1 and i0 == -1:
			i0 = i
			break
	
	if method == 0:
		return d[i0 - 1 ] + (d[i0] - d[i0 - 1])*(r1 - r[i0 - 1])/(r[i0] - r[i0 - 1])
	if method == 1:
		n = i0-1
		A = Continuity_coeffs(r,d,n)
		x = num.linspace(0.,r[n+1]-r[n],20)
		y = A[0] + A[1]*x + A[2]*x**2. + A[3]*x**3.
		return Solve(y,x+r[n],r1)
	else:
		return Get_d1(d,r,r1)

def RemovePoints(tab,i0 = 0, i1 = 2):
	if i0 <= i1 and i0 >= 0 and i1 < size(tab):
		t = num.empty(size(tab))
		k = 0
		for i in range(size(tab)):
			if i >= i0 and i <= i1:
				t[i] = 'nan'
			else:
				t[i] = tab[i]
		return t
	return tab

def Primitive(F_x,x,from_zeros = False,method = 'trapeze'):
	Int = num.zeros(size(x))
	
	if from_zeros == True and x[0] > 0.:
		Int[0] = 0.5*F_x[0]*x[0]
	
	for i in range(size(x)-1):
		Int[i+1] = Int[0] + Integrate(F_x,x,imax = i+1,method = method)
	
	return Int

def Integrate(F_x,x,method = 'trapeze', imax = -1,from_zeros = True):
	I = 0.0
	
	if imax == -1:
		imax = size(x) - 1
	
	if from_zeros == True and x[0]>0.0:
		I+= 0.5*F_x[0]*x[0]
	
	if method == 'trapeze':
		for i in range(imax):
			I+= (F_x[i+1] + F_x[i])*(x[i+1] - x[i])/2.
		
	if method == 'cubic':
		for i in range(imax):
			A = Continuity_coeffs(x,F_x,i)
			if size(A) == 1:
				I += (F_x[i+1] + F_x[i])*(x[i+1] - x[i])/2.
			else:
				Delta = x[i+1] - x[i]
				I += A[0]*Delta + A[1]*Delta**2./2. + A[2]*Delta**3./3. + A[3]*Delta**4./4.
	
	return I

def Mean(y):
	m = 0.
	for i in range(size(y)):
		m += y[i]
	
	return m/size(y)

def Maximum_x_y(x,y):
	xmax = x[0]
	ymax = y[0]
	
	for k in range(x.size):
		if y[k] > ymax:
			ymax = y[k]
			xmax = x[k]
	
	return xmax,ymax

def Minimum_x_y(x,y):
	xmin = x[0]
	ymin = y[0]
	
	for k in range(size(x)):
		if y[k] < ymin:
			ymin = y[k]
			xmin = x[k]
	#print('founded : '+str(xmin)+', '+str(ymin))
	return xmin,ymin

def Get_delta_from_f(r,f,method = 0):
	print('computing delta from f ...')
	
	if method == 0:
		d = num.empty(size(r))
		for i in range(size(r) - 1):
			d[i] = (f[i+1]*r[i+1]**3.-f[i]*r[i]**3.)/(r[i+1]**3. - r[i]**3.)
		d[size(d) - 1] = d[size(d) - 2]
		return d - 1.
	else:
		return r/3.*D(f,r) + f - 1.0

def Get_r1r(r,f):
	d = Get_delta_from_f(r,f)
	return Get_r1(r,d + 1.)

def flip(x):
	z = num.empty(size(x))
	for i in range(size(x)):
		z[i] = x[size(x)-1-i]
	return z

def DetailProfiles(r1,f1,multiplicator = 2):
	r0 = flip(r1)
	f0 = flip(f1)
	n = int(multiplicator)
	if n >= 2:
		r = num.zeros(n*(size(r0)-1.)) 
		f = num.zeros(n*(size(r0)-1.))
		a,b = AB(r0,f0)
		for i in range(size(r0)-1):
			x = num.linspace(0.0,r0[i+1]-r0[i],n+1)
			for j in range(n):
				r[n*i + j] = r0[i] + x[j]
				f[n*i + j] = f0[i] + a[i]*x[j] + b[i]*x[j]**2.
		return flip(r),flip(f)
	else:
		r,f = r1,f1
	return r,f 

def Get_r1(r,f,method = 1):
	i0 = -1
	if f[0] > 1.:
		for i in range(size(r)):
			if f[i] < 1.0 and i0 == -1:
				i0 = i
				break
	else:
		for i in range(size(r)):
			if f[i] > 1.0 and i0 == -1:
				i0 = i
				break
	if i0 == -1:
		print('WARNING : Get_r1 returns 0')
		return 0
	else:
		if method == 0:
			return r[i0-1] + (1. - f[i0-1])*(r[i0] - r[i0 - 1])/(f[i0] - f[i0 - 1])
		if method == 1:
			n = i0 - 1
			A = Continuity_coeffs(r,f,n)
			x = num.linspace(0.,r[n+1]-r[n],20)
			y = A[0] + A[1]*x + A[2]*x**2. + A[3]*x**3.
			r1 = r[n] + Solve(x,y,1.)
			return r1
		else :
			return Get_r1(r,f,method = 0)

def Get_rmin_fmin(r,f):
	i0 = 0
	fmin = f[0]
	if f[0] < 1.:
		for i in range(size(r)):
			if f[i] > fmin:
				i0 = i
				fmin = f[i]
	else:
		for i in range(size(r)):
			if f[i] < fmin:
				i0 = i
				fmin = f[i]
	#on calcule l'approximation au second ordre
	
	
	r_used = num.empty(4)
	f_used = num.empty(4)
	if i0 >= 2 and i0 + 2 < size(r):
		for i in range(4):
			r_used[i] = r[i0 - 2 + i]
			f_used[i] = f[i0 - 2 + i]
		A = PowerLaw_fit(r_used,f_used)
		rm = -A[1]/(2.*A[2])
		return rm,A[0] + rm*A[1] + rm**2.*A[2] 
	
	
	if i0 == size(r) - 1 or i0 == 0:
		return 'nan','nan'
	
	return r[i0],f[i0]


#FONCTIONS DE DESSIN

def RBColor(i,N):
	if N == 1:
		return [1.,0.,0.]
	else:
		return [float(i)/(N-1),0.0,1.0-float(i)/(N-1)]

def errorplot(x,y,dy,Color = 'b'):
	if Color == 'b':
		Color = 'blue'
	if Color == 'g':
		Color = 'green'
	if Color == 'r':
		Color = 'red'
	plot(x,y, linestyle = '-', color = Color)
	fill_between(x,y + dy, y - dy,color = Color, alpha=.3)

