from pylab import *
import numpy as num
import sys
import os

from scipy import integrate as Nintegrate
from scipy.special import erfc
from scipy.special import erf
from scipy.special import hyp2f1
from scipy.special import gamma
from scipy.misc import factorial


#tedious functions

def f(x):
	return (x**3.-3.*x)*(erf(sqrt(5./2.)*x)+erf(sqrt(5./2.)*x/2.))/2. + sqrt(2./(5.*pi))*((31.*x**2./4.+8./5.)*exp(-5.*x**2./8.)+(x**2./2.-8./5.)*exp(-5.*x**2./2.))

def I_dn_dr1(x):
	return -(-(8.*x**4.+2.*x**2.+5.0)/(x**4.*(5.+4*x**2.)**2.)*9.*sqrt(5./(2.*pi)) + 3.*(x**2. - 1.)/(x**5.*sqrt(2.*pi))*(arctan(sqrt(5.)/x) + arctan(sqrt(5.)/(2.*x))))

def Wg(a,b):
	return exp(-(a+b)**2./2.)-exp(-(a-b)**2./2.) + b/2.*sqrt(2.*num.pi)*(erf((a-b)/sqrt(2.)) + erf((a+b)/sqrt(2.)))

def In(y,n,x0 = None):
	#defini par l'integrale entre 0 et +infini, donc pas de (-1)^n	
	if x0 is None:
		x = y 
		return   2.**(n/2.-2.)*gamma((1.+n)/2.)/(5.*sqrt(5.*num.pi))*( 
			 10.*(n+1.)*(5.+x**2.)**(-(3.+n)/2.)
			 - 32.*(5.+x**2.)**(-(1.+n)/2.)
			 + 155.*2.**(3.+n)*(1.+n)*(5.+4.*x**2.)**(-(3.+n)/2.)
			 + 2.**(6.+n)*(5.+4.*x**2.)**(-(1.+n)/2.)
			 - 150.*(1.+n)*x**(-(n+3.))*hyp2f1(1./2.,(3.+n)/2.,3./2.,-5./x**2.)
			 - 75.*(1.+n)*x**(-(n+3.))*hyp2f1(1./2.,(3.+n)/2.,3./2.,-5./(4.*x**2.))
			 + 50.*(1.+n)*(3.+n)*x**(-(n+5.))*hyp2f1(1./2.,(5.+n)/2.,3./2.,-5./x**2.)
			 + 25.*(1.+n)*(3.+n)*x**(-(n+5.))*hyp2f1(1./2.,(5.+n)/2.,3./2.,-5./(4.*x**2.)))
	else:
		x_tab = num.linspace(x0,3.*sqrt(n + 3.)/y,1000)
		F = x_tab**n*f(x_tab)*exp(-(x_tab*y)**2./2.)
		return integrate(x_tab,F)

def Cn(g,n):
	return (-g)**n/(factorial(n)*(1.-g**2.)**n)*In(1./sqrt(1.-g**2.),n)

#files functions

def listdirHidden(dir_path):
	flist = []
	ini = os.listdir(dir_path)
	for i in range(size(ini)):
		if not ini[i].startswith('.'):
			flist.append(ini[i])
	return flist

def filesThatEndAs(dir_path,end,end_char_to_cut = 0):
	print 'looking in '+str(dir_path)+' for files ending as '+end
	files = os.listdir(dir_path)
	nE = len(end)
	file_list = []
	for i in range(size(files)):
		if files[i][-nE:] == end:
			if end_char_to_cut > 0:
				file_list.append(files[i][:-end_char_to_cut])
			else:
				file_list.append(files[i])
		else:
			print 'comparing '+files[i][-nE:]+' to '+end
			
	return list(set(file_list))

#useful functions

def fraction(x1,y1,x2,y2,Npoints = 100):
	x = num.linspace(max(min(x1),min(x2)),min(max(x1),max(x2)),Npoints)
	f1 = SplineFit(x1,y1)
	f2 = SplineFit(x2,y2)
	
	return x,f2(x)/f1(x)
	

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
					return x[i+1] + (y[i+1] - y0)*(x[i] - x[i+1])/(y[i+1] - y[i])
		else:
			for i in range(size(y) - 1):
				if y[i] <= y0 and y[i+1] >= y0:
					return x[i + 1] + (y[i+1] - y0)*(x[i] - x[i+1])/(y[i+1] - y[i])
	return None


def Wth(x):
	W =  3.*(num.sin(x) - x*num.cos(x))/(num.power(x,3))
	if size(x) > 1:
		nan_pos = num.isnan(W)
		W[nan_pos] = 1.0 
	else:
		if x is 0.0:
			return 1.0
	return W
	
def integrate(x,y,x0 = None,x1 = None):
	#x0 = max(i for i in [x0,min(x)] if i is not None)
	#x1 = min(i for i in [x1,max(x)] if i is not None)
	fs = SplineFit(x,y)
	return fs.Integrate(x0,x1)
	"""
	else:
		I = 0.		
		for i in range(size(x) - 1):
			I += 0.5*(max(0.0,min(x[i+1],x1) - max(x[i],x0)))*(y[i] + (y[i+1]-y[i])*(max(x[i],x0)-x[i])/(x[i+1]-x[i]) + y[i] + (y[i+1]-y[i])*(min(x[i+1],x1)-x[i])/(x[i+1]-x[i]))
		return I"""
	
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

def getCICCentralDensity(r,f,R):
	if r[0] > 0.0:
		f3 = f[0:3]
		r3 = r[0:3]
		A = getPowerLawFitCoefficients(r3,f3,2)
		
		#adding 2 points
		r = num.insert(r,0,r[0]/2.)
		r = num.insert(r,0,0.0)
		f = num.insert(f,0,A[0]+A[1]*r[1]+A[2]*r[1]**2.)
		f = num.insert(f,0,A[0])
	
	#defining densities
	d = getDensity(r,f)
	d[0] = f[0]
		
	ds = SplineFit(r,d)
	
	X = num.linspace(0.0,R,20)	
	Y = ds(X)*X*X
	return 3./R**3.*integrate(X,Y)

def CICDensitySmoothing(r,f,R,plotting = False):
	if R == 0.0 :
		return f0
	
	if plotting:
		print 'spherical CIC smoothing on '+str(R)+' [Mpc/h]'	
	
	#getting the f0
	f3 = f[0:3]
	r3 = r[0:3]
	A = getPowerLawFitCoefficients(r3,f3,2)
	
	#adding 2 points
	r = num.insert(r,0,r[0]/2.)
	r = num.insert(r,0,0.0)
	f = num.insert(f,0,A[0]+A[1]*r[1]+A[2]*r[1]**2.)
	f = num.insert(f,0,A[0])
	
	#defining densities
	d = getDensity(r,f)
	d[0] = f[0]
	plot(r,d,'b--')
	dcic = num.zeros(size(d))	
	ds = SplineFit(r,d)
	
	x = num.copy(r)
	for i in range(size(r)):
		if r[i] <= R:
			if r[i] == 0.0:
				X = num.linspace(0.0,R,10)	
				Y = ds(X)*X*X
				dcic[i] = 3./R**3.*integrate(X,Y)
			else:
				X = num.linspace(0.0,R-r[i],10)
				Y = ds(X)*X*X
				dcic[i] = 3./R**3.*integrate(X,Y)
				
				X = num.linspace(R-r[i],R+r[i],10)
				Y = ds(X)*X*(R**2 - (X-r[i])**2.)/(2.*r[i])
				dcic[i] += integrate(X,Y)
		else:
			X = num.linspace(r[i]-R,r[i]+R,10)
			Y = ds(X)*X*(R**2 - (X-r[i])**2.)/(2.*r[i])
			dcic[i] = 3./(2.*R**3.)*integrate(X,Y)
	
	#correcting the last 2 points
	dcic[size(dcic)-1] = d[size(d)-1]
	dcic[size(dcic)-2] = d[size(d)-2]
	
	#computing fcic
	fcic = getMassContrast(r,dcic)
	
	if plotting:
		figure(1)
		grid(True)
		plot(r,f,'b-',label = 'f')
		plot(r,d,'b--', label = 'd')
		plot(r,fcic,linestyle = '-', marker = 'o',color = 'r',label = 'f cic')
		plot(r,dcic,linestyle = '--', marker = '+',color = 'r',label = 'd cic')
		show()
	
	return fcic[2:]
	
	

def gaussianDensitySmoothing(r,f0,R,end_point_to_remove = 5,r1_factor = 2):
	if R == 0.0 :
		return f0
	
	#print 'GAUSSIAN smoothing on '+str(R)+' [Mpc/h]'
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

def getMassContrast(r,d):
	if r[0] > 0.:
		#getting the f0
		d3 = d[0:3]
		r3 = r[0:3]
		A = getPowerLawFitCoefficients(r3,d3,2)
		
		#adding 1 points
		ru = num.insert(r,0,0.0)
		du = num.insert(d,0,A[0])
	else:
		ru = r
		du = d
	
	f = num.zeros(size(ru))
	f[0] = du[0]
	for i in range(size(d) - 1):
		if i == 0:
			f[1] = du[0]/4. + 3.*du[1]/4.
		else:
			x = ru[i]/ru[i+1]
			f[i+1] = f[i] *x**3. + 1./4.*(1-x)*(x**2.*(d[i+1]+3.*d[i]) + 2.*x*(d[i+1]+d[i]) + d[i] + 3.*d[i+1])

	if size(ru) > size(r):
		return f[1:]
	return f

def getDensity(r,f,method = 'simple'):
	d = num.ones(size(r))
	
	if method == 'defaut':
		return f + r/3.*derivative(r,f)
	elif method == 'simple':
		for i in range(size(r)-1):
			d[i] = 1./3.*(f[i+1]*r[i+1]**3. - f[i]*r[i]**3.)/(r[i]**2.*(r[i+1]-r[i]))
	else:
		x=r[size(r)-2]/r[size(r)-1]
		d[size(r)-1] = (f[size(r)-1]-x**3.*f[size(r)-2])/((2.+x**2.)*(1.-x))
		d[size(r)-2] = d[size(r)-1]
		for i in range(size(r)-2):
			j = size(r)-3-i
			x = r[j]/r[j+1]
			d[j] = 4./3.*(f[j+1]-x**3.*f[j])/((1.-x)*(1.+x**2.)) - d[j+1]/3.*(5.+x**2.)/(1.+x**2.)
	
	d[size(r)-1]=d[size(r)-2]
	return d

def Epsilon0(w,Wm0,zcmb):
	if Wm0 == 1.0:
		return 0.0
	return (1.-Wm0)/Wm0*(zcmb+1.)**(3.*w)

def Omega_m_0(Wmtoday,w,zinit):
	x = 1./(zinit+1.0)**(3.*w)
	return Wmtoday*x/(1.-Wmtoday*(1.-x))

def Omega_m(alpha,Wm0,w,zcmb):
	return 1./(1.+Epsilon0(w,Wm0,zcmb)*alpha**(-3.*w))

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
			
	def Integrate(self,x0 = None,x1 = None):
		x0 = max(i for i in [x0,min(self._x)] if i is not None)
		x1 = min(i for i in [x1,max(self._x)] if i is not None)
		
		I = 0.0
		for j in range(size(self._x) - 1):
			X = max((min(x1,self._x[j+1]) - max(x0,self._x[j]))/(self._x[j+1]-self._x[j]),0.0)
			I += (self._x[j+1]-self._x[j])*(self._c[j][0]*X + self._c[j][1]*X**2./2. + self._c[j][2]*X**3./3. + self._c[j][3]*X**4./4.)
		return I
		
	def Primitive(self):
		I = num.zeros(size(self._x))
		for j in range(size(self._x) - 1):
			I[j+1] = I[j] + (self._x[j+1]-self._x[j])*(self._c[j][0] + self._c[j][1]/2. + self._c[j][2]/3. + self._c[j][3]/4.)
		return I

def IncreaseResolution(x,y,factor = 2):
	X = num.linspace(x[0],x[size(x)-1],factor*size(x))
	Sx = SplineFit(x,y)
	return X,Sx(X)

#plotting function

def errorplot(x,y,dy,linestyle = '-',color = 'b',marker = 'o',label = None):
	errorbar(x, y, yerr=dy, xerr=None,fmt= None, ecolor=color)
	plot(x,y,linestyle = linestyle,color = color,marker = marker,label = label)

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

def Plot(x,y):
	figure(0)
	grid(True)
	plot(x,y)
	show()

#miscelanous functions

def selectFile(folder_path,nb_char_to_remove = None):
	folder = os.listdir(folder_path)
	print '\nchoose one file to load :'
	for i in range(size(folder)):
		if nb_char_to_remove is not None:
			print "- "+str(i)+" : "+folder[i][:-nb_char_to_remove]
		else:
			print "- "+str(i)+" : "+folder[i]
	
	num = int(raw_input("\nenter choice : "))
	if num >= 0 and num < size(folder):
		if nb_char_to_remove is not None:
			return folder[num][:-nb_char_to_remove]
		else:
			return folder[num]
	else:
		print 'index out of bounds'
		return None

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

def evolvePsiLinear(f0,a_tab,w,Wm0,zcmb,s0):
	if f0 == 1.:
		if size(a_tab) > 1:
			psi = num.empty(size(a_tab))
			for i in range(size(psi)):
				psi[i] = 1.
			return psi
		else:
			return 1.0
	
	n = Eta(a_tab,w,Wm0,s0,zcmb)
	return 1.0 - (f0 - 1.)*n,0.0

def evolvePsi(f0,a_tab,w,Wm0,zcmb,s0):
	if f0 == 1.:
		if size(a_tab) > 1:
			psi = num.empty(size(a_tab))
			for i in range(size(psi)):
				psi[i] = 1.
			return psi
		else:
			return 1.0
	
	e0 = Epsilon0(w,Wm0,zcmb)
	y0 = [1.,s0*(1. - f0)/3.*sqrt(2.*Wm0)]
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
	return retour,v/sqrt(2.*Omega_m(a,Wm0,w,zcmb))

def evolveProfile(r0,f0,af,w,Wm0,zcmb,dlogD_dloga_cmb):
	print 'evolving parameters : \n\ta = '+str(af)+'\n\tw = '+str(w)+'\n\tWm0 = '+str(Wm0)+'\n\tzcmb = '+str(zcmb)+'\n\tdlogD_dloga_cmb = '+str(dlogD_dloga_cmb)
	s0 = dlogD_dloga_cmb
	if size(af) == 1:		
		r = num.zeros(size(r0))
		f = num.zeros(size(r0))
		v = num.zeros(size(r0))
		for j in range(size(r0)):
			psi, v[j] = evolvePsi(f0[j],af,w,Wm0,zcmb,s0)
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
			psi,v[i] = evolvePsi(f0[i],af,w,Wm0,zcmb,s0)
			for j in range(size(af)):
				r[j][i] = r0[i]*psi[j]
				if psi[j] > 0. and psi[j] < 1.0e+5:
					f[j][i] = f0[i]/(psi[j]**3.)
				else:
					f[j][i] = None
	return r,f,v
