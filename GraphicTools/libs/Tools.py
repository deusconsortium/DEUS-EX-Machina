from pylab import *
import numpy as num

def f(x):
	return (x**3.-3.*x)*(erf(sqrt(5./2.)*x)+erf(sqrt(5./2.)*x/2.))/2. + sqrt(2./(5.*pi))*((31.*x**2./4.+8./5.)*exp(-5.*x**2./8.)+(x**2./2.-8./5.)*exp(-5.*x**2./2.))

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
		return 'nan'

def Wth(x):
	return 3.*(sin(x) - x*cos(x))/(x*x*x)

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



