from pylab import *
import numpy as num
from scipy.special import erfc
from scipy.special import erf

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
		return None

def Wth(x):
	return 3.*(num.sin(x) - x*num.cos(x))/(num.power(x,3.0))

def integrate(x,y,imax = -1,from_zeros = True):
	Int = 0.0
	if imax == -1:
		imax = size(x)

	if from_zeros == True and x[0] > 0.0:
		Int += 0.5*y[0]*x[0]
		
	for i in range(imax - 1):
		Int += 0.5*(y[i+1] + y[i])*(x[i+1] - x[i])
	return Int





