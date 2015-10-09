from DEUStools import*

class DEUSanalytics:
	
	def __init__(self,zinit):		
		self._power_spectrum_path = '/efiler2/bingo_save/Babel/data/'
	
		self._zinit = zinit
		self._dlogD_dloga_cmb = 1.0
		self._P = num.empty([1])
		self._Psmooth = num.empty([1])
		self._k = num.empty([1])
		
		self._s2_0 = num.zeros(3)
		self._s2_r = None
		self._S2_r = None
		self._r = None
	
	def Load(self,cosmo,rtab):
		self._r = rtab
				
		self._s2_r = num.zeros((3,size(self._r)))
		self._S2_r = num.zeros((3,size(self._r)))		
		founded = True
		
		try:
			k,P = num.loadtxt(self._power_spectrum_path + 'pk_' + cosmo + '.dat', usecols=(0, 1), unpack=True)
			a,D,dD = num.loadtxt(self._power_spectrum_path + 'mpgrafic_input_' + cosmo + '.dat', usecols=(0, 2,3), unpack=True)
		except IOError:
			founded = False
			print 'error : no P/k or D/a file in ' + self._power_spectrum_path
		
		if founded:
			self._dlogD_dloga_cmb = solve(dD,a,1./(self._zinit + 1.0))
			if self._dlogD_dloga_cmb is None:
				self._dlogD_dloga_cmb = 1.0
			
			#transcripting P/k in linear scale
			self._k = num.copy(k)
			self._P = num.copy(P)
			self._Psmooth = num.copy(self._P)

		return founded
		
	def getPowerSpectrum(self):
		return self._k,self._P,self._Psmooth
	
	def resetSmoothing(self):
		self._Psmooth = num.copy(self._P)
	
	def globalSmoothSpectrum(self,smoothing_radius,function = 'th'):
		R = smoothing_radius
		if function == 'th':
			print 'global TOP-HAT  smoothing on R = '+str(R)+' [Mpc/h]'
			self._Psmooth = num.copy(self._P*Wth(self._k*R)**2.)
			self._P = self._P*Wth(self._k*R)**2.						
		elif function == 'exp':
			print 'global GAUSSIAN smoothing on R = '+str(R)+' [Mpc/h]'
			self._Psmooth = num.copy(self._P*exp(-(self._k*R)**2.))
			self._P = self._P*exp(-(self._k*R)**2.)
		else:
			print 'unknown smoothing function : '+ function
		self._computeSigmas()
	
	def localSmoothSpectrum(self,smoothing_radius,function = 'th'):
		R = smoothing_radius
		if function == 'th':
			print 'TH  smoothing on R = '+str(R)+' [Mpc/h]'
			self._Psmooth = self._Psmooth*Wth(self._k*R)**2.
			self._computeSigmas(self._Psmooth)
		elif function == 'exp':
			print 'EXP smoothing on R = '+str(R)+' [Mpc/h]'
			self._Psmooth = self._Psmooth*exp(-(self._k*R)**2.)
			self._computeSigmas(self._Psmooth)
		else:
			print 'unknown smoothing function : '+ function
	
	def _computeSigmas(self,P = 'self'):
		if P == 'self':
			P = self._P	
		for i in range(size(self._s2_0)):
			self._s2_0[i] = integrate(self._k,self._k**(2. + 2.*i)*P)
		
		# ! careful to the pi factor in numpy sinc(x) = sin(pi*x)/(pi*x) !
		for i in range(3):
			for r in range(size(self._r)):
				self._s2_r[i][r] = integrate(self._k,self._k**(2. + 2.*i)*P*num.sinc(self._k*self._r[r]/num.pi))
				self._S2_r[i][r] = integrate(self._k,self._k**(2. + 2.*i)*P*Wth(self._k*self._r[r]))
	
	def computeR1Distribution(self,R1,individualR1smoothFactor = 0.5,indiviualSmoothFunction = 'th',normalized = True):			
		Nth = num.zeros(size(R1))		
		Nthd = num.zeros(size(R1))		
		
		for i in range(size(R1)):
			r1 = R1[i]
			if individualR1smoothFactor > 0.0:
				self.resetSmoothing()
				self.localSmoothSpectrum(individualR1smoothFactor*r1,indiviualSmoothFunction)
			
			B2 = self._Beta_r1(r1)
			B2d = self._Beta_r1_d(r1)
			G = self._Gamma()
			x = sqrt((B2**2. + G**2. - 2.*B2*G**2.)/(B2 - G**2.)**2.)
			xd = sqrt((B2d**2. + G**2. - 2.*B2d*G**2.)/(B2d - G**2.)**2.)
			y = self._s2_r1(0,r1)/self._S2_r1(0,r1) - self._s2_r1(1,r1)/self._S2_r1(1,r1)
			yd = self._S2_r1(2,r1)/self._s2_r1(1,r1) - self._S2_r1(1,r1)/self._s2_r1(0,r1)
			
			Nth[i] = 3.*B2/(2.*r1*num.pi**2.*sqrt(2.))*G*sqrt(1. - G**2.)/(G**2. - B2)**2.*y*I_dn_dr1(x)				
			Nthd[i] = r1*B2/(6.*num.pi**2.*sqrt(2.))*G*sqrt(1. - G**2.)/(G**2. - B2d)**2.*yd*I_dn_dr1(xd)				
		
		if normalized:
			fact = integrate(R1,Nth,from_zeros = False)
			Nth /= fact
			factd = integrate(R1,Nthd,from_zeros = False)
			Nthd /= factd
			
		return Nth,Nthd
		
	def computeMeanProfile(self,d0,Wm0,w,r1,af = None):
		f0 = 1. + d0*(self._S2_r[0]- self._S2_r[1]*self._S2_r1(0,r1)/self._S2_r1(1,r1))/(self._s2_0[0]*(1. - self._Beta_r1(r1)))
		
		if af is None:
			return self._r,f0,-(f0 - 1.)*self._dlogD_dloga_cmb/3.
		else:
			return evolveProfile(self._r,f0,af,w,Wm0,self._zinit,self._dlogD_dloga_cmb)
	
	#miscellanous functions
	
	def _Gamma(self):
		return self._s2_0[1]/sqrt(self._s2_0[2]*self._s2_0[0])
	
	def _S2_r1(self,i,r1):
		f = SplineFit(self._r,self._S2_r[i])
		return f(r1)
	
	def _s2_r1(self,i,r1):
		f = SplineFit(self._r,self._s2_r[i])
		return f(r1)
	
	def _Beta_r1(self,r1):		
		return self._S2_r1(0,r1)*self._s2_0[1]/(self._S2_r1(1,r1)*self._s2_0[0])
	
	def _Beta_r1_d(self,r1):		
		return self._s2_r1(0,r1)*self._s2_0[1]/(self._s2_r1(1,r1)*self._s2_0[0])
