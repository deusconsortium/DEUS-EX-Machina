from DEUSTools import*
from DEUSCosmo import*

class DEUSAnalytics(DEUSCosmo):
	
	def __init__(self):
		DEUSCosmo.__init__(self)	
		
		self._s2_0 = num.zeros(3)
		self._s2_r = None
		self._S2_r = None
		self._r = None
		self._s8 = None
	
	def Load(self,cosmo,rtab):
		self._r = rtab
		
		self._s2_r = num.zeros((3,size(self._r)))
		self._S2_r = num.zeros((3,size(self._r)))
		
		if self._LoadSpectrum():
			self._Psmooth = num.copy(self._P0)
			self._computeSigmas()
			
			return True
		return False
		
	def getPowerSpectrum(self):
		return self._k,self._P0,self._Psmooth
	
	def resetSmoothing(self):
		self._Psmooth = num.copy(self._P0)
		self._computeSigmas()
	
	def globalSmoothSpectrum(self,smoothing_radius,function = 'th'):
		R = smoothing_radius
		if function == 'th':
			print 'global TOP-HAT  smoothing on R = '+str(R)+' [Mpc/h]'
			self._Psmooth = num.copy(self._P0*Wth(self._k*R)**2.)
			self._P0 = self._P0*Wth(self._k*R)**2.						
		elif function == 'exp':
			print 'global GAUSSIAN smoothing on R = '+str(R)+' [Mpc/h]'
			self._Psmooth = num.copy(self._P0*exp(-(self._k*R)**2.))
			self._P0 = self._P0*exp(-(self._k*R)**2.)
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
	
	def _computeSigmas0(self,P = 'self'):
		if P == 'self':
			P = self._P0	
			
		prefactor = 1./(2.*(num.pi)**2.)
		for i in range(size(self._s2_0)):
			self._s2_0[i] = prefactor*(integrate(self._k,self._k**(2. + 2.*i)*P))


	def _computeSigmas(self,P = 'self'):
		if P == 'self':
			P = self._P0
			
		prefactor = 1./(2.*(num.pi)**2.)
		for i in range(size(self._s2_0)):
			self._s2_0[i] = prefactor*(integrate(self._k,self._k**(2. + 2.*i)*P))
		
		# ! careful to the pi factor in numpy sinc(x) = sin(pi*x)/(pi*x) !
		for i in range(3):
			for r in range(size(self._r)):
				self._s2_r[i][r] = prefactor*integrate(self._k,self._k**(2. + 2.*i)*P*num.sinc(self._k*self._r[r]/num.pi))
				self._S2_r[i][r] = prefactor*integrate(self._k,self._k**(2. + 2.*i)*P*Wth(self._k*self._r[r]))
	
	def computeR1Distribution(self,R1,v0 = 0.,individualR1smoothFactor = 0.5,indiviualSmoothFunction = 'th',normalized = True):			
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
			F = (G**2. - B2)/(G*(1. - B2))
			x = sqrt((B2**2. + G**2. - 2.*B2*G**2.)/(B2 - G**2.)**2.)
			xd = sqrt((B2d**2. + G**2. - 2.*B2d*G**2.)/(B2d - G**2.)**2.)
			y = self._s2_r1(0,r1)/self._S2_r1(0,r1) - self._s2_r1(1,r1)/self._S2_r1(1,r1)
			yd = self._S2_r1(2,r1)/self._s2_r1(1,r1) - self._S2_r1(1,r1)/self._s2_r1(0,r1)
			
			Nth[i] = 3.*B2/(2.*r1*num.pi**2.*sqrt(2.))*G*sqrt(1. - G**2.)/(G**2. - B2)**2.*y*In(x,1,v0/F)				
			Nthd[i] = r1*B2/(6.*num.pi**2.*sqrt(2.))*G*sqrt(1. - G**2.)/(G**2. - B2d)**2.*yd*In(xd,1,v0/F)				
			
			if individualR1smoothFactor > 0.0:
				self.resetSmoothing()
	
		if normalized:
			fact = integrate(R1,Nth)
			Nth /= fact
			factd = integrate(R1,Nthd)
			Nthd /= factd
			
		return Nth,Nthd
		
	def computeMeanProfile(self,d0,Wm0,w,r1,af = None):
		f0 = 1. + d0*(self._S2_r[0]- self._S2_r[1]*self._S2_r1(0,r1)/self._S2_r1(1,r1))/(self._s2_0[0]*(1. - self._Beta_r1(r1)))
		
		if af is None:
			return self._r,f0,-(f0 - 1.)*self._dlogD_dloga_init/3.
		else:
			return evolveProfile(self._r,f0,af,w,Wm0,self._zinit,self._dlogD_dloga_init)
	
	def computeLinearExtremumDistributionFromHeigh(self,z = None,nmax = 200):
		if z is None:
			z = self._zinit
		
		g = self._Gamma()
		s0 = sqrt(self._getSpectrumFactor(z)*self._s2_0[0])	
		
		print 'founded s0 = '+str(s0)+' and g = '+str(g)+' and s8 = '+str(self._s8)+ ' for z = '+str(z)+' with '+str(nmax)+' recursions'	
		
		v_tab = num.linspace(-5.0,-0.01,100)
		Np = num.zeros(size(v_tab))
		
		for i in range(size(Np)):
			v = v_tab[i]
			CN = 0.0
			for j in range(nmax):
				CN += Cn(g,j)*(v)**j
				
			Np[i] = exp(-v**2./(2.*(1.-g**2.)))/((2.*num.pi)**2.*sqrt(2.*(1.-g**2.)))*CN
		
		return s0*v_tab + 1.0,Np*40.*sqrt(5.)*(num.pi)**(3./2.)/(s0*(29.-6.*sqrt(6.)))
	
	def computeLinearDensityPDF(self,z = None,nmax = 100):			
		s0 = sqrt(self._getSpectrumFactor(z)*self._s2_0[0])
		
		v_tab = num.linspace(-3.0,3.,100)
		
		Np = 1./(num.sqrt(2.*num.pi*s0**2.))*exp(-v_tab**2./2.)
		return s0*v_tab + 1.,Np
	
	def computeZeldovitchExtremumDistributionFromHeigh(self,w,Wm0,z = None,nmax = 100):		
		#computing the initial distribution at zinit
		d0,N0 = self.computeLinearExtremumDistributionFromHeigh(self._zinit,nmax)
		
		if z is None:
			return d0,N0
		else:
			print "evolving N(d) by Zel'Dovitch..."
			
			psi = num.ones(size(d0))
			for i in range(size(psi)):
				psi[i],psip = evolvePsiLinear(d0[i],1./(z+1.),w,Wm0,self._zinit,self._dlogD_dloga_init)
			N = psi**3./(1. - 3.*(d0)/psi*derivative(d0,psi))*N0
			
			return d0/psi**3.,N	
		
	
	def computeEvolvedExtremumDistributionFromHeigh(self,w,Wm0,z = None,nmax = 100,zinit = None):
		if zinit is None:
			zinit = self._zinit
			dlogD = self._dlogD_dloga_init
		else:
			dlogD = self.getDlogD_dloga(zinit)
		
		#computing the initial distribution at zinit
		d0,N0 = self.computeLinearExtremumDistributionFromHeigh(zinit,nmax)
		
		if z is None:
			return d0,N0
		else:
			print 'evolving N(d) ...'
			
			psi = num.ones(size(d0))
			for i in range(size(psi)):
				psi[i],psip = evolvePsi(d0[i],1./(z+1.),w,Wm0,zinit,dlogD)
			N = psi**3./(1. - 3.*(d0)/psi*derivative(d0,psi))*N0
			
			return d0/psi**3.,N	
		
			
	#miscellanous functions
	
	def _getSpectrumFactor(self,z = None):
		if z is None:
			z = self._zinit
		elif z < 0.0:
			return 1.0
		else:
			D0 = solve(self._D_tab,self._a_tab,1.0)
			a = 1./(z+1.)
			Dz = solve(self._D_tab,self._a_tab,a)
			return (Dz/D0)**2.
	
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
