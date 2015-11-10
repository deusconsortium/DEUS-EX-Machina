from DEUSTools import*

class DEUSCosmo:	
	def __init__(self,zinit):
		self._power_spectrum_path = '/data_bingo/Babel/data/'
		self._zinit = zinit
		
		self._a = 0.0
		self._WmToday = 0.0
		self._Wm0 = 0.0
		self._w = -1.0
		self._H0 = 0.0
		
		self._boxlen = 0
		self._npart = 0
		self._cosmo = None
		
		self._dlogD_dloga_init = 1.0
		
		self._P0 = num.empty([1])
		self._k = num.empty([1])
		self._a_tab = num.empty([1])
		self._D_tab = num.empty([1])
		self._dlogD_dloga_tab = num.empty([1])
	
	def _SetSimu(self,boxlen,npart,cosmo):
		self._boxlen = boxlen
		self._npart = npart
		self._cosmo = cosmo
		if cosmo=="lcdmw5":
			self._WmToday = 0.26
			self._H0 = 0.72
			self._w = -1.0
		elif cosmo=="rpcdmw5":
			self._WmToday = 0.23
			self._H0 = 0.72
			self._w = -0.86
		else:
			print "error : unknown cosmo : '"+cosmo+"'"
			self._cosmo = 'lcdmw5'
			self._WmToday = 0.26
			self._H0 = 0.72
			self._w = -1.0		
		self._Wm0 = Omega_m_0(self._WmToday,self._w,self._zinit)
		
		print 'actual Cosmo : \n\tWm0 = '+str(self._Wm0)+'\n\tWm = '+str(self._WmToday)+'\n\tw = '+str(self._w) 
		
	
	def _LoadSpectrum(self):
		if self._cosmo is None:
			print 'error : must set cosmo before calling _LoadSpectrum'
			return False
		else:
			founded = True
			try:
				k,P = num.loadtxt(self._power_spectrum_path + 'pk_' + self._cosmo + '.dat', usecols=(0, 1), unpack=True)
				a,D,dD = num.loadtxt(self._power_spectrum_path + 'mpgrafic_input_' + self._cosmo + '.dat', usecols=(0,2,3), unpack=True)
			except IOError:
				founded = False
				print 'error : no P/k or D/a file in ' + self._power_spectrum_path
			
			if founded:
				self._dlogD_dloga_init = solve(dD,a,1./(self._zinit + 1.0))
				if self._dlogD_dloga_init is None:
					self._dlogD_dloga_init = 1.0
					
				self._a_tab = num.copy(a)
				self._D_tab = num.copy(D)
				self._dlogD_dloga_tab = num.copy(dD)
				
				self._k = num.copy(k)
				self._P0 = num.copy(P)
			
			return founded

	def _GaussianSmoothing(self,R):
		self._P0 = self._P0*exp(-(self._k*R)**2.)
	
	def _NumericalSmoothing(self):
		"""kmin = 1./float(self._boxlen)
		kmax = float(self._npart)*kmin
		
		kused = num.linspace(kmin,kmax,size(self._k))
		Pused = num.zeros(size(kused))
		for i in range(size(kused)):
			Pused = solve(self._P0,self._k,kused[i])
		
		self._k = num.copy(kused)
		self._P0 = num.copy(Pused)
		"""
		
		#R = 0.62*float(self._boxlen)/float(self._npart)
		#self._P0 = self._P0*Wth(self._k*R)**2.
		#self._P0 = self._P0*exp(-(self._k*R)**2.)
