from DEUSgraphics import DEUSgraphics

class DEUSmodel(DEUSgraphics):
	
	def __init__(self,isOverDensity,datapath = ''):
		(DEUSgraphics.__init__(self,isOverDensity,datapath)
		
		self._power_spectrum_path = '/efiler2/bingo_save/Babel/data/'
		
		P_file = "/home/pdefromont/advil/pk_"+cosmo+'.dat'
		D_a_file = "/home/pdefromont/advil/efiler2/bingo_save/Babel/data/mpgrafic_input_"+cosmo+".dat"
		
		self._zcmb = 1100.0
		self._dlogD_dloga_cmb = 1.0
		self._cosmo = ''
		self._P = num.empty([1])
		self._k = num.empty([1])
		
		self._s2_0 = num.zeros(3)
		self._s2_r = num.empty([1])
		self._S2_r = num.empty([1])
	
	def Load(self,cosmo,file_name):
		super.Load(file_name)
		
		self._s2_r = num.zeros(size(self._r))
		
		self._cosmo = cosmo
		founded = True
		
		try:
			k,P = num.loadtxt(self._power_spectrum_path + 'pk_' + self._cosmo + '.dat', usecols=(0, 1), unpack=True)
			a,D,dD = num.loadtxt(self._power_spectrum_path + 'mpgrafic_input_' + self._cosmo + '.dat', usecols=(0, 2,3), unpack=True)
		except IOError:
			founded = False
			print 'error : no P/k or D/a file in ' + self._power_spectrum_path
		
		if founded:
			self._dlogD_dloga_cmb = solve(dD,a,1./(self._zcmb + 1.0))
			if self._dlogD_dloga_cmb is None:
				self._dlogD_dloga_cmb = 1.0
			
			self._k = k
			self._P = P*num.square(Wth(.87*self._k*float(self._boxlen)/float(self._npart)))
			
			self._computeSigmas()
	
	def smoothSpectrum(self,ramses_smoothing_radius,function = 'th'):
		R = ramses_smoothing_radius*float(self.boxlen)
		if function == 'th':
			self._P = self._P*num.square(Wth(self._k*R))
			self._computeSigmas()
		elif function == 'exp':
			self._p = self._P*num.exp(-num.square(self._k*R))
			self._computeSigmas()
		else:
			print 'unknwon smoothing function : '+ function
	
	def _computeSigmas(self):	
		for i in range(size(self._s2_0)):
			self._s2_0[i] = integrate(self._k,num.power(self._k,2 + 2*i)*self._P)
		
		# ! careful to the pi factor in numpy sinc(x) = sin(pi*x)/(pi*x) !
		for r in range(size(self._r)):
			self._s2_r[r] = integrate(self._k,num.power(self._k,2 + 2*i)*self._P*num.sinc(self._k*self._r[r]/num.pi))
			self._S2_r[r] = integrate(self._k,num.power(self._k,2 + 2*i)*self._P*Wth(self._k*self._r[r]))
		
		
		
		
		
