from DEUSgraphics import*

class DEUSmodel(DEUSgraphics):
	
	def __init__(self,datapath = '../ProfileTracer/data/output/'):
		DEUSgraphics.__init__(self,datapath)
		
		self._power_spectrum_path = '/efiler2/bingo_save/Babel/data/'
		
		self._zcmb = 1100.0
		self._dlogD_dloga_cmb = 1.0
		self._cosmo = ''
		self._P = num.empty([1])
		self._k = num.empty([1])
		
		self._s2_0 = num.zeros(3)
		self._s2_r = None
		self._S2_r = None
		
		self._do_plot = False
	
	def Load(self,file_name):
		DEUSgraphics.Load(self,file_name)
		
		self._s2_r = num.zeros((2,size(self._r)))
		self._S2_r = num.zeros((2,size(self._r)))
		
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
		for i in range(2):
			for r in range(size(self._r)):
				self._s2_r[i][r] = integrate(self._k,num.power(self._k,2 + 2*i)*self._P*num.sinc(self._k*self._r[r]/num.pi))
				self._S2_r[i][r] = integrate(self._k,num.power(self._k,2 + 2*i)*self._P*Wth(self._k*self._r[r]))
		
		
	#overiding functions
	def PlotMeanProfile(self,R1value,Dr1 = 'dr'):
		r = self._r
		fig,f,v = DEUSgraphics.PlotMeanProfile(R1value,Dr1)
		r1 = solve(r,f,1.0)
		S10 = solve(self._s2_r[0],self_r,r1)
		S11 = solve(self._s2_r[1],self_r,r1)
		
		D_on_D0 = (self._S2_r[0]/self._s2_0[0] - self._S2_r[1]/self._s2_0[1]*S10/S11)
		
		#getting height of the peak from d1
		d = f - 1.0 - r*derivative(r,f)/3.
		
		fit = SplineFit(r,d)
		d1 = fit(r1)
			
		d10 = d1/(1. + 3.*Eta(a,self._w,self._Wm0,self._dlogD_dloga_cmb)*(1. + d1))
		d0 = D_on_D0 - r*derivative(r,D_on_D0)/3.
		
		fit2 = SplineFit(r,d0)
		d1cmb = fit2(r1)		
		fcmb = 1.0 + d10/d1cmb*D_on_D0
		
		#evolving profile
		r_t,f_t,v_t = evolveProfile(r,fcmb,af,self._w,self._Wm0,self._zcmb,self._dlogD_dloga_cmb)
		
		#plotting
		fig.subplot(211)
		plot(r_t,f_t,linestyle = '-',marker = '+', color = 'r',label = 'theory')
		
		show()
		
