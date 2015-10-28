from DEUSgraphics import*
from DEUSanalytics import*

class DEUSmodel(DEUSgraphics,DEUSanalytics):
	
	def __init__(self,zinit = 1100.0,datapath = '../ProfileTracer/data/output/'):
		DEUSgraphics.__init__(self,datapath)
		DEUSanalytics.__init__(self,zinit)
		
		self._zinit = zinit	
		self._Wm0 = None			
		self._do_plot = False
	
	def Load(self,file_name = None):
		"""
			Loading DEUSmodel 
			:param arg1: file_name
			:type arg1: string, None by defaut
		"""
		
		DEUSgraphics.Load(self,file_name)
		founded = DEUSanalytics.Load(self,self._cosmo,self._r)
		
		self._Wm0 = Omega_m_0(self._WmToday,self._w,self._zinit)
		
		print 'cosmological model : \n\tWm0 = '+str(self._Wm0)+'\n\tWm = '+str(self._WmToday)+'\n\tw = '+str(self._w) 
		
		if founded:
			"""
			k = num.copy(self._k)
			P = num.copy(self._P)
			kmin = self._k[0]
			kmax = float(self._npart)*kmin/2.
			Kmax = k[size(k)-1]
			self._k = num.linspace(kmin,Kmax,2*self._npart)
			self._P = num.zeros(size(self._k))
			
			for i in range(size(self._k)):
				self._P[i] = solve(P,k,self._k[i])
			"""
			
			self._Psmooth = num.copy(self._P)
			self._P0 = num.copy(self._P)
			
			#lissage du spectre pour prendre en compte les effets numeriques			
			self.numericalSmoothing(self._boxlen,self._npart)
	
	#overiding functions
	def PlotMeanProfile(self,R1value,Dr1 = 'dr',Rsmooth = 0.0):
		f,v = DEUSgraphics.PlotMeanProfile(self,R1value,Dr1,Rsmooth)
		
		print 'computing initial profile ...'
		r1 = solve(self._r,f,1.0)
		
		if r1 is not None:
			self.localSmoothSpectrum(0.5*r1,'th')
			
			if Rsmooth > 0.0:
				self.localSmoothSpectrum(Rsmooth,'exp')		
			
			#getting height of the peak from d1
			d0 = self._getInitialHeight(self._r,f,r1)
			
			#computing profile
			R,F,V = self.computeMeanProfile(d0,self._Wm0,self._w,r1,self._a)
			
			#plotting
			subplot(211)
			plot(R,F,linestyle = '',marker = '+', color = 'r',label = 'theory')
			legend()
			
			subplot(212)
			#plot(r_t,v_t,linestyle = '-',marker = '+', color = 'r',label = 'theory')
			legend()
			
			self.resetSmoothing()
		show()

	def PlotStatistics(self,Npoints = 0,normalized = True):
		R1,Nr1,Nr1d = DEUSgraphics.PlotStatistics(self,Npoints,normalized)
		
		#R1 = num.add(R1, 0.5*(R1[2]-R1[1]))
		Nth,Nthd = self.computeR1Distribution(R1,normalized = normalized)
		
		subplot(211)
		plot(R1,Nth,linestyle = '-',marker = '+', color = 'r',label = 'theory')
		legend()
		
		subplot(212)
		plot(R1,Nthd,linestyle = '-',marker = '+', color = 'r',label = 'theory')
		legend()
		
		show()
		
	#proper functions
	
	def PlotSpectrum(self):
		figure(1)
		grid(True)
		yscale('log')
		xscale('log')
		plot(self._k,self._P0,linestyle = '-', color = 'b',label = 'linear spectrum')
		plot(self._k,self._P,linestyle = '--', color = 'r',label = 'smoothed spectrum')
		xlabel('$k$ in $h^{-1}.Mpc$')
		ylabel('$P(k)$')
		legend()
		show()
		
	def _getInitialHeight(self,r,f,r1):
		d1init = (self._s2_r1(0,r1)- self._s2_r1(1,r1)*self._S2_r1(0,r1)/self._S2_r1(1,r1))/(self._s2_0[0]*(1. - self._Beta_r1(r1)))
		d = getDensity(r,f) - 1.
		fit = SplineFit(r,d)
		d1 = fit(r1)			
		d10 = d1/(1. + 3.*Eta(self._a,self._w,self._Wm0,self._dlogD_dloga_cmb,self._zinit)*(1. + d1))
		return d10/d1init
