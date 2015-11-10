from DEUSGraphics import*
from DEUSAnalytics import*

class DEUSProfile(DEUSGraphics,DEUSAnalytics):
	
	def __init__(self,zinit = 50.0,datapath = '../ProfileTracer/data/output/'):
		DEUSGraphics.__init__(self,zinit,datapath)
		DEUSAnalytics.__init__(self,zinit)
			
		self._do_plot = False
	
	def Load(self,file_name = None):
		"""
			Loading DEUSmodel 
			:param arg1: file_name
			:type arg1: string, None by defaut
		"""
		
		DEUSGraphics.Load(self,file_name)		
		if DEUSAnalytics.Load(self,self._cosmo,self._r):
			self._NumericalSmoothing()
			self._computeSigmas()
	
	#overiding functions
	def PlotMeanProfile(self,R1value,Dr1 = 'dr',Rsmooth = 0.0):
		f,v = DEUSGraphics.PlotMeanProfile(self,R1value,Dr1,Rsmooth)
		
		print 'computing initial profile ...'
		r1 = solve(self._r,f,1.0)
		
		if r1 is not None:
			self.localSmoothSpectrum(0.5*r1,'th')
			#self.localSmoothSpectrum(r1,'th')
			
			if Rsmooth > 0.0:
				self.localSmoothSpectrum(Rsmooth,'exp')	
				
			#smoothing with CIC kernel
			R0 = 0.62035*float(self._boxlen)/float(self._npart)
			fcic = CICDensitySmoothing(self._r,f,R0)
			subplot(211)
			plot(self._r,fcic,linestyle = '--',marker = '', color = 'g',label = 'CIC')
				
			
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

	def PlotRadiusStatistics(self,v0 = 0.0,Npoints = None,normalized = True):
		R1,Nr1,Nr1d = DEUSGraphics.PlotRadiusStatistics(self,Npoints,normalized)
		Nth,Nthd = self.computeR1Distribution(R1,v0,normalized = normalized)
		
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
		xlabel('$k$ in $h^{-1}.Mpc$')
		ylabel('$P(k)$')
		legend()
		show()
		
	def _getInitialHeight(self,r,f,r1):
		d1init = (self._s2_r1(0,r1)- self._s2_r1(1,r1)*self._S2_r1(0,r1)/self._S2_r1(1,r1))/(self._s2_0[0]*(1. - self._Beta_r1(r1)))
		d = getDensity(r,f) - 1.
		fit = SplineFit(r,d)
		d1 = fit(r1)			
		d10 = d1/(1. + 3.*Eta(self._a,self._w,self._Wm0,self._dlogD_dloga_init,self._zinit)*(1. + d1))
		return d10/d1init
