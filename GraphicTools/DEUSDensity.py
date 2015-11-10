from DEUSDensityGraphics import*
from DEUSAnalytics import*
from DEUSCosmo import*
from DEUSTools import*

class DEUSDensity(DEUSDensityGraphics,DEUSAnalytics):
	def __init__(self,zinit = 100.):
		DEUSDensityGraphics.__init__(self)
		DEUSAnalytics.__init__(self,zinit)
		
		self._Rg = 0.0
		self._Z = None
		self._S = None
		
		self._P = None
		
		self._dtot = []
		self._dPDF = []
		self._dmin = []
		self._minPDF = []
	
	def _extractInfoFromName(self,file_name):
		arguments = file_name.split('_')
		self._boxlen = int(arguments[0][6:])
		self._npart = int(arguments[1][1:])
		self._cosmo = arguments[2]
		self._Z = arguments[3][1:]
		self._a = 1./(1. + float(self._Z))	
		self._S = arguments[4][1:]
		self._Rg = float(self._boxlen)/float(self._npart)*float(self._S)
		
	def Load(self,file_name = None):
		self._extractInfoFromName(file_name)
		self._LoadSpectrum()
		self._NumericalSmoothing()
		self._GaussianSmoothing(self._Rg)
		self._SetSimu(self._boxlen,self._npart,self._cosmo)
		
		self._computeSigmas0(self._P0)
		
		simu = 'boxlen'+str(self._boxlen)+'_n'+str(self._npart)+'_'+self._cosmo
		self._dtot,self._dPDF = self.loadGlob(simu,self._Z,self._S)
		self._dmin,self._minPDF = self.loadMin(simu,self._Z,self._S)
		self._normalizeData()
	
	def _normalizeData(self):
		self._dPDF /= integrate(self._dtot,self._dPDF)
		self._minPDF /= integrate(self._dmin,self._minPDF)
		
	def PlotPDF(self):
		figure(1)
		grid(True)
		xlabel('$\\delta$')
		ylabel('PDF')
		dtot,Ptot = self.computeLinearDensityPDF(1./(self._a) - 1.) 
		#bar(self._dtot,self._dPDF,width=self._dtot[1]-self._dtot[0],color='b')
		plot(self._dtot,self._dPDF,linestyle = '',marker = 'o',color = 'b',label = 'data')
		plot(dtot-1.,Ptot,linestyle = '-',marker = '+',color = 'g',label = 'gaussian prediction')
		legend()
		
		figure(2)
		grid(True)
		xlabel('$\\delta$')
		ylabel('PDF')
		dtheo_lin, Ndtheo_lin = DEUSAnalytics.computeLinearExtremumDistributionFromHeigh(self,1./self._a -1.)
		dtheo_evo, Ndtheo_evo = DEUSAnalytics.computeEvolvedExtremumDistributionFromHeigh(self,self._w,self._Wm0,1./self._a - 1.)
		dtheo_zeldo, Ndtheo_zeldo = DEUSAnalytics.computeZeldovitchExtremumDistributionFromHeigh(self,self._w,self._Wm0,1./self._a - 1.)

		Ndtheo_lin /= integrate(dtheo_lin - 1.,Ndtheo_lin)
		Ndtheo_evo /= integrate(dtheo_evo - 1.,Ndtheo_evo)
		Ndtheo_zeldo /= integrate(dtheo_zeldo - 1.,Ndtheo_zeldo)
		
		#bar(self._dmin,self._minPDF,width=self._dmin[1]-self._dmin[0],color='b')
		plot(self._dmin,self._minPDF,linestyle = '',marker = 'o',color = 'b',label = 'data')
		plot(dtheo_lin-1.,Ndtheo_lin,linestyle = '-',marker = '+',color = 'g',label = 'gaussian prediction')
		plot(dtheo_evo-1.,Ndtheo_evo,linestyle = '-',marker = '+',color = 'r',label = 'evolved prediction')
		plot(dtheo_zeldo-1.,Ndtheo_zeldo,linestyle = '--',marker = '',color = 'k',label = 'Zeldovitch prediction')
		legend(loc=2)
		show()
		
