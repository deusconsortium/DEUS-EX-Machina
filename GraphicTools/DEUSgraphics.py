from pylab import *
import numpy as num

#~ 
#~ 
#~ 
#~ A toolkit class to analyse data from the ProfileTracer software
#~ 
#~ 
#~ 
#~ 


class DEUSgraphics :
	def __init__(self,boxlen,npart,isOverDensity):
		self._boxlen = boxlen
		self._npart = npart
		self._isOverDensity = isOverDensity
		
		self._dataPath = ''
		
		self._Nprofile = 0
		self._Nradius = 0
		self._r = []
		self._r1 = []
		self._f_tab = None
		self._v_tab = None
	
	def Load(self,file_name = ''):
		
		fileName = self._dataPath + '/' + file_name + '.DEUSprofile'
		File = open(fileName, mode='rb')
		if File is None:
			print 'Error : no file ' + fileName
		
		else
			data = File.read()
			
			#file heading
			(self._Nprofile,self._Nradius,R0,DR) = struct.unpack("iiii", data[:20])
			
			#radius reading
			self._r = struct.unpack("f" * (self._Nradius), fileContent[20:20 + 4*self._Nradius])
			
			#conversion in Mpc/h
			self._r = self._r*boxlen/npart
			
			#initialisation
			self._f_tab = num.zeros((self._Nprofile,size(self._r))
			self._v_tab = num.zeros((self._Nprofile,size(self._r)))
			self._r1 = num.zeros(self._Nprofile)
			
			#reading
			cursor = 20 + 4*self._Nradius
			for i in range(self._Nprofile):
				self._f_tab[i] = struct.unpack("f" * (self._Nradius), fileContent[cursor:cursor + 4*self._Nradius])
				cursor = cursor + cursor + 4*self.Nradius
				self._v_tab[i] = struct.unpack("f" * (self._Nradius), fileContent[cursor:cursor + 4*self._Nradius])
				
				#r1 mass
				self._r1[i] = _getR1(self._r,self._f_tab[i])
	
	## setters and getters
	
	def setDataPath(self,datapath):
		self._dataPath = datapath
	
	def getRadiusTab(self):
		return self._r
	
	def getProfileNumber(self):
		return self._Nprofile
	
	def getMassProfile(self,index):
		if index == -1:
			return self._f_tab
		elif index >= 0 && index < self._Nprofile:
			return self._f_tab[index]
	
	def getSpeedProfile(self,index):
		if index == -1:
			return self._v_tab
		elif index >= 0 && index < self._Nprofile:
			return self._v_tab[index]
	
	def getR1(self,index):
		if index == -1:
			return self._r1
		elif index >= 0 && index < self._Nprofile:
			return self._r1[index]
	
	##plotting functions
	
	def PlotSingleProfile(self,index):
		if index >= 0 and index < self._Nprofile:
			figure(1)
			f = self._f_tab[index]
			v = self._v_tab[index]
			
			subplot(121)
			gris(True)
			xlabel('$r$ in $[Mpc/h]$')
			ylabel('$f(r)$')
			if self._isOverDensity:
				yscale('log')
			
			plot(r,f,linestyle = '--', marker = 'o', color = 'b')
			
			subplot(122)
			gris(True)
			xlabel('$r$ in $[Mpc/h]$')
			ylabel('$v_p(r)$ in ?')
			plot(r,v,linestyle = '--', marker = 'o', color = 'b')
			show()
	
	def PlotMeanProfile(self,R1value,Dr1 = 'dr'):
		if Dr1 == 'dr':
			Dr1 = self._r[2] - self._r[1]
		mask = _mask(self._r1,[R1value,R1value + Dr1])
		mf,sf = getMeanAndSigma(self._f_tab,mask)
		mv,sv = getMeanAndSigma(self._v_tab,mask)
		
		figure(1)
		subplot(121)
		gris(True)
		xlabel('$r$ in $[Mpc/h]$')
		ylabel('$f(r)$')
		if self._isOverDensity:
			yscale('log')
		
		plot(self._r,mf,linestyle = '-', color = 'b')
		fill_between(self._r,mf + sf, mf - sf,color = 'b', alpha=.3)
		
		subplot(122)
		gris(True)
		xlabel('$r$ in $[Mpc/h]$')
		ylabel('$v_p(r)$ in ?')
		plot(self._r,mv,linestyle = '-', color = 'b')
		fill_between(self._r,mv + sv, mv - sv,color = 'b', alpha=.3)
		
		show()
	
	def PlotStatistics(self):
		Nr1 = num.zeros(size(self._r) - 1)
		r1 = num.zeros(size(self._r) - 1)
		for i in range(size(Nr1)):
			Nr1[i] = num.count_nonzero(_mask(self._r1,self._r[i]))
			r1 = 
		
		Nr /= self._r[2] - self._r[1]
		
		figure(1)
		gris(True)
		xlabel('$r_1$ in $[Mpc/h]$')
		ylabel('$\partial n/\partial r_1$')
		
		plot(self)
	
	#private functions
	
	def _mask(masking_tab,limit_values):
		mask_m = (masking_tab >= limit_values[0])
		mask_p = (masking_tab <= limit_values[1])
		return = num.logical_and(mask_m,mask_p)
	
	def getMeanAndSigma(full_tab,mask):
		if shape(full)[0] == size(mask):
			Nr = shape(full)[1]
			mean = num.zeros(size(Nr))
			sigma2 = num.zeros(size(Nr))
			index = []
			for i in range(size(mask)):
				if mask[i]:
					mean += full[i]
					index.append(i)
			mean /= size(index)
			
			for j in range(size(index)):
				sigma2 +=  (full[index[j]] - mean)*(full[index[j]] - mean)
			return mean,sqrt(sigma2)
		else:
			print 'input arrays have not the same lenght : ' + shape(full) + ' vs ' + shape(mask)
			return 0.0,0.0
	
	def _getR1(r,f):
		if _isOverDensity && f[0] >= 1.0:
			for i in range(size(r) -1):
				if f[i + 1] < 1.0 and f[i] >= 1.0
					return r[i + 1] + (f[i + 1] - 1.0)*(r[i] - r[i + 1])/(f[i + 1] - f[i])
		elif !_isOverDensity && f[0] <= 1.0:
			for i in range(size(r) -1):
				if f[i + 1] > 1.0 and f[i] <= 1.0
					return r[i + 1] + (f[i + 1] - 1.0)*(r[i] - r[i + 1])/(f[i + 1] - f[i])
		else:
			return -1.
