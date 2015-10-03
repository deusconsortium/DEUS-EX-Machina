from pylab import *
import numpy as num
import struct
import os

#~ 
#~ 
#~ 
#~ A toolkit class to analyse data from the ProfileTracer software
#~ 
#~ 
#~ 
#~ 


class DEUSgraphics :
	def __init__(self,isOverDensity):
		self._boxlen = 0
		self._npart = 0
		self._isOverDensity = isOverDensity
		
		self._dataPath = ''
		
		self._Nprofile = 0
		self._Nradius = 0
		self._r = []
		self._r1 = []
		self._f_tab = None
		self._v_tab = None
	
	def __init__(self,isOverDensity,data_path):
		self._boxlen = 0
		self._npart = 0
		self._isOverDensity = isOverDensity
		
		self._dataPath = data_path
		
		self._Nprofile = 0
		self._Nradius = 0
		self._r = []
		self._r1 = []
		self._f_tab = None
		self._v_tab = None
	
	def ListOutput(self):
		folder = os.listdir(self._dataPath)
		print '\nchoose one output to load :'
		for i in range(size(folder)):
			print "- "+str(i)+" : "+folder[i][:-12]
		
		num = int(raw_input("\nenter choice : "))
		if num >= 0 and num < size(folder):
			self.Load(folder[num][:-12])
		else:
			print 'index out of bounds'
	
	def Load(self,file_name):
		fileName = self._dataPath + '/' + file_name + '.DEUSprofile'
		File = open(fileName, mode='rb')
		if File is None:
			print 'Error : no file ' + fileName
		
		else:
			data = File.read()			
			
			#file heading
			(self._boxlen,self._npart,self._Nprofile,self._Nradius,R0,DR) = struct.unpack("iiiiff", data[:24])
			
			#radius reading
			self._r = num.asarray(struct.unpack("f" * (self._Nradius), data[24:24 + 4*self._Nradius]))
			
			#conversion in Mpc/h
			self._r *= float(self._boxlen)
			
			#initialisation
			self._f_tab = num.zeros((self._Nprofile,size(self._r)))
			self._v_tab = num.zeros((self._Nprofile,size(self._r)))
			self._r1 = num.zeros(self._Nprofile)
			
			print 'extracting '+str(self._Nprofile)+' profiles ...'
			
			#reading
			cursor0 = 24 + 4*self._Nradius
			dc = 4*self._Nradius
			for i in range(self._Nprofile):
				#print str(cursor0 + 2*i*dc) + ' , ' + str(cursor0 + (2*i+1)*dc)
				self._f_tab[i] = num.asarray(struct.unpack("f" * (self._Nradius), data[cursor0 + 2*i*dc : cursor0 + (2*i+1)*dc]))
				#print str(cursor0 + (2*i+1)*dc) + ' , ' + str(cursor0 + (2*i+2)*dc)
				self._v_tab[i] = num.asarray(struct.unpack("f" * (self._Nradius), data[cursor0 + (2*i+1)*dc : cursor0 + 2*(i+1)*dc]))
				
				#r1 mass
				self._r1[i] = self._getR1(self._r,self._f_tab[i])
	
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
		elif index >= 0 and index < self._Nprofile:
			return self._f_tab[index]
	
	def getSpeedProfile(self,index):
		if index == -1:
			return self._v_tab
		elif index >= 0 and index < self._Nprofile:
			return self._v_tab[index]
	
	def getR1tab(self,index):
		if index == -1:
			return self._r1
		elif index >= 0 and index < self._Nprofile:
			return self._r1[index]
	
	##plotting functions
	
	def PlotSingleProfile(self,index):
		if index >= 0 and index < self._Nprofile:
			figure(1)
			f = self._f_tab[index]
			v = self._v_tab[index]
			r = self._r
			d = f - r*self._derivative(r,f)/3.
			
			subplot(121)
			grid(True)
			xlabel('$r$ in $[Mpc/h]$')
			if self._isOverDensity:
				yscale('log')
			
			plot(r,f,linestyle = '-', color = 'b',label = '$f(r)$')
			plot(r,d,linestyle = '--', color = 'b',label = '$\\delta(r) + 1$')
			legend()
			
			subplot(122)
			grid(True)
			xlabel('$r$ in $[Mpc/h]$')
			ylabel('$v_p(r)$ in ?')
			plot(r,v,linestyle = '--', marker = 'o', color = 'b')
			show()
	
	def PlotMeanProfile(self,R1value,Dr1 = 'dr'):
		if Dr1 == 'dr':
			Dr1 = self._r[2] - self._r[1]
		mask = self._mask(self._r1,[R1value,R1value + Dr1])
		mf,sf,N = self._getMeanAndSigma(self._f_tab,mask)
		mv,sv,N = self._getMeanAndSigma(self._v_tab,mask)
		
		print str(N) + ' profiles have been selected'
		
		figure(1)
		subplot(121)
		grid(True)
		xlabel('$r$ in $[Mpc/h]$')
		ylabel('$f(r)$')
		if self._isOverDensity:
			yscale('log')
		
		plot(self._r,mf,linestyle = '-', color = 'b')
		fill_between(self._r,mf - 1.96*sf/sqrt(N), mf + 1.96*sf/sqrt(N),color = 'b', alpha=.3)
		
		subplot(122)
		grid(True)
		xlabel('$r$ in $[Mpc/h]$')
		ylabel('$v_p(r)$ in ?')
		plot(self._r,mv,linestyle = '-', color = 'b')
		fill_between(self._r,mv - 1.96*sv/sqrt(N), mv + 1.96*sv/sqrt(N),color = 'b', alpha=.3)
		
		show()
	
	def PlotStatistics(self,Npoints = 0):
		if Npoints is 0:
			Npoints = size(self._r) - 1
		Nr1 = num.zeros(Npoints)
		r1 = num.zeros(Npoints)
		dr = (self._r[size(self._r) - 2] - self._r[0])/(Npoints - 1)
		for i in range(Npoints):
			r1[i] = self._r[0] + i*dr
			Nr1[i] = num.count_nonzero(self._mask(self._r1,[r1[i],r1[i] + dr]))
		
		Nr1 /= dr*(self._boxlen**3.)
		
		figure(1)
		subplot(121)
		grid(True)
		xlabel('$r_1$ in $[Mpc/h]$',fontsize = 20)
		ylabel('$\partial n/\partial r_1$ in $[Mpc/h]^{-4}$',fontsize = 20)
		bar(r1,Nr1,width=dr,color='b')
		subplot(122)
		grid(True)
		xlabel('$r_1$ in $[Mpc/h]$',fontsize = 20)
		ylabel('log scaled',fontsize = 20)
		bar(r1,Nr1,width=dr,color='b',log = True)
		show()		
	
	#private functions
	
	def _derivative(self,x,y):
		dy = num.zeros(size(y))
		
		for i in range(size(x) - 2):
			dy[i + 1] = (y[i + 2] - y[i])/(x[i + 2] - x[i])
		
		dy[0] = (y[1] - y[0])/(x[1]-x[0])
		dy[size(dy)-1] = (y[size(y) - 1] -  y[size(y) - 2])/(x[size(y) - 1] - x[size(y) - 2])
		return dy
 	
	def _mask(self,masking_tab,limit_values):
		mask_m = (masking_tab >= limit_values[0])
		mask_p = (masking_tab <= limit_values[1])
		return num.logical_and(mask_m,mask_p)
	
	def _getMeanAndSigma(self,full_tab,mask):
		if num.shape(full_tab)[0] == size(mask):
			Nr = num.shape(full_tab)[1]
			mean = num.zeros(size(Nr))
			sigma2 = num.zeros(size(Nr))
			index = []
			for i in range(size(mask)):
				if mask[i]:
					mean = num.add(mean,full_tab[i])
					index.append(i)
			mean /= size(index)
			
			for j in range(size(index)):
				sigma2 =  num.add(sigma2,(full_tab[index[j]] - mean)*(full_tab[index[j]] - mean))
			return mean,sqrt(sigma2/size(index)),size(index)
		else:
			print 'input arrays have not the same lenght : ' + shape(full) + ' vs ' + shape(mask)
			return 0.0,0.0,0.0
	
	def _getR1(self,r,f):
		if self._isOverDensity and f[0] >= 1.0:
			for i in range(size(r) -1):
				if f[i + 1] < 1.0 and f[i] >= 1.0:
					return r[i + 1] + (f[i + 1] - 1.0)*(r[i] - r[i + 1])/(f[i + 1] - f[i])
		elif not self._isOverDensity and f[0] <= 1.0:
			for i in range(size(r) -1):
				if f[i + 1] > 1.0 and f[i] <= 1.0:
					return r[i + 1] + (f[i + 1] - 1.0)*(r[i] - r[i + 1])/(f[i + 1] - f[i])
		return -1.
