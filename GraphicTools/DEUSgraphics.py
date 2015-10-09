from DEUStools import*
import struct
from matplotlib.colors import LogNorm
from random import randint

#~ 
#~ 
#~ 
#~ A toolkit class to analyse data from the ProfileTracer software
#~ 
#~ 
#~ 
#~ 

#DEUSgraphics class

class DEUSgraphics :
	def __init__(self,data_path = '../ProfileTracer/data/output/'):
		
		self._H0 = 0.0
		self._a = 0.0
		self._Wm0 = 0.0
		self._w = -1.0
		self._boxlen = 0
		self._npart = 0
		self._cosmo = None
		self._isOverDensity = True
		
		self._dataPath = data_path
		
		self._Nprofile = 0
		self._Nradius = 0
		self._r = []
		self._r1 = []
		self._r1_d = []
		self._r1_full = []
		self._f_tab = None
		self._v_tab = None
		
		self._do_plot = True
	
	def ListOutput(self):
		folder = os.listdir(self._dataPath)
		print '\nchoose one output to load :'
		for i in range(size(folder)):
			print "- "+str(i)+" : "+folder[i][:-12]
		
		num = int(raw_input("\nenter choice : "))
		if num >= 0 and num < size(folder):
			DEUSgraphics.Load(self,folder[num][:-12])
		else:
			print 'index out of bounds'
	
	def Load(self,file_name = None):
		if file_name is None:
			self.ListOutput()
		else:
			fileName = self._dataPath + '/' + file_name + '.DEUSprofile'
			File = open(fileName, mode='rb')
			if File is None:
				print 'Error : no file ' + fileName
			
			else:
				data = File.read()			
				
				#file heading
				(self._H0,self._a,self._Wm0,Nc) = struct.unpack("fffi", data[:16])
				self._cosmo = ''.join(num.asarray(struct.unpack("c" * Nc, data[16:16+Nc])))
				(self._boxlen,self._npart,self._isOverDensity,self._Nprofile,self._Nradius,R0,DR) = struct.unpack("<ii?iiff", data[16 + Nc:41 + Nc])
				
				#radius reading
				self._r = num.asarray(struct.unpack("f" * (self._Nradius), data[41 + Nc:41 + Nc + 4*self._Nradius]))
				
				#conversion in Mpc/h
				self._r *= float(self._boxlen)
				
				#initialisation
				self._f_tab = num.zeros((self._Nprofile,size(self._r)))
				self._v_tab = num.zeros((self._Nprofile,size(self._r)))
				self._r1_full = num.zeros(self._Nprofile)
				self._r1 = num.empty([1])
				self._r1_d = num.empty([1])
				
				print 'extracting '+str(self._Nprofile)+' profiles ...'
				
				#reading
				cursor0 = 41 + Nc + 4*self._Nradius
				dc = 4*self._Nradius
				for i in range(self._Nprofile):
					self._f_tab[i] = num.asarray(struct.unpack("f" * (self._Nradius), data[cursor0 + 2*i*dc : cursor0 + (2*i+1)*dc]))
					self._v_tab[i] = num.asarray(struct.unpack("f" * (self._Nradius), data[cursor0 + (2*i+1)*dc : cursor0 + 2*(i+1)*dc]))
					
					#r1 mass
					r1 = solve(self._r,self._f_tab[i],1.)
					d = getDensity(self._r,self._f_tab[i])
					r1_d = solve(self._r,d,1.0)
					self._r1_full[i] = r1
					
					if (r1 is not None) and (r1_d is not None):
						self._r1 = num.append(self._r1,r1)
						self._r1_d = num.append(self._r1_d,r1_d)
	
	def Save(self):
		return True
	
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
	
	def PlotSingleProfile(self,index = None):
		if index is None:
			self.PlotSingleProfile(randint(0,self._Nprofile - 1))
		elif index >= 0 and index < self._Nprofile:
			fig = figure(1)
			
			f = self._f_tab[index]
			v = self._v_tab[index]
			r = self._r
			d = f + r*derivative(r,f)/3.
			
			subplot(211)
			grid(True)
			xlabel('$r$ in $[Mpc/h]$')
			if self._isOverDensity:
				yscale('log')
			
			plot(r,f,linestyle = '-', color = 'b',label = '$f(r)$')
			plot(r,d,linestyle = '--', color = 'b',label = '$\\delta(r) + 1$')
			legend()
			
			subplot(212)
			grid(True)
			xlabel('$r$ in $[Mpc/h]$')
			ylabel('$v_p(r)$ in ?')
			plot(r,v,linestyle = '--', marker = 'o', color = 'b')
			
			show()
			
			return fig
	
	def PlotMeanProfile(self,R1value,Dr1 = 'dr',Rsmooth = 0.0):
		if Dr1 == 'dr':
			Dr1 = self._r[2] - self._r[1]
		mask = ma.masked_outside(self._r1_full,R1value, R1value + Dr1).mask
		mf,sf,N = self._getMeanAndSigma(self._f_tab,mask)
		mv,sv,N = self._getMeanAndSigma(self._v_tab,mask)
		
		print str(N) + ' profiles have been selected'
		
		if Rsmooth > 0.0:
			Mf = gaussianDensitySmoothing(self._r,mf,Rsmooth)
			Mfp = gaussianDensitySmoothing(self._r,mf+1.96*sf/sqrt(N),Rsmooth)
			Mfm = gaussianDensitySmoothing(self._r,mf-1.96*sf/sqrt(N),Rsmooth)
			mf = Mf
			sf = 0.5*num.fabs(Mfp - Mfm)/1.96*sqrt(N)
		
		if N > 0:
			fig = figure(1)
			
			subplot(211)
			grid(True)
			xlabel('$r$ in $[Mpc/h]$')
			ylabel('$f(r)$')
			if self._isOverDensity:
				yscale('log')
			
			plot(self._r,mf,linestyle = '-', color = 'b',label = 'data')
			fill_between(self._r,mf - 1.96*sf/sqrt(N), mf + 1.96*sf/sqrt(N),color = 'b', alpha=.3)
			
			subplot(212)
			grid(True)
			xlabel('$r$ in $[Mpc/h]$')
			ylabel('$v_p(r)$ in ?')
			plot(self._r,mv,linestyle = '-', color = 'b',label = 'data')
			fill_between(self._r,mv - 1.96*sv/sqrt(N), mv + 1.96*sv/sqrt(N),color = 'b', alpha=.3)
			
			if self._do_plot:
				show()
			
			return mf,mv
	
	def PlotStatistics(self,Npoints = 0,normalized = False):
		if Npoints is 0:
			Npoints = size(self._r) - 1
		Nr1 = num.zeros(Npoints)
		Nr1d = num.zeros(Npoints)
		r1 = num.zeros(Npoints)
		dr = (self._r[size(self._r) - 2] - self._r[0])/(Npoints - 1)
		for i in range(Npoints):
			r1[i] = self._r[0] + i*dr
			Nr1[i] = num.count_nonzero(self._mask(self._r1,[r1[i],r1[i] + dr]))
			Nr1d[i] = num.count_nonzero(self._mask(self._r1_d,[r1[i],r1[i] + dr]))
		
		Nr1 /= dr*(self._boxlen**3.)
		Nr1d /= dr*(self._boxlen**3.)
		
		if normalized:
			fact = 0.0
			factd = 0.0
			for i in range(size(Nr1)):
				fact += Nr1[i]*dr
			for i in range(size(Nr1)):
				factd += Nr1d[i]*dr
			
			Nr1 /= fact
			Nr1d /= factd
		
		fig = figure(1)
		
		subplot(211)
		grid(True)
		xlabel('$r_1$ in $[Mpc/h]$',fontsize = 20)
		ylabel('$\partial n/\partial r_1$ in $[Mpc/h]^{-4}$',fontsize = 20)
		bar(r1,Nr1,width=dr,color='b')
		subplot(212)
		grid(True)
		xlabel('$r_1^\\delta$ in $[Mpc/h]$',fontsize = 20)
		ylabel('$\partial n/\partial r_1^\\delta$ in $[Mpc/h]^{-4}$',fontsize = 20)
		bar(r1,Nr1d,width=dr,color='b')
		
		#~ xlabel('$r_1$ in $[Mpc/h]$',fontsize = 20)
		#~ ylabel('log scaled',fontsize = 20)
		#~ bar(r1,Nr1,width=dr,color='b',log = True)
		#~ xscale('log')
		
		if self._do_plot:
			show()
		
		return r1,Nr1,Nr1d
	
	def PlotRadiusDistribution(self):
		figure(1)
		grid(True)	
		xlabel('$r_1$ in $[Mpc/h]$')
		ylabel('$r_1^\\delta$ in $[Mpc/h]$')
		
		xerr = num.zeros(size(self._r)-1)
		yerr = num.zeros(size(xerr))
		r1m = num.zeros(size(xerr))
		r1dm = num.zeros(size(xerr))
		for i in range(size(xerr)):
			rM = ma.masked_outside(self._r1,self._r[i], self._r[i+1])
			r1m[i] = rM.mean()
			xerr[i] = rM.std()
			rD = ma.masked_array(self._r1_d, rM.mask)
			r1dm[i] = rD.mean()
			yerr[i] = rD.std()
		
		errorbar(r1m, r1dm, xerr = xerr,yerr=yerr, fmt='o', ecolor='g')
		show()
	
	#private functions
	
	def _doPlot(self,doplot):
		self._do_plot = doplot
 	
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
			print 'input arrays have not the same lenght : ' + str(shape(full_tab)) + ' vs ' + str(shape(mask))
			return 0.0,0.0,0.0
