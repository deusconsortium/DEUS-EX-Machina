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
		self._WmToday = 0.0
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
		self._d0 = []
		self._r1_d = []
		self._r1_full = []
		self._selected_profile = []
		self._f_tab = None
		self._v_tab = None
		self._pos_tab = None
		
		#eventually Extremums files
		
		#self._FOFextremaPath = "../../../../jpasdeloup/DEUS-EX-Machina/compute_extrema/data/"
		self._FOFextremaPath = "../compute_extrema/data/"
		self._nFOFextrema = 0
		self._density = num.array([])
		self._seuil = num.array([])
		self._mean_seuil = num.array([])
		
		self._do_plot = True
	
	def ListOutput(self):
		folder = os.listdir(self._dataPath)
		print '\nchoose one output to load :'
		for i in range(size(folder)):
			print "- "+str(i)+" : "+folder[i][:-12]
		
		num = int(raw_input("\nenter choice : "))
		if num >= 0 and num < size(folder):
			DEUSgraphics.LoadProfiles(self,folder[num][:-12])
		else:
			print 'index out of bounds'
			
	def ListFOFextrema(self):
		folder = listdirHidden(self._FOFextremaPath)
		print '\nchoose one FOFextrema file  to load :'
		for i in range(size(folder)):
			print "- "+str(i)+" : "+folder[i]
		
		num = int(raw_input("\nenter choice : "))
		if num >= 0 and num < size(folder):
			possibilities = filesThatEndAs(self._FOFextremaPath + folder[num]+ "/",".deus_extrema",18)
			if size(possibilities) > 1:
				for i in range(size(possibilities)):
					print "\t- "+str(i)+" : "+str(possibilities[i])
				simu = int(raw_input("\nenter choice : "))
				if simu >= 0 and simu < size(possibilities):
					DEUSgraphics.LoadFOFextrema(self,folder[num] + "/" + possibilities[simu])
			else:
				DEUSgraphics.LoadFOFextrema(self,folder[num] + "/" + possibilities[0])
		else:
			print 'index out of bounds'
	
	def _LoadSingleFOFextrema(self,filename):
		try:
			File = open(filename, mode='rb')
		except IOError :
			return None,None,None,None
			
		data = File.read()
		
		#ignore first octet
		Nextr = (struct.unpack(">i",data[4:8]))[0]
		x = (struct.unpack(">i",data[8:12]))[0]
		x = (struct.unpack(">i",data[12:16]))[0]
		density = num.zeros(Nextr)
		seuil = num.zeros(Nextr)
		mean_seuil = num.zeros(Nextr)
		
		print 'extracting '+str(Nextr)+' extremas'
		
		for i in range(Nextr):
			(x,y,z,d,s,ms) = struct.unpack(">ffffff",data[16 + 24*i:16 + 24*(i+1)])
			density[i] = d
			seuil[i] = s
			mean_seuil[i] = ms
			k = randint(0,1000)
			if k >= 999:
				print '-'
				print x,y,z,d,s,ms
		
		return Nextr,density,seuil,mean_seuil
				
	
	def LoadFOFextrema(self,file_name = None,nProc = None):
		if file_name is None:
			self.ListFOFextrema()
		else:
			root_name = self._FOFextremaPath  + file_name
			
			#getting file name
			i = 0
			while i is not None:
				name = root_name + str(i).zfill(5) + ".deus_extrema"
				n,d,s,ms = self._LoadSingleFOFextrema(name)				
				if n is not None:
					self._nFOFextrema += n
					self._density = num.concatenate((self._density,d))
					self._seuil = num.concatenate((self._seuil,s))
					self._mean_seuil = num.concatenate((self._mean_seuil,ms))
					i += 1
					print 'file '+str(i)+' done'
				else:
					i = None
	
	def LoadProfiles(self,file_name = None):
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
				(self._H0,self._a,self._WmToday,Nc) = struct.unpack("fffi", data[:16])
				self._cosmo = ''.join(num.asarray(struct.unpack("c" * Nc, data[16:16+Nc])))
				
				(self._boxlen,self._npart,self._isOverDensity,self._Nprofile,self._Nradius) = struct.unpack("<ii?ii", data[16 + Nc:33 + Nc])
				
				#radius reading
				self._r = num.asarray(struct.unpack("f" * (self._Nradius), data[33 + Nc:33 + Nc + 4*self._Nradius]))
				
				#conversion in Mpc/h
				self._r *= float(self._boxlen)
				
				#initialisation
				self._f_tab = num.zeros((self._Nprofile,size(self._r)))
				self._v_tab = num.zeros((self._Nprofile,size(self._r)))
				self._pos_tab = num.zeros((self._Nprofile,3))
				self._r1_full = num.zeros(self._Nprofile)
				self._d0 = num.zeros(self._Nprofile)
				self._selected_profile = np.ones(self._Nprofile)
				self._r1 = []
				self._r1_d = []
				
				print 'extracting '+str(self._Nprofile)+' profiles ...'
				#print 'CIC smoothing for each profile on R = '+str(0.62035*float(self._boxlen)/float(self._npart))+' [Mpc/h]'
				sys.stdout.flush()
				
				#reading
				cursor0 = 33 + Nc + 4*self._Nradius
				dc = 16 + 8*self._Nradius
				for i in range(self._Nprofile):
					self._pos_tab[i] = num.asarray(struct.unpack("fff", data[cursor0 + i*dc : cursor0 + i*dc + 12]))
					self._d0[i] = (struct.unpack("f",data[cursor0 + i*dc + 12:cursor0 + i*dc + 16]))[0]
					self._f_tab[i] = num.asarray(struct.unpack("f" * (self._Nradius), data[cursor0 + i*dc + 16 : cursor0 + i*dc + 16 + 4*self._Nradius]))
					self._v_tab[i] = num.asarray(struct.unpack("f" * (self._Nradius), data[cursor0 + i*dc + 16 + 4*self._Nradius : cursor0 + i*dc + 16 + 8*self._Nradius]))
					
					#CIC smoothing
					#self._f_tab[i] = CICDensitySmoothing(self._r,self._f_tab[i],0.62035*float(self._boxlen)/float(self._npart))
					
					#getting d0 for each profile
					Rcic = 0.62035*float(self._boxlen)/float(self._npart)
					#self._d0[i] = getCICCentralDensity(self._r[:10],self._f_tab[i][:10],Rcic)
					
					#r1 mass
					r1 = solve(self._r,self._f_tab[i],1.)
					d = getDensity(self._r,self._f_tab[i])
					r1_d = solve(self._r,d,1.0)
					self._r1_full[i] = r1
					
					self._selected_profile[i] = False
					
					if (r1 is not None) and (r1_d is not None):
						self._r1.append(r1)
						self._r1_d.append(r1_d)
						self._selected_profile[i] = True
				
				self._r1 = num.asarray(self._r1)
				self._r1_d = num.asarray(self._r1_d)
				
				print 'Warning : w = -1'
	
	#save compact data
	def Save(self,filename = ''):
		if filename is not None:
			name = filename
		else:
			name = 'box'+str(self._boxlen)+'_n'+str(self._npart)+self._cosmo
		
		print 'saving file '+name+'.txt ...'
		sys.stdout.flush()
		my_file = open(name+".txt", "w")
		k = 0
		for i in range(self._Nprofile):
			if self._selected_profile[i]:
				my_file.write(str(self._pos_tab[i][0])+"\t")
				my_file.write(str(self._pos_tab[i][1])+"\t")
				my_file.write(str(self._pos_tab[i][2])+"\t")
				my_file.write(str(self._r1[k])+"\t")				
				my_file.write(str(self._r1_d[k])+"\t")
				
				d = getDensity(self._r,self._f_tab[i])
				dFit = SplineFit(self._r,d)
				my_file.write(str(dFit(self._r1[k]))+"\t")
				
				fFit = SplineFit(self._r,self._f_tab[i])
				my_file.write(str(fFit(self._r1_d[k]))+"\n")
				k = k+1
		
		my_file.close()
		print 'saving done'
		
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
			d = getDensity(r,f)
			
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
	
	def PlotMeanProfile(self,R1value,Dr1 = 'dr',Rsmooth = 0.0,Precision = 2):
		if Dr1 == 'dr':
			Dr1 = self._r[2] - self._r[1]
		mask = ma.masked_inside(self._r1_full,R1value, R1value + Dr1).mask
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
			
			r,f = IncreaseResolution(self._r,mf,Precision)
			r,df = IncreaseResolution(self._r,sf,Precision)
			r,v = IncreaseResolution(self._r,mv,Precision)
			r,dv = IncreaseResolution(self._r,sv,Precision)
			
			subplot(211)
			grid(True)
			xlabel('$r$ in $[Mpc/h]$')
			ylabel('$f(r)$')
			if self._isOverDensity:
				yscale('log')
			
			plot(r,f,linestyle = '-', color = 'b',label = 'data')
			fill_between(r,f - 1.96*df/sqrt(N), f + 1.96*df/sqrt(N),color = 'b', alpha=.3)
			plot(r,getDensity(r,f),linestyle = '--', color = 'b',label = 'density')
			
			subplot(212)
			grid(True)
			xlabel('$r$ in $[Mpc/h]$')
			ylabel('$v_p(r)$ in ?')
			plot(r,v,linestyle = '-', color = 'b',label = 'data')
			fill_between(r,v - 1.96*dv/sqrt(N), v + 1.96*dv/sqrt(N),color = 'b', alpha=.3)
			
			if self._do_plot:
				show()
			else:
				return mf,mv
	
	def PlotHeightStatistics(self,Npoints = 100,normalized = False):
		if self._nFOFextrema is 0:
			print 'non FOFextrema file loaded. Call LoadFOFextrema first'
			return None
		else:
			mask = self._mask(self._seuil,[0.0,None])
			if min(self._density < 0.0):
				density = self._density[mask] + 1.
			else:
				density = self._density[mask]
			
			print 'selected '+str(num.count_nonzero(mask))+' minimums : '+str(100.*num.count_nonzero(mask)/float(self._nFOFextrema))+' % of initial densities'
			
			d0_tab = num.linspace(min(density),2.0,Npoints)
			Dd = 0.5*(d0_tab[1] - d0_tab[0])
			Nd = num.zeros(Npoints)
			
			for i in range(Npoints):
				Nd[i] = num.count_nonzero(self._mask(density,[d0_tab[i] - Dd,d0_tab[i] + Dd]))
			
			if normalized:
				Nd /= num.count_nonzero(mask)*2.*Dd
			
			figure(1)
			grid(True)
			xlabel('$\\delta\\rho$',fontsize = 20)
			
			bar(d0_tab,Nd,width=2.*Dd,color='b')
			
			if self._do_plot:
				show()
			
			else:
				return d0_tab,Nd
	
	def PlotRadiusStatistics(self,Npoints = None,normalized = False):
		if Npoints is None:
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
		
		N = size(self._r) - 1
		xerr = num.zeros(N)
		yerr = num.zeros(N)
		r1m = num.zeros(N)
		r1dm = num.zeros(N)
		
		for i in range(N):
			rM = ma.masked_outside(self._r1,self._r[i], self._r[i+1])
			r1m[i] = rM.mean()
			xerr[i] = 1.96*rM.std()/float(sqrt(rM.count()))
			rD = ma.masked_array(self._r1_d, rM.mask)
			r1dm[i] = rD.mean()
			yerr[i] = 1.96*rD.std()/float(sqrt(rM.count()))
		
		errorbar(r1m, r1dm, xerr = xerr,yerr=yerr, fmt='o', ecolor='g')
		plot(r1m,2./3.*r1m,linestyle = '--', color = 'r',label = '$r_1^\\delta =2/3r_1$')
		legend()
		
		figure(2)
		grid(True)	
		xlabel('$r_1$ in $[Mpc/h]$')
		ylabel('$r_1^\\delta$ in $[Mpc/h]$')
		hist2d(self._r1, self._r1_d, bins = size(self._r))
		colorbar()
		show()
	
	#private functions
	
	def _doPlot(self,doplot):
		self._do_plot = doplot
 	
	def _mask(self,masking_tab,limit_values):
		if limit_values[0] is None:
			return (masking_tab <= limit_values[1])
		elif limit_values[1] is None:
			return (masking_tab >= limit_values[0])
		else:
			mask_m = (masking_tab >= min(limit_values))
			mask_p = (masking_tab <= max(limit_values))
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
