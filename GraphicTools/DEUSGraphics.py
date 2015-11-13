import struct
from matplotlib.colors import LogNorm
from random import randint

from DEUSCosmo import*
from DEUSTools import*

#~ 
#~ 
#~ 
#~ A toolkit class to analyse data from the ProfileTracer software
#~ 
#~ 
#~ 
#~ 

#DEUSgraphics class

class DEUSGraphics(DEUSCosmo) :
	def __init__(self,zinit,data_path = '../ProfileTracer/data/output/'):
		DEUSCosmo.__init__(self,zinit)
		
		self._isOverDensity = True		
		self._dataPath = data_path
		
		self._Nprofile = 0
		self._Nradius = 0
		self._r = []
		self._r1 = []
		self._d0 = []
		self._r1_d = []
		self._r1_full = []
		self._r1d_full = []
		self._selected_profile = []
		self._f_tab = None
		self._v_tab = None
		self._pos_tab = None
		
		#eventually Extremums files
		
		self._do_plot = True
	
	def setDoPlot(self,doPlot = False):
		self._do_plot = doPlot
	
	def Load(self,file_name = None):
		if file_name is None:
			DEUSGraphics.Load(self,selectFile(self._dataPath,12))
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
				self._r1d_full = num.zeros(self._Nprofile)
				self._d0 = num.zeros(self._Nprofile)
				self._selected_profile = np.ones(self._Nprofile)
				self._r1 = []
				self._r1_d = []
				
				#cosmo settings
				self._SetSimu(self._boxlen,self._npart,self._cosmo)
				
				print 'extracting '+str(self._Nprofile)+' profiles ...'
				sys.stdout.flush()
				
				#reading
				cursor0 = 33 + Nc + 4*self._Nradius
				dc = 16 + 8*self._Nradius
				for i in range(self._Nprofile):
					self._pos_tab[i] = num.asarray(struct.unpack("fff", data[cursor0 + i*dc : cursor0 + i*dc + 12]))
					self._d0[i] = (struct.unpack("f",data[cursor0 + i*dc + 12:cursor0 + i*dc + 16]))[0]
					self._f_tab[i] = num.asarray(struct.unpack("f" * (self._Nradius), data[cursor0 + i*dc + 16 : cursor0 + i*dc + 16 + 4*self._Nradius]))
					self._v_tab[i] = num.asarray(struct.unpack("f" * (self._Nradius), data[cursor0 + i*dc + 16 + 4*self._Nradius : cursor0 + i*dc + 16 + 8*self._Nradius]))
					
					#R = float(self._boxlen)/float(self._npart)
					#self._f_tab[i] = gaussianDensitySmoothing(self._r,self._f_tab[i],R)
					
					#r1 mass
					r1 = solve(self._r,self._f_tab[i],1.0)
					d = getDensity(self._r,self._f_tab[i])
					r1_d = solve(self._r,d,1.0)
					self._r1_full[i] = r1
					self._r1d_full[i] = r1_d
					
					self._selected_profile[i] = False
					
					if (r1 is not None) and (r1_d is not None):
						self._r1.append(r1)
						self._r1_d.append(r1_d)
						self._selected_profile[i] = True
				
				self._r1 = num.asarray(self._r1)
				self._r1_d = num.asarray(self._r1_d)
	
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
	
	def getMeanProfile(self,R1value,Dr1 = 'dr',select_density = False,r1d_factor = 0.5):
		if Dr1 == 'dr':
			Dr1 = 0.5*(self._r[1] - self._r[0])
		
		mNone = num.logical_not(ma.masked_invalid(self._r1_full).mask)
		dNone = num.logical_not(ma.masked_invalid(self._r1d_full).mask)
		
		mask_m = num.logical_and(ma.masked_inside(self._r1_full,R1value, R1value + Dr1).mask,mNone)
		
		if select_density:
			r1dm = ma.masked_array(self._r1d_full,num.logical_not(num.logical_and(mask_m,dNone))).mean()
		
			mask_d = num.logical_and(ma.masked_inside(self._r1d_full,max(r1dm-r1d_factor*Dr1,0.0), r1dm + r1d_factor*Dr1).mask,mNone)
		
			mask = num.logical_and(mask_m,mask_d)
		else:
			mask = mask_m
		
		mf,sf,N = self._getMeanAndSigma(self._f_tab,mask)
		mv,sv,N = self._getMeanAndSigma(self._v_tab,mask)
		
		return self._r,mf,sf,mv,sv
	
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
			if self._do_plot:
				fig = figure(1)
			
			"""
			r,f = IncreaseResolution(self._r,mf,Precision)
			r,df = IncreaseResolution(self._r,sf,Precision)
			r,v = IncreaseResolution(self._r,mv,Precision)
			r,dv = IncreaseResolution(self._r,sv,Precision)
			"""
			r,f,df,v,dv = self._r,mf,sf,mv,sv
			
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
	
	def PlotRadiusDistribution(self,N=100):
		figure(1)
		grid(True)	
		xlabel('$r_1$ in $[Mpc/h]$')
		ylabel('$r_1^\\delta$ in $[Mpc/h]$')
		
		xerr = num.zeros(N)
		yerr = num.zeros(N)
		r1m = num.zeros(N)
		r1dm = num.zeros(N)
		r = num.linspace(min(self._r1_d),max(self._r1_d),N+1)
		
		for i in range(N):
			rM = ma.masked_outside(self._r1,r[i],r[i+1])
			r1m[i] = rM.mean()
			xerr[i] = 1.96*rM.std()/float(sqrt(rM.count()))
			rD = ma.masked_array(self._r1_d, rM.mask)
			r1dm[i] = rD.mean()
			yerr[i] = 1.96*rD.std()/float(sqrt(rD.count()))
		
		errorbar(r1m, r1dm, xerr = xerr,yerr=yerr, fmt='o', ecolor='g')
		plot(r1m,2./3.*r1m,linestyle = '--', color = 'r',label = '$r_1^\\delta =2/3r_1$')
		legend()
		
		figure(2)
		grid(True)	
		xlabel('$r_1$ in $[Mpc/h]$')
		ylabel('$r_1^\\delta$ in $[Mpc/h]$')
		hist2d(self._r1, self._r1_d, bins = size(self._r),norm=LogNorm())
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
			if size(index)>0:
				mean /= size(index)
				
				for j in range(size(index)):
					sigma2 =  num.add(sigma2,(full_tab[index[j]] - mean)*(full_tab[index[j]] - mean))
				return mean,sqrt(sigma2/size(index)),size(index)
			else:
				return None,None,0
		else:
			print 'input arrays have not the same lenght : ' + str(shape(full_tab)) + ' vs ' + str(shape(mask))
			return 0.0,0.0,0.0
