__author__ = 'jpasdeloup'
import numpy as np
import matplotlib.pyplot as plt

class DEUSExtremaGraphics :

    def __init__(self):
        self.dir = "../compute_extrema/data/"
        self.nb_histo = 1000
        self.glob_start = -1.0
        self.glob_end = 2.0
        self.min_start = -1.0
        self.min_end = 0.0

    def loadGlob(self, simu, Z, S):
        file_root = self._fileRoot(simu, Z, S)
        data = np.loadtxt(self.dir + file_root + '/' + file_root + '_all.deus_histo.txt')
        data = np.delete(data, self.nb_histo-1)        
        densityscale = np.linspace(self.glob_start,self.glob_end,self.nb_histo,0)
        densityscale = np.delete(densityscale, self.nb_histo-1)
        delta = densityscale[2] - densityscale[1]
        return densityscale, data, delta

    def loadMin(self, simu, Z, S):
        file_root = self._fileRoot(simu, Z, S)
        data = np.loadtxt(self.dir + file_root + '/' + file_root + '_min.deus_histo.txt')        
        densityscale = np.linspace(self.min_start,self.min_end,self.nb_histo,0)        
        delta = densityscale[2] - densityscale[1]
        return densityscale, data, delta

    @staticmethod
    def _fileRoot(simu, Z, S):
        return simu + "_Z" + Z + "_S" + S

    def showDensity(self):
        plt.figure(1)

        simu = ['boxlen648_n1024_lcdmw5'] # ['boxlen2592_n1024_lcdmw5','boxlen648_n1024_rpcdmw5','boxlen648_n1024_lcdmw5']
        S = ["1","2"] #["2"] # ["1","2","3"]
        Z = ["0"] #["56","18","8","4","2","1","0"] # ["81","19","9","4","2","1","0"] ["93","36","17","9","4","2","1","0"] 
        
        # GLOBAL
        plt.subplot(211)
        plt.grid(True)
        
        for mySimu in simu:
            for myZ in Z:
                for myS in S:
                    
                    myLegend =  " "+ mySimu + " Z="+ myZ #" S="+ myS
                    myLabel = "S=" + myS
                    #myLegend =  " Z="+ myZ + " S="+ myS
                    #myLabel = mySimu
                    
                    # Load data
                    densityscale,data,delta = self.loadGlob(mySimu,myZ,myS)                                        
                    
                    # If histo
                    #plt.bar(densityscale, data, delta, label = myLabel)
                    
                    # If points
                    densityscale += delta/2  # move point to middle of area
                    plt.plot(densityscale, data, label = myLabel, marker = '.')                                      
        
                    # Legend
                    plt.legend(title="Total density" + myLegend, loc='upper right', ncol=3) 

                    # Limits
                    plt.xlim([-1.0,2.0])
                    #plt.ylim([0.0,0.55])

        # MINS
        plt.subplot(212)
        plt.grid(True)
        
        for mySimu in simu:
            for myZ in Z:
                for myS in S:
                    
                    myLegend =  " "+ mySimu + " Z="+ myZ #" S="+ myS
                    myLabel = "S=" + myS
                    
                    # Load data
                    densityscale,data,delta = self.loadMin(mySimu,myZ,myS)                                        
                    
                    # If histo
                    #plt.bar(densityscale, data, delta, label = myLabel)
                    
                    # If points
                    densityscale += delta/2  # move point to middle of area
                    plt.plot(densityscale, data, label = myLabel, marker = '.')                                      
        
                    # Legend
                    plt.legend(title="Minima density" + myLegend, loc='upper right', ncol=3)

                    # Limits
                    plt.xlim([-1.0,0.0])
                    #plt.ylim([0.0,0.55])

        plt.show()
        plt.clf()


