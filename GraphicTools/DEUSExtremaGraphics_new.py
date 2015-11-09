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
        return densityscale, data

    def loadMin(self, simu, Z, S):
        file_root = self._fileRoot(simu, Z, S)
        data = np.loadtxt(self.dir + file_root + '/' + file_root + '_min.deus_histo.txt')        
        densityscale = np.linspace(self.min_start,self.min_end,self.nb_histo,0)        
        return densityscale, data

    @staticmethod
    def _fileRoot(simu, Z, S):
        return simu + "_Z" + Z + "_S" + S + "b"

    def showDensity(self):
        plt.figure(1)

        plt.subplot(211)
        plt.grid(True)

        simu = 'boxlen648_n1024_lcdmw5'
        S = "2" # ["1","2","3"]
        Z = "93" #["81","19","9","4","1","0"]

        densityscale,data = self.loadGlob(simu,Z,S)
        plt.bar(densityscale, data, 3.0/self.nb_histo, label = 'All cells')
        plt.legend(title="Total density S="+S, loc='upper right', ncol=3)

        # for myZ in Z:
        #     densityscale,data = self.loadGlob(simu,myZ,S)
        #     plt.plot(densityscale, data, label = 'Z='+myZ, marker = '.')
        # plt.legend(title="Total density S="+S, loc='upper center', ncol=3)

        # for myS in S:
        #     densityscale,data = self.loadGlob(simu,Z,myS)
        #     plt.plot(densityscale, data, label = 'S='+myS, marker = '.')
        # plt.legend(title="Total density Z="+Z, loc='upper right', ncol=3)

        plt.xlim([-0.3,0.3])
        #plt.ylim([0.0,0.55])

        plt.subplot(212)
        plt.grid(True)

        densityscale,data = self.loadMin(simu,Z,S)
        plt.bar(densityscale, data, 1.0/self.nb_histo, label = 'Minima')
        #plt.plot(densityscale, data, label = 'Avant', marker = '.')
        
        plt.xlim([-0.3,0.3])
        
        plt.legend(title="Total density S="+S, loc='upper right', ncol=3)

        plt.show()
        plt.clf()


