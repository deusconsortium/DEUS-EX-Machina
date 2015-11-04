__author__ = 'jpasdeloup'
import numpy as np
import matplotlib.pyplot as plt

class DEUSExtremaGraphics :

    def __init__(self):
        self.dir = "../compute_extrema/data/"

    def loadGlob(self, simu, Z, S):
        file_root = self._fileRoot(simu, Z, S)
        data = np.loadtxt(self.dir + file_root + '/' + file_root + '_all.deus_histo.txt')
        data = np.delete(data, 100)
        densityscale = np.linspace(-1.0,2.0,100)
        return densityscale, data

    def loadMin(self, simu, Z, S):
        file_root = self._fileRoot(simu, Z, S)
        data = np.loadtxt(self.dir + file_root + '/' + file_root + '_min.deus_histo.txt')
        densityscale = np.linspace(-1.0,0.0,101)
        return densityscale, data

    @staticmethod
    def _fileRoot(simu, Z, S):
        return simu + "_Z" + Z + "_S" + S

    def showDensity(self):
        plt.figure(1)

        plt.subplot(211)
        plt.grid(True)

        simu = 'boxlen648_n1024_lcdmw5'
        S = ["1","2","3"]
        Z = "93" #["93","36","17","9","4","2","1","0"]

        # for myZ in Z:
        #     data = self.loadGlob(simu,myZ,S)
        #     plt.plot(densityscale, data, label = 'Z='+myZ, marker = '.')
        # plt.legend(title="Total density S="+S, loc='upper center', ncol=3)

        for myS in S:
            densityscale,data = self.loadGlob(simu,Z,myS)
            plt.plot(densityscale, data, label = 'S='+myS, marker = '.')
        plt.legend(title="Total density Z="+Z, loc='upper center', ncol=3)

        plt.xlim([-1.0,0.0])
        #plt.ylim([0.0,0.35])

        plt.subplot(212)
        plt.grid(True)

        # for myZ in Z:
        #     data = self.loadMin(simu,myZ,S)
        #     plt.plot(densityscale, data, label = 'Z='+myZ, marker = '.')
        # plt.legend(title="Minima density S="+S, loc='upper center', ncol=3)

        for myS in S:
            densityscale,data = self.loadMin(simu,Z,myS)
            plt.plot(densityscale, data, label = 'S='+myS, marker = '.')
        plt.legend(title="Minima density Z="+Z, loc='upper center', ncol=3)

        plt.xlim([-1.0,0.0])
        #plt.ylim([0.0,0.25])

        # plt.subplot(212)
        # plt.grid(True)
        # densityscale = np.linspace(-1.0,0.0,100)
        # plt.plot(densityscale, self.mins0, color = 'r', label = 'Z= 0 (414 330 minima)', marker = '.')
        # plt.plot(densityscale, self.mins36, color = 'g', label = 'Z=36 (559 930 minima)', marker = '.')
        # plt.plot(densityscale, self.mins93, color = 'b', label = 'Z=93 (568 781 minima)', marker = '.')
        # plt.legend(title="Minima density", loc='upper center')

        plt.show()
        plt.clf()


