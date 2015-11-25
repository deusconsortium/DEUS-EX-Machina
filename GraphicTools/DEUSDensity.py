from DEUSExtremaGraphics_new import *
from DEUSAnalytics import *
from DEUSCosmo import *
from DEUSTools import *


class DEUSDensity(DEUSExtremaGraphics, DEUSAnalytics):
    def __init__(self):
        DEUSExtremaGraphics.__init__(self)
        DEUSAnalytics.__init__(self)

        self._data_path = "../compute_extrema/data"

        self._Rg = 0.0
        self._Z = None
        self._S = None

        self._P = None

        self._dtot = []
        self._dPDF = []
        self._dmin = []
        self._minPDF = []

    def _extractInfoFromName(self, file_name):
        arguments = file_name.split('_')
        self._boxlen = int(arguments[0][6:])
        self._npart = int(arguments[1][1:])
        self._cosmo = arguments[2]
        self._Z = arguments[3][1:]
        self._a = 1. / (1. + float(self._Z))
        self._S = arguments[4][1:]
        self._Rg = float(self._boxlen) / float(self._npart) * float(self._S)

    def load(self, file_name=None):
        if file_name is None:
            self.load(selectFile(self._data_path))
        else:
            self._extractInfoFromName(file_name)
            self._load_spectrum()
            self._numerical_smoothing()
            self._gaussian_smoothing(self._Rg)
            self._set_simulation(self._boxlen, self._npart, self._cosmo)

            self._compute_sigmas0(self._P0)

            simu = 'boxlen' + str(self._boxlen) + '_n' + str(self._npart) + '_' + self._cosmo
            self._dtot, self._dPDF = self.loadGlob(simu, self._Z, self._S)
            self._dmin, self._minPDF = self.loadMin(simu, self._Z, self._S)
            self._normalizeData()

    def _normalizeData(self):
        self._dPDF /= integrate(self._dtot, self._dPDF)
        self._minPDF /= integrate(self._dmin, self._minPDF)

    def getNumPDF(self):
        return self._dtot, self._dPDF, self._dmin, self._minPDF

    def PlotPDF(self, theory=True):
        if theory:
            figure(1)
            grid(True)
            xlabel('$\\delta$')
            ylabel('PDF')
            dtot, Ptot = self.compute_linear_density_PDF(1. / (self._a) - 1.)
            # bar(self._dtot,self._dPDF,width=self._dtot[1]-self._dtot[0],color='b')
            plot(self._dtot, self._dPDF, linestyle='', marker='o', color='b', label='data')
            plot(dtot - 1., Ptot, linestyle='-', marker='+', color='g', label='gaussian prediction')
            title("density PDF at z=" + str(1. / (self._a) - 1.) + " for boxlen" + str(self._boxlen) + "_n" + str(
                self._npart) + " Rs = " + str(self._Rg) + " [Mpc/h]")
            legend()

            figure(2)
            grid(True)
            xlabel('$\\delta$')
            ylabel('PDF')
            title("extremum PDF at z=" + str(1. / (self._a) - 1.) + " for boxlen" + str(self._boxlen) + "_n" + str(
                self._npart) + " Rs = " + str(self._Rg) + " [Mpc/h]")
            dtheo_lin, Ndtheo_lin = DEUSAnalytics.compute_linear_extremum_distribution_from_heigh(self, 1. / self._a - 1.)
            dtheo_evo, Ndtheo_evo = DEUSAnalytics.compute_evolved_extremum_distribution_from_heigh(self, self._w, self._Wm0,
                                                                                                   1. / self._a - 1.)
            dtheo_evo_cmb, Ndtheo_evo_cmb = DEUSAnalytics.compute_evolved_extremum_distribution_from_heigh(self, self._w,
                                                                                                           self._Wm0,
                                                                                                           1. / self._a - 1.,
                                                                                                           zinit=1100.)
            dtheo_zeldo, Ndtheo_zeldo = DEUSAnalytics.compute_zeldovitch_extremum_distribution_from_heigh(self, self._w,
                                                                                                          self._Wm0,
                                                                                                          1. / self._a - 1.)

            Ndtheo_lin /= integrate(dtheo_lin - 1., Ndtheo_lin)
            Ndtheo_evo /= integrate(dtheo_evo - 1., Ndtheo_evo)
            Ndtheo_evo_cmb /= integrate(dtheo_evo_cmb - 1., Ndtheo_evo_cmb)
            Ndtheo_zeldo /= integrate(dtheo_zeldo - 1., Ndtheo_zeldo)

            # bar(self._dmin,self._minPDF,width=self._dmin[1]-self._dmin[0],color='b')
            plot(self._dmin, self._minPDF, linestyle='', marker='o', color='b', label='data')
            plot(dtheo_lin - 1., Ndtheo_lin, linestyle='-', marker='+', color='g', label='gaussian prediction')
            plot(dtheo_evo - 1., Ndtheo_evo, linestyle='-', marker='+', color='r', label='evolved prediction')
            plot(dtheo_evo_cmb - 1., Ndtheo_evo_cmb, linestyle='--', marker='*', color='r', label='evolved from z=1100')
            plot(dtheo_zeldo - 1., Ndtheo_zeldo, linestyle='--', marker='', color='k', label='Zeldovitch prediction')
            legend(loc=2)
            show()
        else:
            figure(1)
            grid(True)
            xlabel('$\\delta$')
            ylabel('PDF')
            plot(self._dtot, self._dPDF, linestyle='', marker='o', color='b', label='data')
            title("density PDF at z=" + str(1. / (self._a) - 1.) + " for boxlen" + str(self._boxlen) + "_n" + str(
                self._npart) + " Rs = " + str(self._Rg) + " [Mpc/h]")
            legend()

            figure(2)
            grid(True)
            xlabel('$\\delta$')
            ylabel('PDF')
            title("extremum PDF at z=" + str(1. / (self._a) - 1.) + " for boxlen" + str(self._boxlen) + "_n" + str(
                self._npart) + " Rs = " + str(self._Rg) + " [Mpc/h]")
            plot(self._dmin, self._minPDF, linestyle='', marker='o', color='b', label='data')
            legend(loc=2)
            show()
