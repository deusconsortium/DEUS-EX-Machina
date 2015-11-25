from DEUSGraphics import *
from DEUSAnalytics import *


class DEUSProfile(DEUSGraphics, DEUSAnalytics):
    def __init__(self, datapath='../ProfileTracer/data/output/'):
        DEUSGraphics.__init__(self, datapath)
        DEUSAnalytics.__init__(self)

        self._do_plot = False

    def load(self, file_name=None):
        """
            Loading DEUSmodel
            :param arg1: file_name
            :type arg1: string, None by defaut
        """

        DEUSGraphics.Load(self, file_name)
        if DEUSAnalytics.load(self, self._cosmo, self._r):
            self._numerical_smoothing()
            self._compute_sigmas()

    # overiding functions
    def PlotMeanProfile(self, R1value, Dr1='dr', Rsmooth=0.0, velocity = False):
        figure(1)
        f, v = DEUSGraphics.PlotMeanProfile(self, R1value, Dr1, Rsmooth, velocity=velocity)
        R0 = 0.62035 * float(self._boxlen) / float(self._npart)
        fcic = CICDensitySmoothing(self._r, f, R0)

        print 'computing initial profile ...'
        r1 = solve(self._r, f, 1.0)

        if r1 is not None:
            self.local_smooth_spectrum(0.5 * r1, 'th')
            # self.localSmoothSpectrum(r1/(2.0*sqrt(5.0)), 'exp')
            # self.localSmoothSpectrum(10.0,'th')
            # self.localSmoothSpectrum(r1,'th')
            # self.localSmoothSpectrum(R0, 'exp')

            if Rsmooth > 0.0:
                self.local_smooth_spectrum(Rsmooth, 'exp')

            # getting height of the peak from d1 and getting profile
            d0 = self._getInitialHeight(self._r, f, r1)
            R, F, V = self.compute_mean_profile(d0, self._Wm0, self._w, r1, self._a)

            # smoothing with CIC kernel


            if velocity:
                subplot(211)

            plot(self._r, fcic, linestyle='--', marker='', color='g', label='CIC')
            plot(R, F, linestyle='', marker='+', color='r', label='theory')
            legend()

            if velocity:
                subplot(212)
                plot(R, V, linestyle='-', marker='+', color='r', label='theory')
                legend()

            self.reset_smoothing()
        show()

    def PlotRadiusStatistics(self, v0=0.0, Npoints=None, normalized=True):
        R1, Nr1, Nr1d = DEUSGraphics.PlotRadiusStatistics(self, Npoints, normalized)
        Nth, Nthd = self.compute_R1_distribution(R1, v0, normalized=normalized)

        subplot(211)
        plot(R1, Nth, linestyle='-', marker='+', color='r', label='theory')
        legend()

        subplot(212)
        plot(R1, Nthd, linestyle='-', marker='+', color='r', label='theory')
        legend()

        show()

    # proper functions

    def PlotSpectrum(self):
        figure(1)
        grid(True)
        yscale('log')
        xscale('log')
        plot(self._k, self._P0, linestyle='-', color='b', label='linear spectrum')
        xlabel('$k$ in $h^{-1}.Mpc$')
        ylabel('$P(k)$')
        legend()
        show()

    def _getInitialHeight(self, r, f, r1):
        d1init = (self._s2_r1(0, r1) - self._s2_r1(1, r1) * self._S2_r1(0, r1) / self._S2_r1(1, r1)) / (
            self._s2_0[0] * (1. - self._Beta_r1(r1)))
        d = get_density(r, f) - 1.
        fit = SplineFit(r, d)
        d1 = fit(r1)
        d10 = d1 / (1. + 3. * physics_tools.eta(self._a, self._w, self._Wm0, self._dlogD_dloga_init, self._zinit) * (1. + d1))
        return d10 / d1init
