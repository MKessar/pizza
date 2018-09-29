# -*- coding: utf-8 -*-
import os, re
import matplotlib.pyplot as plt
import numpy as np
from .log import PizzaSetup
from .libpizza import scanDir, fast_read
from .npfile import *


class PizzaSpectrum(PizzaSetup):
    """
    This class can be used to read and display the spectra 'spec_#.TAG'
    or 'spec_avg.TAG' and the vorticity balance spectra 'vort_terms_avg.TAG'

    >>> # display the content of spec_1.TAG
    >>> # where TAG is the most recent file in the current directory
    >>> sp = PizzaSpectrum(ispec=1)
    """

    def __init__(self, field='spec', datadir='.', iplot=True, ispec=None,
                 ave=True, tag=None):
        """
        :param field: a string to decide which spectra one wants to load and plot
        :type field: str
        :param iplot: display the output plot when set to True (default is True)
        :type iplot: bool
        :param ispec: the number of the spectrum you want to plot
        :type ispec: int
        :param tag: file suffix (tag), if not specified the most recent one in
                    the current directory is chosen
        :type tag: str
        :param ave: plot a time-averaged spectrum when set to True
        :type ave: bool
        :param datadir: current working directory
        :type datadir: str
        """

        if ispec is not None:
            ave = False

        self.ave = ave

        if field in ('vort_terms', 'forces', 'forcebal', 'force_bal'):
            self.name = 'vort_terms_avg'
            self.ave = True
            ispec = None
        else:
            if self.ave:
                self.name = 'spec_avg'
            else:
                self.name = 'spec_'

        if tag is not None:
            if ispec is not None:
                file = '%s%i.%s' % (self.name, ispec, tag)
                filename = os.path.join(datadir, file)
            else:
                pattern = os.path.join(datadir, '%s*%s' % (self.name, tag))
                files = scanDir(pattern)
                if len(files) != 0:
                    filename = files[-1]
                else:
                    print('No such tag... try again')
                    return

            if os.path.exists(os.path.join(datadir, 'log.%s' % tag)):
                PizzaSetup.__init__(self, datadir=datadir, quiet=True,
                                    nml='log.%s' % tag)
        else:
            if ispec is not None:
                pattern = os.path.join(datadir, '%s%i*' % (self.name, ispec))
                files = scanDir(pattern)
                filename = files[-1]
            else:
                pattern = os.path.join(datadir, '%s*' % self.name)
                files = scanDir(pattern)
                filename = files[-1]
            # Determine the setup
            mask = re.compile(r'.*\.(.*)')
            ending = mask.search(files[-1]).groups(0)[0]
            if os.path.exists(os.path.join(datadir, 'log.%s' % ending)):
                PizzaSetup.__init__(self, datadir=datadir, quiet=True,
                                    nml='log.%s' % ending)

        if not os.path.exists(filename):
            print('No such file')
            return

        data = fast_read(filename)

        self.index = data[:, 0]


        if self.name == 'vort_terms_avg':
            self.buo_mean = data[:, 1]
            self.buo_std = data[:, 2]
            self.cor_mean = data[:, 3]
            self.cor_std = data[:, 4]
            self.adv_mean = data[:, 5]
            self.adv_std = data[:, 6]
            self.domdt_mean = data[:, 7]
            self.domdt_std = data[:, 8]
            self.visc_mean = data[:, 9]
            self.visc_std = data[:, 10]
            self.pump_mean = data[:, 11]
            self.pump_std = data[:, 12]
            self.thwind_mean = data[:, 13]
            self.thwind_std = data[:, 14]
            self.iner_mean = data[:, 15]
            self.iner_std = data[:, 16]
            self.cia_mean = data[:, 17]
            self.cia_std = data[:, 18]
        else:
            if not self.ave:
                self.us2 = data[:, 1]
                self.up2 = data[:, 2]
                self.enst = data[:, 3]
            else:
                self.us2_mean = data[:, 1]
                self.us2_std = data[:, 2]
                self.up2_mean = data[:, 3]
                self.up2_std = data[:, 4]
                self.enst_mean = data[:, 5]
                self.enst_std = data[:, 6]

        if iplot:
            self.plot()

    def plot(self):
        """
        Plotting function
        """
        if self.name == 'vort_terms_avg':
            fig = plt.figure()
            ax = fig.add_subplot(111)
            ax.fill_between(self.index, self.buo_mean-self.buo_std, \
                            self.buo_mean+self.buo_std, alpha=0.1)
            ax.plot(self.index, self.buo_mean, label='Buoyancy')

            ax.fill_between(self.index, self.cor_mean-self.cor_std, \
                            self.cor_mean+self.cor_std, alpha=0.1)
            ax.plot(self.index, self.cor_mean, label='Coriolis')

            ax.fill_between(self.index, self.iner_mean-self.iner_std, \
                            self.iner_mean+self.iner_std, alpha=0.1)
            ax.plot(self.index, self.iner_mean, label='Inertia')

            ax.fill_between(self.index, self.visc_mean-self.visc_std, \
                            self.visc_mean+self.visc_std, alpha=0.1)
            ax.plot(self.index, self.visc_mean, label='Viscosity')

            ax.fill_between(self.index, self.pump_mean-self.pump_std, \
                            self.pump_mean+self.pump_std, alpha=0.1)
            ax.plot(self.index, self.pump_mean, label='Ekman pumping')

            ax.set_yscale('log')
            ax.set_xscale('log')
            ax.set_xlabel('m')
            ax.set_ylabel('Forces')
            ax.set_xlim(1, self.index[-1])
            ax.legend(loc='best', frameon=False)
            fig.tight_layout()

        else:
            if not self.ave:
                fig = plt.figure()
                ax = fig.add_subplot(111)
                if abs(self.up2[0]) > 0.:
                    ax.loglog(self.index[1:]+1, self.us2[1:], label='us**2')
                    ax.loglog(self.index+1, self.up2, label='up**2')
                    ax.loglog(self.index+1, self.enst, label='omega**2')
                else:
                    ax.loglog(self.index[1:]+1, self.us2[1:], label='us**2')
                    ax.loglog(self.index[1:]+1, self.up2[1:], label='up**2')
                    ax.loglog(self.index[1:]+1, self.enst[1:], label='omega**2')
                ax.set_xlabel('m+1')
                ax.set_xlim(1, self.index[-1])
                ax.legend(loc='best', frameon=False)
                fig.tight_layout()
            else:
                fig = plt.figure()
                ax = fig.add_subplot(111)

                if abs(self.up2_mean[0]) > 0.:
                    ax.fill_between(self.index[1:]+1, self.us2_mean[1:]-self.us2_std[1:], \
                                    self.us2_mean[1:]+self.us2_std[1:], alpha=0.1)
                    ax.plot(self.index[1:]+1, self.us2_mean[1:], label='us**2')

                    ax.fill_between(self.index+1, self.up2_mean-self.up2_std, \
                                    self.up2_mean+self.up2_std, alpha=0.1)
                    ax.plot(self.index+1, self.up2_mean, label='up**2')

                    ax.fill_between(self.index+1, self.enst_mean-self.enst_std, \
                                    self.enst_mean+self.enst_std, alpha=0.1)
                    ax.plot(self.index+1, self.enst_mean, label='omega**2')
                else:
                    ax.fill_between(self.index[1:]+1, self.us2_mean[1:]-self.us2_std[1:], \
                                    self.us2_mean[1:]+self.us2_std[1:], alpha=0.1)
                    ax.plot(self.index[1:]+1, self.us2_mean[1:], label='us**2')

                    ax.fill_between(self.index[1:]+1, self.up2_mean[1:]-self.up2_std[1:], \
                                    self.up2_mean[1:]+self.up2_std[1:], alpha=0.1)
                    ax.plot(self.index[1:]+1, self.up2_mean[1:], label='up**2')

                    ax.fill_between(self.index[1:]+1, self.enst_mean[1:]-self.enst_std[1:], \
                                    self.enst_mean[1:]+self.enst_std[1:], alpha=0.1)
                    ax.plot(self.index[1:]+1, self.enst_mean[1:], label='omega**2')

                ax.set_yscale('log')
                ax.set_xscale('log')
                ax.set_xlabel('m+1')
                ax.set_xlim(1, self.index[-1])
                ax.legend(loc='best', frameon=False)
                fig.tight_layout()


class Pizza2DSpectrum(PizzaSetup):
    """
    This class can be used to read and display the 2D spectra '2D_spec_avg.TAG'

    >>> # display the content of 2D_spec_avg.TAG
    >>> # where TAG is the most recent file in the current directory
    >>> sp = Pizza2DSpectrum()
    """

    def __init__(self,  datadir='.', iplot=True, tag=None, endian='l'):
        """
        :param iplot: display the output plot when set to True (default is True)
        :type iplot: bool
        :param tag: file suffix (tag), if not specified the most recent one in
                    the current directory is chosen
        :type tag: str
        :param datadir: current working directory
        :type datadir: str
        :param endian: control the endianness of the file
        :type endian: str
        """


        self.name = '2D_spec_avg'

        if tag is not None:
            pattern = os.path.join(datadir, '%s*%s' % (self.name, tag))
            files = scanDir(pattern)
            if len(files) != 0:
                filename = files[-1]
            else:
                print('No such tag... try again')
                return

            if os.path.exists(os.path.join(datadir, 'log.%s' % tag)):
                PizzaSetup.__init__(self, datadir=datadir, quiet=True,
                                    nml='log.%s' % tag)
        else:
            pattern = os.path.join(datadir, '%s*' % self.name)
            files = scanDir(pattern)
            filename = files[-1]
            # Determine the setup
            mask = re.compile(r'.*\.(.*)')
            ending = mask.search(files[-1]).groups(0)[0]
            if os.path.exists(os.path.join(datadir, 'log.%s' % ending)):
                PizzaSetup.__init__(self, datadir=datadir, quiet=True,
                                    nml='log.%s' % ending)

        if not os.path.exists(filename):
            print('No such file')
            return

        self.read(filename, endian)

        if iplot:
            self.plot()

    def read(self, filename, endian='l'):
        """
        :param filename: name of the input file
        :type filename: str
        :param endian: endianness of the binary file
        :type endian: str
        """
        file = open(filename, 'rb')
        dt = np.dtype("i4, 6Float64")
        self.version, params = np.fromfile(file, dtype=dt, count=1)[0]
        self.ra, self.pr, self.raxi, self.sc, self.ek, self.radratio = params
        dt = np.dtype("4i4")
        self.n_r_max, self.n_m_max, self.m_max, self.minc = np.fromfile(file, dtype=dt, count=1)[0]

        dt = np.dtype("%iFloat64" % self.n_r_max)
        self.radius = np.fromfile(file, dtype=dt, count=1)[0]
        self.idx2m = np.zeros(self.n_m_max)
        for i in range(self.n_m_max):
            self.idx2m[i] = i*self.minc

        dt = np.dtype("(%i,%i)Float64" % (self.n_r_max,self.n_m_max))
        self.us2 = np.fromfile(file, dtype=dt, count=1)[0, :, :]
        self.us2 = self.us2.T
        self.up2 = np.fromfile(file, dtype=dt, count=1)[0, :, :]
        self.up2 = self.up2.T
        self.enst = np.fromfile(file, dtype=dt, count=1)[0, :, :]
        self.enst = self.enst.T

        self.ekin = self.us2+self.up2

        file.close()

    def plot(self, levels=65, cm='magma', cut=1., solid_contour=False, 
             log_yscale=True):
        """
        Plotting function

		:param cm: name of the colormap ('jet', 'seismic', 'RdYlBu_r', etc.)
		:type cm: str
		:param levels: the number of levels used in the contour plot
		:type levels: int
		:param  cut: a coefficient to change the dynamics of the contour levels
		:param cut: float
        :param solid_contour: a boolean to decide whether one also wants the solid
                              contour lines
        :type solid_contour: bool
        :param log_yscale: a boolean to decide whether one wants a logarithmic
                           y-axis
        :type log_yscale: bool
        """

        vmax = cut*np.log10(self.ekin[1:,:]+1e-34).max()
        vmin = cut*np.log10(self.ekin[self.ekin>1e-15]).min()
        levs = np.linspace(vmin, vmax, levels)

        fig = plt.figure(figsize=(8,4))
        ax1 = fig.add_subplot(121)
        im = ax1.contourf(self.radius[1:-1], self.idx2m+1,
                          np.log10(self.ekin[:,1:-1]+1e-20),
                         levs, extend='both', cmap=plt.get_cmap(cm))
        if solid_contour:
            ax1.contour(self.radius[1:-1], self.idx2m+1,
                        np.log10(self.ekin[:,1:-1]+1e-20), levs,
                        extend='both', linestyles=['-'], colors=['k'],
                        linewidths=[0.5])
        fig.colorbar(im)
        ax1.set_title('Kinetic energy')
        if log_yscale: ax1.set_yscale('log')
        ax1.set_xlabel('Radius')
        ax1.set_ylabel('m+1')
        
        vmax = cut*np.log10(self.enst[1:,:]+1e-34).max()
        vmin = cut*np.log10(self.enst[self.enst>1e-15]).min()
        levs = np.linspace(vmin, vmax, levels)
        ax2 = fig.add_subplot(122, sharey=ax1, sharex=ax1)
        im = ax2.contourf(self.radius, self.idx2m+1, np.log10(self.enst+1e-20),
                         levs, extend='both', cmap=plt.get_cmap(cm))
        if solid_contour:
            ax2.contour(self.radius, self.idx2m+1, np.log10(self.enst+1e-20),
                       levs, extend='both', linestyles=['-'], colors=['k'],
                       linewidths=[0.5])
        fig.colorbar(im)
        ax2.set_title(r'Enstrophy')
        if log_yscale: ax2.set_yscale('log')
        ax2.set_xlabel('Radius')

        plt.setp(ax2.get_yticklabels(), visible=False)
        ax1.set_xlim(self.radius[-1], self.radius[0])

        fig.tight_layout()
