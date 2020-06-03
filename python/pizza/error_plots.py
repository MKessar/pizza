# -*- coding: utf-8 -*-
import copy
import os
import re
from .npfile import npfile
import matplotlib.pyplot as plt
import numpy as np
from .log import PizzaSetup
from .libpizza import scanDir, fast_read,cc2real,chebgrid,intcheb,costf,spec_spat,spat_spec
from .frame import Frame
import scipy.interpolate as inp



class PizzaErrorInstant(PizzaSetup):
    """
    This class can be used to read and display the spectra 'spec_#.TAG'
    or 'spec_avg.TAG' and the vorticity balance spectra 'vort_terms_avg.TAG'

    >>> # display the content of spec_1.TAG
    >>> # where TAG is the most recent file in the current directory
    >>> sp = PizzaSpectrum(ispec=1)
    >>> # stack the content of spec_ave.test_a to spec_ave.test_c
    >>> sp = PizzaSpectrum(tag='test_[a-c]', iplot=False)
    """
    """
    This module is used to read the frames, transform them back on the
    physical grid and display them.
    """

    def __init__(self, ivar=None, datadir='.', filename1=None,filename2=None, endian='l',
                 verbose=False,r_icb=0.1, r_cmb=0.9 ):

        """
        :param ivar: the number of the snapshot file
        :type ivar: int
        :param datadir: working directory
        :type datadir: str
        :param filename1: name of the snapshot file
        :type filename1: str
        :param filename2: name of the snapshot file
        :type filename2: str
        :param endian: endianness of the file ('B' or 'l')
        :type endian: str
        :param verbose: a boolean to display some informations
        :type verbose: bool
        :param r_icb: inner radius
        :type r_icb: float
        :param r_cmb: outer radius
        :type r_cmb: float
        """

        #filename = self.get_filename('frame_us', ivar, datadir, tag, verbose)
        f1 = Frame(filename1, endian=endian)
        f2 = Frame(filename2, endian=endian)
        #print "spec_inst a shape (f.field_m)=", np.shape(f.field_m)

        #self.f_m = f.field_m
        self.n_r_max1 = f1.n_r_max
        self.n_m_max1 = f1.n_m_max
        print "n_r_max1",self.n_r_max1
        print "n_m_max1",self.n_m_max1
        self.n_phi1= (self.n_m_max1-1)*2
        self.dphi1= 2.0*np.pi/(self.n_phi1)
        print "self.n_phi1=", self.n_phi1


        self.n_r_max2 = f2.n_r_max
        self.n_m_max2 = f2.n_m_max
        print "n_r_max2",self.n_r_max2
        print "n_m_max2",self.n_m_max2
        self.n_phi2 = (self.n_m_max2-1)*2
        self.dphi2  = 2.0*np.pi/(self.n_phi2)
        print "self.n_phi2=", self.n_phi2
        
        #self.f_m = costf(f.field_m)
        rr1 = chebgrid(self.n_r_max1-1, r_icb, r_cmb)
        rr2 = chebgrid(self.n_r_max2-1, r_icb, r_cmb)

        # assuming f1 is larger resolution than f2
        
        
        self.n_m_max_all=np.maximum(self.n_m_max1,self.n_m_max2)
        self.n_r_max_all=np.maximum(self.n_r_max1,self.n_r_max2)
        self.n_m_max=self.n_m_max_all
        self.n_r_max=self.n_r_max_all
        #for i in range(0,self.n_r_max2):
            #for k in range(0,self.n_r_max2):
                #print "i=",i,"k=", k, "f1.field_m[i,k]=", f1.field_m[i,k]    , "f2.field_m[i,k]=", f2.field_m[i,k]     
        
        #f1_cheb = costf(f1.field_m,fac=False)
        #f2_cheb = costf(f2.field_m,fac=False)
        f1_cheb = costf(f1.field_m)
        f2_cheb = costf(f2.field_m)
        #print "shape (self.f2_cheb)=", np.shape(self.f2_cheb)
        #print "shape (self.f2_cheb)=", np.shape(self.f2_cheb)

        #for i in range(0,5):
            #for k in range(0,self.n_r_max2):
                #print "i=",i,"k=", k, "self.f1_cheb[i,k]=", f1_cheb[i,k]    , "self.f2_cheb[i,k]=", f2_cheb[i,k]     
        f2_cheb_test=2.0*f2_cheb
        for i in range(0,self.n_r_max2):
            for k in range(0,self.n_r_max2):
                #f2_cheb_test[i,k]=(((f2_cheb.shape[-1]-1)/(f1_cheb.shape[-1]-1)))*f2_cheb[i,k]
                #f2_cheb_test[i,k]=(((f2_cheb.shape[-1]-1)))*f2_cheb[i,k]
                f2_cheb_test[i,k]=np.sqrt((1.0*(f1_cheb.shape[-1]-1))/(1.0*(f2_cheb.shape[-1]-1)))*f2_cheb[i,k]
                #f2_cheb_test[i,k]=(f1_cheb.shape[-1]-1)*f2_cheb[i,k]

        #print "f1_cheb.shape[-1]-1", f1_cheb.shape[-1]-1
        #print "f2_cheb.shape[-1]-1", f2_cheb.shape[-1]-1
        #print "f2_cheb.shape[-1]-1/f2_cheb.shape[-1]-1 test1", (((f2_cheb.shape[-1]-1))/((f1_cheb.shape[-1]-1)))
        #print "f2_cheb.shape[-1]-1/f1_cheb.shape[-1]-1 test2", ((1.0*(f2_cheb.shape[-1]-1))/(1.0*(f1_cheb.shape[-1]-1)))
        #print "f1_cheb.shape[-1]-1/f2_cheb.shape[-1]-1 test2", ((1.0*(f1_cheb.shape[-1]-1))/(1.0*(f2_cheb.shape[-1]-1)))
        #print "sqrt f2_cheb.shape[-1]-1/f1_cheb.shape[-1]-1 test2", np.sqrt((1.0*(f2_cheb.shape[-1]-1))/(1.0*(f1_cheb.shape[-1]-1)))
        #print "sqrt f1_cheb.shape[-1]-1/f2_cheb.shape[-1]-1 test2", np.sqrt((1.0*(f1_cheb.shape[-1]-1))/(1.0*(f2_cheb.shape[-1]-1)))

        #f1_cheb_backreal = costf(f1_cheb)
        #f2_cheb_backreal = costf(f2_cheb)
        #for i in range(0,5):
            #for k in range(0,self.n_r_max2):
                #print "i=",i,"k=", k, "self.f1_cheb[i,k]=", f1_cheb[i,k] , "self.f2_cheb[i,k]=", f2_cheb[i,k], "self.f2_cheb_test[i,k]=", f2_cheb_test[i,k]


        #for i in range(0,self.n_r_max2):
            #for k in range(0,self.n_r_max2):
                #print "i=",i,"k=", k, "f1.field_m[i,k]=", f1.field_m[i,k]    , "f2.field_m[i,k]=", f2.field_m[i,k]     
                #print "i=",i,"k=", k, "self.f1_cheb_backreal[i,k]=", f1_cheb_backreal[i,k]    , "self.f2_cheb_backreal[i,k]=", f2_cheb_backreal[i,k]     


        #for i in range(0,self.n_r_max2):
            #for k in range(0,self.n_r_max2):
                #print "i=",i,"k=", k, "self.f1_cheb_backreal[i,k]=", f1_cheb_backreal[i,k]    , "self.f2_cheb_backreal[i,k]=", f2_cheb_backreal[i,k]     

                
        f2_bis_cheb=2.0*f1_cheb
        
        
        for i in range(0,self.n_m_max_all):
            for k in range(0,self.n_r_max_all):
                #self.f1_bis[i,k] = 0.000000000000000
                f2_bis_cheb[i,k] = 0.000000000000000
        #print "shape (f2_bis_cheb)=", np.shape(f2_bis_cheb)

        #for i in range(0,self.n_m_max):
            #for k in range(0,self.n_r_max):
                #print "i=",i,"k=", k, "f2_bis_cheb[i,k]=", f2_bis_cheb[i,k],"f1_cheb[i,k]=", f1_cheb[i,k]

        #for i in range(0,self.n_m_max1):
            #for k in range(0,self.n_r_max1):
                #self.f1_bis_cheb[i,k] = self.f1_cheb[i,k]

        for i in range(0,self.n_m_max2):
            for k in range(0,self.n_r_max2):
                #f2_bis_cheb[i,k] = f2_cheb[i,k]
                #f2_bis_cheb[i,k] = (self.n_r_max1-1)*f2_cheb[i,k]/(self.n_r_max2-1)
                #f2_bis_cheb[i,k] = 2.0*(np.sqrt(0.5/(self.n_r_max1-1)))*f2_cheb[i,k]/(np.sqrt(0.5/(self.n_r_max2-1)))
                #f2_bis_cheb[i,k] = 4.0*(np.sqrt(0.5/(f2_bis_cheb.shape[-1]-1)))*f2_cheb[i,k]/(np.sqrt(0.5/(f2_cheb.shape[-1]-1)))
                f2_bis_cheb[i,k]=np.sqrt((1.0*(f1_cheb.shape[-1]-1))/(1.0*(f2_cheb.shape[-1]-1)))*f2_cheb[i,k]
                
        #test_shape=f2_cheb.shape[-1]
        #print "shape (f2_bis_cheb)=", np.shape(f2_bis_cheb)
        #print "shape (f2_bis_cheb)=", (f2_bis_cheb.shape[-1])
        #print "shape (f2_cheb)=", (f2_cheb.shape[-1])
        #print "shape (f2_cheb)=", test_shape

        #for i in range(0,self.n_m_max2):
            #for k in range(0,self.n_r_max2):
                #print "i=",i,"k=", k, "self.f2_bis_cheb[i,k]=", f2_bis_cheb[i,k],"self.f1_cheb[i,k]=", f1_cheb[i,k]

        #f1_cheb = costf(f1.field_m)
        f2_large_res = costf(f2_bis_cheb)
        for i in range(0,5):
            for k in range(0,self.n_r_max2):
                print "i=",i,"k=", k, "f2_large_res[i,k]=", f2_large_res[i,k],"f1.field_m[i,k]=", f1.field_m[i,k]

        f1_cheb_abs= np.sqrt(f1_cheb*np.conjugate(f1_cheb))
        f1_cheb_abs_real= f1_cheb_abs.real

        f2_bis_cheb_abs= np.sqrt(f2_bis_cheb*np.conjugate(f2_bis_cheb))
        f2_bis_cheb_abs_real= f2_bis_cheb_abs.real
        
        f_error_spec_cheb=abs(f1_cheb_abs_real - f2_bis_cheb_abs_real)
        #/(self.f1_cheb_abs_real)
        #print "shape (self.f_error_spec_cheb)=", np.shape(self.f_error_spec_cheb)

        #for i in range(0,self.n_m_max):
            #for k in range(0,self.n_r_max):
                #print "i=",i,"k=", k, "self.f_error_spec_cheb[i,k]=", f_error_spec_cheb[i,k]
        #for i in range(0,self.n_m_max2):
            #for k in range(0,self.n_r_max2):
                #print "i=",i,"k=", k, "self.f2_bis_cheb[i,k]=", f2_bis_cheb[i,k],"self.f1_cheb[i,k]=", f1_cheb[i,k], "self.f_error_spec_cheb[i,k]=", f_error_spec_cheb[i,k]
                
        
        for i in range(0,self.n_m_max_all):
            for k in range(0,self.n_r_max_all):
                #self.f1_bis[i,k] = 0.000000000000000
                if (abs(f1_cheb_abs_real[i,k]))>0.000:
                   f_error_spec_cheb[i,k] = f_error_spec_cheb[i,k]/f1_cheb_abs_real[i,k]
                if (abs(f1_cheb_abs_real[i,k]))==0.000:
                    #print "i=",i,"k=", k,"f1_cheb_abs_real=0.0"
                    f_error_spec_cheb[i,k]=0.0
                    
        #for i in range(0,self.n_m_max2):
            #for k in range(0,self.n_r_max2):
                ##self.f1_bis[i,k] = 0.000000000000000
                #if (abs(f1_cheb_abs_real[i,k]))>(1.0*10**(-10.0)):
                   #print "i=",i,"k=", k,"f_error_spec_cheb=", f_error_spec_cheb[i,k] 
                     
                    
        #for i in range(0,self.n_m_max2):
            #for k in range(0,self.n_r_max2):
                #print "i=",i,"k=", k, "self.f_error_spec_cheb[i,k]=", f_error_spec_cheb[i,k]
                
        self.Rspec_f_r_2D_spec=f_error_spec_cheb
        
        #print "shape (self.Rspec_f_r_2D_spec)=", np.shape(self.Rspec_f_r_2D_spec)

        #for i in range(0,self.n_m_max):
            #for k in range(0,self.n_r_max):
                #print "i=",i,"k=", k, "self.Rspec_f_r_2D_spec[i,k]=", self.Rspec_f_r_2D_spec[i,k]

        
        #self.f_real= spec_spat(f.field_m, self.n_phi)
  
        
        
        #self.f_cheb = self.f_real

        #self.f_m = spat_spec(self.f_cheb, self.n_m_max)
     
        
        
        #self.Rspec_f=np.sqrt( self.f_m*np.conjugate( self.f_m))


        
        #for k in range(0,self.n_r_max):
            #self.Rspec_f[:,k] = np.pi*self.Rspec_f[:,k]*rr[k]
            #self.f_cheb = costf(self.Rspec_f)

        self.index = range(self.n_m_max_all)

        #print "shape (self.index)=", np.shape(self.index)
        #self.Rspec_f_r = costf(self.Rspec_f_r)

        #self.Rspec_f_r = self.Rspec_f.real
        #self.Rspec_f_r_2D_spec = costf(self.Rspec_f_r)

        #print "shape (self.Rspec_f_r)=", np.shape(self.Rspec_f_r)

      


# Rspec_f_r_2D_spec is now the 2d spectra
        
        #self.Rspec_fm   = intcheb(self.Rspec_f)
        #self.Rspec_fm_r = self.Rspec_fm.real
        #self.Rspec_fm   = intcheb(self.Rspec_f_r)
        #self.Rspec_fm_r = self.Rspec_fm.real
# Rspec_fm is now the power spectra as a function of m+1
# 

# need to plot those
        
        
        
    def readinst(self, filename, endian='l'):
        """
        :param filename: name of the input file
        :type filename: str
        :param endian: endianness of the binary file
        :type endian: str
        """
        file = open(filename, 'rb')
        dt = np.dtype("i4, 6f8")
        self.version, params = np.fromfile(file, dtype=dt, count=1)[0]
        self.ra, self.pr, self.raxi, self.sc, self.ek, self.radratio = params
        dt = np.dtype("4i4")
        self.n_r_max, self.n_m_max, self.m_max, self.minc = \
            np.fromfile(file, dtype=dt, count=1)[0]

        dt = np.dtype("%if8" % self.n_r_max)
        self.radius = np.fromfile(file, dtype=dt, count=1)[0]
        self.idx2m = np.zeros(self.n_m_max)
        for i in range(self.n_m_max):
            self.idx2m[i] = i*self.minc

        dt = np.dtype("(%i,%i)f8" % (self.n_r_max, self.n_m_max))
        if self.version == 1:
            data = np.fromfile(file, dtype=dt, count=3)
        elif self.version == 2:
            data = np.fromfile(file, dtype=dt, count=6)

        file.close()

        return data        
        
    def plot2D_inst(self,r_icb,r_cmb, levels=17, cm='magma', cut=1., solid_contour=True,
             log_yscale=True):        
        """
        Plotting function
        :type r_icb: float
        :param r_cmb: outer radius
        :type r_cmb: float   
        :param cm: name of the colormap ('jet', 'seismic', 'RdYlBu_r', etc.)
        :type cm: str
        :param levels: the number of levels used in the contour plot
        :type levels: int
        :param  cut: a coefficient to change the dynamics of the contour levels
        :param cut: float
        :param solid_contour: a boolean to decide whether one also wants the
                              solid contour lines
        :type solid_contour: bool
        :param log_yscale: a boolean to decide whether one wants a logarithmic
                           y-axis
        :type log_yscale: bool
        :param r_icb: inner radius
     
        """      
        
        # data=PizzaSpectrumInstant.readinst( filename='spec_avg.test', endian='l')
        #self.radius =   chebgrid(self.n_r_max-1, r_icb, r_cmb)
        #print "self.radius=",self.radius
        #f=open("radius.test", "r")
        #centents = f.read()
        #print "contents=",contents
        
        self.radius = np.zeros(self.n_r_max)
        for i in range(1,self.n_r_max):
            self.radius[i] = i*1.0        
            
        self.idx2m = np.zeros(self.n_m_max)
        for i in range(self.n_m_max):
            self.idx2m[i] = i*1.0
        print "np.shape(self.radius)=", np.shape(self.radius), "np.shape(idx2m)=", np.shape(self.idx2m)
        print "np.shape(self.Rspec_f_r_2D_spec)=", np.shape(self.Rspec_f_r_2D_spec)
        #print "idx2m=", self.idx2m
        #print "idx2m=", n_m_max
        ind1 = self.Rspec_f_r_2D_spec[:, 1:-1].argmax(axis=1)
        #print "ind1", ind1
        
        ind = self.Rspec_f_r_2D_spec[:, 1:-1].argmax(axis=0)
        #print "ind", ind
        self.peaks = np.zeros_like(ind)
        for k, idx in enumerate(ind):
            self.peaks[k] = self.idx2m[idx]
        
        vmax = cut*np.log10(self.Rspec_f_r_2D_spec[1:, :]+1e-34).max()
        vmin = cut*np.log10(self.Rspec_f_r_2D_spec[self.Rspec_f_r_2D_spec > 1e-15]).min()
        #vmax = cut*(self.Rspec_f_r_2D_spec[1:, :]+1e-34).max()
        #vmin = cut*(self.Rspec_f_r_2D_spec[self.Rspec_f_r_2D_spec > 1e-15]).min()
        levs = np.linspace(vmin, vmax, levels)

        fig = plt.figure()
        ax1 = fig.add_subplot(111)
        #im = ax1.contourf(self.radius, self.idx2m,
                          #(self.Rspec_f_r_2D_spec+1e-20),
                          #levs, extend='both', cmap=plt.get_cmap(cm))
        im = ax1.contourf(self.radius, self.idx2m,
                          np.log10(self.Rspec_f_r_2D_spec+1e-20),
                          levs, extend='both', cmap=plt.get_cmap(cm))
        if solid_contour:
            #ax1.contour(self.radius[1:-1], self.idx2m[1:],
                        #(self.Rspec_f_r_2D_spec[1:, 1:-1]+1e-20), levs,
                        #extend='both', linestyles=['-'], colors=['k'],
                        #linewidths=[0.5])
            ax1.contour(self.radius[1:-1], self.idx2m[1:],
                        np.log10(self.Rspec_f_r_2D_spec[1:, 1:-1]+1e-20), levs,
                        extend='both', linestyles=['-'], colors=['k'],
                        linewidths=[0.5])
        ax1.plot(self.radius[1:-1], self.peaks, ls='--')
        ax1.set_title('logT(r,m+1)')
        #if log_yscale:
            #ax1.set_yscale('log')
        ax1.set_xlabel('Tn')
        ax1.set_ylabel('m')
        
        ax1.legend(loc='best', frameon=False)
        fig.colorbar(im)



        #fichier = open("spec2D_m.txt", "w")
        #fichier.write(self.index, self.Rspec_fm_r)


        #for k in range(0,self.n_r_max):
            #for i in range(0,self.n_m_max):
              ##fichier.write( str(self.Rspec_f_r_2D_spec[i,k])+"    ")
                #print "i=",i,"k=", k, "self.Rspec_f_r_2D_spec[i,k]=", self.Rspec_f_r_2D_spec[i,k]
          ##fichier.write( "\n")
        ##fichier.close() 
        
        
        
    def plot_inst(self,r_icb,r_cmb, levels=17, cm='magma', cut=1., solid_contour=True,
             log_yscale=True):        
        """
        Plotting function
        :type r_icb: float
        :param r_cmb: outer radius
        :type r_cmb: float   
        :param cm: name of the colormap ('jet', 'seismic', 'RdYlBu_r', etc.)
        :type cm: str
        :param levels: the number of levels used in the contour plot
        :type levels: int
        :param  cut: a coefficient to change the dynamics of the contour levels
        :param cut: float
        :param solid_contour: a boolean to decide whether one also wants the
                              solid contour lines
        :type solid_contour: bool
        :param log_yscale: a boolean to decide whether one wants a logarithmic
                           y-axis
        :type log_yscale: bool
        :param r_icb: inner radius
     
        """              
        
        fig = plt.figure()
        ax = fig.add_subplot(111)
        sd = self.Rspec_fm_r/np.sqrt(self.Rspec_fm_r)/2.
        #ax.fill_between(self.index, np.sqrt(self.Rspec_fm_r)-sd,
                          #np.sqrt(self.Rspec_fm_r)+sd, alpha=0.1)
        ax.plot(self.index, self.Rspec_fm_r, label='T(m+1)')


        ax.set_yscale('log')
        ax.set_xscale('log')
        ax.set_xlabel('m+1')
        ax.set_ylabel('T(m+1)')
        ax.set_xlim(1, self.index[-1])
        ax.legend(loc='best', frameon=False)
        fig.tight_layout()        
        
        fichier = open("spec1D_m.txt", "w")
        #fichier.write(self.index, self.Rspec_fm_r)

        for k in range(1,self.n_m_max):

          fichier.write( str(self.index[k])+"    ")
          fichier.write( str(self.Rspec_fm_r[k]))
          fichier.write( "\n")
        fichier.close() 
