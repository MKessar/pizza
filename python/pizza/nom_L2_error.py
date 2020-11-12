import copy
import os
import re
from .npfile import npfile
import matplotlib.pyplot as plt
import numpy as np
from .log import PizzaSetup
from .libpizza import scanDir, fast_read,cc2real,chebgrid,intcheb,costf,spec_spat,spat_spec,symmetrize, scanDir
from .frame import Frame
import scipy.interpolate as inp
from .plotlib import equatContour

class norm_L2_error_instdef(PizzaSetup):
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
                 verbose=False ):

#def maxerror_instdef(self, ivar=None, datadir='.', filename1=None,filename2=None, endian='l',
                 #verbose=False,r_icb=0.1, r_cmb=0.9 ):
 


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

        """


    def norm_L2_err_rel(self,filename1=None,filename2=None, endian='l'):

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

                
      
        f1 = Frame(filename1, endian=endian)
        f2 = Frame(filename2, endian=endian)

        self.n_r_max1 = f1.n_r_max
        self.n_m_max1 = f1.n_m_max
        self.n_phi1= (self.n_m_max1-1)*2
        self.dphi1= 2.0*np.pi/(self.n_phi1)


        self.n_r_max2 = f2.n_r_max
        self.n_m_max2 = f2.n_m_max
        self.n_phi2 = (self.n_m_max2-1)*2
        
        f1_r = spec_spat(f1.field_m, self.n_phi1)
        f2_r = spec_spat(f2.field_m, self.n_phi2)

        #for i in range(0,self.n_phi1):
            #for j in range(0,self.n_r_max1):
                #f1_r[i,j]=0.0
                #f2_r[i,j]=0.0
        #f2_r[4,4]=5.0
        #f1_r[4,4]=1.0

        #f_error_real_abs=abs(f1_r - f2_r)
        #f_diff=np.zero(self.n_phi1,self.n_r_max1)
        f_diff=np.zeros((f1_r.shape[0],f1_r.shape[1]), dtype=f1_r.dtype)
        f_diff = (f1_r-f2_r)**2.0
        f1_r_square=f1_r**2.0
                
        
#cm=, levels=65, deminc=True,
              #normed=True, vmax=None, vmin=None, normRad=False, stream=False,
              #streamNorm='vel', streamDensity=1.5, cbar=True, label=None,
              #streamColor='k'
        f_diff_intr      = intcheb(f_diff)
        f1_r_square_intr = intcheb(f1_r_square)
        
        #print "f_diff_intr=", f_diff_intr
        #print "f2_r_square_intr=", f2_r_square_intr
        
        
        f_diff_intr_intm      = 0.000000000000000
        f1_r_square_intr_intm = 0.000000000000000
        
        for i in range(0,self.n_phi1):
            f_diff_intr_intm      = f_diff_intr_intm      + f_diff_intr[i]
            f1_r_square_intr_intm = f1_r_square_intr_intm + f1_r_square_intr[i]
            
        L2_norm_rel_error=np.sqrt(f_diff_intr_intm/f1_r_square_intr_intm)
        print "L2_norm_rel_error", L2_norm_rel_error


        #data = symmetrize(np.sqrt(f_diff), ms=True)
        ##file = npfile.npfile(filename1, endian=endian)

        
        #self.fig, xx, yy = equatContour(data, f1.radius, minc=True,
                                        #levels=365, cm='seismic', deminc=True,
                                        #normed=True, vmax=None, vmin=None,
                                        #normRad=False, cbar=True,
                                        #label=None)
        #plt.savefig('convergence_error_l2.pdf')


        return L2_norm_rel_error        
