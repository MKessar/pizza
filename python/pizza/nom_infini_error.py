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


class maxerror_instdef(PizzaSetup):
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

        ##filename = self.get_filename('frame_us', ivar, datadir, tag, verbose)
        #f1 = Frame(filename1, endian=endian)
        #f2 = Frame(filename2, endian=endian)
        ##print "spec_inst a shape (f.field_m)=", np.shape(f.field_m)

        ##self.f_m = f.field_m
        #self.n_r_max1 = f1.n_r_max
        #self.n_m_max1 = f1.n_m_max
        #print "n_r_max1",self.n_r_max1
        #print "n_m_max1",self.n_m_max1
        #self.n_phi1= (self.n_m_max1-1)*2
        #self.dphi1= 2.0*np.pi/(self.n_phi1)
        #print "self.n_phi1=", self.n_phi1


        #self.n_r_max2 = f2.n_r_max
        #self.n_m_max2 = f2.n_m_max
        #print "n_r_max2",self.n_r_max2
        #print "n_m_max2",self.n_m_max2
        #self.n_phi2 = (self.n_m_max2-1)*2
        #self.dphi2  = 2.0*np.pi/(self.n_phi2)
        #print "self.n_phi2=", self.n_phi2
        
        ##self.uphi_m = f1.field_m
        #f1_r = spec_spat(f1.field_m, self.n_phi1)
        #f2_r = spec_spat(f2.field_m, self.n_phi2)

        #print "made it here"

        ##f_error_real_abs=abs(f1_r - f2_r)
        #f_error_real_abs=2.0*f2_r
        #for i in range(0,self.n_phi1):
            #for k in range(0,self.n_r_max1):
                ##if (abs(f_error_real_abs[i,k]))>max_f_error_real_abs:
                   #f_error_real_abs[i,k] = abs(f1_r[i,k] - f2_r[i,k])
                   

        #print "made it here a"

        #max_f_error_real_abs=0.000000000000000
        
        #for i in range(0,self.n_phi1):
            #for k in range(0,self.n_r_max1):
                #if (abs(f_error_real_abs[i,k]))>max_f_error_real_abs:
                   #max_f_error_real_abs = f_error_real_abs[i,k]

        #print "max_f_error_real_abs", max_f_error_real_abs

                
        #f_error_real_rel=2.0*f2_r
        #for i in range(0,self.n_phi1):
            #for k in range(0,self.n_r_max1):
                ##self.f1_bis[i,k] = 0.000000000000000
                #if (abs(f1_r[i,k]))>(0.00000000001):
                   #f_error_real_rel[i,k] = f_error_real_abs[i,k]/f1_r[i,k]
                #if (abs(f1_r[i,k]))<=(0.00000000001):
                    ##print "i=",i,"k=", k,"f1_cheb_abs_real=0.0"
                    #f_error_real_rel[i,k]=0.0
        #max_f_error_real_rel=0.000000000000000

        #for i in range(0,self.n_phi1):
            #for k in range(0,self.n_r_max1):
                ##print "f_error_real_rel(i=",i,",k=",k,")=",f_error_real_rel[i,k]
                #if (abs(f_error_real_rel[i,k]))>max_f_error_real_rel:
                   #max_f_error_real_rel = f_error_real_rel[i,k]

        ##max_f_error_real_rel=max(abs(f1_r - f2_r))
                    
        #print "max_f_error_real_rel", max_f_error_real_rel
        #return max_f_error_real_rel

    def max_err_rel(self,filename1=None,filename2=None, endian='l'):

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
        #print "n_r_max1",self.n_r_max1
        #print "n_m_max1",self.n_m_max1
        self.n_phi1= (self.n_m_max1-1)*2
        self.dphi1= 2.0*np.pi/(self.n_phi1)
        #print "self.n_phi1=", self.n_phi1


        self.n_r_max2 = f2.n_r_max
        self.n_m_max2 = f2.n_m_max
        #print "n_r_max2",self.n_r_max2
        #print "n_m_max2",self.n_m_max2
        self.n_phi2 = (self.n_m_max2-1)*2
        self.dphi2  = 2.0*np.pi/(self.n_phi2)
        #print "self.n_phi2=", self.n_phi2
        
        #self.uphi_m = f1.field_m
        f1_r = spec_spat(f1.field_m, self.n_phi1)
        f2_r = spec_spat(f2.field_m, self.n_phi2)

        #print "made it here"

        #f_error_real_abs=abs(f1_r - f2_r)
        f_error_real_abs=2.0*f2_r
        for i in range(0,self.n_phi1):
            for k in range(0,self.n_r_max1):
                #if (abs(f_error_real_abs[i,k]))>max_f_error_real_abs:
                   f_error_real_abs[i,k] = abs(f1_r[i,k] - f2_r[i,k])
                   #print "f_error_real_abs",i,k, f_error_real_abs[i,k]
                   #print "f1_r",i,k, f1_r[i,k]
                   #print "f2_r",i,k, f2_r[i,k]


        #print "made it here a"

        max_f_error_real_abs=0.000000000000000
        
        for i in range(0,self.n_phi1):
            for k in range(0,self.n_r_max1):
                if (abs(f_error_real_abs[i,k]))>max_f_error_real_abs:
                   max_f_error_real_abs = f_error_real_abs[i,k]

        print "max_f_error_real_abs", max_f_error_real_abs

                
        f_error_real_rel=2.0*f2_r
        for i in range(0,self.n_phi1):
            for k in range(0,self.n_r_max1):
                #self.f1_bis[i,k] = 0.000000000000000
                if (abs(f1_r[i,k]))>(0.000000000001):
                   f_error_real_rel[i,k] = f_error_real_abs[i,k]/f1_r[i,k]
                if (abs(f1_r[i,k]))<=(0.000000000001):
                    #print "i=",i,"k=", k,"f1_cheb_abs_real=0.0"
                    f_error_real_rel[i,k]=0.0
        max_f_error_real_rel=0.000000000000000

        for i in range(0,self.n_phi1):
            for k in range(0,self.n_r_max1):
                #print "f_error_real_rel(i=",i,",k=",k,")=",f_error_real_rel[i,k]
                if (abs(f_error_real_rel[i,k]))>max_f_error_real_rel:
                   max_f_error_real_rel = f_error_real_rel[i,k]

        #max_f_error_real_rel=max(abs(f1_r - f2_r))
                    
        print "max_f_error_real_rel", max_f_error_real_rel
        return max_f_error_real_rel        
