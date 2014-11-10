from numpy import *
from scipy.interpolate import interp1d
import time
import pickle

#IMPORTS DATA TABLES AND PROVIDES FUNCTIONAL FORMS FOR COMPUTATIONS

class BesselGen:
    """
    Generates interpolated functional forms from data tables for computations.
    """
    def __init__(self,files):
        """
        Initializes class instance and imports data files from list <files>.
        Returns interpolated forms.
        """
        datafile = open('{0}'.format(files[0]),"rb")
        data = pickle.load(datafile)
        self.ainter = interp1d(data[:,0],data[:,1])
        data = loadtxt('{0}'.format(files[1]))
        self.xinter = interp1d(data[:,0],data[:,1])
        data = loadtxt('{0}'.format(files[2]))
        self.binter = interp1d(data[:,0],data[:,1])
        data = loadtxt('{0}'.format(files[3]))
        self.minter = interp1d(data[:,0],data[:,1])
        
    def alpham(self,m):
        """
        m - index of a zero of the J0 bessel function
        Returns mth zero of J0 bessel function (with approximations for m>1 million).
        """
        if isinstance(m,(ndarray,list))==True:
            alpha1 = self.ainter(m[(m<10**6)])
            alpha2 = pi*(m[(m>=10**6)] - 0.25)
            return concatenate((alpha1,alpha2))
        elif isinstance(m,(float,int))==True:
            if m<10**6:
                return self.ainter(m)
            elif m>=10**6:
                return pi*(m-0.25)
        
    def xi(self,q):
        """
        q - ratio of theta_dyn to theta_lc squared
        Returns piecewise definition of xi(q).
        """
        qimin = -10.
        qimax = 1.
        dqi = 0.03
        xichange = 10**(floor((qimax-qimin)/dqi)*dqi + qimin)
        if isinstance(q,(ndarray,list))==True:
            xi1 = sqrt(4*q[(q<10**qimin)]/pi)
            xi2 = self.xinter(q[(q>=10**qimin)&(q<=xichange)])
            xi3 = q[(q>xichange)]/q[(q>xichange)]
            return concatenate((xi3,xi2,xi1))
        elif isinstance(q,(float,int))==True:
            if q < 10**qimin:
                return array(sqrt(4*q/pi))
            elif q >= 10**qimin and q < xichange:
                return self.xinter(q)
            elif q >= xichange:
                return array(1.)

    def besselfin(self,m,u):
        """
        m - index of a zero of J0 bessel function
        u - square root of pericentre radius over tidal disruption radius
        Returns J0(alpham*u).
        """
        zimin = 10**-2
        zimax = 10**2
        dzi = 10**-2
        z = self.alpham(m)*u
        if isinstance(z,(ndarray,list))==True:
            besselfin1 = z[(z<zimin)]/z[(z<zimin)]
            besselfin2 = self.binter(z[(z>=zimin)&(z<=zimax)])
            besselfin3 = sqrt(2./((m[(z>zimax)]-0.25)*u*pi**2))*cos(pi*(u*(m[(z>zimax)]-0.25)-0.25))
            return concatenate((besselfin1,besselfin2,besselfin3))
        elif isinstance(z,(float,int))==True:
            if z<zimin:
                return 1
            elif z>=zimin and z<=zimax:
                return self.binter(z)
            elif z>zimax:
                return sqrt(2./((m-0.25)*u*pi**2))*cos(pi*(u*(m-0.25)-0.25))

    def mpiece(self,m):
        """
        m - index of a zero of J0 bessel function
        Returns J1(alpham).
        """
        mlim = 200.
        if isinstance(m,(ndarray,list))==True:
            mpiece1 = self.minter(m[(m<=mlim)])
            mpiece2 = (-1)**(m[m>mlim]-1)*sqrt(2*m[m>mlim] - 0.5)
            return concatenate((mpiece1,mpiece2))
        elif isinstance(m,(float,int))==True:
            if m<=mlim:
                return self.minter(m)
            elif m>mlim:
                return (-1)**(m-1)*sqrt(2*m - 0.5)
