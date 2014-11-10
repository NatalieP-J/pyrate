
from numpy import *
import math
from subprocess import call
import pickle
from scipy.interpolate import interp1d
from suppressor import RedirectStdStreams
import os
from matplotlib.backends.backend_pdf import PdfPages
import time

mos = {1:'Jan',2:'Feb',3:'Mar',4:'Apr',5:'May',6:'Jun',7:'Jul',8:'Aug',9:'Sep',10:'Oct',11:'Nov',12:'Dec'}

def gettime():
    t = str(time.gmtime(time.time())).split(' ')
    cyr = t[0][25:-1]
    cmo = mos[int(t[1][7:-1])]
    cda = t[2][8:-1]
    hr = t[3][8:-1]
    mi = t[4][7:-1]
    se = t[5][7:-1]
    dtimestr = '{0}-{1}-{2}.{3}.{4}.{5}'.format(cyr,cmo,cda,hr,mi,se)
    return dtimestr

devnull = open(os.devnull,'w')

Gconst = 6.67259e-8
realMsun = 1.989e33
Rsun = 6.9599e10
pc = 3.1e18
km = 10**5
yr = 365*24*3600
MsunV = 4.83
call(['mkdir','NukerRhoGals'],stdout = devnull,stderr = devnull)
call(['mkdir','SersicRhoGals'],stdout = devnull,stderr = devnull)

from construction import *

class NukerModeldIdR:
    """
    Call <modelname>.getrho() after model has been initialized as <modelname>
    """
    #initialize variables that constitute our model
    def __init__(self,model_name,alpha,beta,gamma,r0pc,rho0,MBH_Msun,generate=False,memo = False):
        #name of model
        self.name = model_name
        #Nuker fit parameters
        self.a = alpha
        self.b = beta
        self.g = gamma
        self.G = self.g-1
        self.B = self.b-1
        #starting radius
        self.r0 = r0pc
        #starting density
        self.rho0 = rho0
        #black holes mass in units of Msun
        self.MBH = MBH_Msun
        #Coulomb logarithm
        self.Lam = self.MBH*0.4
        #black hole mass normalized to galaxy density and radius
        self.Mnorm = self.MBH/(self.rho0*(self.r0)**3)
        #tidal disruption radius
        self.rT = Rsun*(self.MBH)**(1./3)
        #number of tidal radii to span galaxy
        self.r0_rT=(self.r0*pc)/self.rT
        #dynamical timescale (currently unused)
        self.tdyn0 = ((Gconst*self.rho0*realMsun)/pc**3)**(-1./2)
        #prefactor in I
        self.factor = 2**((self.b-self.g)/self.a)
        #start a new directory?
        self.generate = generate
        #directory name
        self.directory = 'NukerRhoGals/{0}_dIdR_a{1}_b{2}_g{3}_MBH{4}'.format(self.name,self.a,self.b,self.g,self.MBH)
        call(['mkdir','{0}'.format(self.directory)],stdout = devnull,stderr = devnull)
        timestr = gettime()
        self.statfile = open('{0}/stats_{1}.dat'.format(self.directory,timestr),'wb')
        self.pdfdump = PdfPages('{0}/{1}_master.pdf'.format(self.directory,self.name))
        self.memo = memo
        self.p1f = {}
        self.p2f = {}
        self.p3f = {}
        self.p1bG = {}
        self.p2bG = {}
        self.p3bG = {}
    
    #compute luminosity density
    def I(self,r):
        return self.factor*(r**-self.G)*(1+(r**self.a))**(-(self.B-self.G)/self.a)
    #and its first
    def dIdR(self,r):
        return -(self.I(r)/(r*(1+r**self.a)))*(self.G + self.B*r**self.a)

    #second
    def d2IdR2(self,r):
        return (self.I(r)/(r*(1+r**self.a))**2)*((r**(2*self.a))*self.B*(1+self.B) + self.G + (self.G**2) + (r**self.a)*(self.B-self.a*self.B + self.G + 2*self.B*self.G))
    
    #and third derivatives
    def d3IdR3(self,r):
        return (self.I(r)/(r*(1+r**self.a))**3)*((-r**(3*self.a))*self.B*(1+self.B)*(2+self.B) - self.G*(1 + self.G)*(2+self.G) + (r**(2*self.a))*((-1+self.a)*self.B*(4+self.a+3*self.B) - (2+(self.a**2) + 3*self.a*(1+self.B) + 3*self.B*(2+self.B))*self.G) + (r**self.a)*((1+self.a)*(-4+self.a-3*self.G)*self.G - self.B*(2+(self.a**2) - 3*self.a*(1+self.G) + 3*self.G*(2+self.G))))
    
    #use luminosity density to compute rho
    def rhointerior(self,theta,r):
        return self.dIdR(r/cos(theta))/cos(theta)

    def rho(self,r,verbose = False):
        try:
            r.shape
            rargs = [tuple((i,)) for i in r]
            lows = zeros(len(r))
            highs = zeros(len(r)) + pi/2.   
        except (AttributeError,TypeError):
            rargs = tuple((r,))
            lows = 0
            highs = pi/2.
        return -(1./pi)*integrator(r,[self.rhointerior,'rho'],lows,highs,args = rargs,fileobj=self.statfile)[0]
    
    #and its first
    def drhodrinterior(self,theta,r):
        return self.d2IdR2_2(r/cos(theta))/cos(theta)**2

    def drhodr(self,r,verbose = False):
        try:
            r.shape
            rargs = [tuple((i,)) for i in r]
            lows = zeros(len(r))
            highs = zeros(len(r)) + pi/2. 
        except (AttributeError,TypeError):
            rargs = tuple((r,))
            lows = 0
            highs = pi/2.
        return (1./pi)*integrator(r,[self.drhodrinterior,'drhodr'],lows,highs,args = rargs,fileobj=self.statfile)[0]
    
    #and second derivative
    def d2rhodr2interior(self,theta,r):
        return self.d3IdR3_2(r/cos(theta))/cos(theta)**3

    def d2rhodr2(self,r,verbose = False):
        try:
            r.shape
            rargs = [tuple((i,)) for i in r]
            lows = zeros(len(r))
            highs = zeros(len(r)) + pi/2. 
        except (AttributeError,TypeError):
            rargs = tuple((r,))
            lows = 0
            highs = pi/2.
        return (1./pi)*integrator(r,[self.d2rhodr2interior,'d2rhodr2'],lows,highs,args = rargs,fileobj = self.statfile)[0]


class NukerModelGenRho:
    """
    Call <modelname>.getrho() after model has been initialized as <modelname>
    """
    #initialize variables that constitute our model
    def __init__(self,model_name,alpha,beta,gamma,r0pc,rho0,MBH_Msun,generate=False,memo = False):
        #name of model
        self.name = model_name
        #Nuker fit parameters
        self.a = alpha
        self.b = beta
        self.g = gamma
        self.B = self.b-1
        self.G = self.g-1
        #starting radius
        self.r0 = r0pc
        #starting density
        self.rho0 = rho0
        #black holes mass in units of Msun
        self.MBH = MBH_Msun
        #Coulomb logarithm
        self.Lam = self.MBH*0.4
        #black hole mass normalized to galaxy density and radius
        self.Mnorm = self.MBH/(self.rho0*(self.r0)**3)
        #tidal disruption radius
        self.rT = Rsun*(self.MBH)**(1./3)
        #number of tidal radii to span galaxy
        self.r0_rT=(self.r0*pc)/self.rT
        #dynamical timescale (currently unused)
        self.tdyn0 = ((Gconst*self.rho0*realMsun)/pc**3)**(-1./2)
        #prefactor in I
        self.factor = 2**((self.b-self.g)/self.a)
        #start a new directory?
        self.generate = generate
        #directory name
        self.directory = 'NukerRhoGals/{0}_GenRho_a{1}_b{2}_g{3}_MBH{4}'.format(self.name,self.a,self.b,self.g,self.MBH)
        call(['mkdir','{0}'.format(self.directory)],stdout = devnull,stderr = devnull)
        timestr = gettime()
        self.statfile = open('{0}/stats_{1}.dat'.format(self.directory,timestr),'wb')
        self.pdfdump = PdfPages('{0}/{1}_master.pdf'.format(self.directory,self.name))
        self.memo = memo
        self.p1f = {}
        self.p2f = {}
        self.p3f = {}
        self.p1bG = {}
        self.p2bG = {}
        self.p3bG = {}
                               
    
        #compute luminosity density
    def I(self,r):
        return self.factor*(r**-self.G)*(1+(r**self.a))**(-(self.B-self.G)/self.a)
    #and its first
    def dIdR(self,r):
        return -(self.I(r)/(r*(1+r**self.a)))*(self.G + self.B*r**self.a)

    #second
    def d2IdR2(self,r):
        return (self.I(r)/(r*(1+r**self.a))**2)*((r**(2*self.a))*self.B*(1+self.B) + self.G + (self.G**2) + (r**self.a)*(self.B-self.a*self.B + self.G + 2*self.B*self.G))
    
    #and third derivatives
    def d3IdR3(self,r):
        return (self.I(r)/(r*(1+r**self.a))**3)*((-r**(3*self.a))*self.B*(1+self.B)*(2+self.B) - self.G*(1 + self.G)*(2+self.G) + (r**(2*self.a))*((-1+self.a)*self.B*(4+self.a+3*self.B) - (2+(self.a**2) + 3*self.a*(1+self.B) + 3*self.B*(2+self.B))*self.G) + (r**self.a)*((1+self.a)*(-4+self.a-3*self.G)*self.G - self.B*(2+(self.a**2) - 3*self.a*(1+self.G) + 3*self.G*(2+self.G))))
    
    #use luminosity density to compute rho
    def rhointerior(self,theta,r):
        return self.dIdR(r/cos(theta))/cos(theta)

    def funcrho(self,r,verbose = False):
        try:
            r.shape
            rargs = [tuple((i,)) for i in r]
            lows = zeros(len(r))
            highs = zeros(len(r)) + pi/2.   
        except (AttributeError,TypeError):
            rargs = tuple((r,))
            lows = 0
            highs = pi/2.
        return -(1./pi)*integrator(r,[self.rhointerior,'rho'],lows,highs,args = rargs,fileobj = self.statfile)[0]
    
    #and its first
    def drhodrinterior(self,theta,r):
        return self.d2IdR2(r/cos(theta))/cos(theta)**2

    def funcdrhodr(self,r,verbose = False):
        try:
            r.shape
            rargs = [tuple((i,)) for i in r]
            lows = zeros(len(r))
            highs = zeros(len(r)) + pi/2. 
        except (AttributeError,TypeError):
            rargs = tuple((r,))
            lows = 0
            highs = pi/2.
        return (1./pi)*integrator(r,[self.drhodrinterior,'drhodr'],lows,highs,args = rargs,fileobj = self.statfile)[0]
    
    #and second derivative
    def d2rhodr2interior(self,theta,r):
        return self.d3IdR3(r/cos(theta))/cos(theta)**3

    def funcd2rhodr2(self,r,verbose = False):
        try:
            r.shape
            rargs = [tuple((i,)) for i in r]
            lows = zeros(len(r))
            highs = zeros(len(r)) + pi/2. 
        except (AttributeError,TypeError):
            rargs = tuple((r,))
            lows = 0
            highs = pi/2.
        return (1./pi)*integrator(r,[self.d2rhodr2interior,'d2rhodr2'],lows,highs,args = rargs,fileobj = self.statfile)[0]
    
    def getrho(self,rtest = 0):
        if isinstance(rtest,(ndarray,)) == False:
            rtest = arange (-7,7,0.01)
            rtest = 10**rtest
            rtest1 = append(rtest,1e40)
            rtest1 = insert(rtest1,0,1e-40)
        if self.generate == True:
            tab1 = self.funcrho(rtest)
            tab2 = abs(self.funcdrhodr(rtest))
            tab3 = abs(self.funcd2rhodr2(rtest))
            inter0 = interp1d(log10(rtest),log10(tab1))
            inter1 = interp1d(log10(rtest),log10(tab2))
            inter2 = interp1d(log10(rtest),log10(tab3))
            piecerho = piecewise2(rtest1,inter0,tab1[0],tab1[len(tab1)-1],rtest[0],rtest[len(rtest)-1],-self.g,-self.b)
            self.inter0 = interp1d(log10(rtest1),log10(piecerho))
            piecedrhodr = piecewise2(rtest1,inter1,tab2[0],tab2[len(tab2)-1],rtest[0],rtest[len(rtest)-1],-self.g-1,-self.b-1)
            self.inter1 = interp1d(log10(rtest1),log10(piecedrhodr))
            pieced2rhodr2 = piecewise2(rtest1,inter2,tab3[0],tab3[len(tab3)-1],rtest[0],rtest[len(rtest)-1],-self.g-2,-self.b-2)
            self.inter2 = interp1d(log10(rtest1),log10(pieced2rhodr2))
            pklwrite('{0}/{1}.pkl'.format(self.directory,'rho'),column_stack((rtest1,piecerho)))
            pklwrite('{0}/{1}.pkl'.format(self.directory,'drhodr'),column_stack((rtest1,piecedrhodr)))
            pklwrite('{0}/{1}.pkl'.format(self.directory,'d2rhodr2'),column_stack((rtest1,pieced2rhodr2)))
            print 'Saved densities and interpolated'
            
        elif self.generate != True:
            rhovals = pklread('{0}/{1}.pkl'.format(self.directory,'rho'))
            rtest1 = rhovals[:,0]
            piecerho = rhovals[:,1]
            drhovals = pklread('{0}/{1}.pkl'.format(self.directory,'drhodr'))
            piecedrhodr = drhovals[:,1]
            d2rhovals = pklread('{0}/{1}.pkl'.format(self.directory,'d2rhodr2'))
            pieced2rhodr2 = d2rhovals[:,1]
            self.inter0 = interp1d(log10(rtest1),log10(piecerho))
            self.inter1 = interp1d(log10(rtest1),log10(piecedrhodr))
            self.inter2 = interp1d(log10(rtest1),log10(pieced2rhodr2))
            print 'Loaded densities and interpolated'
        
    def rho(self,r):
        return 10**self.inter0(log10(r))
    
    def drhodr(self,r):
        return -10**self.inter1(log10(r))
    
    def d2rhodr2(self,r):
        return 10**self.inter2(log10(r))


class NukerModelRho:
    #initialize variables that constitute our model
    def __init__(self,model_name,alpha,beta,gamma,r0pc,rho0,MBH_Msun,generate=False,memo = False):
        #name of model
        self.name = model_name
        #Nuker fit parameters
        self.a = alpha
        self.b = beta
        self.g = gamma
        #starting radius
        self.r0 = r0pc
        #starting density
        self.rho0 = rho0
        #black holes mass in units of Msun
        self.MBH = MBH_Msun
        #Coulomb logarithm
        self.Lam = self.MBH*0.4
        #black hole mass normalized to galaxy density and radius
        self.Mnorm = self.MBH/(self.rho0*(self.r0)**3)
        #tidal disruption radius
        self.rT = Rsun*(self.MBH)**(1./3)
        #number of tidal radii to span galaxy
        self.r0_rT=(self.r0*pc)/self.rT
        #dynamical timescale
        self.tdyn0 = ((Gconst*self.rho0*realMsun)/pc**3)**(-1./2)
        #start a new directory?
        self.generate = generate
        #directory name
        self.directory = 'NukerRhoGals/{0}_Rho_a{1}_b{2}_g{3}_MBH{4}'.format(self.name,self.a,self.b,self.g,self.MBH)
        call(['mkdir','{0}'.format(self.directory)],stdout = devnull,stderr = devnull)
        timestr = gettime()
        self.statfile = open('{0}/stats_{1}.dat'.format(self.directory,timestr),'wb')
        self.pdfdump = PdfPages('{0}/{1}_master.pdf'.format(self.directory,self.name))
        self.memo = memo
        self.p1f = {}
        self.p2f = {}
        self.p3f = {}
        self.p1bG = {}
        self.p2bG = {}
        self.p3bG = {}

    #compute density
    def rho(self,r):
        return (r**-self.g)*(1+r**self.a)**((self.g-self.b)/self.a)
    #and its first
    def drhodr(self,r):
        #return (-r**(-1-self.g))*((1+r**self.a)**((self.g-self.a-self.b)/self.a))*(self.g+self.b*r**self.a)
        return -self.rho(r)*((self.b*(r**self.a)+self.g)/(r+r**(1+self.a)))
    #and second derivatives
    def d2rhodr2(self,r):
        #part1 = r**(-2-self.g)
        #part2 = (1+r**self.a)**((self.g-self.b-2*self.a)/self.a)
        #part3a = self.b*(1+self.b)*r**(2*self.a)
        #part3b = self.g + self.g**2
        #part3c = (self.b - (self.a*self.b) + self.g + (self.a*self.g) + (2*self.b*self.g))*r**self.a
        #part3 = part3a + part3b + part3c
        #return part1*part2*part3
        nume = (1+self.b)*self.b*(r**(2*self.a)) + self.g + self.g**2 + (self.b-self.a*self.b + self.g + self.a*self.g + 2*self.b*self.g)*r**self.a
        deno = (r**2)*(1+r**self.a)**2
        return self.rho(r)*(nume/deno)

class SersicModeldIdR:
    #initialize variables that constitute the model
    def __init__(self,model_name,n,rho0,Re,M2L,I0,MBH_Msun,r0,generate,memo = False):
        #model name
        self.name = model_name
        #Sersic index
        self.n = n
        #Effective radius (half-light radius)
        self.Re = Re
        #Scale radius
        self.r0 = r0
        #Mass to light ratio
        self.M2L = M2L
        #Intensity at r0
        self.I0 = I0
        #Mass of the black hole in units of Msun
        self.MBH = MBH_Msun
        #determine other parameters based on n
        if self.n <10. and self.n > 0.6:
            self.b = 2*self.n - (1./3.) + 0.009876/self.n
            self.p = 1. - 0.6097/self.n + 0.05563/self.n**2
        #density at r0
        self.rho0 = M2L*I0*(self.b**(self.n*(1-self.p)))*(math.gamma(2*self.n)/(2*Re*math.gamma(self.n*(3-self.p))))
        #Coulomb logarithm
        self.Lam = self.MBH*0.4
        #black hole mass normalized to galaxy density and radius
        self.Mnorm = self.MBH/(self.rho0*(self.r0)**3)
        #tidal disruption radius
        self.rT = Rsun*(self.MBH)**(1./3)
        #number of tidal radii to span galaxy
        self.r0_rT=(self.r0*pc)/self.rT
        #dynamical timescale
        self.tdyn0 = ((Gconst*self.rho0*realMsun)/pc**3)**(-1./2)
        #start a new directory?
        self.generate = generate
        #directory name
        self.directory = 'SersicRhoGals/{0}_dIdR_n{1}_MBH{2}'.format(self.name,self.n.self.MBH)
        call(['mkdir','{0}'.format(self.directory)],stdout = devnull,stderr = devnull)
        timestr = gettime()
        self.statfile = open('{0}/stats_{1}.dat'.format(self.directory,timestr),'wb')
        self.pdfdump = PdfPages('{0}/{1}_master.pdf'.format(self.directory,self.name))
        self.memo = memo
        if self.memo == True:
            self.p1bG = {}
            self.p2bG = {}
            self.p3bG = {}

    def I(self,r):
        return exp(-self.b*(r**(1./self.n)))

    def dIdR(self,r):
        return -self.I(r)*(self.b/self.n)*r**(-1.+(1./n))

    def dI2dR2(self,r):
        return self.I(r)*(self.b/self.n**2)*(-1+self.n+self.b*r**(1./n))*r**(-2.+(1./n))
    
    def dI3dR3(self,r):
        return self.I(r)*(self.b/self.n**3)*(-1+3*self.n-2*self.n**2-3*self.b*(-1+n)*r**(1./n)-k**2*r**(2./n))*r**(-3+(1./n))

    #use luminosity density to compute rho
    def rhointerior(self,theta,r):
        return self.dIdR(r/cos(theta))/cos(theta)

    def rho(self,r,verbose = False):
        try:
            r.shape
            rargs = [tuple((i,)) for i in r]
            lows = zeros(len(r))
            highs = zeros(len(r)) + pi/2.   
        except (AttributeError,TypeError):
            rargs = tuple((r,))
            lows = 0
            highs = pi/2.
        return -(1./pi)*integrator(r,[self.rhointerior,'rho'],lows,highs,args = rargs,fileobj = self.statfile)[0]
    
    #and its first
    def drhodrinterior(self,theta,r):
        return self.d2IdR2_2(r/cos(theta))/cos(theta)**2

    def drhodr(self,r,verbose = False):
        try:
            r.shape
            rargs = [tuple((i,)) for i in r]
            lows = zeros(len(r))
            highs = zeros(len(r)) + pi/2. 
        except (AttributeError,TypeError):
            rargs = tuple((r,))
            lows = 0
            highs = pi/2.
        return (1./pi)*integrator(r,[self.drhodrinterior,'drhodr'],lows,highs,args = rargs,fileobj = self.statfile)[0]
    
    #and second derivative
    def d2rhodr2interior(self,theta,r):
        return self.d3IdR3_2(r/cos(theta))/cos(theta)**3

    def d2rhodr2(self,r,verbose = False):
        try:
            r.shape
            rargs = [tuple((i,)) for i in r]
            lows = zeros(len(r))
            highs = zeros(len(r)) + pi/2. 
        except (AttributeError,TypeError):
            rargs = tuple((r,))
            lows = 0
            highs = pi/2.
        return (1./pi)*integrator(r,[self.d2rhodr2interior,'d2rhodr2'],lows,highs,args = rargs,fileobj = self.statfile)[0]

class SersicModelGenRho:
    #initialize variables that constitute the model
    def __init__(self,model_name,n,rho0,Re,M2L,I0,MBH_Msun,r0,generate,memo = False):
        #model name
        self.name = model_name
        #Sersic index
        self.n = n
        #Effective radius (half-light radius)
        self.Re = Re
        #Scale radius
        self.r0 = r0
        #Mass to light ratio
        self.M2L = M2L
        #Intensity at r0
        self.I0 = I0
        #Mass of the black hole in units of Msun
        self.MBH = MBH_Msun
        #determine other parameters based on n
        if self.n <10. and self.n > 0.6:
            self.b = 2*self.n - (1./3.) + 0.009876/self.n
            self.p = 1. - 0.6097/self.n + 0.05563/self.n**2
        #density at r0
        self.rho0 = M2L*I0*(self.b**(self.n*(1-self.p)))*(math.gamma(2*self.n)/(2*Re*math.gamma(self.n*(3-self.p))))
        #Coulomb logarithm
        self.Lam = self.MBH*0.4
        #black hole mass normalized to galaxy density and radius
        self.Mnorm = self.MBH/(self.rho0*(self.r0)**3)
        #tidal disruption radius
        self.rT = Rsun*(self.MBH)**(1./3)
        #number of tidal radii to span galaxy
        self.r0_rT=(self.r0*pc)/self.rT
        #dynamical timescale
        self.tdyn0 = ((Gconst*self.rho0*realMsun)/pc**3)**(-1./2)
        #start a new directory?
        self.generate = generate
        #directory name
        self.directory = 'SersicRhoGals/{0}_GenRho_n{1}_MBH{2}'.format(self.name,self.n.self.MBH)
        call(['mkdir','{0}'.format(self.directory)],stdout = devnull,stderr = devnull)
        timestr = gettime()
        self.statfile = open('{0}/stats_{1}.dat'.format(self.directory,timestr),'wb')
        self.pdfdump = PdfPages('{0}/{1}_master.pdf'.format(self.directory,self.name))
        self.memo = memo
        if self.memo == True:
            self.p1bG = {}
            self.p2bG = {}
            self.p3bG = {}

    def I(self,r):
        return exp(-self.b*(r**(1./self.n)))

    def dIdR(self,r):
        return -self.I(r)*(self.b/self.n)*r**(-1.+(1./n))

    def dI2dR2(self,r):
        return self.I(r)*(self.b/self.n**2)*(-1+self.n+self.b*r**(1./n))*r**(-2.+(1./n))
    
    def dI3dR3(self,r):
        return self.I(r)*(self.b/self.n**3)*(-1+3*self.n-2*self.n**2-3*self.b*(-1+n)*r**(1./n)-k**2*r**(2./n))*r**(-3+(1./n))

    #use luminosity density to compute rho
    def rhointerior(self,theta,r):
        return self.dIdR(r/cos(theta))/cos(theta)

    def funcrho(self,r,verbose = False):
        try:
            r.shape
            rargs = [tuple((i,)) for i in r]
            lows = zeros(len(r))
            highs = zeros(len(r)) + pi/2.   
        except (AttributeError,TypeError):
            rargs = tuple((r,))
            lows = 0
            highs = pi/2.
        return -(1./pi)*integrator(r,[self.rhointerior,'rho'],lows,highs,args = rargs,fileobj = self.statfile)[0]
    
    #and its first
    def drhodrinterior(self,theta,r):
        return self.d2IdR2_2(r/cos(theta))/cos(theta)**2

    def funcdrhodr(self,r,verbose = False):
        try:
            r.shape
            rargs = [tuple((i,)) for i in r]
            lows = zeros(len(r))
            highs = zeros(len(r)) + pi/2. 
        except (AttributeError,TypeError):
            rargs = tuple((r,))
            lows = 0
            highs = pi/2.
        return (1./pi)*integrator(r,[self.drhodrinterior,'drhodr'],lows,highs,args = rargs,fileobj = self.statfile)[0]
    
    #and second derivative
    def d2rhodr2interior(self,theta,r):
        return self.d3IdR3_2(r/cos(theta))/cos(theta)**3

    def funcd2rhodr2(self,r,verbose = False):
        try:
            r.shape
            rargs = [tuple((i,)) for i in r]
            lows = zeros(len(r))
            highs = zeros(len(r)) + pi/2. 
        except (AttributeError,TypeError):
            rargs = tuple((r,))
            lows = 0
            highs = pi/2.
        return (1./pi)*integrator(r,[self.d2rhodr2interior,'d2rhodr2'],lows,highs,args = rargs,fileobj = self.statfile)[0]
    
    def getrho(self):
        rtest = arange(-7,7,0.01)
        rtest = append(rtest,40)
        rtest = insert(rtest,0,-40)
        rtest = 10**rtest
        if self.generate == True:
            call(['mkdir','{0}'.format(self.directory)])
            pklrfile = open('{0}/{1}.pkl'.format(self.directory,'rrho'),"wb")
            pickle.dump(rtest,pklrfile)
            pklrfile.close()
            tab1 = self.funcrho(rtest)
            pklrfile = open('{0}/{1}.pkl'.format(self.directory,'rho'),"wb")
            pickle.dump(tab1,pklrfile)
            pklrfile.close()
            tab2 = abs(self.funcdrhodr(rtest))
            pklrfile = open('{0}/{1}.pkl'.format(self.directory,'drhodr'),"wb")
            pickle.dump(tab2,pklrfile)
            pklrfile.close()
            tab3 = abs(self.funcd2rhodr2(rtest))
            pklrfile = open('{0}/{1}.pkl'.format(self.directory,'d2rhodr2'),"wb")
            pickle.dump(tab3,pklrfile)
            pklrfile.close()
            self.inter0 = interp1d(log10(rtest),log10(tab1))
            self.inter1 = interp1d(log10(rtest),log10(tab2))
            self.inter2 = interp1d(log10(rtest),log10(tab3))
            print 'Saved densities and interpolated'
            
        elif self.generate != True:
            pklrfile = open('{0}/{1}.pkl'.format(self.directory,'rrho'),"rb")
            rtest = pickle.load(pklrfile)
            pklrfile.close()
            pklrfile = open('{0}/{1}.pkl'.format(self.directory,'rho'),"rb")
            tab1 = pickle.load(pklrfile)
            pklrfile.close()
            pklrfile = open('{0}/{1}.pkl'.format(self.directory,'drhodr'),"rb")
            tab2 = pickle.load(pklrfile)
            pklrfile.close()
            pklrfile = open('{0}/{1}.pkl'.format(self.directory,'d2rhodr2'),"rb")
            tab3 = pickle.load(pklrfile)
            pklrfile.close()
            self.inter0 = interp1d(log10(rtest),log10(tab1))
            self.inter1 = interp1d(log10(rtest),log10(tab2))
            self.inter2 = interp1d(log10(rtest),log10(tab3))
            print 'Loaded densities and interpolated'
        
    def rho(self,r):
        return 10**self.inter0(log10(r))
    
    def drhodr(self,r):
        return -10**self.inter1(log10(r))
    
    def d2rhodr2(self,r):
        return 10**self.inter2(log10(r))

class SersicModelRho:
    #initialize variables that constitute the model
    def __init__(self,model_name,n,rho0,Re,M2L,I0,MBH_Msun,r0,generate,memo = False):
        #model name
        self.name = model_name
        #Sersic index
        self.n = n
        #Effective radius (half-light radius)
        self.Re = Re
        #Scale radius
        self.r0 = r0
        #Mass to light ratio
        self.M2L = M2L
        #Intensity at r0
        self.I0 = I0
        #Mass of the black hole in units of Msun
        self.MBH = MBH_Msun
        #determine other parameters based on n
        if self.n <10. and self.n > 0.6:
            self.b = 2*self.n - (1./3.) + 0.009876/self.n
            self.p = 1. - 0.6097/self.n + 0.05563/self.n**2
        #density at r0
        self.rho0 = M2L*I0*(self.b**(self.n*(1-self.p)))*(math.gamma(2*self.n)/(2*Re*math.gamma(self.n*(3-self.p))))
        #Coulomb logarithm
        self.Lam = self.MBH*0.4
        #black hole mass normalized to galaxy density and radius
        self.Mnorm = self.MBH/(self.rho0*(self.r0)**3)
        #tidal disruption radius
        self.rT = Rsun*(self.MBH)**(1./3)
        #number of tidal radii to span galaxy
        self.r0_rT=(self.r0*pc)/self.rT
        #dynamical timescale
        self.tdyn0 = ((Gconst*self.rho0*realMsun)/pc**3)**(-1./2)
        #start a new directory?
        self.generate = generate
        #directory name
        self.directory = 'SersicRhoGals/{0}_Rho_n{1}_MBH{2}'.format(self.name,self.n.self.MBH)
        call(['mkdir','{0}'.format(self.directory)],stdout = devnull,stderr = devnull)
        timestr = gettime()
        self.statfile = open('{0}/stats_{1}.dat'.format(self.directory,timestr),'wb')
        self.pdfdump = PdfPages('{0}/{1}_master.pdf'.format(self.directory,self.name))
        self.memo = memo
        if self.memo == True:
            self.p1bG = {}
            self.p2bG = {}
            self.p3bG = {}
    
    #Compute density
    def rho(self,r):
        return (r**-self.p)*exp(-self.b*(r**(1./self.n)))
    #and its first
    def drhodr(self,r):
        pre = (r**-self.p)*exp(-self.b*(r**(1./self.n)))*((self.n*r)**-1)
        post = (self.n*self.p)+self.b*(r**(1./n))
        return pre*post
    #and second derivatives
    def d2rhodr2(self,r):
        pre = (r**-p)*exp(-b*(r**(1./self.n)))*((self.n*r)**-2)
        post = (self.p*(1+self.p)*self.n**2) + self.b*(-1 + self.n + 2*self.n*self.p)*(r**(1./n)) + (b**2)*(r**-p)
        return pre*post
