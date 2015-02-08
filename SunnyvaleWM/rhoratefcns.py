import numpy as np
nsum = np.sum
from numpy import *
import scipy.integrate as intg
from scipy.optimize import root
import matplotlib.pyplot as plt
from besselgen import BesselGen
from construction import *
from subprocess import call
from suppressor import RedirectStdStreams

def findrho0(rb,M2L,mub):
    MsunV = 4.83
    return (1./rb)*(1./(10)**2)*(206265**2)*M2L*10**((MsunV-mub)/2.5) 

########******************* ENCLOSED MASS *******************########

def Minterior(r,prereqs):
    """
    interior of the Menc integral
    """
    model, = prereqs
    return model.rho(r)*r**2

def funcMenc(r,prereqs):
    """
    functional form of Menc
    relies on Minterior
    returns Menc(r)
    """
    model, = prereqs
    if isinstance(r,(list,ndarray))==True:
        dls = zeros(len(r))
        uls = r
        plist = [tuple((prereqs,))]*len(r)
        pre = [4*pi]*len(r)
    elif isinstance(r,(int,float))==True:
        dls = 0
        uls = r
        plist = tuple((prereqs,))
        pre = 4*pi
    return integrator(r,[Minterior,'Menc'],dls,uls,args = plist,fileobj = model.statfile,prefactor = pre)

########******************* RADIUS OF INFLUENCE *******************######## 

def rHimplicit(r,prereqs):
    """
    equation that has its minimum when r = rH
    """
    model, = prereqs
    return abs(model.Mnorm-funcMenc(abs(r),prereqs)[0])

def rH(prereqs):
    """
    finds root of rHimplicit
    """
    model, = prereqs
    rresult=root(rHimplicit,1e-4,args = prereqs)
    if rresult.success == True:
        return abs(rresult.x)
    elif rresult.success == False:
        model.statfile.write('\nFailed to evaluate rH\n')
        model.statfile.write('{0}\n'.format(rresult.message))
        return abs(rresult.x)

########******************* GRID CREATION  *******************######## 

def stdgrid(prereqs,upstep=5,downstep=-5,step=0.03):
    """
    constructs a grid in radius and adds one point at each extreme (up and down)
    returns 10**grid
    """
    model, = prereqs
    rarray = arange(downstep,upstep,step)
    rarray = 10**rarray
    rchange = rarray[len(rarray)-1]
    rstart = rarray[0]
    return rarray,rchange,rstart

def rgrid(prereqs,upstep=5,downstep=-5,step=0.03):
    """
    constructs a grid in radius and adds one point at each extreme (up and down)
    returns 10**grid
    """
    model, = prereqs
    rmin = min([rH([model]),[1.]])
    rmax = max([rH([model]),[1.]])
    rimin = log10(rmin) + downstep
    rimax = log10(rmax) + upstep
    dri = step
    rarray = arange(rimin,rimax,dri)
    rarray = 10**rarray
    rchange = rarray[len(rarray)-1]
    rstart = rarray[0]
    return rarray,rchange,rstart

def Egrid(prereqs,upstep=5,downstep=-3,step=0.1):
    """
    constructs a grid in energy and adds one point at each extreme (up and down)
    returns 10**grid
    """
    model, = prereqs
    rmin = min([rH([model]),[1.]])[0]
    rmax = max([rH([model]),[1.]])[0]
    eimin = log10(funcMenc(rmax,prereqs)[0]/rmax) + downstep
    eimax = log10(model.Mnorm/rmin) + upstep
    dei = step
    Earray = arange(eimin,eimax,dei)
    Earray = 10**Earray
    Echange = Earray[len(Earray)-1]
    Estart = Earray[0]
    return Earray,Echange,Estart

########******************* POTENTIAL *******************######## 
        
def psi2interior(r,prereqs):
    """
    interior of psi part 2 integral
    """
    model,Mencgood = prereqs
    return model.rho(r)*r

tol = 1e-3
def psi2(r,prereqs):
    """
    functional form of psi part 2
    relies on psi2interior
    returns psi2(r)
    """
    model,Mencgood= prereqs
    if isinstance(r,(list,ndarray))==True:
        dls = r
        uls = zeros(len(r)) + inf
        plist = [tuple((prereqs,))]*len(r)
        pre = [4*pi]*len(r)
    elif isinstance(r,(int,float))==True:
        dls = r
        uls = inf
        plist = tuple((prereqs,))
        pre = 4*pi
    return integrator(r,[psi2interior,'psi2'],dls,uls,args = plist,fileobj = model.statfile,prefactor = pre) 

def funcpsi(r,prereqs):
    """
    returns potential as a function of r
    """
    model,Mencgood = prereqs
    part1 = (model.Mnorm/r)
    part2 = 10**Mencgood(log10(r))
    part3 =  psi2(r,prereqs)
    problems = array(part3[1])
    #problems = array([])
    #if part2[1] != []:
     #   problems = concatenate((problems,array(part2[1])))
    #if part3[1] != []:
    #    problems = concatenate((problems,array([val for val in part3[1] if val not in part2[1]])))
    return part1 + (part2/r) + part3[0],problems

########******************* APOCENTER RADIUS *******************######## 

def rapoimplicit(r,E,prereqs):
    """
    function with a minimum at r=rapo
    """
    psigood, = prereqs
    return abs(10**psigood(log10(abs(r)))-E)

def rapo(E,prereqs):
    """
    finds root of rapoimplicit
    """
    if E**-1 > 0.2:
        rguess = 10*E**-1
    elif E**-1 < 0.2:
        rguess = 0.01*E**-1
    rresult = root(rapoimplicit,rguess,args=(E,prereqs))
    if rresult.success == True:
        return abs(rresult.x)
    elif rresult.success == False:
        print 'Failed to evaluate rapo'
        print rresult.message
        return abs(rresult.x)

########******************* CIRCULAR ANGULAR MOMENTUM *******************######## 

def Jc2implicit(r,E,prereqs):
    """
    """
    model,Mencgood,psigood = prereqs
    return abs(10**psigood(log10(abs(r)))-E-((10**Mencgood(log10(abs(r)))+model.Mnorm)/(2*r)))

def funcJc2(E,prereqs):
    """
    see Jc2implicit
    """
    model,Mencgood,psigood = prereqs
    try:
        t = E.shape
        Jcs = []
        problems = []
        for i in range(len(E)):
            if E[i]**-1 > 0.4:
                rguess = 10*E[i]**-1
            elif E[i]**-1 < 0.4:
                rguess = 0.01*E[i]**-1
            rresult = root(Jc2implicit,rguess,args=(E[i],prereqs))
            rresult.x = abs(rresult.x)
            if rresult.success == True:
                Jcs.append(((10**Mencgood(log10(rresult.x))+model.Mnorm)*rresult.x)[0])
            elif rresult.success == False:
                model.statfile.write('\nFailed to evaluate Jc2')
                model.statfile.write('\n{0}'.format(rresult.message))
                Jcs.append(((10**Mencgood(log10(rresult.x))+model.Mnorm)*rresult.x)[0])
                problems.append(i)
        return array(Jcs),problems
    except AttributeError:
        rresult = root(Jc2implicit,0.01*E**-1,args=(E,prereqs))
        problem = []
        rresult.x = abs(rresult.x)
        if rresult.success == True:
            Jc = ((10**Mencgood(log10(rresult.x))+model.Mnorm)*rresult.x)[0]
        elif rresult.success == False:
            model.statfile.write('\nFailed to evaluate Jc2')
            model.statfile.write('\n{0}'.format(rresult.message))
            Jc = ((10**Mencgood(log10(rresult.x))+model.Mnorm)*rresult.x)[0]
            problem = [E]
        return Jc, problem

########******************* g *******************######## 

def lginterior(r,E,prereqs):
    """
    interior of g integral
    """
    model,psigood = prereqs
    return (model.drhodr(1./r))*(r**-2)*((sqrt(abs(E-10**psigood(log10(1./r)))))**-1)
    
def funclg(E,prereqs):
    """
    functional form of g
    relies on ginterior
    returns g(E)
    """
    model,psigood = prereqs
    if isinstance(E,(list,ndarray))==True:
        dls = zeros(len(E))
        uls = []
        plist = []
        for index in range(len(E)):
            uls.append(1./rapo(E[index],[psigood]))
            plist.append(tuple((E[index],prereqs)))
        pre = [-pi]*len(E)
    elif isinstance(E,(int,float))==True:
        dls = 0
        uls = 1./rapo(E,[psigood])
        plist = tuple((E,prereqs))
        pre = -pi
    return integrator(E,[lginterior,'g'],dls,uls,args = plist,fileobj = model.statfile,prefactor = pre) 

########******************* mathcalG *******************########

def bGinterior(theta,r,E,prereqs):
    """
    interior of G integral
    """
    model,psigood,ggood = prereqs
    if model.memo == True:
        psir = 10**psigood(log10(r))
        if not r in model.p1bG:
            model.p1bG[r] = (r**2)/sqrt(psir-E)
        part1 = model.p1bG[r]
        if not log10(psir*(1-theta) + E*theta) in model.p2bG:
            model.p2bG[log10(psir*(1-theta) + E*theta)] = 10**ggood(log10(psir*(1-theta) + E*theta))
        part2 = model.p2bG[log10(psir*(1-theta) + E*theta)]
        if not theta in model.p3bG:
            model.p3bG[theta]= (1./sqrt(theta))-sqrt(theta)
        part3 = model.p3bG[theta]
    if model.memo == False:
        psir = 10**psigood(log10(r))
        part1 = (r**2)/sqrt(psir-E)
        part2 = 10**ggood(log10(psir*(1-theta) + E*theta))
        part3 =(1./sqrt(theta))-sqrt(theta) 
    return part1*part2*part3

def funcbG(E,prereqs):
    """
    functional form of mathcalG
    relies on Ginterior
    returns mathcalG(E)
    """
    model,psigood,ggood = prereqs
    try:
        t = E.shape
        Gans = []
        problems = []
        for i in range(len(E)):
            print i+1, 'of', len(E)
            rapoval = rapo(E[i],[psigood])
            try:
                temp = intg.dblquad(bGinterior,0,rapoval,lambda r: 1e-4, lambda r: 1,args = (E[i],prereqs))
            except UserWarning as e:
                if verbose == True:
                    print 'G, E = ', E[i], 'message = ', e
                problems.append(i)
            Gans.append(temp[0])
        return array(Gans),problems
    except AttributeError:
        rapoval = rapo(E,psigood)
        problem = []
        try:
            temp = intg.dblquad(bGinterior,0,rapoval,lambda r: 0, lambda r: 1,args = (E,prereqs))
        except UserWarning as e:
            if verbose == True:
                print 'G, E = ', E, 'message = ', temp[3]
            problem = [E]
        return temp[0],problem

def funcbG2(E,prereqs):
    model,psigood,ggood = prereqs
    if isinstance(E,(list,ndarray))==True:
        dls = []
        dls.append(zeros(len(E)))
        dls.append(array([lambda r:1e-4]*len(E)))
        uls = []
        uls1 = array([])
        uls2 = array([lambda r:1]*len(E))
        plist = []
        for i in range(len(E)):
            rapoval = rapo(E[i],[psigood])
            uls1 = append(uls1,rapoval)
            plist.append(tuple((E[i],prereqs)))
        uls.append(uls1)
        uls.append(uls2)
    elif isinstance(E,(int,float))==True:
        dls = [0,lambda r:1e-4]
        uls = [rapo(E[i],[psigood]),lambda r:1]
        plist = tuple((E,prereqs))
    return dblintegrator(E,[bGinterior,'G'],dls,uls,args = plist,fileobj = model.statfile)
'''
def funcbG3(E,prereqs):
    model,psigood,ggood = prereqs
    if isinstance(E,(list,ndarray))==True:
        dls1 = zeros(len(E))
        dls2 = array([lambda r:1e-4]*len(E))
        uls1 = array([])
        uls2 = array([lambda r:1]*len(E))
        plist = []
        for i in range(len(E)):
            rapoval = rapo(E[i],[psigood])
            uls1 = append(uls1,rapoval)
            plist.append(tuple((E[i],prereqs)))
    elif isinstance(E,(int,float))==True:
        dls = [0,lambda r:1e-4]
        uls = [rapo(E[i],[psigood]),lambda r:1]
        plist = tuple((E,prereqs))
    return dblintegrator(E,[bGinterior,'G'],dls,uls,args = plist,fileobj = model.statfile)
'''
########******************* DISTRIBUTION FUNCTION *******************######## 

def finterior(r,E,rapoval,prereqs):
    model,Mencgood,psigood = prereqs
    var = float(rapoval/r)
    psi = (10**psigood(log10(var)))
    Mencvar = (10**Mencgood(log10(var)))
    Mencrap = (10**Mencgood(log10(rapoval)))
    if model.memo == True:
        if not var in model.p1f:
            model.p1f[var] = (var**3)*(1./sqrt(abs(E-psi)))*model.d2rhodr2(var) 
        result1 = model.p1f[var]
        if not var in model.p2f:
            model.p2f[var] = (var**2)*(1./sqrt(abs(E-psi)))*model.drhodr(var) 
        result2 = model.p2f[var]
        if not var in model.p3f:
            model.p3f[var] = -(var**2)*(1./sqrt(abs(E-psi)))*(1./(2*rapoval))*model.drhodr(var)*((model.Mnorm*(r-1) + r*Mencvar - Mencrap)/abs(E-psi)) 
        result3 = model.p3f[var]
    if model.memo == False:
        result1 = (var**3)*(1./sqrt(abs(E-psi)))*model.d2rhodr2(var) 
        result2 = (var**2)*(1./sqrt(abs(E-psi)))*model.drhodr(var) 
        result3 = -(var**2)*(1./sqrt(abs(E-psi)))*(1./(2*rapoval))*model.drhodr(var)*((model.Mnorm*(r-1) + r*Mencvar - Mencrap)/abs(E-psi))  
    return result1+result2+result3

def funcf(E,prereqs):
    """
    functional form of f
    relies on finterior1,finterior2,finterior3
    returns f(E)
    """
    model,Mencgood,psigood = prereqs
    if isinstance(E,(list,ndarray))==True:
        dls = zeros(len(E))
        uls = zeros(len(E))+1.
        plist = []
        pre = []
        for index in range(len(E)):
            rapoval = rapo(E[index],[psigood])
            pre.append(1./(sqrt(8)*pi**2*(model.Mnorm + 10**float(Mencgood(log10(rapoval))))))
            plist.append(tuple((E[index],rapoval,prereqs)))
    elif isinstance(E,(int,float))==True:
        dls = 0
        uls = 1.
        rapoval = rapo(E,[psigood])
        plist = tuple((E,rapoval,prereqs))
        pre = (1./(sqrt(8)*pi**2*(model.Mnorm + 10**Mencgood(log10(rapoval)))))
    return integrator(E,[finterior,'f'],dls,uls,args = plist,fileobj = model.statfile,prefactor = pre) 

########******************* ADDITIONAL FUNCTIONS *******************######## 

def funcq(r,prereqs):
    model,Ggood = prereqs
    return (4./pi)*log(model.Lam)*(model.r0_rT/model.MBH)*10**Ggood(log10(r))

def Rlc(r,prereqs):
    model,Jc2good = prereqs
    interior = 2*(model.Mnorm/model.r0_rT)*(1./10**Jc2good(log10(r)))
    return -log(interior)

########******************* IMPORT DATA TABLES *******************########

bessel = BesselGen(['alpham_table.pkl','xi_table.txt','Bessel_table.txt','mpiece_table.txt'])

########******************* RATE *******************########

def dgdlnrpinterior(E,prereqs,u,qmin):
    """
    interior of the integral used to calculate rate as function of pericentre
    """
    model,Jc2good,Ggood,fgood = prereqs
    sumlim = max([200,2*qmin**-0.5])
    ms = arange(1,sumlim,1.)
    alphas = bessel.alpham(ms)
    qval = funcq(E,[model,Ggood])
    fval = 10**fgood(log10(E))
    xival = bessel.xi(qval)
    bfin = bessel.besselfin(ms,u)
    mpiece = bessel.mpiece(ms)
    part1 = array(fval/(1+(qval**-1)*(xival)*Rlc(E,[model,Jc2good])))
    part2list = exp(-array(matrix(alphas**2).T*matrix(qval/4)))
    part2list = array([(bfin/mpiece)[i]*part2list[i] for i in range(len(alphas))])
    part2 = 1-2*nsum(part2list,axis = 0)
    return part1*part2

def funcdgdlnrp(u,prereqs,Emin = 0.01,Emax=100,verbose = False):
    """
    rp - pericentre radius
    Emin, Emax - bounds of the integral
    verbose = True - print error messages from integration
    returns the rate for given rp
    """
    model,Jc2good,Ggood,fgood = prereqs
    #u = sqrt(rp/model.rT)
    prefactor = (8*pi**2)*model.MBH*(model.r0_rT**-1)*((model.tdyn0/(3600*24*365))**-1)
    #print 'prefactor = ',prefactor
    qmin = funcq(Emax,[model,Ggood])
    qmax = funcq(Emin,[model,Ggood])
    if isinstance(u,(list,ndarray)) == True:
        dls = array([Emin]*len(u))
        uls = array([Emax]*len(u))
        plist = []
        pre = []
        for i in range(len(u)):
            plist.append(tuple((prereqs,u[i],qmin)))
            pre.append((8*pi**2)*model.MBH*(model.r0_rT**-1)*((model.tdyn0/(3600*24*365))**-1)*u[i]**2)
    elif isinstance(u,(int,float)) == True:
        dls = Emin
        uls = Emax
        plist = tuple((prereqs,u,qmin))
        pre = (8*pi**2)*model.MBH*(model.r0_rT**-1)*((model.tdyn0/(3600*24*365))**-1)*u**2
    return integrator(u,[dgdlnrpinterior,'rprate'],dls,uls,tol = 1e-3,args = plist,fileobj = model.statfile,prefactor = pre)#,div = 20)

def Ndotinterior(r):
    u = sqrt(r/model.rT)
    return r*dgdlnrp(u)

def Ndot():
    return intg.quad(Ndotinterior,0,model.rT,full_output=1)[0]


def ginterior(E):
    qval = funcq(E)
    return (10**fgood(log10(E))*qval)/((qval/bessel.xi(qval)) + Rlc(E))

def gdirect(Emin = 0.01,Emax = 100,verbose = False):
    prefactor = (8*pi**2)*model.MBH*(model.r0_rT**-1)*(model.tdyn0**-1)
    qmin = funcq(Emax)
    qmax = funcq(Emin)
    result = intg.quad(ginterior,Emin,Emax,full_output = 1)
    t = result[0]
    try:
        if result[3] != '':
            if verbose == True:
                print 'direct rate message = ',result[3]
    except (IndexError,TypeError):
        pass
    return prefactor*t
