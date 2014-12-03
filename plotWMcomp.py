from loadWM import getWM
from rateget import *
from numpy import *
from rhoratefcns import *
from rhomodels import NukerModelRho
import matplotlib.pyplot as plt
plt.ion()

Gconst = 6.67259e-8
realMsun = 1.989e33
pc = 3.1e18
yr = 365*24*3600

def period(E,prereq):
    return sqrt((4*(pi**2)*rapo(E,prereq)**3)/Gconst)

def Rlcint(r,prereq):
    model,Jc2good = prereq
    interior = 2*(model.Mnorm/model.r0_rT)*(1./10**Jc2good(log10(r)))
    return interior

def logR0(E,prereq1,prereq2):
    q = funcq(E,prereq1)
    if isinstance(E,ndarray) == True:
        part1 = -(log(Rlc(E[q > 1],prereq2))-q[q > 1])
        part2 = -(log(Rlc(E[q < 1],prereq2))-0.186*q[q < 1] - 0.824*sqrt(q[q < 1]))
        return concatenate((part1,part2))
    elif isinstance(E,ndarray) == False:
        if q > 1:
            return Rlc(E,prereq2)*exp(-q)
        elif q < 1:
            return Rlc(E,prereq2)*exp(-0.186*q - 0.824*sqrt(q))

def F(E,prereq1,prereq2,prereq3):
    Jc2good,fgood = prereq3
    Jconst = Gconst*(model.rho0*realMsun*((model.r0*pc)**-1)*((10*pc)**-2)*(model.r0*pc)**4)
    qtest = funcq(E,prereq1)
    vals = np.where(isnan(qtest)==False)
    E = E[vals]
    return [E,(100000**2)*4*(pi**2)*(Jconst*10**Jc2good(log10(E)))*funcq(E,prereq1)*Rlcint(E,prereq2)*(fconst*10**fgood(log10(E)))/logR0(E,prereq1,prereq2)]
    

plotarrays = [arange(0.9,2.1,0.005),arange(0.9,2.1,0.005),arange(0.9,4,0.005),arange(0.9,4,0.005),arange(0.9,4,0.005)]

i1 = 17
i2 = 30

masses = arange(4,14,2)

M,names,dists,rbs,mubs,alphas,betas,gammas,M2Ls,MBH1s,MBH2s = getWM()

GENERATE = False

for i in [i1]:#,i2]:
    rapos = []
    fs = []
    qs = []
    Ps = []
    Rlcs = []
    Fs = []
    name = names[i]
    alpha = alphas[i]
    beta = betas[i]
    gamma = gammas[i]
    M2L = M2Ls[i]
    rb = rbs[i]
    mub = mubs[i]
    rho0 = findrho0(rb,M2L,mub)
    dist = dists[i]
    Earrays =[]
    for j in range(len(masses)):
        MBH_Msun = 10**masses[j]
        model = NukerModelRho('{0}'.format(name),alpha,beta,gamma,rb,rho0,MBH_Msun,GENERATE,memo = False)
        partial =[ {'Menc':'OFF','psi':'OFF','Jc2':'OFF','lg':'OFF','bG':'OFF','f':'OFF','dgdlnrp':'FAIL'},{'Menc':False,'psi':False,'Jc2':False,'lg':False,'bG':False,'f':False,'dgdlnrp':False}]
        M,psi,Jc2,g,G,f,rate = getrate(model,partial)
        r0real = model.r0*pc
        rho0real = model.rho0*(realMsun/pc**3)
        Econst = Gconst*rho0real*r0real**2/1e14
        Earrays.append(Econst*(10**plotarrays[j]))
        fconst = (sqrt((Gconst**3)*model.rho0)*realMsun*((r0real)**3))**-1
        rapoval = rapo(10**plotarrays[j],[psi])
        q = funcq(10**plotarrays[j],[model,G])
        R = Rlcint(10**plotarrays[j],[model,Jc2])
        P = period(10**plotarrays[j],[psi])
        Fval = F(10**plotarrays[j],[model,G],[model,Jc2],[Jc2,f])
        rapos.append(rapoval)
        fs.append(f(plotarrays[j]))
        qs.append(q)
        Ps.append(P)
        Rlcs.append(R)
        Fs.append(Fval)
    plt.figure()
    plt.suptitle(name)
    plt.subplot(326,adjustable = 'box',aspect = 0.1)
    for k in range(len(rapos)):
        plt.plot(log(Fs[k][0]),log(Fs[k][1]*yr),label = 'MBH {0}'.format(masses[k]))
    plt.xlabel('E')
    plt.ylabel(r'F')
    plt.subplot(322,aspect = 1)
    for k in range(len(rapos)):
        plt.loglog(10**plotarrays[k],Ps[k],label = 'log10(mass) = {0}'.format(masses[k]))
    plt.xlabel('E')
    plt.ylabel(r'P')
    plt.subplot(321,aspect = 1)
    for k in range(len(rapos)):
        plt.loglog(10**plotarrays[k],rapos[k]/(dist*1000),label = 'log10(mass) = {0}'.format(masses[k]))
    plt.xlabel('E')
    plt.ylabel(r'$r_{apo}$["]')
    plt.subplot(323,aspect = 1)
    for k in range(len(fs)):
        plt.loglog(Earrays[k],10**fs[k],label = 'log10(mass) = {0}'.format(masses[k]))
    plt.xlabel('E')
    plt.ylabel('f')
    plt.subplot(325,aspect = 1)
    for k in range(len(qs)):
        plt.loglog(10**plotarrays[k],qs[k],label = 'log10(mass) = {0}'.format(masses[k]))
    plt.xlabel('E')
    plt.ylabel('q')
    plt.subplot(324,aspect = 1)
    for k in range(len(Rlcs)):
        plt.loglog(10**plotarrays[k],Rlcs[k],label = 'log10(mass) = {0}'.format(masses[k]))
    plt.xlabel('E')
    plt.ylabel(r'$R_{lc}$')
    plt.show()
    plt.savefig('WM{0}allplot.png'.format(name))
