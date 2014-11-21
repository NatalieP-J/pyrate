from loadWM import getWM
from rateget import *
from numpy import *
from rhoratefcns import *
from rhomodels import NukerModelRho
import matplotlib.pyplot as plt

Gconst = 6.67259e-8

def period(E,prereq):
    return sqrt((4*(pi**2)*rapo(E,prereq)**3)/Gconst)

def Rlcint(r,prereq):
    model,Jc2good = prereq
    interior = 2*(model.Mnorm/model.r0_rT)*(1./10**Jc2good(log10(r)))
    return interior

plotarrays = [arange(0.9,2.1,0.005),arange(0.9,2.1,0.005),arange(0.9,4,0.005),arange(0.9,4,0.005),arange(0.9,4,0.005)]

i1 = 17
#i2 = ?

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
    for j in range(len(masses)):
        MBH_Msun = 10**masses[j]
        model = NukerModelRho('{0}'.format(name),alpha,beta,gamma,rb,rho0,MBH_Msun,GENERATE,memo = False)
        M,psi,Jc2,g,G,f,rate = getrate(model)
        rapoval = rapo(plotarrays[j],[psi])
        q = funcq(plotarrays[j],[model,G])
        R = Rlcint(plotarrays[j],[model,Jc2])
        P = period(plotarrays[j],[psi])
        rapos.append(rapoval)
        fs.append(f(plotarrays[j]))
        qs.append(q)
        Ps.append(P)
        Rlcs.append(R)
    plt.figure()
    plt.suptitle(name)
    plt.subplot(326)
    #for i in range(len(rapos)):
    #    plt.loglog(Fs[j][0],Fs[j][1]*yr,label = 'MBH {0}'.format(masses[i]))
    plt.xlabel('E')
    plt.ylabel(r'F')
    plt.subplot(322)
    for k in range(len(rapos)):
        plt.loglog(10**plotarrays[k],Ps[k],label = 'log10(mass) = {0}'.format(masses[k]))
    plt.xlabel('E')
    plt.ylabel(r'P')
    plt.subplot(321)
    for k in range(len(rapos)):
        plt.loglog(10**plotarrays[k],rapos[k]/(dist*1000),label = 'log10(mass) = {0}'.format(masses[k]))
    plt.xlabel('E')
    plt.ylabel(r'$r_{apo}$["]')
    plt.subplot(323)
    for k in range(len(fs)):
        plt.loglog(10**plotarrays[k],10**fs[k],label = 'log10(mass) = {0}'.format(masses[k]))
    plt.xlabel('E')
    plt.ylabel('f')
    plt.subplot(325)
    for k in range(len(qs)):
        plt.loglog(10**plotarrays[k],qs[k],label = 'log10(mass) = {0}'.format(masses[k]))
    plt.xlabel('E')
    plt.ylabel('q')
    plt.subplot(324)
    for i in range(len(Rlcs)):
        plt.loglog(10**plotarrays[k],Rlcs[k],label = 'log10(mass) = {0}'.format(masses[k]))
    plt.xlabel('E')
    plt.ylabel(r'$R_{lc}$')
    plt.show()
    plt.savefig('WM{0}allplot.png'.format(name))
