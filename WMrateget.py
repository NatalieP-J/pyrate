from numpy import *
from subprocess import call
from matplotlib.backends.backend_pdf import PdfPages
from rateget import *
from construction import fromfileplot
from rhomodels import NukerModelGenRho,NukerModelRho
import sys
from loadWM import getWM

i = int(sys.argv[1])
m = int(sys.argv[2])

WM,names,dists,rbs,mubs,alphas,betas,gammas,M2Ls,MBH1s,MBH2s = getWM()

GENERATE = True

name = names[i]
alpha = alphas[i]
beta = betas[i]
gamma = gammas[i]
M2L = M2Ls[i]
MBH_Msun = 10**m
rb = rbs[i]
mub = mubs[i]
rho0 = findrho0(rb,M2L,mub)
model1 = NukerModelGenRho('{0}_gen'.format(name),alpha,beta,gamma,rb,rho0,MBH_Msun,GENERATE,memo = False)
model1.getrho()
model2 = NukerModelRho('{0}'.format(name),alpha,beta,gamma,rb,rho0,MBH_Msun,GENERATE,memo = False)
result = getrate(model2)





