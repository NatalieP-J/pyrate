from numpy import *

def LoadDataTab(fname):
    f=open(fname,'r')
    data=[]
    for line in f.readlines():
        data.append(line.replace('\n','').split('\t'))
    f.close()
    return data

def getWM():
	WM = array(LoadDataTab('WM04.dat'))[:,][:-10]
	names = WM[:,0]
	dists = array([float(i) for i in WM[:,2]])
	rbs = 10**array([float(i) for i in WM[:,3]])
	mubs = array([float(i) for i in WM[:,4]])
	alphas = array([float(i) for i in WM[:,5]])
	betas = array([float(i) for i in WM[:,6]]) + 1
	gammas = array([float(i) for i in WM[:,7]]) + 1
	M2Ls = array([float(i) for i in WM[:,8]])
	MBH1s = 10**array([float(i) for i in WM[:,10]])
	MBH2s = 10**array([float(i) for i in WM[:,12]])
	return WM,names,dists,rbs,mubs,alphas,betas,gammas,M2Ls,MBH1s,MBH2s

if __name__ == '__main__':
	from rhomodels import NukerModelRho
	from rateget import getrate
	from rhoratefcns import findrho0
	WM,names,dists,rbs,mubs,alphas,betas,gammas,M2Ls,MBH1s,MBH2s = getWM()
	rho0s = findrho0(rbs,M2Ls,mubs)
	ilist = array([ 0,  1,  2,  3,  4,  5,  6,  7,  8,  9, 10, 11, 12, 13, 14, 15, 16,
       17, 18, 19, 20, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31, 32, 33, 34,
       35, 36, 37, 38, 39, 40]) #skip 21
	for i in ilist:
		print 'Working on ',names[i],' galaxy ',i+1,' of ',len(ilist)
		alpha = alphas[i]
		beta = betas[i]
		gamma = gammas[i]
		r0pc = rbs[i]
		rho0 = rho0s[i]
		MBH_Msun1 = MBH1s[i]
		MBH_Msun2 = MBH2s[i]
		name1 = '{0}_1'.format(names[i])
		name2 = '{0}_2'.format(names[i])
		GENERATE = False
		model1 = NukerModelRho(name1,alpha,beta,gamma,r0pc,rho0,MBH_Msun1,GENERATE,memo = False)
		result1 = getrate(model1)
		model2 = NukerModelRho(name2,alpha,beta,gamma,r0pc,rho0,MBH_Msun2,GENERATE,memo = False)
		result2 = getrate(model2)
