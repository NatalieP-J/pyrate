from numpy import *
from scipy.optimize import leastsq
import matplotlib.pyplot as plt
plt.ion()

def nuker(x,p):
    a,b,g,R = p
    xr = x/R
    return xr**(-g)*(1 + xr**a)**(-(b-g)/a)

def lognuker(x,p):
    a,b,g,R = p
    xr = x/R
    return log10(xr**(-g)*(1 + xr**a)**(-(b-g)/a))
    
def sersic(x,n):
    p = 1. - 0.6097/n + 0.05563/n**2
    b = 2*n - (1./3.) + 0.009876/n
    return x**(-p)*exp(-b*x**(1./n))

def sersic2(x,n,Re):
	xnew = x/Re
	p = 1. - 0.6097/n + 0.05563/n**2
	b = 2*n - (1./3.) + 0.009876/n
	return xnew**(-p)*exp(-b*xnew**(1./n))

def logsersic(x,n):
    p = 1. - 0.6097/n + 0.05563/n**2
    b = 2*n - (1./3.) + 0.009876/n
    return log10(x**(-p)*exp(-b*x**(1./n)))

def drhodr(r,n):
	p = 1. - 0.6097/n + 0.05563/n**2
	b = 2*n - (1./3.) + 0.009876/n
	pre = (r**-p)*exp(-b*(r**(1./n)))*((n*r)**-1)
	post = (n*p)+b*(r**(1./n))
	return pre*post

def powerlaw(x,p):
	m,b = p
	return 10**(m*log10(x) + b)

def logpowerlaw(x,p):
	m,b = p
	return m*log10(x) + b

def residuals(p,y,x,func):
    return y - func(x,p)

r = arange(-7,7,0.01)
r = 10**r

n = 4.

y1 = sersic(r,n)
y2 = logsersic(r,n)

p03 = [0.85,6]
p04 = [-20.,50]

params3 = leastsq(residuals, p03, args = (y1[:200],r[:200],powerlaw), full_output=True)
params4 = leastsq(residuals, p04, args = (y2[-50:],r[-50:],logpowerlaw), full_output=True)


plt.figure()
plt.loglog(r,y1)
plt.loglog(r,powerlaw(r,params3[0]))
plt.loglog(r[-300:],10**logpowerlaw(r[-300:],params4[0]))

slope = (log10(y1[1:]) - log10(y1[0:-1]))/(log10(r[1:]) - log10(r[0:-1]))

