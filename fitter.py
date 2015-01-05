from numpy import *
from scipy.optimize import leastsq
import matplotlib.pyplot as plt

def nuker(x,p):
    a,b,g,R = p
    xr = x/R
    return xr**(-g)*(1 + xr**a)**(-(b-g)/a)
    
def sersic(x,n):
    p = 1. - 0.6097/n + 0.05563/n**2
    b = 2*n - (1./3.) + 0.009876/n
    return x**(-p)*exp(-b*x**(1./n))

def residuals(p,y,x,func):
    return y - func(x,p)

r = arange(-7,7,0.01)
r = 10**r

n = 4.

y = sersic(r,n)

p0 = [0.5,100.,1.,100.]

params = leastsq(residuals, p0, args = (y,r,nuker), full_output=True)

plt.figure()
plt.loglog(r,y)
plt.loglog(r,nuker(r,params[0]))
