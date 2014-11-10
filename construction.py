from numpy import *
import matplotlib.pyplot as plt
import pickle
from scipy.interpolate import interp1d
import scipy.integrate as intg
import time
import datetime
import os
from matplotlib.backends.backend_pdf import PdfPages
from suppressor import RedirectStdStreams

# dictionaries used in the fromfileplot and fromfileplotall functions
# used to search for files based on a variety of keywords 

#accepted keys for enclosed mass
funcnames = dict.fromkeys(['Menc','Mencgood','menc','funcM','mencgood','M','m',
                           'Mgood','mgood','mass','Mass'],'Menc')

#accepted keys for gravitational potential
funcnames.update(dict.fromkeys(['psi','Psi','psigood','Psigood','potential',
                                'Potential','funcp','P','p','U','u',
                                'potential energy','Potential Energy',
                                'potential Energy'],'psi'))

#accepted keys for circular angular momentum squared
funcnames.update(dict.fromkeys(['Jc2','Jc2good','Jc','jc2','jc2good','jc',
                                'Jcgood','jcgood','J','j','jgood','Jgood',
                                'Angular momentum','Angular Momentum',
                                'angular momentum'],'Jc2'))

#accepted keys for little g
funcnames.update(dict.fromkeys(['g','ggood','lg','little g','Little g'],'lg'))

#accepted keys for mathcalG
funcnames.update(dict.fromkeys(['G','Ggood','mathcalG','mathcalGgood'],'bG'))

#accepted keys for the distribution function
funcnames.update(dict.fromkeys(['f','DF','df','fgood','distribution',
                                'distribution function','Distribution function',
                                'Distribution','Distribution Function',
                                'distribution Function','F'],'f'))

#accepted keys for density
funcnames.update(dict.fromkeys(['rho','density'],'rho'))

#accepted keys for density's first derivative
funcnames.update(dict.fromkeys(['drhodr'],'drhodr'))

#accepted keys for density's second derivative
funcnames.update(dict.fromkeys(['d2rhodr2'],'d2rhodr2'))

#accepted keys for the rate as a function of pericentre
funcnames.update(dict.fromkeys(['dgdlnrp','rate'],'dgdlnrp'))

#dictionary linking dependent variable with appropriate dependent variable
indeps = {'Menc':'r','psi':'r','Jc2':'E','lg':'E','bG':'E','f':'E','rho':'r',
          'drhodr':'r','d2rhodr2':'r','dgdlnrp':'u^2'}
                
#create null file object
devnull = open(os.devnull,'w')

def loaddata(fname):
    """
    Opens a file with name fname with any sort of content and splits on 
    the newline.
    Returns an array with file contents, with length the number of lines 
    in the file.
    """
    f=open(fname,'r')
    data=[]
    for line in f.readlines():
        data.append(line.replace('\n',''))
    f.close()
    return data

def pklread(fname):
    """
    Opens a pickled file with name fname.
    Returns data array stored in pickled file.
    """
    pklffile = open(fname,"rb")
    dat = pickle.load(pklffile)
    pklffile.close()
    return dat

def pklwrite(fname,dat):
    """
    Pickles dat and writes it to file with name fname.
    """
    pklffile = open(fname,"wb")
    pickle.dump(dat,pklffile)
    pklffile.close()

def goodcheck(vals):
    """
    Checks if vals is composed entirely of nans or negative numbers. 
    If it does, return False. If not, return True.
    """
    negs = vals[where(vals<0)]
    nans = vals[where(isnan(vals)==True)]
    if len(negs) == len(vals):
        return False
    elif len(nans) == len(vals):
        return False
    elif (len(nans) + len(negs)) == len(vals):
        return False
    else:
        return True

def piecewise(r,inter,start,end,lim1,lim2,smallrexp,largerexp):
    """
    r - independent variable
    inter - interpolated object
    start - first computed value of function
    end - last computed value of function
    lim1 - first piecewise break
    lim2 - second piecewise break
    smallrexp - log slope at small r or large E
    largerexp - log slope at large r or small E

    Returns value(s) for the function at r.
    """
    set1 = r[(r<lim1)]
    set2 = r[(r>=lim1)&(r<=lim2)]
    set3 = r[(r>lim2)]
    piece1 = start*(set1/lim1)**smallrexp
    piece2 = 10**(inter(log10(set2)))
    piece3 = end*(set3/lim2)**largerexp
    return concatenate((piece1,piece2,piece3))

def makegood(prereqs,func,r,size,grid,smallrexp,largerexp,plotting):
    """
    prereqs - array containing model class instance as first element
    func - function to be evaluated
    r - independent variable array
    size - size of generated independent variable array with format 
    	   [log10(max),log10(min),stepsize]
    grid - choice of grid generator function
    smallrexp - log slope at small r or large E
    largerexp - log slope at large r or small E
    plotting - if False, do not plot. 
               if not False, must be array with ['<xlabel>','<ylabel>']
    
    Returns an interpolated object version of the function based 
    computed values.
    """
    model = prereqs[0]
    #generate independent array grid
    rarray,rchange,rstart = grid([model],size[0],size[1],size[2])
    #compute value of function for grid points
    tab,problems = func(rarray,prereqs)
    frac = float(len(problems))/float(len(tab))
    #report the fraction of problem points to console and file
    print 'fraction reporting a message: {0}'.format(frac)
    model.statfile.write('\nmesg frac = {0}\n'.format(frac))
    #check for problem points not caught in integration process
    gcheck = goodcheck(tab)
    #interpolate in log10 space
    inter = interp1d(log10(rarray),log10(tab))
    #generate array to further extremes using powerlaw behaviour
    m = piecewise(r,inter,tab[0],tab[len(rarray)-1],rstart,rchange,smallrexp,largerexp)
    #interpolate extended array in log10 space
    inter2 = interp1d(log10(r),log10(m))
    #save values used to interpolate to file (NOT in log10 space)
    saver = column_stack((r,m))
    funcname = str(func).split(' ')[1][4:]
    pklwrite('{0}/{1}.pkl'.format(model.directory,funcname),saver)
    #if plotting is possible and the array doesn't consist entirely of problems
    #add plot to pdf and return interpolate functional form
    if plotting != False and gcheck == True:
        xaxis,yaxis = plotting
        plt.figure()
        plt.loglog(r[1:-1],m[1:-1],'c',linewidth = 5)
        plt.loglog(rarray,tab,'.',color = 'DarkOrange')
        plt.ylabel(r'{0}'.format(yaxis))
        plt.xlabel('{0}'.format(xaxis))
        plt.xlim(min(r[1:-1]),max(r[1:-1]))
        plt.ylim(min(m[1:-1]),max(m[1:-1]))
        plt.title(model.name)
        model.pdfdump.savefig()
        plt.close()
        return inter2
    #if plotting isn't possible but array doesn't consist entirely of problems
    #return interpolated functional form
    elif plotting == False and gcheck == True:
        return inter2
    #if computation failed, return 0
    #this signals the rest of the program that computation failed here
    elif gcheck == False:
        return 0

def compute(dependencies,function,rtest,size,grid,exps,plotdat,create):
    """
    dependencies - other functions needed to compute this one, 
                   format [model,"model",func1, "func1",func2,"func2",...]
    function - name of the functional form
    rtest - independent variable array
    size - size of generated independent variable array 
    	   with format [log10(max),log10(min),stepsize]
    grid - choice of grid generator function (rgrid or Egrid)
    exps - extreme r or E power law behaviour,format [smallrexp,largerexp]
    plotdat - if False, do not plot. 
              if not False, must be array with ['<xlabel>','<ylabel>']
    create - if 'ON' generate interpolated functional form with makegood. 
             if 'OFF' load functional form from file.
             if 'FAIL' return zero
    
    Finds interpolated form based on conditions set by create and pickles 
    it or unpickles intepolated form.
    Returns interpolated form.
    """
    model = dependencies[0]
    prereqs = dependencies[0::2]
    strname = str(function).split(' ')[1][4:]
    #if function is set to made from scratch
    if create == "ON":
        try: 
            #check if all dependencies are available
            if len(dependencies) > 2:
                i = 2
                while i<len(dependencies):
                    dependencies[i](1)
                    i+=2
            #if dependencies are available, make functional form with makegood
            #time the process and write it to console and file
            smallrexp,largerexp = exps
            tic = time.clock()
            good = makegood(prereqs,function,rtest,size,grid,smallrexp,largerexp,plotting = plotdat)
            toc = time.clock()
            delt = toc-tic
            print '{0}good ran in \t {1}'.format(strname,str(datetime.timedelta(seconds=delt)))
            model.statfile.write('{0}good ran in \t {1}\n'.format(strname,str(datetime.timedelta(seconds=delt))))
            return good
        except TypeError:
            #if dependencies are unavailable, write error message to console and file
            print 'To compute {0}, please turn {1} ON'.format(strname,dependencies[i+1])
            model.statfile.write('To compute {0}, please turn {1} ON\n'.format(strname,dependencies[i+1]))
    #if function is already saved, load it from file
    elif create == "OFF":
        #if file exists, read it in
        try:
            dat = pklread('{0}/{1}.pkl'.format(model.directory,strname))
            rarray = dat[:,0]
            tab = dat[:,1]
            gcheck = goodcheck(tab)
            #if plotting is available, add plot to pdf
            if plotdat != False and gcheck == True:
                xaxis,yaxis = plotdat
                plt.loglog(rarray[1:-1],tab[1:-1],'c',linewidth = 5)
                plt.ylabel(r'{0}'.format(yaxis))
                plt.xlabel('{0}'.format(xaxis))
                plt.xlim(min(rarray[1:-1]),max(rarray[1:-1]))
                plt.ylim(min(tab[1:-1]),max(tab[1:-1]))
                plt.title(model.name)
                model.pdfdump.savefig()
                plt.close()
            good =  interp1d(log10(rarray),log10(tab))
            print '{0}good loaded in'.format(strname)
            model.statfile.write('{0}good loaded in'.format(strname))
            return good
        #if file does not yet exist, write error to console and file
        #additionally, signal evaluation has failed here
        except IOError:
            print '{0} not yet generated, please turn it ON\n'.format(strname)
            model.statfile.write('{0} not yet generated, please turn it ON\n'.format(strname))
            return 0
    #if function creation failed, signal failed evaluation here
    elif create == 'FAIL':
        return 0
        


def integrator(vals,fcn,downlim,uplim,args,tol=1.49e-7,fileobj=devnull,prefactor = 1):
    """
    vals - list of values for which to compute the integral
    fcn - interior function for the integrator
    downlim - lower limit (or array thereof) for integrator
    uplim - upper limit (or array thereof) for integrator
    args - additional arguments fed to integrator as a tuple or array of tuples
    tol - tolerance, defaults to 1.49e-7
    fileobj - location to write full_output, defaults to devnull
    prefactor - constant by which to multiply the integral as a float or
                an array, defaults to 1

    Returns array of computed integrals and an array of problem points.
    """
    #if vals is an array of values, cycle through them, integrating for each
    if isinstance(vals,(list,ndarray))==True:
        problems = []
        results = []
        for i in range(len(vals)):
            #set up different integrals depending on whether there are 
            #additional arguments
            if args != []:
                temp = intg.quad(fcn[0],downlim[i],uplim[i],args = args[i],epsabs = tol,full_output = 1)
            elif args == []:
                temp = intg.quad(fcn[0],downlim[i],uplim[i],epsabs = tol,full_output = 1)
            #check for messages from the integrator and if they exist
            #write them to file
            try:
                if temp[3] != '':
                    problems.append(i)
                    fileobj.write('{0},\t i = {1},\t message = {2}\n'.format(fcn[1],i,temp[3]))
            except IndexError:
                pass
            #check whether to multiply by prefactor or not
            if prefactor != 1:
                results.append(prefactor[i]*temp[0])
            elif prefactor == 1:
                results.append(temp[0])
        return array(results),problems
    elif isinstance(vals,(int,float))==True:
        problems = []
        #set up different integrals depending on whether there are 
        #additional arguments
        if args != []:
            temp = intg.quad(fcn[0],downlim,uplim,args = args,epsabs = tol,full_output = 1)
        elif args == []:
            temp = intg.quad(fcn[0],downlim,uplim,epsabs = tol,full_output = 1)
        #check for messages from the integrator and if they exist
        #write them to file
        try:
            if temp[3] != '':
                problems.append(vals)
                fileobj.write('\n{0},\t val = {1},\t message = {2}'.format(fcn[1],vals,temp[3]))
        except IndexError:
            pass
        return prefactor*temp[0],problems


def rintegrator(vals,fcn,downlim,uplim,tol=1.49e-7,args = [],fileobj=devnull,prefactor = 1,div=50):
    """
    vals - list of values for which to compute the integral
    fcn - interior function for the integrator
    downlim - lower limit (or array thereof) for integrator
    uplim - upper limit (or array thereof) for integrator
    tol - tolerance
    args - additional arguments fed to integrator as a tuple or array of tuples
    fileobj - location to write full_output
    prefactor - constant by which to multiply the integral as a float or an array
    div - maximum number of divisions

    return array of romberg-computed integrals and an array of problem points
    """
    #if vals is an array of values, cycle through them, integrating for each
    if isinstance(vals,(list,ndarray))==True:
        problems = []
        results = []
        for i in range(len(vals)):
            #write all output to chosen file
            with RedirectStdStreams(stdout = fileobj,stderr = fileobj):
                #set up different integrals depending on whether there are 
                #additional arguments
                if args != []:
                    temp = intg.romberg(fcn[0],downlim[i],uplim[i],args = args[i],tol = tol,divmax = div)
                elif args == []:
                    temp = intg.romberg(fcn[0],downlim[i],uplim[i],tol = tol,divmax = div)
            #check whether to multiply by prefactor or not
            if prefactor != 1:
                results.append(prefactor[i]*temp[0])
            elif prefactor == 1:
                results.append(temp[0])
        return array(results),problems
    elif isinstance(vals,(int,float))==True:
        problems = []
        #write all output to chosen file
        with RedirectStdStreams(stdout = fileobj,stderr = fileobj):
            #set up different integrals depending on whether there are 
            #additional arguments
            if args != []:
                temp = intg.romberg(fcn[0],downlim,uplim,args = args,tol = tol,divmax = div)
            elif args == []:
                temp = intg.romberg(fcn[0],downlim,uplim,tol = tol,divmax = div)
        return prefactor*temp[0],problems


def dblintegrator(vals,fcn,downlim,uplim,tol=1.49e-7,args = [],fileobj=devnull,prefactor = 1):
    """
    vals - list of values for which to compute the integral
    fcn - interior function for the integrator
    downlim - length two array containing two sets of lower limits for the double integrator
    uplim - length two array containing two sets of upper limits for the double integrator
    tol - tolerance
    args - additional arguments fed to integrator as a tuple or array of tuples
    fileobj - location to write full_output
    prefactor - constant by which to multiply the integral as a float or an array

    return array of computed integrals and an array of problem points
    """
    #if vals is an array of values, cycle through them, integrating for each
    if isinstance(vals,(list,ndarray))==True:
        problems = []
        results = []
        for i in range(len(vals)):
            #write all output to chosen file
            with RedirectStdStreams(stdout = fileobj,stderr = fileobj):
                #set up different integrals depending on whether there are 
                #additional arguments
                if args != []:
                    temp = intg.dblquad(fcn[0],downlim[0][i],uplim[0][i],downlim[1][i],uplim[1][i],args = args[i],epsabs = tol)
                elif args == []:
                    temp = intg.dblquad(fcn[0],downlim[0][i],uplim[0][i],downlim[1][i],uplim[1][i],epsabs = to1)
            if prefactor != 1:
                results.append(prefactor[i]*temp)
            elif prefactor == 1:
                results.append(temp)
        return array(results),problems
    elif isinstance(vals,(int,float))==True:
        problems = []
        #write all output to chosen file
        with RedirectStdStreams(stdout = fileobj,stderr = fileobj):
            #set up different integrals depending on whether there are 
            #additional arguments
            if args != []:
                temp = intg.dblquad(fcn[0],downlim[0],uplim[0],downlim[1],uplim[1],args = args,epsabs = tol)
            elif args == []:
                temp = intg.dblquad(fcn[0],downlim[0],uplim[0],downlim[1],uplim[1],epsabs = tol)
        return prefactor*temp,problems

def fromfileplot(galname,funcname,up,down):
    """
    
    Searches for directories containing galname one step down 
    from current directory.
    Plots funcname for chosen directory between limits up and down.

    """
    #set up independent variable array
    r = arange(down,up,0.01)
    r = 10**r
    #search for appropriate directories and choose one
    success = os.system('ls -d */*{0}* > templist.dat'.format(galname))
    if success == 0:
        avails = loaddata('templist.dat')
        os.system('rm -f templist.dat')
        if len(avails) > 1:
            print 'Multiple directories with that name, please choose from the following list'
            i = 0
            a = 'n'
            while a == 'n' or a == '':
                i+=1
                i = i%len(avails)
                a = raw_input('Is this your galaxy [y/n]?\n{0} '.format(avails[i]))
            direc = avails[i]
        elif len(avails) == 1:
            direc = avails[0]
        #find actual function name
        try:
            fname = funcnames[funcname]
            path = '{0}/{1}.pkl'.format(direc,fname)
            dat = pklread(path)
            rarray = dat[:,0]
            tab = dat[:,1]
            good = interp1d(log10(rarray),log10(tab))
            plt.figure()
            plt.loglog(r,10**good(log10(r)))
            plt.ylabel(funcname)
            plt.xlabel(indeps[fname])
            plt.title(galname)
            plt.show()
        except KeyError:
            print 'I can\'t find that function.\nTry one of these: "Menc","psi","Jc2","lg","bG","f","dgdlnrp"'
    elif success != 0:
        print 'There do not appear to be any directories with that galaxy name, terminating plot'
        

def fromfileplotall(galname):
    """
    
    Searches for directories containing galname one step down 
    from current directory.
    Plots all six functions for chosen directory.

    """
    #make standard independent variable array
    r = arange(-4,4,0.01)
    r = 10**r
    #search for directories and choose one
    success = os.system('ls -d */*{0}* > templist.dat'.format(galname))
    if success == 0:
        avails = loaddata('templist.dat')
        os.system('rm -f templist.dat')
        if len(avails) > 1:
            print 'Multiple directories with that name, please choose from the following list'
            i = 0
            a = 'n'
            while a == 'n' or a == '':
                i+=1
                i = i%len(avails)
                a = raw_input('Is this your galaxy [y/n]?\n{0} '.format(avails[i]))
            direc = avails[i]
        elif len(avails) == 1:
            direc = avails[0]
        #load in all functions and plot each
        dat = pklread('{0}/{1}.pkl'.format(direc,'Menc'))
        rarray = dat[:,0]
        tab = dat[:,1]
        Mgood =  interp1d(log10(rarray),log10(tab))
        plt.figure()
        plt.loglog(r,10**Mgood(log10(r)))
        plt.xlabel(r'$r$')
        plt.ylabel(r'$M_{enc}$')
        plt.title(name)
        dat = pklread('{0}/{1}.pkl'.format(direc,'psi'))
        rarray = dat[:,0]
        tab = dat[:,1]
        Mgood =  interp1d(log10(rarray),log10(tab))
        plt.figure()
        plt.loglog(r,10**Mgood(log10(r)))
        plt.xlabel(r'$r$')
        plt.ylabel(r'$\psi$')
        plt.title(name)
        dat = pklread('{0}/{1}.pkl'.format(direc,'Jc2'))
        rarray = dat[:,0]
        tab = dat[:,1]
        Mgood =  interp1d(log10(rarray),log10(tab))
        plt.figure()
        plt.loglog(r,10**Mgood(log10(r)))
        plt.xlabel(r'$E$')
        plt.ylabel(r'$J_c^2$')
        plt.title(name)
        dat = pklread('{0}/{1}.pkl'.format(direc,'lg'))
        rarray = dat[:,0]
        tab = dat[:,1]
        Mgood =  interp1d(log10(rarray),log10(tab))
        plt.figure()
        plt.loglog(r,10**Mgood(log10(r)))
        plt.xlabel(r'$E$')
        plt.ylabel(r'$g$')
        plt.title(name)
        dat = pklread('{0}/{1}.pkl'.format(direc,'bG'))
        rarray = dat[:,0]
        tab = dat[:,1]
        Mgood =  interp1d(log10(rarray),log10(tab))
        plt.figure()
        plt.loglog(r,10**Mgood(log10(r)))
        plt.xlabel(r'$E$')
        plt.ylabel(r'$G$')
        plt.title(name)
        dat = pklread('{0}/{1}.pkl'.format(direc,'f'))
        rarray = dat[:,0]
        tab = dat[:,1]
        Mgood =  interp1d(log10(rarray),log10(tab))
        plt.figure()
        plt.loglog(r,10**Mgood(log10(r)))
        plt.xlabel(r'$E$')
        plt.ylabel(r'$f$')
        plt.title(name)
        plt.show()
        dat = pklread('{0}/{1}.pkl'.format(direc,'dgdlnrp'))
        rarray = dat[:,0]
        tab = dat[:,1]
        Mgood =  interp1d(log10(rarray),log10(tab))
        u = arange(-1,0,0.01)
        u = 10**u
        plt.figure()
        plt.loglog(u**2,10**Mgood(log10(u)))
        plt.xlabel(r'$r_p/r_T$')
        plt.ylabel(r'$d\gamma/d\ln{r_p}$')
        plt.title(name)
        plt.show()
    elif success != 0:
        print 'There do not appear to be any directories with that galaxy name, terminating plot'
