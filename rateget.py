from numpy import *
from rhoratefcns import *

import os
from construction import loaddata
from construction import pklread

#standard independent variable array
rtest = arange(-7,7,0.01)
rtest = append(rtest,40)
rtest = insert(rtest,0,-40)
rtest = 10**rtest

#independent variable array for dgdlnrp
utest = arange(-7,0,0.01)
utest = insert(utest,0,-40)
utest = 10**utest

def displaycheck():
    """
    Check if a display is available.
    Return True if it is, return False if it isn't.
    """
    os.system('echo $DISPLAY > tempdisplay')
    displays = loaddata('tempdisplay')
    os.system('rm -f tempdisplay')
    display = displays[0]
    if display == '':
        return False
    elif display != '':
        return True

def existcheck(directory,dcheck):
    """
    Check which functions have already been created and generate a dictionary. 
     This dictionary tells you which functions are available for loading 
    or plotting.
    """
    #tells compute whether or not to generate the function
    seton = {}
    #tells compute whether or not to plot the function
    plottinglist = {}
    #holds True or False for each function to ensure that a failure in an early
    #function means later ones will not be computed
    gvals = {}
    strnames = ['Menc','psi','Jc2','g','G','f','dgdlnrp']
    #prerequesite functions for each of the six in order of strnames
    prereqs = [[],['Menc'],['Menc','psi'],['psi'],['psi','g'],['Menc','psi'],
               ['Jc2','G','f']]
    #array holding plot values (a list of ['<xlabel>','<ylabel>']
    pvals = [['r','M'],['r',r'$\psi$'],['E',r'$J_c^2$'],['E','g'],['E','G'],
             ['E','f'],[r'$u^2$',r'$\frac{dg}{dlnr_p}$']]
    #for each function, perform the following checks
    for i in range(len(strnames)):
        prechecks = array([])
        #check each prerequesite function exists and has not failed to evaluate
        for j in range(len(prereqs[i])):
            prechecks = append(prechecks,gvals[prereqs[i][j]])
        #prepass is zero iff all prerequesite functions exist
        prepass = len(prechecks) - len(prechecks[where(prechecks == True)])
        #try to load in the function from file
        try:
            vals = pklread('{0}/{1}.pkl'.format(directory,strnames[i]))
            gcheck = goodcheck(vals[:,1])
            gvals[strnames[i]] = gcheck
            #if this function and all of its prerequesites passed, 
            #don't reevaluate, and plot if possible
            if gcheck == True and prepass == 0:
                seton[strnames[i]] = 'OFF'
                if dcheck == True:
                    plottinglist[strnames[i]] = pvals[i]
                if dcheck == False:
                    plottinglist[strnames[i]] = False
            #if this function or one of its prerequesites failed,
            #don't attempt any sort of evaluation
            if gcheck != True or prepass != 0:
                seton[strnames[i]] = 'FAIL'
                plottinglist[strnames[i]] = False
        #if function file doesn't exist:
        except IOError:
            #add to dictionary that this function 'failed' through nonexistence
            gvals[strnames[i]] = False
            #if prerequesites all exist and pass, set this function to evaluate
            #and plot if possible
            if prepass == 0:
                seton[strnames[i]] = 'ON'
                if dcheck == True:
                    plottinglist[strnames[i]] = pvals[i]
                if dcheck == False:
                    plottinglist[strnames[i]] = False
            #if prerequesites fail, fail this function as well
            if prepass != 0:
                seton[strnames[i]] = 'FAIL'
                plottinglist[strnames[i]] = False
    return seton,plottinglist
            
def getrate(model,partial = False):
    """
    partial - option to manually set which functions to create
    Generate rate function for <model> according to seton, and plot if possible.
    """
    #check if display is on
    dcheck = displaycheck()
    #if functions to create have not been manually selected:
    if partial == False:
        
        #if not asked to generate, automatically decide what to evaluate
        if model.generate == False:
            seton,plottinglist = existcheck(model.directory,dcheck)
        
        #if asked to generate, create everything and plot if possible
        elif model.generate == True:
            seton = {'Menc':"ON",'psi':"ON",'Jc2':"ON",'g':"ON",'G':"ON",
                     'f':"ON",'dgdlnrp':"ON"}
            if dcheck == False:
                plottinglist = {'Menc':False,'psi':False,'Jc2':False,'g':False,
                                'G':False,'f':False,'dgdlnrp':False}
            if dcheck == True:
                plottinglist = {'Menc':['r','M'],'psi':['r',r'$\psi$'],
                                'Jc2':['E',r'$J_c^2$'],'g':['E','g'],
                                'G':['E','G'],'f':['E','f'],
                                'dgdlnrp':[r'$u^2$',r'$\frac{dg}{dlnr_p}$']}
    
    #if functions have been chosen, do not plot any
    elif partial != False:
        seton = partial 
        plottinglist = {'Menc':False,'psi':False,'Jc2':False,'g':False,
                        'G':False,'f':False,'dgdlnrp':False}   
    try:
        #dictionary of power law behaviour
        exps = {'Menc':[3-model.g,0],'psi':[-1,-1],'Jc2':[-1,-1],
                'g':[model.b-0.5,model.g-0.5],'G':[model.b-4,model.g-4],
                'f':[model.b-1.5,model.g-1.5],'dgdlnrp':[2,0]}

        sh = {'Menc':[4,-6,0.03],'psi':[4.3,-6,0.03],'Jc2':[3,-4,0.01],
              'g':[3,-3,0.1],'G':[3,-3,0.1],'f':[5,-3,0.03],'rate':[0,-4,0.04]}

        #begin output
        model.statfile.write('GALAXY: {0}\n'.format(model.name))

        #evaluate Menc
        model.statfile.write('Menc:\n')
        up,down,step = sh[Menc]
        rarray,rchange,rstart = rgrid([model],up,down,step)
        Mencgood = compute([model],funcMenc,rtest,sh['Menc'],rgrid,exps['Menc'],
                           plottinglist['Menc'],seton['Menc'])

        #if failed at Menc, abort computation of further functions
        if Mencgood == 0:
            model.statfile.write('Failed to evaluate Menc')
            seton = {'Menc':"FAIL",'psi':"FAIL",'Jc2':"FAIL",'g':"FAIL",
                     'G':"FAIL",'f':"FAIL",'dgdlnrp':"FAIL"}
          
        #evaluate psi
        pprereqs = [model,'Model',Mencgood,'Menc']
        model.statfile.write('\npsi:\n')
        psigood = compute(pprereqs,funcpsi,rtest,sh['psi'],rgrid,exps['psi'],
                          plottinglist['psi'],seton['psi'])
        
        #if failed at psi, abort computation of further functions
        if psigood == 0:
            model.statfile.write('Failed to evaluate psi')
            seton = {'Menc':"FAIL",'psi':"FAIL",'Jc2':"FAIL",'g':"FAIL",
                     'G':"FAIL",'f':"FAIL",'dgdlnrp':"FAIL"}
        
        #evaluate Jc2
        Jprereqs = [model,'Model',Mencgood,"Menc",psigood,"psi"]
        model.statfile.write('\nJc2:\n')
        Jc2good = compute(Jprereqs,funcJc2,rtest,sh['Jc2'],Egrid,exps['Jc2'],
                          plottinglist['Jc2'],seton['Jc2'])
        
        #if failed at Jc2, abort computation of further functions
        if Jc2good == 0:
            model.statfile.write('Failed to evaluate Jc2')
            seton = {'Menc':"FAIL",'psi':"FAIL",'Jc2':"FAIL",'g':"FAIL",
                     'G':"FAIL",'f':"FAIL",'dgdlnrp':"FAIL"}

        #evaluate g
        lgprereqs = [model,'Model',psigood,"psi"]
        model.statfile.write('\ng:\n')
        ggood = compute(lgprereqs,funclg,rtest,sh['g'],Egrid,exps['g'],
                        plottinglist['g'],seton['g'])
        
        #if failed at g, abort computation of further functions
        if ggood == 0:
            model.statfile.write('Failed to evaluate g')
            seton = {'Menc':"FAIL",'psi':"FAIL",'Jc2':"FAIL",'g':"FAIL",
                     'G':"FAIL",'f':"FAIL",'dgdlnrp':"FAIL"}
                
        #evaluate mathcalG
        bGprereqs = [model,'Model',psigood, "psi",ggood,"g"]
        
        Gtest = arange(-5,5,0.01)
        Gtest = append(Gtest,40)
        Gtest = insert(Gtest,0,-40)
        Gtest = 10**Gtest
                        
        model.statfile.write('\nG:\n')
        Ggood = compute(bGprereqs,funcbG,Gtest,sh['G'],Egrid,exps['G'],
                        plottinglist['G'],seton['G'])
    
        if model.memo == True:
            model.p1bG = {}
            model.p2bG = {}
            model.p3bG = {}
        
        #if failed at mathcalG, abort computation of further functions
        if Ggood == 0:
            model.statfile.write('Failed to evaluate G')
            seton = {'Menc':"FAIL",'psi':"FAIL",'Jc2':"FAIL",'g':"FAIL",
                     'G':"FAIL",'f':"FAIL",'dgdlnrp':"FAIL"}

        #evaluate f
        fprereqs = [model,'Model',Mencgood,"Menc",psigood,"psi"]
            
        ftest = arange(-3,5,0.01)
        ftest = append(ftest,40)
        ftest = insert(ftest,0,-40)
        ftest = 10**ftest
        
        model.statfile.write('\nf:\n')
        fgood = compute(fprereqs,funcf,ftest,sh['f'],Egrid,exps['f'],
                        plottinglist['f'],seton['f'])
        
        #if failed at f, abort computation of further functions
        if fgood == 0:
            model.statfile.write('Failed to evaluate f')
            seton = {'Menc':"FAIL",'psi':"FAIL",'Jc2':"FAIL",'g':"FAIL",
                     'G':"FAIL",'f':"FAIL",'dgdlnrp':"FAIL"}

        #evaluate dgdlnrp
        rprereqs = [model,'Model',Jc2good,'Jc2',Ggood,'Ggood',fgood,'fgood']
        model.statfile.write('\nrate:\n')
        rategood = compute(rprereqs,funcdgdlnrp,utest,sh['dgdlnrp'],stdgrid,
                           exps['dgdlnrp'],plottinglist['dgdlnrp'],
                           seton['dgdlnrp'])
        
        if rategood == 0:
            model.statfile.write('Failed to evaluate dgdlnrp')
        
        #close plot pdf and statusfile
        model.statfile.close()
        model.pdfdump.close()
        #if all functions were plotted, change name of pdf
        if plottinglist == {Menc:['r','M'],psi:['r',r'$\psi$'],
                            Jc2:['E',r'$J_c^2$'],g:['E','g'],G:['E','G'],
                            f:['E','f'],rate:[r'$u^2$',r'$\frac{dg}{dlnr_p}$']}:
            oldname = '{0}/{1}_master.pdf'.format(model.directory,model.name)
            newname = '{0}/{1}_complete.pdf'.format(model.directory,model.name)
            os.system('mv {0} {1}'.format(oldname,newname))
        print('\a')
        return Mencgood,psigood,Jc2good,ggood,Ggood,fgood,rategood
                
    #if process interrupted, close files safely exit
    except KeyboardInterrupt:
        model.statfile.write('\n\nFunction creation cancelled')
        model.statfile.close()
        model.pdfdump.close()
        raise

#a sample of how to run        
if __name__ == '__main__':

    alpha = 1.0
    beta = 4.0
    gamma = 1.5
    r0pc = 1.0
    rho0 = 1e5
    MBH_Msun = 1e3
    name = 'testform'
    GENERATE = False
    from rhomodels import NukerModelRho
    model = NukerModelRho(name,alpha,beta,gamma,r0pc,rho0,MBH_Msun,GENERATE,memo = False)
    Mencgood,psigood,Jc2good,ggood,Ggood,fgood = getrateplot(model)





