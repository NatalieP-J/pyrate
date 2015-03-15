from numpy import *
from rhoratefcns import *
from construction import loaddata,pklread,displaycheck,existcheck,integrator

#standard independent variable array
rtest = arange(-7,7,0.01)
rtest = append(rtest,40)
rtest = insert(rtest,0,-40)
rtest = 10**rtest

#independent variable array for dgdlnrp
utest1 = arange(-7,-4,0.01)
utest2 = arange(-4,1e-4,1e-4)
utest = concatenate((utest1,utest2))
utest = insert(utest,0,-40)
utest = 10**utest
            
def getrate(model,partial = False):
    """
    partial - option to manually set which functions to create
    Generate rate function for <model> and plot if possible.
    """
    #check if display is on
    dcheck = displaycheck()
    #if functions to create have not been manually selected:
    if partial == False:
        
        #if not asked to generate, automatically decide what to evaluate
        if model.generate == False:
            seton,plottinglist = existcheck(model.directory,dcheck)
            #print 'seton = ', seton
        
        #if asked to generate, create everything and plot if possible
        elif model.generate == True:
            seton = {'Menc':"ON",'psi':"ON",'Jc2':"ON",'lg':"ON",'bG':"ON",
                     'f':"ON",'dgdlnrp':"ON"}
            if dcheck == False:
                plottinglist = {'Menc':False,'psi':False,'Jc2':False,'lg':False,
                                'bG':False,'f':False,'dgdlnrp':False}
            if dcheck == True:
                plottinglist = {'Menc':['r','M'],'psi':['r',r'$\psi$'],
                                'Jc2':['E',r'$J_c^2$'],'lg':['E','g'],
                                'bG':['E','G'],'f':['E','f'],
                                'dgdlnrp':[r'$u^2$',r'$\frac{dg}{dlnr_p}$']}
    
    #if functions have been chosen, do not plot any
    elif partial != False:
        seton = partial[0] 
        plottinglist = partial[1] 

    if model.g < 0 or model.b < 2:
        seton = {'Menc':"OFF",'psi':"OFF",'Jc2':"OFF",'lg':"OFF",'bG':"OFF",
                     'f':"OFF",'dgdlnrp':"OFF"}
        plottinglist = {'Menc':False,'psi':False,'Jc2':False,'lg':False,
                                'bG':False,'f':False,'dgdlnrp':False}

    try:
        #print seton,plottinglist
        #dictionary of power law behaviour
        exps = {'dgdlnrp':[2,0]}
        
        if model.b > 3:
            exps['Menc'] = [3-model.g,0]
            exps['psi'] = [-1,-1]
            exps['Jc2'] = [-1,-1][::-1]
            exps['lg'] = [model.g - 0.5, model.b - 0.5][::-1]
            exps['f'] = [model.g - 1.5, model.b - 1.5][::-1]
            exps['bG'] = [model.g - 4, model.b - 4][::-1]
        if model.b < 3:
            exps['Menc'] = [3-model.g, 3-model.b]
            exps['psi'] = [-1,2-model.b]
            exps['Jc2'] = [-1,((4-model.b)/(2-model.b))][::-1]
            exps['lg'] = [model.g - 0.5,(((model.b - 1) + 1)/((model.b - 1) - 1)) - 0.5][::-1]
            exps['f'] = [model.g - 1.5, (((model.b - 1) + 1)/((model.b - 1) - 1)) - 1.5][::-1]
            exps['bG'] = [model.g - 4, (((model.b - 1) - 2)/((model.b - 1) - 1)) - 1][::-1]

        sh = {'Menc':[6,-6,0.03],'psi':[4,-4,0.03],'Jc2':[4,-1,0.01],
              'lg':[4,-1,0.01],'bG':[3,-2,0.1],'f':[5,-1,0.03],'dgdlnrp':[0,-4,0.04]}

        #begin output
        model.statfile.write('GALAXY: {0}\n'.format(model.name))

        #evaluate Menc
        model.statfile.write('Menc:\n')
        up,down,step = sh['Menc']
        rarray,rchange,rstart = rgrid([model],up,down,step)
        #print 'Start Mencgood'
        Mencgood = compute([model],funcMenc,rtest,sh['Menc'],rgrid,exps['Menc'],
                           plottinglist['Menc'],seton['Menc'])
        #print 'End Mencgood'
        #if failed at Menc, abort computation of further functions
        if Mencgood == 0:
            model.statfile.write('Failed to evaluate Menc')
            seton = {'Menc':"FAIL",'psi':"FAIL",'Jc2':"FAIL",'lg':"FAIL",
                     'bG':"FAIL",'f':"FAIL",'dgdlnrp':"FAIL"}
          
        #evaluate psi
        pprereqs = [model,'Model',Mencgood,'Menc']
        model.statfile.write('\npsi:\n')
        #print 'Start psigood'
        psigood = compute(pprereqs,funcpsi,rtest,sh['psi'],rgrid,exps['psi'],
                          plottinglist['psi'],seton['psi'])
        #print 'End psigood'
        #if failed at psi, abort computation of further functions
        if psigood == 0:
            model.statfile.write('Failed to evaluate psi')
            seton = {'Menc':"FAIL",'psi':"FAIL",'Jc2':"FAIL",'lg':"FAIL",
                     'bG':"FAIL",'f':"FAIL",'dgdlnrp':"FAIL"}
        
        #evaluate Jc2
        Jprereqs = [model,'Model',Mencgood,"Menc",psigood,"psi"]
        model.statfile.write('\nJc2:\n')
        #print 'Start Jc2good'
        Jc2good = compute(Jprereqs,funcJc2,rtest,sh['Jc2'],Egrid,exps['Jc2'],
                          plottinglist['Jc2'],seton['Jc2'])
        #print 'End Jc2good'
        #if failed at Jc2, abort computation of further functions
        if Jc2good == 0:
            model.statfile.write('Failed to evaluate Jc2')
            seton = {'Menc':"FAIL",'psi':"FAIL",'Jc2':"FAIL",'lg':"FAIL",
                     'bG':"FAIL",'f':"FAIL",'dgdlnrp':"FAIL"}

        #evaluate g
        lgprereqs = [model,'Model',psigood,"psi"]
        model.statfile.write('\ng:\n')
        #print 'Start ggood'
        ggood = compute(lgprereqs,funclg,rtest,sh['lg'],Egrid,exps['lg'],
                        plottinglist['lg'],seton['lg'])
        #print 'End ggood'
        #if failed at g, abort computation of further functions
        if ggood == 0:
            model.statfile.write('Failed to evaluate g')
            seton = {'Menc':"FAIL",'psi':"FAIL",'Jc2':"FAIL",'lg':"FAIL",
                     'bG':"FAIL",'f':"FAIL",'dgdlnrp':"FAIL"}
                
        #evaluate mathcalG
        bGprereqs = [model,'Model',psigood, "psi",ggood,"lg"]
        
        Gtest = arange(-5,5,0.01)
        Gtest = append(Gtest,40)
        Gtest = insert(Gtest,0,-40)
        Gtest = 10**Gtest
                        
        model.statfile.write('\nG:\n')
        #print 'Start Gggood'
        Ggood = compute(bGprereqs,funcbG,Gtest,sh['bG'],Egrid,exps['bG'],
                        plottinglist['bG'],seton['bG'])
        #print 'End Ggood'
        if model.memo == True:
            model.p1bG = {}
            model.p2bG = {}
            model.p3bG = {}
        
        #if failed at mathcalG, abort computation of further functions
        if Ggood == 0:
            model.statfile.write('Failed to evaluate G')
            seton = {'Menc':"FAIL",'psi':"FAIL",'Jc2':"FAIL",'lg':"FAIL",
                     'bG':"FAIL",'f':"FAIL",'dgdlnrp':"FAIL"}

        #evaluate f
        fprereqs = [model,'Model',Mencgood,"Menc",psigood,"psi"]
            
        ftest = arange(-3,5,0.01)
        ftest = append(ftest,40)
        ftest = insert(ftest,0,-40)
        ftest = 10**ftest
        
        model.statfile.write('\nf:\n')
        #print 'Start fgood'
        fgood = compute(fprereqs,funcf,ftest,sh['f'],Egrid,exps['f'],
                        plottinglist['f'],seton['f'])
        #print 'End fgood'
        #if failed at f, abort computation of further functions
        if fgood == 0:
            model.statfile.write('Failed to evaluate f')
            seton = {'Menc':"FAIL",'psi':"FAIL",'Jc2':"FAIL",'lg':"FAIL",
                     'bG':"FAIL",'f':"FAIL",'dgdlnrp':"FAIL"}

        #evaluate dgdlnrp
        rprereqs = [model,'Model',Jc2good,'Jc2',Ggood,'Ggood',fgood,'fgood']
        model.statfile.write('\nrate:\n')
        #print 'Start rate'
        rategood = compute(rprereqs,funcdgdlnrp,utest,sh['dgdlnrp'],stdgrid,
                           exps['dgdlnrp'],plottinglist['dgdlnrp'],
                           seton['dgdlnrp'])
        
        if rategood != 0:
            ratetot = integrator(1,rategood,-40,0)
        #print 'End rate'
        if rategood == 0:
            ratetot = 0
            model.statfile.write('Failed to evaluate dgdlnrp')
        
        #close plot pdf and statusfile
        model.statfile.close()
        model.pdfdump.close()
        #if all functions were plotted, change name of pdf
        if plottinglist == {'Menc':['r','M'],'psi':['r',r'$\psi$'],
                                'Jc2':['E',r'$J_c^2$'],'lg':['E','g'],
                                'bG':['E','G'],'f':['E','f'],
                                'dgdlnrp':[r'$u^2$',r'$\frac{dg}{dlnr_p}$']}:
            oldname = '{0}/{1}_master.pdf'.format(model.directory,model.name)
            newname = '{0}/{1}_complete.pdf'.format(model.directory,model.name)
            os.system('mv {0} {1}'.format(oldname,newname))
        print('\a')
        return Mencgood,psigood,Jc2good,ggood,Ggood,fgood,rategood,ratetot
                
    #if process interrupted, close files safely exit
    except KeyboardInterrupt:
        model.statfile.write('\n\nFunction creation cancelled')
        model.statfile.close()
        model.pdfdump.close()
        raise

#a sample of how to run        
if __name__ == '__main__':

    alpha = 1.0
    beta = 3.5
    gamma = 1.5
    r0pc = 1.0
    rho0 = 1e5
    MBH_Msun = 1e3
    name = 'testform'
    GENERATE = False
    from rhomodels import NukerModelRho
    model = NukerModelRho(name,alpha,beta,gamma,r0pc,rho0,MBH_Msun,GENERATE,memo = False)
    #partial = [{'Menc':"ON",'psi':"ON",'Jc2':"FAIL",'lg':"FAIL",'bG':"FAIL",'f':"FAIL",'dgdlnrp':"FAIL"},{'Menc':['r','M'],'psi':['r',r'$\psi$'],'Jc2':False,'lg':False,'bG':False,'f':False,'dgdlnrp':False}]
    result = getrate(model)#,partial)





