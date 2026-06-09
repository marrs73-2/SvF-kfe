# -*- coding: UTF-8 -*-
from __future__ import division
#from  numpy import *
import numpy as np
import time 
#from  time import *
import os

import COMMON as co
from ssop_session import *
import ssop_config
#from Ssop import *

# import Lego_symFun


#from Lego    import *
from CVSets  import *
from GaKru   import *
import SurMin
from SurMin  import SurMinMethod, SurMinOld
from Pars    import *
from Tools   import *
from Task    import Var_to_Grd, FillNaNAll, setUse_var
from ModelFiles import to_logOut


from SolverTools import *

#from pyomo.environ import *
from pyomo.opt import SolverFactory


#from PyomoEverestEnv import *
#==========================


#MU_LABAL = 2147483647

buf = ""

#Mng = ''
# import sys
# import io
# sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding='utf-8', errors='replace') #kfe_added

def SvFstart19 ( Task ) :
    """
    Main method for SvF, envokes optimization.
    It's only called from the end of written StartModel.py file for a specific Task.
    """

    full_start = time.time()
    print ('\n\n\nStart SvFstart')

    # Switch gr to grd bacause numerical optimization is about to start
    setUse_var(False)  # 25.10

    print(f"co.resF = {co.resF}, co.Penalty = {co.Penalty}")
    # If co.resF is None, Penalty is not read from the .res file. Otherwise, it is.
    if co.resF is None: 
        if co.ResAux != '':
            # Set custom name of .res file for writing in later
            co.resF = co.ResAux + '.res'
        else:
            # If co.resF isn't custom set name of .res file based on .mng file name for writing in later. 
            co.resF = co.mngF[:co.mngF.rfind('.')] + '.res'
    else :
        # Set name of .res file based on .mng file name if co.resF isn't custom. 
        if co.resF == '' :  co.resF = co.mngF[:co.mngF.rfind('.')]+'.res'

        # Try getting Penalty from .res file
        try:
                with (open(co.resF,'r') as f):
                    s = f.readline().strip().replace(',', ' ').replace('   ', ' ').replace('  ', ' ').replace(' ', ',')  
                    Penal = readListFloat19 (s)
                    print ('read PENALTY:', Penal)\
                    
                    # save read initinal regularization coefficients back to configuration file
                    co.Penalty = Penal
        except IOError as e:
                print ("******* Can''t open RES file: ", co.resF)


    # Idk what's that for
    co.Penalty = [0.03 if x is None else x for x in co.Penalty]

    # Make solver for low-level of surrogate optimization
    co.optFact = Factory(SvF.SolverConfigFileNames[SvF.SolverNameLow], co.SolverNameLow)

    # optStep is basically Penalty, but a little smaller. Idk what was it for...
    if type(co.OptStep) is str :
        co.OptStep = [float(co.OptStep) * p for p in co.Penalty[:]]

    if co.CVNumOfIter < 0: pass
    # If iterations number is set to zero then just calculate CVError of initial Penalty without optimization
    elif co.CVNumOfIter == 0:
            get_sigCV(co.Penalty, -1)
    else :
        # Start surrogate optimization to find better regularization parameters if iterations are > 0.
        # Two ways of optimization are available: 1) new based on spotoptim module 2) old self-written
        if co.Use_spotoptim == True:
            points = SurMinMethod ( co.CVNumOfIter, co.OptStep, co.Penalty, get_sigCV_for_spotoptim, Task )
        else:
            points, step = SurMinOld ( co.CVNumOfIter, co.OptStep, co.ExitStep, co.Penalty, get_sigCV, Task )

        # Write all points to a .res file
        with open(co.resF,'a') as f:     
            f.write( 'Points:' )
            for p in points :  f.write( 'Num '+str(p.Num) + ' Val ' + str(p.Val) + ' Arg ' + str(p.Arg) + '\n')

    Task.ReadSols('')
    Gr = Task.Gr

    if co.CVNumOfIter > 0 :
      if  co.CVNoBorder  :   #  на границах не учитываем
        for s in co.notTrainingSets[ 0] : Gr.mu0[s]=0
        for s in co.notTrainingSets[-1] : Gr.mu0[s]=0
        NoRnoB = sum ( Gr.mu0[s]() for s in Gr.F[0].sR )
        print ('***** MSD_NoBorder', np.sqrt(Gr.F[0].NoR* Task.defMSDVal ( Gr, 0 ) /NoRnoB)*Gr.F[0].V.sigma)
#        print '***** MSD_NoBorder', np.sqrt(Gr.F[0].NoR*Gr.F[0].MSD()()/NoRnoB)*Gr.F[0].V.sigma
        for s in Gr.F[0].sR : Gr.mu0[s]=1

    for f in Task.Funs :
      if not f.param :
#        f.SaveTbl('')
        if co.SavePoints : f.SavePoints()
#        if co.SaveDeriv and f.type != 'p' :  f.SaveDeriv ( "" )
#        if co.SaveGrid=='Y' and f.dim == 2 and f.type != 'p':     f.SaveGrid ( co.TranspGrid, '' )

    # Draw optimization graph (currently only for 1D) with all points evaluated in optimization
    # and other evenly spaced additionaly evalueted points 
    if co.DrawOpt == True:
        # If optimization didn't happen fill the OptPoints array with initial value
        if Task.OptPoints == []: Task.OptPoints = [co.Penalty]

        plotOptim (Task)

    print ("EofCulc")
    print ('TIME: ', time.time() - full_start)
    return


##################################################################


#optEstim = sys.float_info.max

def printMSD () :
    printS ("sol SD%: |")
    for ifu, fu in enumerate(co.Task.Funs) :
            if  fu.type == 'tensor' : continue
            if fu.V.dat is None or fu.param: continue       #  23.11
#            if fu.mu is None:  continue                    #  23.11
            if not co.Task.DeltaVal is None: fu.MSDv = co.Task.defMSDVal ( Gr, ifu )
            else:
                if fu.MSDmode == 'MSDrel':  fu.MSDv = fu.MSDrel(fu.measurement_accur)          # 21.02.2023
                else :                      fu.MSDv = fu.MSDnan()
#            if co.Task.defMSDVal is None : fu.MSDv = fu.MSDnan()  ### ()
 #           else                         : fu.MSDv = co.Task.defMSDVal ( Gr, ifu )
            printS (fu.V.name, np.sqrt(fu.MSDv)*100.,' |')
    print ('')


def testEstim (Gr, k) :  # k - ValidationSets
    Var_to_Grd()
    for ifu, fu in enumerate(co.Task.Funs):
#        print ("BBBMMMMMMMMMMMMMMMMMMMMMMMMMMMMM", ifu, fu.name)
        if fu.type == 'tensor' : continue
        if fu.mu is None: continue                     #  2024.01
        if fu.V.dat is None or fu.param: continue       #  23.11
 #       print ("MMMMMMMMMMMMMMMMMMMMMMMMMMMMM", fu.name)
#        if fu.NoR > co.CV_NoRs[0] : continue               # 25/04
        spart = 0
        npart = 0
        if fu.CVerr is None: fu.CVerr = np.zeros(fu.NoR, np.float64)   # 04.2023

        for s in fu.ValidationSets[k]:
            if s >= fu.NoR        : continue                    #  25/04
            if np.isnan(fu.V.dat[s]) : continue

            if not co.Task.DeltaVal is None: err = Task.DeltaVal(Gr, ifu, fu.V.dat, s) #** 2
            elif  fu.MSDmode == 'MSDrel':    err = fu.delta_rel(s) #** 2          # 21.02.2023
            else:                            err = fu.delta(s) #** 2
   #         if s==0 :
    #            print ('err=', err)
     #           1/0
            spart += err ** 2
            npart += 1
            fu.CVerr[s] = err                                   # 04.2023
#            print ('CVerr', s, fu.CVerr[s])
 #       print (fu.CVerr)
        if npart == 0:  print (spart, npart, 'NoVal', 'NoVal',)
        else:
            if fu.MSDmode == 'MSDrel':
                print ('  ', k, fu.name, spart, npart, np.sqrt(spart / npart), np.sqrt(spart / npart) * 100.,)
            else :
                print ('  ', k, fu.name, spart, npart, np.sqrt(spart / npart), np.sqrt(spart / npart) / fu.V.sigma * 100.,)
    # OLTCHEV
    #                NDT = Gr.F[d].NDT
    #                sumTbl  = sum((Gr.F[d].tbl[s, Gr.F[d].V.num] != NDT) * Gr.F[d].tbl[s,Gr.F[d].V.num] for s in ValidationSets[k])
    #               sumFTbl = sum((Gr.F[d].tbl[s, Gr.F[d].V.num] != NDT) * Gr.F[d].Ftbl ( s )() for s in ValidationSets[k])
    #              sumDelta = sum((Gr.F[d].tbl[s, Gr.F[d].V.num] != NDT) * Gr.F[d].delta ( s )() for s in ValidationSets[k])
    #             print sumTbl,sumFTbl, sumDelta/ npart,

####    if co.CVNoBorder:  # на границах не учитываем
####        if k == 0 or k == len(fu.ValidationSets) - 1:   continue
        fu.CVresult.append([spart, npart])


def getEstimCV(Gr) :
        printS ("\nParts |")      #    printS ("\n???????Parts |")
        Estim = 0
        NumOfFuns = 0
        for ifu, fu in enumerate (co.Task.Funs ):
            if fu.type == 'tensor': continue
            if fu.mu is None : continue             #  20.01
            if (fu.V.dat is None) or fu.param:  continue
            sCrVa = 0;  nCrVa = 0
            for CVr in fu.CVresult :
                sCrVa += CVr[0]
                nCrVa += CVr[1]
            if nCrVa > 0 :
                sCrVa = np.sqrt(sCrVa / nCrVa)
            if fu.MSDmode == 'MSDrel':  print (sCrVa, sCrVa * 100.,)
            else:                       print (sCrVa, sCrVa / fu.V.sigma * 100.,)
            fu.sCrVa = sCrVa
            if fu.MSDmode == 'MSDrel':  Estim += sCrVa                  # 21.02.2023
            else :                      Estim += sCrVa / fu.V.sigma
            NumOfFuns += 1
        Estim = Estim / NumOfFuns * 100
        return Estim


def get_sigCV( Penal, itera):
    co.CV_Iter = itera
    Task = co.Task
#    reload (Model)
    Task.ReadSols('')
#    SvF.Penalty = Penal
    if not (SvF.feasibleSol is None) : SvF.feasibleSol(Penal)

    print ('for Penal: ', Penal ) #, end=' ')

    if SvF.OptMode == 'SurMinOpt' :
  #      co.Use_var = True       # 29
        setUse_var(True)       # 25.10

        Gr = Task.createGr(Task, Penal)   # обновлем на каждой итерации
        Grd_to_Var()
        #co.Use_var = False       # 29
        setUse_var(False)  # 25.10

        resultss = solveProblemsNl(Gr, '', co.RunMode[0])               #  tmp
        Gr.solutions.load_from(resultss[0])
        Var_to_Grd()
        Task.SaveSols('.tmp')

        print ('OBJ',Gr.OBJ())
        Estim = Gr.OBJ()

    elif SvF.OptMode == 'SvF':
        
        FillNaNAll ()
     #   co.Use_var = True       # 29
        setUse_var(True)  # 25.10

        Gr = Task.createGr(Task, Penal)   # обновлем на каждой итерации
        Grd_to_Var()
        #        co.Use_var = False       # 29
        setUse_var(False)  # 25.10

        for fu in Task.Funs :  fu.CVresult = []

        if itera <= 0 : printS (' Load on Start ');   printMSD()

        resultss = solveProblemsNl(Gr, '', co.RunMode[0])               #  tmp
        Gr.solutions.load_from(resultss[0])
        Var_to_Grd()

        Task.SaveSols('.tmp')

        print ('OBJ',Gr.OBJ())
        printMSD()


        if not Task.OBJ_U is None :
            Estim = Task.OBJ_U(Task)()
            print ('**************KK=', Estim)

        elif co.CVNumOfIter != 0 :
            star = time.time()

            if co.RunMode[2] != 'L':
    #            resultss = solveProblemsNl ( Gr, co.notTrainingSets, co.RunMode[2] )
                resultss = solveProblemsNl ( Gr, '*', co.RunMode[2] )              # All tests
                for nres, res in enumerate (resultss) :
                    Gr.solutions.load_from(res)
                    testEstim(Gr, nres)
            else:       ## co.RunMode[2] == 'L' :
                res_num = 0
 #               print(co.CV_NoSets, "OOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOO")
                for k in range(co.CV_NoSets):                                  # LOAD RES,  culculation
                    #               if co.NotCulcBorder :  #  границ не считаем
                    #                  if k == 0 or k == len(ValidationSets) - 1:   1/0;  continue  #########  ???????????????????
                    #             printS (str(k)+' |')
                    #                results = solveProblemsNl(Gr, [co.notTrainingSets[k]], co.RunMode[2])[0]  #!! РАБОТАЕТ ТОЛЬКО ДЛЯ ОДНОГО resultss
                    results = solveProblemsNl(Gr, k, co.RunMode[2])[0]  #!! РАБОТАЕТ ТОЛЬКО ДЛЯ ОДНОГО resultss
                    res_num += 1
                    Gr.solutions.load_from(results)
                    testEstim(Gr, k)

            Estim = getEstimCV(Gr)
            print ('\tEstim% '+str(Estim)+"\tTime  "+ str(time.time() - star))

        else : Estim = -1

    elif SvF.OptMode == 'SurMin':
       Estim = SvF.ObjectiveFun (Penal)

    if Estim < co.optEstim :
            co.optEstim = Estim
            Task.RenameSols( '.tmp', '.sol' )

            setMuToTeach_k('')     #  mu = 1   ##############  26/04/24
            Task.ReadSols()                    ##############  26/04/24
            Grd_to_Var()
            with open(co.resF,'w') as f:      #  RES filewrite
#                print >> f, [p for  p in Penal]
                f.write (str(Penal))
                for fu in Task.Funs :           # v21
#                      if fu.mu is None: continue                     #  2023.11
                      if fu.type == 'tensor': continue
                      if fu.V.dat is None or fu.param: continue  # 23.11
#                      str_wr = '\n'+fu.nameFun()+' '
                      str_wr = fu.nameFun()+' '
                      if co.CVNumOfIter > 0:
                          if fu.MSDmode == 'MSDrel':   str_wr += " CV% " + str(fu.sCrVa*100)
                          else:                     str_wr += " CV% " + str(fu.sCrVa/fu.V.sigma*100)
                      str_wr += ' SD% ' + str(np.sqrt(fu.MSDv)*100) + " CV " + str(fu.sCrVa) \
                                +' SD ' + str(np.sqrt(fu.MSDv)*fu.V.sigma) + ' sig '+str(fu.V.sigma)

                      f.write ( '\n' + str_wr )
                      print (str_wr)
                      to_logOut ( str_wr )
                f.write( '\n'+'Estim ' + str(Estim))
                to_logOut ( 'Estim ' + str(Estim) )
                to_logOut ( 'OBJ ' + str(Gr.OBJ()) )

#                print >> f, 'Estim',Estim
                if Task.print_res != None :
                    Task.print_res(Task, Penal, f)
    return Estim

def get_sigCV_kfe( Penal):
    print_here = False
    Task = co.Task
    Task.ReadSols('')

    print ('for Penal: ', Penal ) #, end=' ')

    FillNaNAll ()
    setUse_var(True)  # 25.10

    Gr = Task.createGr(Task, Penal)   # обновлем на каждой итерации
    Grd_to_Var()
    setUse_var(False)  # 25.10

    for fu in Task.Funs :  fu.CVresult = []

    resultss = solveProblemsNl(Gr, '', co.RunMode[0])               #  tmp
    Gr.solutions.load_from(resultss[0])
    Var_to_Grd()

    #Task.SaveSols('.tmp')

    print ('OBJ',Gr.OBJ())
    printMSD()


    if co.CVNumOfIter != 0 :
        star = time.time()

        if co.RunMode[2] != 'L':
#            resultss = solveProblemsNl ( Gr, co.notTrainingSets, co.RunMode[2] )
            resultss = solveProblemsNl ( Gr, '*', co.RunMode[2] )              # All tests
            for nres, res in enumerate (resultss) :
                Gr.solutions.load_from(res)
                testEstim(Gr, nres)
        else:       ## co.RunMode[2] == 'L' :
            res_num = 0
#               print(co.CV_NoSets, "OOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOO")
            for k in range(co.CV_NoSets):                                  # LOAD RES,  culculation
                #               if co.NotCulcBorder :  #  границ не считаем
                #                  if k == 0 or k == len(ValidationSets) - 1:   1/0;  continue  #########  ???????????????????
                #             printS (str(k)+' |')
                #                results = solveProblemsNl(Gr, [co.notTrainingSets[k]], co.RunMode[2])[0]  #!! РАБОТАЕТ ТОЛЬКО ДЛЯ ОДНОГО resultss
                results = solveProblemsNl(Gr, k, co.RunMode[2])[0]  #!! РАБОТАЕТ ТОЛЬКО ДЛЯ ОДНОГО resultss
                res_num += 1
                Gr.solutions.load_from(results)
                testEstim(Gr, k)

        Estim = getEstimCV(Gr)
        print ('\tEstim% '+str(Estim)+"\tTime  "+ str(time.time() - star))

    else : Estim = -1
    return Estim
import sys
from io import StringIO

def suppress_print(func):
    """Декоратор для подавления print внутри функции"""
    def wrapper(*args, **kwargs):
        original_stdout = sys.stdout
        sys.stdout = StringIO()
        try:
            result = func(*args, **kwargs)
        finally:
            sys.stdout = original_stdout
        return result
    return wrapper

#@suppress_print
def get_sigCV_for_spotoptim(Penal):
    co.CV_Iter += 1
    #Penal = [[0.9]]
    print(f"PRINTING init_x_to_y_dict: on {co.CV_Iter} = {SurMin.init_x_to_y_dict}, \n Penal = {Penal}")
    Penal = np.atleast_2d(np.array(Penal))
    n_samples = Penal.shape[0]
    if n_samples > 1:
        print("В get_sigCV передано несколько точек - возвращаю значения для точек инициализации")
        print(f"Массив возвращаемых значений: {np.array(list(SurMin.init_x_to_y_dict.values()))}")
        return np.array(list(SurMin.init_x_to_y_dict.values()))
    
    if Penal.ndim == 2 and Penal.shape[0] == 1:
        Penal = Penal[0]
    
    print(f"Penal after 'squeeze': {Penal}")

    Task = co.Task
#    reload (Model)
    Task.ReadSols('')
#    SvF.Penalty = Penal
    if not (co.feasibleSol is None) : co.feasibleSol(Penal)

    if co.OptMode == 'SvF':

        setUse_var(True)  # 25.10

        Gr = Task.createGr(Task, Penal)   # обновлем на каждой итерации
        Grd_to_Var()
        setUse_var(False)  # 25.10

        for fu in Task.Funs :  fu.CVresult = []

        if co.CV_Iter <= 0 : printS (' Load on Start ');   printMSD()

        resultss = solveProblemsNl(Gr, '', co.RunMode[0])               #  tmp
        Gr.solutions.load_from(resultss[0])
        Var_to_Grd()

        Task.SaveSols('.tmp')
        print ('OBJ',Gr.OBJ())
        printMSD()


        if co.CVNumOfIter != 0 :
            star = time.time()

            if co.RunMode[2] != 'L':
    #            resultss = solveProblemsNl ( Gr, co.notTrainingSets, co.RunMode[2] )
                resultss = solveProblemsNl ( Gr, '*', co.RunMode[2] )              # All tests
                for nres, res in enumerate (resultss) :
                    Gr.solutions.load_from(res)
                    testEstim(Gr, nres)
            else:       ## co.RunMode[2] == 'L' :
                res_num = 0
 #               print(co.CV_NoSets, "OOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOO")
                for k in range(co.CV_NoSets):                                  # LOAD RES,  culculation
                    #               if co.NotCulcBorder :  #  границ не считаем
                    #                  if k == 0 or k == len(ValidationSets) - 1:   1/0;  continue  #########  ???????????????????
                    #             printS (str(k)+' |')
                    #                results = solveProblemsNl(Gr, [co.notTrainingSets[k]], co.RunMode[2])[0]  #!! РАБОТАЕТ ТОЛЬКО ДЛЯ ОДНОГО resultss
                    results = solveProblemsNl(Gr, k, co.RunMode[2])[0]  #!! РАБОТАЕТ ТОЛЬКО ДЛЯ ОДНОГО resultss
                    res_num += 1
                    Gr.solutions.load_from(results)
                    testEstim(Gr, k)

            Estim = getEstimCV(Gr)
            print ('\tEstim% '+str(Estim)+"\tTime  "+ str(time.time() - star))

        else : Estim = -1

    if Estim < co.optEstim :
            co.optEstim = Estim
            Task.RenameSols( '.tmp', '.sol' )

            setMuToTeach_k('')     #  mu = 1   ##############  26/04/24
            Task.ReadSols()                    ##############  26/04/24
            Grd_to_Var()
            with open(co.resF,'w') as f:      #  RES filewrite
#                print >> f, [p for  p in Penal]
                f.write (str(Penal))
                for fu in Task.Funs :           # v21
#                      if fu.mu is None: continue                     #  2023.11
                      if fu.type == 'tensor': continue
                      if fu.V.dat is None or fu.param: continue  # 23.11
#                      str_wr = '\n'+fu.nameFun()+' '
                      str_wr = fu.nameFun()+' '
                      if co.CVNumOfIter > 0:
                          if fu.MSDmode == 'MSDrel':   str_wr += " CV% " + str(fu.sCrVa*100)
                          else:                     str_wr += " CV% " + str(fu.sCrVa/fu.V.sigma*100)
                      str_wr += ' SD% ' + str(np.sqrt(fu.MSDv)*100) + " CV " + str(fu.sCrVa) \
                                +' SD ' + str(np.sqrt(fu.MSDv)*fu.V.sigma) + ' sig '+str(fu.V.sigma)

                      f.write ( '\n' + str_wr )
                      print (str_wr)
                      to_logOut ( str_wr )
                f.write( '\n'+'Estim ' + str(Estim))
                to_logOut ( 'Estim ' + str(Estim) )
                to_logOut ( 'OBJ ' + str(Gr.OBJ()) )

#                print >> f, 'Estim',Estim
                if Task.print_res != None :
                    Task.print_res(Task, Penal, f)
    return Estim


def plotOptim (task): #kfe_added - отображение "оптимизации"
    import matplotlib.pyplot as plt 
    
    for i, func in enumerate(task.Funs):
        if func.type == 'tensor': continue
        dim = func.dim
        points = co.DrawOptPoints
        print(f"Начало отрисовки коэффициентов регуляризации для функции {func.name} с dim {func.dim}")

        if dim == 1:
            if co.DrawOptMode == "Relative":
                if len(task.OptPoints) > 1:
                    center = task.OptPoints[0].Arg
                else:
                    center = task.OptPoints[0]

                width = center / co.DrawOptWidthCoef  
                x = np.linspace(center - width, center + width, points)
            elif co.DrawOptMode == "Absolute":
                x = np.linspace(co.DrawOptSegment[0], co.DrawOptSegment[1], points)

            y = np.array([get_sigCV_kfe([xi]) for xi in x])
            #y = np.array([get_sigCV([xi], 1) for xi in x])
            
            plt.figure(figsize=(10, 6))
            plt.plot(x, y, 'b-o', linewidth=2, markersize=6)
            plt.grid(True, alpha=0.3)
            plt.title('График оптимизации: CV% от Penal')
            plt.xlabel('Коэффициент регуляризации')
            plt.ylabel('CV% - ошибка кросс-валидации')
            
        elif dim == 2:
            if co.DrawOptMode == "Relative":
                if len(task.OptPoints) > 1:
                    center = [task.OptPoints[0].Arg, task.OptPoints[1].Arg]
                else:
                    center = [task.OptPoints[0], task.OptPoints[1]]
                width = [center[0] / co.DrawOptWidthCoef, center[1] / co.DrawOptWidthCoef]
                
                x = np.linspace(center[0] - width[0]/2, center[0] + width[0]/2, points)
                y = np.linspace(center[1] - width[1]/2, center[1] + width[1]/2, points)
            elif co.DrawOptMode == "Absolute":
                x = np.linspace(co.DrawOptSegment[0][0], co.DrawOptSegment[0][1], points)
                y = np.linspace(co.DrawOptSegment[1][0], co.DrawOptSegment[1][1], points)

            X, Y = np.meshgrid(x, y)
            Z = get_sigCV(X, Y)
            
            plt.figure(figsize=(10, 8))
            plt.contourf(X, Y, Z, levels=20, cmap='viridis')
            plt.colorbar(label='f(x, y)')
            plt.scatter(center[0], center[1], color='red', s=100, marker='*')
            plt.title('2D функция на сетке')
            plt.xlabel('X')
            plt.ylabel('Y')
            plt.grid(True, alpha=0.3)

        if co.CVNumOfIter >= 0:
            x_coef = [p.Arg for p in task.OptPoints]
            y_sigma = [p.Val for p in task.OptPoints]
            plt.scatter( x_coef, y_sigma, linewidth=2, color='red', label='Точки Соколова')
            labels = [f"{num}" for num in range(len(task.OptPoints))]

            # подписи к каждой точке
            for _, (xi, yi, label) in enumerate(zip(x_coef, y_sigma, labels)):
                plt.annotate(label, (xi, yi), 
                            xytext=(5, 5),  # смещение текста
                            textcoords='offset points',
                            fontsize=10,
                            bbox=dict(boxstyle='round,pad=0.3', facecolor='yellow', alpha=0.7))
            
        plt.savefig(f'optPlot{i}.png')

