# -*- coding: cp1251 -*-
from __future__ import division
#from  numpy import *
import numpy as np
from copy   import *
#from pyomo.environ import *
import pyomo.environ as py

from SolverTools import Factory #*
from CVSets  import *
from GaKru   import *
from Pars    import *
from Tools   import *
#from Task    import Grd_to_Var
from Task    import Var_to_Grd, FillNaNAll, setUse_var
from ModelFiles import to_logOut

#from StartModel import Model
#import Model as Model


from SolverTools import *

import COMMON as co


#   for index, string in enumerate(strings):
#   np.linspace(0, 2, 9)  # 9 ����� �� 0 �� 2 ������������

# ������������� ��������� 2 ������� pol ( C, B1, B2) 
# 1) ��������� ��������� - ���������� ��� �����  B1i, B2i, SIGi
# 2) ��������� �������� ���������������� ������������ 2-��� �������:
#    - ����� ����������� � ������� ������ (��� ������ ����� �� �������� ��������, ��� ������ ���)
#    - ����������� ����� �� �������� ����������� (�� �������� ������������� ��� ������������ ������),
#      �������� �������� ����������� ���, ����� ������� �� �������� � ���������� ����� ��� ������ ��������.
# 3) ����������� ����� �������� B1 � B2, ��� ��� ����������� SIG. ����������� �� ��� 2.
# ����� �� ���-�� �������� ��� �� �������� ���� (���������� �����)
# ����� �   (sum ( (pol ( C, B1i, B2i) - SIGi))**2/ri**(2+q) i in range(I)) + ??? (C20+C02+C11)*reg => min
#  ��� ri = sqrt ()
#

coprintL = 0    # ����� ���������� �� ������
#from pyomo.opt import SolverFactory

#def LittleFactory (optFile , py_max_iter, py_tol ): # , warm_start_bound_push      =1e-6,
 #   opt = SolverFactory('ipopt')
#
 #   opt.options["print_level"] = 4 #4 #6
  #  opt.options['warm_start_init_point']      = 'yes'
   # opt.options['max_iter']             = py_max_iter
    #opt.options["tol"]                  = py_tol
#    #opt.options['acceptable_tol']       = 1e-10
#    opt.options['print_user_options']      = 'yes'

#    if  optFile != None :
 #       if co.RunMode[0] != 'L' or co.RunMode[2] != 'L' :
  #          makeSolverOptionsFile(co.tmpFileDir+'/'+optFile, "ipopt", opt.options)
 #   return opt


def CrePolyPow2 ( dim, maxP ) :            # np.arrays
        polyPow = []
        for i in range((maxP+1)**dim) :
            powers = np.zeros(dim)
            sumP = 0
            tmpP = i
            for d1 in range(dim) :
                d = dim-d1-1
                powers[d] = int( tmpP/((maxP+1)**d) )
                sumP += powers[d]
                tmpP -= powers[d] * ((maxP+1)**d)
            if sumP <= maxP : polyPow.append(powers)
        sortPow = []
        for pw in range (maxP+1) :
            for p in polyPow :
                if sum ( p ) == pw :
                    sortPow.append ( p )
        return sortPow


class Poly :                 # 
  def __init__ ( self, coef, pows ) :     # 
      self.coef  = coef
      self.pows  = pows
      self.NoPo  = len(pows)

  def Culc ( self, arg ) :         # sum i   coef[i] * arg[0]**pows[i][0] * .....
    val = 0
    for p in range ( self.NoPo ) :
        term = 1
        for d in range (len(arg)) :
            if self.pows[p][d] != 0 : term *= arg[d]**self.pows[p][d]
        val += term * self.coef[p]
    return val

  def CulcNNN ( self, arg ) :    # sum i   coef[i] * arg[0]**pows[i][0] * .....   �� ������� ��� pyoma  :(
      return sum ( self.coef[p] * prod ( arg**self.pows[p] ) for p in range (self.NoPo) )

  def CulcShift ( self, arg, point ) :
      return point.Val + self.Culc (arg-point.Arg)


                       
def cosFi ( v1, v2 ) :
#            dim = len (v1)
 #           scalProd = sum ( v1[c]*v2[c] for c in range (dim) ) 
  #          v1_l2    = sum ( v1[c]**2    for c in range (dim) ) 
   #         v2_l2    = sum ( v2[c]**2    for c in range (dim) )
    #        return  scalProd / sqrt (v1_l2*v2_l2)                                         # cos ��������
            return   sum(v1*v2) / ( norma(v1) * norma(v2) )
        
def angleV1V2 ( v1, v2 ) :
            try:

                return  np.arccos ( cosFi ( v1, v2 ) )    * 180 / np.pi                      # ���� ��������
            except :
                print ('cosFi', cosFi ( v1, v2 ))
                exit(-1)

def int_angleV1V2 ( v1, v2 ) :
            return  int ( round ( angleV1V2 ( v1, v2 ) ) )                             # ���� �������� int


def norma ( v ) :
        return np.sqrt ( sum ( v**2 ) )


def distance ( p1, p2 ) :
        return norma ( p1.Arg-p2.Arg )


class Point :                 # 
    def __init__ ( self, Arg, Val=0, Num=0 ) :     # 
        self.Arg  = Arg.copy()
        self.Val  = Val
        self.Num  = Num
        self.wieght = 0
      
    def prin(self) :
        print ('Num', self.Num, 'Val', self.Val, 'Arg', self.Arg)


def AddPoint ( points, Arg, Val ) :
    poi = Point ( Arg, Val, len(points) )
    grad = (Val - points[-1].Val) / distance (poi, points[-1])
    points.append ( poi )
    if points[-2].Val < points[-1].Val :    # swap
        poin = deepcopy (points[-2])
        points[-2] = points[-1] 
        points[-1] = poin
        ismin = '***'
    else :
        ismin = 'min'
    return points[-1].Arg, ismin, grad


def wieght ( npp, points, farWieght ) :
                dim = len (points[-1].Arg)
                if npp < len(points)-1-dim : return np.exp ( - distance (points[npp],points[-1]) * farWieght )
                else                    : return 1
#                return exp ( - distance (points[np],points[-1]) * farWieght )


def wieghtCulc ( points, farWieght ) : 
    for ip, p in enumerate(points) :
        p.wieght = wieght ( ip, points, farWieght )


        
def CulcCoef ( opt, points, curvPenal, farWieght ) :
                dim = len(points[0].Arg)
                wieghtCulc ( points, farWieght )
                lastWeight = sum ( points[ip].wieght  for ip in range (len(points)-dim-1) )
                normDist   = sum ( points[ip].wieght  for ip in range (len(points)-1    ) )
                partWeight = lastWeight / normDist
#                print 'normDist', normDist, lastWeight, partWeight

                CM = py.ConcreteModel()

                polPow = CrePolyPow2 ( dim, 2 )
                polLen = len(polPow)
                CM.Cpoly  = py.Var (range(polLen), domain=py.Reals, initialize = 1 )
                pol = Poly (CM.Cpoly, polPow)
                
                CM.Cpoly[0].value = 0;    CM.Cpoly[0].fixed = True
                
                def obj_expression(CM):
                    return ( 1 / normDist
                              * sum ( (pol.Culc(points[p].Arg-points[-1].Arg)-(points[p].Val-points[-1].Val))**2
                                  * points[p].wieght for p in range (len(points)-1) )
                           + curvPenal
                              * sum ( CM.Cpoly[c]**2 for c in range (dim+1, polLen) )
                           )
                CM.OBJ = py.Objective(rule=obj_expression)
                
                results = opt.solve(CM)                        
                CM.solutions.load_from(results)
                if str(results.solver.termination_condition) != 'optimal':
                    print ("Stst:", results.solver.termination_condition, '\n')
                    
                obj = CM.OBJ()
#                print 'CM.OBJ', obj
                if polLen > dim and obj > 0 :
                    partPen = curvPenal * sum (CM.Cpoly[c]()**2 for c in range (dim+1, polLen)) / obj
                else :    partPen = 0

                pol.coef = np.array([ CM.Cpoly[c]() for c in range(polLen) ])
                return pol, partWeight, partPen



def Prognose ( opt, pol, point, step ) :
                dim = len(point.Arg)
                CMP = py.ConcreteModel()

                norm = np.sqrt ( sum ( pol.coef**2 ) )
                def ini_circ (CMP, c): return - pol.coef[c+1] / norm * step * 0.999999
                CMP.nInc  = py.Var ( range(dim), domain=py.Reals, initialize = ini_circ )

                def nInc_ge(CMP,a) :
                    if  step < 0.7*point.Arg[a] : return py.Constraint.Skip                            #  ��� �� �������� ������ ��������
                    else                        : return ( CMP.nInc[a] >= - 0.7*point.Arg[a] )      #  for Arg >= 0
                CMP.cnInc_ge = py.Constraint( range (dim), rule=nInc_ge )

                def circArg(CMP) :  return ( sum ((CMP.nInc[a])**2 for a in range(dim)) / step**2 <= 1 )
#                def circArg(CMP) :  return ( sum ((CMP.nInc[a])**2 for a in range(dim)) <= step**2 )
                CMP.cnInc = py.Constraint( rule=circArg )

                def obj_expression(CMP):  return ( pol.Culc(CMP.nInc) )
                CMP.OBJ = py.Objective(rule=obj_expression)

                try:
                    results = opt.solve(CMP)                        
                    CMP.solutions.load_from(results)
                    if str(results.solver.termination_condition) != 'optimal':
                        print ("Stst:", results.solver.termination_condition, '\n')
                        print ("************��������� �������� ������� (�������)")
                        for a in range(dim) :
                            CMP.nInc[a].value = - pol.coef[a+1] / norm * step * 1 
                except:
                    print ("EXCEPT ************��������� �������� ������� (�������)")
                    for a in range(dim) :
                        CMP.nInc[a].value = - pol.coef[a+1] / norm * step * 1 
                    

                print ('Incr:    ', [ CMP.nInc[a]() for a in range(dim) ], '    step=', step )
                prognVal = CMP.OBJ()
                if  sum ( (CMP.nInc[a]())**2 for a in range(dim))/step**2  < 0.94 : Constr = 'Inside'
                else :   Constr = 'Const'
#                print prognVal, point.Val, pol.Culc(CMP.nInc)(), \
 #                     'ang', int_angleV1V2 ( [ CMP.nInc[a]() for a in range(dim) ],pol.coef[1:dim+1] )
                Incr = np.array([ CMP.nInc[a]() for a in range(dim) ])
                return  prognVal+point.Val, Incr+point.Arg, Constr, Incr
#                return  prognVal+point.Val, np.array([ CMP.nInc[a]()+point.Arg[a] for a in range(dim) ]), \
 #                                   Constr, np.array([ CMP.nInc[a]()              for a in range(dim) ])

    
#maxpartPen = .5
maxpartPen = .001
#maxpartPen = 1e-4 #0.3 #0.1
minpartPen = 1e-5

def arrange_farWieght_curvPenal (opt, points, curvPenal, farWieght, Val, nArg) :
 #       print "sAR", '  ', farWieght, curvPenal
        dim = len (nArg)
        for n in range(100) :                    # �������� partPen � �������
            pol, partWeight, partPen = CulcCoef (opt, points, curvPenal, farWieght)
            if   partWeight > maxpartPen :  farWieght *= 1.07
            elif partWeight < minpartPen :  farWieght /= 1.07
            else                         : break
#            print 'QQQQQQQQQQQQQQ', n, farWieght, 'pW', partWeight
        prognVal = pol.CulcShift (nArg, points[-1])
        print ("+AR", '**', farWieght, 'pW', partWeight, curvPenal, partPen, (Val-prognVal), )
        
        farWieghtN = farWieght * 1.001
        pol, partWeight, partPen = CulcCoef (opt, points, curvPenal, farWieghtN)
        prognValN = pol.CulcShift (nArg, points[-1])
        dW = (prognValN-prognVal)/0.001

        curvPenalN = curvPenal * 1.001
        pol, partWeight, partPen = CulcCoef (opt, points, curvPenalN, farWieght)
        prognValN = pol.CulcShift (nArg, points[-1])
        dP = (prognValN-prognVal)/0.001

        if dW==0 :  return farWieght, curvPenal

        delteVal = Val-prognVal    #  > 0 ���� ������� ����
#        ad = 0.03
        ad = 0.05
        if abs (dP) <= abs(dW) :  addW = ad; addP = ad* abs(dP/dW)
        else                   :  addP = ad; addW = ad* abs(dW/dP)
        if dW * delteVal > 0 :  multW=1+addW
        else                 :  multW=1-addW
        if dP * delteVal > 0  : multP=1+addP
        elif dP == 0          : multP=1
        else                  : multP=1-addP

        print ('dW', dW, multW, 'dP', dP, multP, delteVal)

#        for n in range (100) :
#        for n in range (50) :
        for n in range (25) :
            farWieghtN = farWieght * multW
            curvPenalN = curvPenal * multP
            pol, partWeight, partPen = CulcCoef (opt, points, curvPenalN, farWieghtN)
            if multW > 1 and partWeight < minpartPen : break                    #  ����� �� ������� �������
            if multW < 1 and partWeight > maxpartPen : break                    #  ����� �� ����������� �������
            prognValN = pol.CulcShift (nArg, points[-1])
            if delteVal > 0 :
                if prognValN >= Val : break         #  �������������
                elif prognValN <= prognVal : break  #  ������� ���������
                else :
                    prognVal  = prognValN
                    curvPenal = curvPenalN
                    farWieght = farWieghtN
            else :
                if prognValN <= Val : break         #  �������������
                elif prognValN >= prognVal : break  #  ������� ���������
                else :
                    prognVal  = prognValN
                    curvPenal = curvPenalN
                    farWieght = farWieghtN
        print ("+AR", n, farWieght, 'pW', partWeight, curvPenal, partPen, (Val-prognVal))
        return farWieght, curvPenal


def condition ( points, farWieght, opt ) :
    dim = len(points[0].Arg)
    M = py.ConcreteModel()
    M.plane  = py.Var ( range(dim+1), domain=py.Reals, initialize = 1 )              # ���������� ��������� ��������� :
    def Eq1 (M) : return ( sum ( M.plane[i]**2 for i in range(dim) ) == 1 )    #     ����� ��������� ����� �������
    M.cEq1 = py.Constraint( rule=Eq1 )
    def Eq2 (M) : return ( M.plane[dim] <= 0 )                                 #     ��������� .. <=0
    M.cEq2 = py.Constraint( rule=Eq2 )

    normDist = sum ( p.wieght for p in points )
    
    def point_plane2 ( point ) :
        return  ( sum (M.plane[i]*point.Arg[i] for i in range(dim))+M.plane[dim] )**2
                
    def obj_expression(M):
        return  sum ( point_plane2( p ) * p.wieght for p in points ) / normDist
    M.OBJ = py.Objective(rule=obj_expression)

    results = opt.solve(M)                        
    M.solutions.load_from(results)
    if str(results.solver.termination_condition) != 'optimal':
                    print ("Stst:", results.solver.termination_condition, '\n')
    obus = np.sqrt ( M.OBJ() )
 #   for cp in range(dim+1): print cp, M.plane[cp]() 
#    print 'obus', obus
#    return obus, [M.plane[c]() for c in range(dim+1)]
    return obus, np.array([M.plane[c]() for c in range(dim)])

def SetAllArgs (Arg, InArg, stepsIN ) :
    a_in = 0
    for iss,s in enumerate (stepsIN) :
        if s != 0:
            InArg[iss] = Arg[a_in]
            a_in += 1
    return InArg


init_x_to_y_dict = None

def SurMinMethod ( CVNumOfIter, stepsIN, ExitStep, InArg, getVal, Task=None) :
    old_points = []
    old_cCos = []
    old_stepMult = 1
    par_der = []   # deriv
    firstDerec = True
#    opt = LittleFactory ( None, 10000, 1e-9 )  # ����� ����������
    global init_x_to_y_dict
    init_x_to_y_dict = dict()
    opt = Factory (None, co.SolverNameHigh)  # ����� ����������

    Arg = []
    steps = []
    for ia, a in enumerate(InArg) :
        if stepsIN[ia] != 0 :
            Arg.append(a)
            steps.append(stepsIN[ia])
    Arg = np.array (Arg)
    dim = len (Arg)

    curvPenal = 0.001
    print ("len", len (steps) )
    if len (steps)==0 : farWieght = 1   #  ��� ����������
    else:               farWieght = 2/dim/abs(steps[0])
    #  exp ( - farWieght * dist...
    print ('\nstart  farWieght', farWieght)
    step = 1e37
    if CVNumOfIter == 0:
            AllArgs = SetAllArgs(Arg, InArg, stepsIN)
            Val = getVal (AllArgs)
            points = [Point (Arg, Val, 0)]
            for p in points:  p.prin()
            return points, nan

    for itera in range (dim+1) :
        print(f"Started iteration {itera} of SurMinMethod with CVNumOfIter={CVNumOfIter}, dim={dim}")
        if itera == 0 :
            AllArgs = SetAllArgs(Arg, InArg, stepsIN)
            Val = getVal (AllArgs)
#            Val = getVal ( Arg, itera )
            print ('\nITER', itera, 'start', Val, 'st', np.nan,  Arg, '\n')
            points = [Point (Arg, Val, 0)]

        elif itera <= dim :
            init_x_to_y_dict.clear()
            for p in points:
                print(f"p.Arg = {p.Arg}")
                print(f"p.Val = {p.Val}")
            init_x_to_y_dict.update({
                tuple(np.round(np.asarray(p.Arg).flatten(), 8)): float(p.Val)
                for p in points
            })
            nArg = Arg.copy()
#            if stepIN > 0 : stepp = stepIN
 #           else          : stepp = Arg[itera-1] * (-stepIN)
            step_i = steps[itera-1]
            Nattempt = 5
            for attempt in range (Nattempt) :
                nArg[itera-1] += step_i
                AllArgs = SetAllArgs(nArg, InArg, stepsIN)
                Val = getVal(AllArgs)
#                Val = getVal ( nArg, itera )
                Arg, mi, grad = AddPoint ( points, nArg, Val )
                print ('\nITER', itera, mi, Val, grad, 'st', step_i, nArg, '\n')
                if mi != '***' : break
                if attempt != Nattempt-1 :
                    nArg[itera - 1] -= step_i
                    if attempt == 0 :  step_i = - step_i
                    else            :  step_i *= -0.05
            steps[itera - 1] = step_i                           # 24.05
            par_der.append(grad*np.sign(step_i))
            if (itera == dim) :
                min_step = min([abs(steps[ipd] / pd) for ipd, pd in enumerate(par_der)])
                step = np.sqrt(sum((pd * min_step) ** 2 for pd in par_der))
                print('par_der', par_der, 'min_step', min_step, step)
#           step = min (step, abs(step_i))
    
    print(f"dim = {dim}, CVNumOfIter = {CVNumOfIter}, len(points) = {len(points)}")
    if CVNumOfIter > dim:
        init_x_to_y_dict.clear()
        init_x_to_y_dict.update({
            tuple(np.round(np.asarray(p.Arg).flatten(), 8)): float(p.Val)
            for p in points
        })
    

        from spotoptim import SpotOptim
        from spotoptim.plot.visualization import plot_progress, plot_surrogate

        print(f"Args начальных данных: ", [p.Arg for p in points])
        print(f"Vals начальных данных: ", [p.Val for p in points])
        bounds =  [[0.0, 0.8] for _ in range(dim)]
        X_init = np.array([np.asarray(point.Arg).flatten() for point in points]) 

        print("* Creation of SpotOPtim optimizer")
        opt = SpotOptim(
            fun=getVal,
            bounds=bounds,
            max_iter=CVNumOfIter - dim - 1 + X_init.shape[0],
            n_initial=0,
            selection_method='distant',
            n_jobs=1,
            verbose=True
        )

        # print("Начинается заполнение начальных точек в оптимизатор")
        # opt.X_ = np.atleast_2d(X_init)
        # opt.y_ = np.asarray(Y_init).flatten()
        # opt.counter = len(X_init)

        # idx_best = np.argmin(opt.y_)  
        # opt.best_x_ = opt.X_[idx_best].copy()
        # opt.best_y_ = float(opt.y_[idx_best])
        # print("Заполнение начальных точек в оптимизатор закончено. Начинается оптимизация")

        # Run optimization
        result = opt.optimize(X0=X_init)
        print(f"points = {points}, result.X = {result.X}, result.y = {result.y}")
        new_points = points + [Point(result.X[i], result.y[i], i+X_init.shape[0]) for i in range(result.X.shape[0])]

        #for p in new_points :  p.prin()
        print(f"Args финальных данных: ", [p.Arg for p in new_points])
        print(f"Vals финальных данных: ", [p.Val for p in new_points])
        #print("Added points: ", points)
        Task.OptPoints = new_points #kfe_added

        plot_progress(opt)
        if dim == 2:
            plot_surrogate(opt)
        points = new_points

#   STEP тут не настоящий! Он нужен для того, чтобы не менять интерфейс функции SurMin, а также для того, чтобы можно было использовать его в других местах, где он нужен. В данном случае, он не используется, так как оптимизация выполняется с помощью библиотеки spotoptim, которая сама определяет шаги оптимизации.
    return points, step
